library(EnsDb.Hsapiens.v75)
library(data.table)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

# read Harmonizome LINCS L1000 data
harmonizome <- fread("gunzip -c ../dat/harmonizome_L1000_gene_attribute_edges.txt.gz", header = TRUE, skip = 1, select = c(1, 3, 4, 7))
# read LINCS to PubChem mappings
lincs.mappings <- fread("../dat/LINCS_LSM_Pubchem_ID_mappings.tsv")

# map Entrez IDs to Ensembl
harmonizome.genes <- harmonizome[, as.character(unique(GeneID))]
anno <- as.data.table(select(EnsDb.Hsapiens.v75, keys = harmonizome.genes, keytype = "ENTREZID", columns = c("GENEID", "ENTREZID")))
harmonizome <- merge(anno[, ENTREZID := as.numeric(ENTREZID)], harmonizome, by.x = "ENTREZID", by.y = "GeneID", allow.cartesian = TRUE)

# map LINCS IDs to ChEMBL IDs
unichem.url <- "https://www.ebi.ac.uk/unichem/rest/src_compound_id/"
unichem.res <- foreach(i = seq(lincs.mappings[, pubchem_cid]), .combine = rbind) %dopar% {
    lincs.id <- lincs.mappings[i, pert_id]
    pubchem.id <- lincs.mappings[i, pubchem_cid]
    chembl.id <- as.character(fromJSON(content(GET(paste0(unichem.url, lincs.mappings[i, pubchem_cid], "/22/1")), as = "text", encoding = "UTF-8")))
    if (length(chembl.id > 0) && startsWith(chembl.id, "CHEMBL")) {
        tmp <- data.table(lincs.id, pubchem.id, chembl.id)
    }
}

# annotate and tidy
harmonizome[, lincs.id := substr(`Perturbation ID_Perturbagen_Cell Line_Time_Time Unit_Dose_Dose Unit`, 1, 13)]
harmonizome <- merge(harmonizome, unichem.res, by = "lincs.id")
harmonizome <- harmonizome[, .(ensembl.id = GENEID, gene.symbol = GeneSym, lincs.id, pubchem.id, chembl.id, perturbation = `Perturbation ID_Perturbagen_Cell Line_Time_Time Unit_Dose_Dose Unit`, direction = weight)]
fwrite(harmonizome, "../dat/harmonizome.tsv")
