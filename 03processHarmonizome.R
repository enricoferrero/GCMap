library(data.table)
library(EnsDb.Hsapiens.v75)

# read Harmonizome LINCS L1000 data
harmonizome <- fread("gunzip -c ../dat/harmonizome_L1000_gene_attribute_edges.txt.gz", header = TRUE, skip = 1, select = c(1, 3, 4, 7))

# map Entrez IDs to Ensembl
harmonizome.genes <- harmonizome[, as.character(unique(GeneID))]
anno <- as.data.table(select(EnsDb.Hsapiens.v75, keys = harmonizome.genes, keytype = "ENTREZID", columns = c("GENEID", "ENTREZID")))
harmonizome <- merge(anno[, ENTREZID := as.numeric(ENTREZID)], harmonizome, by.x = "ENTREZID", by.y = "GeneID", allow.cartesian = TRUE)

# tidy
harmonizome <- harmonizome[, .(ensembl.id = GENEID, gene.symbol = GeneSym, perturbation = `Perturbation ID_Perturbagen_Cell Line_Time_Time Unit_Dose_Dose Unit`, direction = weight)]
fwrite(harmonizome, "../dat/harmonizome.tsv")
