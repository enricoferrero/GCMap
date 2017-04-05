library(data.table)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

# L1000CDS2 API configuration
lincs.url <- "http://amp.pharm.mssm.edu/L1000CDS2/query"
lincs.config <- list(aggravate = FALSE, searchMethod = "CD", share = FALSE, combination = FALSE, `db-version` = "latest")

# UniChem API configuration
unichem.url <- "https://www.ebi.ac.uk/unichem/rest/src_compound_id/"

# read STOPGAP and split by diease/trait
stopgap <- fread("../dat/stopgap.tsv")
stopgap <- split(stopgap, stopgap$efo.id)

# read Open Targets and split by disease/trait
opentargets <- fread("../dat/opentargets.tsv")
opentargets <- split(opentargets, opentargets$efo.id)

# loop through diseases
lincs <- foreach(i = seq(names(stopgap)), .combine = rbind) %dopar% {
#for (i in seq(names(stopgap))) {
    
    # create genetic signature
    genes <- stopgap[[i]][, gene.symbol]
    vals <- stopgap[[i]][, score]

    # query the L1000CDS2 API
    lincs.body <- list(data = list(genes = genes, vals = vals), config = lincs.config)
    lincs.res <- POST(url = lincs.url, body = lincs.body, encode = "json")
    lincs.res <- fromJSON(content(lincs.res, as = "text"))

    if (class(lincs.res$topNeta) == "data.frame") {

        lincs.res <- as.data.table(flatten(lincs.res$topMeta))
        lincs.res <- na.omit(lincs.res[, .(pubchem_id, pert_id, pert_desc, score, pert_dose, pert_time, cell_id)])

        # loop through compounds
        unichem.res <- foreach(j = seq(lincs.res[, pubchem_id]), .combine = rbind) %dopar% {
            # query the UniChem API
            pubchem_id <- lincs.res[j, pubchem_id]
            chembl.id <- as.character(fromJSON(content(GET(paste0(unichem.url, lincs.res[j, pubchem_id], "/22/1")), as = "text", encoding = "UTF-8")))
            if (length(chembl.id) > 0) {
                id <- data.table(pubchem_id, chembl.id)
            }
        }

        # map L1000CDS2 output to ChEMBL and EFO
        lincs.res <- merge(lincs.res, unichem.res)
        lincs.res[, efo.id := names(stopgap)[i]]
    }
}

# export LINCS and split by disease/trait
fwrite(lincs, "../dat/lincs.tsv", sep ="\t")
lincs <- split(lincs, lincs$efo.id)
