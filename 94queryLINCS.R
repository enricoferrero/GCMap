library(data.table)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

# L1000CDS2 API
lincs.url <- "http://amp.pharm.mssm.edu/L1000CDS2/query"

# UniChem API
unichem.url <- "https://www.ebi.ac.uk/unichem/rest/src_compound_id/"

# read STOPGAP and split by diease/trait
stopgap <- fread("../dat/stopgap.tsv")
stopgap.list <- split(stopgap, stopgap$efo.id)

# reverse and mimic
for (search.mode in c("reverse", "mimic")) {

    if (search.mode == "reverse") {
        aggravate.value <- FALSE
    } else if (search.mode == "mimic") {
        aggravate.value <- TRUE
    }

    # loop through diseases
    lincs <- foreach(i = seq(names(stopgap.list)), .combine = rbind) %dopar% {
        
        # create genetic signature
        genes <- stopgap.list[[i]][, gene.symbol]
        vals <- stopgap.list[[i]][, stopgap.score]

        # query the L1000CDS2 API
        lincs.body <- list(data = list(genes = genes, vals = vals), config = list(aggravate = aggravate.value, searchMethod = "CD", share = FALSE, combination = FALSE, `db-version` = "latest"))
        lincs.res <- POST(url = lincs.url, body = lincs.body, encode = "json")
        lincs.res <- fromJSON(content(lincs.res, as = "text", encoding = "UTF-8"))

        if (class(lincs.res$topMeta) == "data.frame") {

            lincs.res <- as.data.table(flatten(lincs.res$topMeta))
            lincs.res <- na.omit(lincs.res[, .(pubchem.id = pubchem_id, lincs.id = pert_id, lincs.name = pert_desc, lincs.cell = cell_id, lincs.dose = pert_dose, lincs.dose.unit = pert_dose_unit, lincs.time = pert_time, lincs.time.unit = pert_time_unit, lincs.score = score)])

            if (nrow(lincs.res) > 0) {

                # loop through compounds
                unichem.res <- foreach(j = seq(lincs.res[, pubchem.id]), .combine = rbind) %dopar% {
                    # query the UniChem API to map to ChEMBL IDs
                    pubchem.id <- lincs.res[j, pubchem.id]
                    chembl.id <- as.character(fromJSON(content(GET(paste0(unichem.url, lincs.res[j, pubchem.id], "/22/1")), as = "text", encoding = "UTF-8")))
                    if (length(chembl.id > 0) && startsWith(chembl.id, "CHEMBL")) {
                        id <- data.table(pubchem.id, chembl.id)
                    }
                }

                # map L1000CDS2 output to ChEMBL and EFO
                lincs.res <- merge(unique(lincs.res), unique(unichem.res))
                lincs.res[, efo.id := names(stopgap.list)[i]]

            }

        }

    }

    # export
    fwrite(lincs, paste0("../dat/lincs.", search.mode, ".tsv"), sep ="\t")

}

