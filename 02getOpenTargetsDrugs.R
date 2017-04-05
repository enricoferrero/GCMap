library(data.table)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

# read STOPGAP data to get EFO IDs for diseases
stopgap <- fread("../dat/stopgap.tsv")

# query the Open Targets API to get all known drugs for each disease
efo.ids <- unique(stopgap[, efo.id])
opentargets  <- foreach(i = seq(efo.ids), .combine = rbind) %dopar% {
    tmp <- fromJSON(content(GET(paste0("https://www.targetvalidation.org/api/latest/public/evidence/filter?disease=", efo.ids[i], "&datatype=known_drug&scorevalue_min=0.2&size=10000")), type = "text", encoding = "UTF-8"))
    if (class(tmp$data$drug) == "data.frame") {
        tmp <- unique(as.data.table(flatten(tmp$data$drug)))
        tmp[, c("efo.id", "chembl.id") := .(efo.ids[i], sub("http://identifiers.org/chembl.compound/", "", id))]
    }
}
opentargets <- opentargets[!sapply(opentargets, is.null)]
fwrite(opentargets, "../dat/opentargets.tsv", sep = "\t")
