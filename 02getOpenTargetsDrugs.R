library(data.table)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

# read STOPGAP data to get EFO IDs for diseases
stopgap <- fread("../dat/stopgap.tsv")

# request token for Open Targets API
token <- content(GET("https://www.targetvalidation.org/api/latest/public/auth/request_token", query = list(app_name = "personal", secret = "8r88SQNM2A2nayCd47rSWVTh7x5DSw47")))$token

# query the Open Targets API to get all known drugs for each disease
efo.ids <- unique(stopgap[, efo.id])
opentargets  <- foreach(i = seq(efo.ids), .combine = rbind) %dopar% {
    tmp <- fromJSON(content(GET("https://www.targetvalidation.org/api/latest/public/evidence/filter", query = list(disease = efo.ids[i], datatype = "known_drug", scorevalue_min = 0.2, size = 10000), add_headers(`Auth-Token` = token)), type = "text", encoding = "UTF-8"))
    if (class(tmp$data$drug) == "data.frame") {
        tmp <- unique(as.data.table(flatten(tmp$data$drug)))
        tmp[, c("efo.id", "chembl.id") := .(efo.ids[i], sub("http://identifiers.org/chembl.compound/", "", id))]
    }
}

# tidy
opentargets <- opentargets[, .(chembl.id, chembl.name = molecule_name, chembl.type = molecule_type, phase = max_phase_for_all_diseases.numeric_index, efo.id)]
fwrite(opentargets, "../dat/opentargets.tsv", sep = "\t")
