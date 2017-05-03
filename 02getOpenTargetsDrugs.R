library(data.table)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

# read STOPGAP data to get EFO IDs for diseases
stopgap <- fread("../dat/stopgap.tsv")

## query the Open Targets API to get all known drugs for each disease
# get token
token <- content(GET("https://www.targetvalidation.org/api/latest/public/auth/request_token", query = list(app_name = "personal", secret = "8r88SQNM2A2nayCd47rSWVTh7x5DSw47")))$token

# loop thorugh diseases
efo.ids <- unique(stopgap[, efo.id])
opentargets  <- foreach(i = seq(efo.ids), .combine = rbind) %dopar% {

    # query API
    tmp <- GET("https://www.targetvalidation.org/api/latest/public/evidence/filter", query = list(disease = efo.ids[i], datatype = "known_drug", size = 10000), add_headers(`Auth-Token` = token))

    # check response and request new token if necessary
    if (tmp$status_code == 419) {
        token <- content(GET("https://www.targetvalidation.org/api/latest/public/auth/request_token", query = list(app_name = "personal", secret = "8r88SQNM2A2nayCd47rSWVTh7x5DSw47")))$token
        tmp <- GET("https://www.targetvalidation.org/api/latest/public/evidence/filter", query = list(disease = efo.ids[i], datatype = "known_drug", size = 10000), add_headers(`Auth-Token` = token))
    }

    # extract and check data
    tmp <- fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
    if (length(tmp$data) > 0) {
        tmp <- unique(as.data.table(flatten(tmp$data$drug)))
        tmp[, c("efo.id", "chembl.id") := .(efo.ids[i], sub("http://identifiers.org/chembl.compound/", "", id))]
    }
}

# tidy
opentargets <- opentargets[, .(chembl.id, chembl.name = molecule_name, chembl.type = molecule_type, phase = max_phase_for_all_diseases.numeric_index, efo.id)]
fwrite(opentargets, "../dat/opentargets.tsv", sep = "\t")
