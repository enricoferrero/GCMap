library(data.table)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

# options
min.drugs.score <- 0.5
min.degs.score <- 0.2

# read STOPGAP data to get EFO IDs for diseases
stopgap <- fread("../dat/stopgap.tsv")

## query the Open Targets API to get all known drugs for each disease
# get token
token <- content(GET("https://www.targetvalidation.org/api/latest/public/auth/request_token", query = list(app_name = "personal", secret = "8r88SQNM2A2nayCd47rSWVTh7x5DSw47")))$token

# diseases to process
efo.ids <- unique(stopgap[, efo.id])

## drugs
# loop thorugh diseases
opentargets.drugs  <- foreach(i = seq(efo.ids), .combine = rbind) %dopar% {

    # query API
    tmp <- GET("https://www.targetvalidation.org/api/latest/public/evidence/filter", query = list(disease = efo.ids[i], datatype = "known_drug", scorevalue_min = min.drugs.score, size = 10000), add_headers(`Auth-Token` = token))

    # check response and request new token if necessary
    if (tmp$status_code == 419) {
        token <- content(GET("https://www.targetvalidation.org/api/latest/public/auth/request_token", query = list(app_name = "personal", secret = "8r88SQNM2A2nayCd47rSWVTh7x5DSw47")))$token
        tmp <- GET("https://www.targetvalidation.org/api/latest/public/evidence/filter", query = list(disease = efo.ids[i], datatype = "known_drug", scorevalue_min = min.drugs.score, size = 10000), add_headers(`Auth-Token` = token))
    }

    # extract and check data
    tmp <- fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
    if (length(tmp$data) > 0) {
        tmp <- unique(as.data.table(flatten(tmp$data)))
        tmp[, chembl.id := sub("http://identifiers.org/chembl.compound/", "", drug.id)]
        tmp <- tmp[, .(efo.id = disease.id, efo.term = disease.efo_info.label, chembl.id, chembl.name = drug.molecule_name, chembl.type = drug.molecule_type, phase = drug.max_phase_for_all_diseases.numeric_index)]
    }
}

# export
fwrite(opentargets.drugs, "../dat/opentargets.drugs.tsv", sep = "\t")

## DEGs
# loop thorugh diseases
opentargets.degs <- foreach(i = seq(efo.ids), .combine = rbind) %dopar% {
    
    # query API
    tmp <- GET("https://www.targetvalidation.org/api/latest/public/evidence/filter", query = list(disease = efo.ids[i], datatype = "rna_expression", scorevalue_min = min.degs.score, size = 10000), add_headers(`Auth-Token` = token))

    # check response and request new token if necessary
    if (tmp$status_code == 419) {
        token <- content(GET("https://www.targetvalidation.org/api/latest/public/auth/request_token", query = list(app_name = "personal", secret = "8r88SQNM2A2nayCd47rSWVTh7x5DSw47")))$token
        tmp <- GET("https://www.targetvalidation.org/api/latest/public/evidence/filter", query = list(disease = efo.ids[i], datatype = "rna_expression", scorevalue_min = min.degs.score, size = 10000), add_headers(`Auth-Token` = token))
    }

    # extract and check data
    tmp <- fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
    if (length(tmp$data) > 0) {
        tmp <- unique(as.data.table(flatten(tmp$data)))
        tmp <- tmp[, .(efo.id = disease.id, efo.term = disease.efo_info.label, ensembl.id = target.id, gene.symbol = target.gene_info.symbol, comparison = evidence.comparison_name, lfc = evidence.log2_fold_change.value, pval = evidence.resource_score.value)]
    }
}

# export
fwrite(opentargets.degs, "../dat/opentargets.degs.tsv", sep = "\t")
