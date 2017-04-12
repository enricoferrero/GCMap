library(data.table)
library(EnsDb.Hsapiens.v75)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

# read data in
stopgap <- fread("../dat/stopgap.gene.mesh.txt")
# keep relevant columns
stopgap <- unique(stopgap[, .(gene.v19, msh, pvalue, gene.score, gene.rank.min)])
# replace p-values of zero with arbitrarily low p-value
stopgap[pvalue == 0, pvalue := 3e-324]
# create score
stopgap[, score := -log10(pvalue) * gene.score / gene.rank.min]

# map gene symbols to Ensembl
stopgap.genes <- unique(stopgap[, gene.v19])
anno <- as.data.table(select(EnsDb.Hsapiens.v75, keys = stopgap.genes, keytype = "SYMBOL", columns = c("GENEID", "SYMBOL")))
stopgap <- merge(stopgap, anno, by.x = "gene.v19", by.y = "SYMBOL", all = FALSE)

# map MeSH terms to EFO using Zooma 
mesh.terms <- gsub(" ", "+", unique(stopgap[, msh]))
zooma.results <- foreach(i = seq(mesh.terms), .combine = rbind) %dopar% {
    tmp <- fromJSON(content(GET(paste0("http://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate?propertyValue=", mesh.terms[i])), as = "text"))
    if (length(tmp) > 0) {
        tmp <- data.table(mesh.term = gsub("\\+", " ", mesh.terms[i]), efo.term = tmp$annotatedProperty$propertyValue[1], efo.url = tmp$semanticTags[[1]][1], confidence = tmp$confidence[1])
    }
}

# manual curation of Zooma results
fwrite(zooma.results, "../dat/zooma.results.tsv", sep = "\t")
zooma <- fread("../dat/zooma.cleaned.results.tsv")
# merge
stopgap <- merge(stopgap, zooma, by.x = "msh", by.y = "mesh.term", all = FALSE)

# tidy
stopgap[, efo.id := sub(".+\\/([A-Za-z]+_[0-9]+)", "\\1", efo.url)]
stopgap <- stopgap[, .(ensembl.id = GENEID, gene.symbol = gene.v19, efo.id, efo.term, stopgap.score = score, stopgap.pvalue = pvalue, stopgap.gene.score = gene.score, stopgap.gene.rank = gene.rank.min)]
fwrite(stopgap, "../dat/stopgap.tsv", sep = "\t")
