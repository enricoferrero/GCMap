library(data.table)
library(EnsDb.Hsapiens.v75)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

# options
min.gene.score <- 1
max.gene.rank.min <- 3

# read data in
stopgap <- fread("../dat/stopgap.gene.mesh.txt")
# exclude rare disease (OMIM) evidence and keep relevant columns
stopgap <- unique(stopgap[source %in% c("nhgri", "gwasdb", "grasp", "paper"), .(gene.v19, msh, pvalue, gene.score, gene.rank.min)])
# replace p-values of zero with arbitrarily low p-value
stopgap[pvalue == 0, pvalue := 3e-324]

# map gene symbols to Ensembl
stopgap.genes <- unique(stopgap[, gene.v19])
anno <- as.data.table(select(EnsDb.Hsapiens.v75, keys = stopgap.genes, keytype = "SYMBOL", columns = c("GENEID", "SYMBOL")))
stopgap <- merge(stopgap, anno, by.x = "gene.v19", by.y = "SYMBOL", all = FALSE)

# map MeSH terms to EFO using Zooma 
mesh.terms <- gsub(" ", "+", unique(stopgap[, msh]))
zooma <- foreach(i = seq(mesh.terms), .combine = rbind) %dopar% {
    tmp <- fromJSON(content(GET(paste0("http://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate?propertyValue=", mesh.terms[i])), as = "text"))
    if (length(tmp) > 0) {
        tmp <- data.table(mesh.term = gsub("\\+", " ", mesh.terms[i]), efo.term = tmp$annotatedProperty$propertyValue[1], efo.url = tmp$semanticTags[[1]][1], confidence = tmp$confidence[1])
    }
}

# merge
stopgap <- merge(stopgap, zooma, by.x = "msh", by.y = "mesh.term", all = FALSE)
# filter
stopgap <- stopgap[gene.rank.min <= max.gene.rank.min & gene.score > min.gene.score]
# tidy
stopgap[, efo.id := sub(".+\\/([A-Za-z]+_[0-9]+)", "\\1", efo.url)]
stopgap <- stopgap[, .(ensembl.id = GENEID, gene.symbol = gene.v19, efo.id, efo.term, stopgap.pvalue = pvalue, stopgap.gene.score = gene.score, stopgap.gene.rank = gene.rank.min)]
fwrite(stopgap, "../dat/stopgap.tsv", sep = "\t")
