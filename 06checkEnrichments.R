library(EnsDb.Hsapiens.v75)
library(data.table)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)
library(ggplot2)

# function to calculate enrichment
checkOverlapSignificance <- function(set1, set2, universe) {
    set1 <- unique(set1)
    set2 <- unique(set2)
    a <- sum(set1 %in% set2)
    b <- length(set1) - a
    c <- length(set2) - a
    d <- length(universe) - a - b - c
    fisher <- fisher.test(matrix(c(a,b,c,d), nrow=2, ncol=2, byrow=TRUE), alternative="greater")
    data.table(set1 = length(set1), set2 = length(set2), overlap = a, universe = length(universe), odds.ratio = fisher$estimate, p.value = fisher$p.value)
}

# read datasets
stopgap <- fread("../dat/stopgap.tsv")
opentargets <- fread("../dat/opentargets.tsv")
harmonizome <- fread("../dat/harmonizome.tsv")

# get token for Open Targets API
token <- content(GET("https://www.targetvalidation.org/api/latest/public/auth/request_token", query = list(app_name = "personal", secret = "8r88SQNM2A2nayCd47rSWVTh7x5DSw47")))$token

### question 1: are genes differentially expressed in disease X enriched for genes genetically associated with disease X?
# split by disease
stopgap.list <- split(stopgap, stopgap$efo.id)

# loop thorugh diseases
efo.ids <- names(stopgap.list)
opentargets <- foreach(i = seq(efo.ids), .combine = rbind) %dopar% {
    
    # query API
    tmp <- GET("https://www.targetvalidation.org/api/latest/public/evidence/filter", query = list(disease = efo.ids[i], datatype = "rna_expression", size = 10000), add_headers(`Auth-Token` = token))

    # check response and request new token if necessary
    if (tmp$status_code == 419) {
        token <- content(GET("https://www.targetvalidation.org/api/latest/public/auth/request_token", query = list(app_name = "personal", secret = "8r88SQNM2A2nayCd47rSWVTh7x5DSw47")))$token
        tmp <- GET("https://www.targetvalidation.org/api/latest/public/evidence/filter", query = list(disease = efo.ids[i], datatype = "rna_expression", size = 10000), add_headers(`Auth-Token` = token))
    }

    # extract and check data
    tmp <- fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
    if (length(tmp$data) > 0) {
        tmp <- unique(as.data.table(flatten(tmp$data)))
        tmp <- tmp[, .(efo.id = disease.id, efo.term = disease.efo_info.label, ensembl.id = target.id, gene.symbol = target.gene_info.symbol, comparison = evidence.comparison_name, lfc = evidence.log2_fold_change.value, pval = evidence.resource_score.value)]
    }
}

# split by disease
opentargets.list <- split(opentargets, opentargets$efo.id)

# create gene universe
ensembl.ids <- keys(EnsDb.Hsapiens.v75)

## perform Fisher's test
# loop through diseases
fisher.genes <- foreach (i = seq(opentargets.list), .combine = rbind, .errorhandling = "remove") %dopar% {
    efo.id <- names(opentargets.list)[i]
    tmp <- checkOverlapSignificance(opentargets.list[[efo.id]]$ensembl.id, stopgap.list[[efo.id]]$ensembl.id, ensembl.ids)
    # annotate
    tmp <- cbind(efo.id, tmp)
    merge(unique(opentargets[, .(efo.id, efo.term)]), tmp, all.x = FALSE, all.y = TRUE)
}

# correct p-values
fisher.genes <- fisher.genes[order(p.value), ]
fisher.genes[, padj := p.adjust(p.value, method = "fdr")]

# check results 
fisher.genes[, sum(padj < 0.05) / .N]
fisher.genes[padj < 0.05, ]

# export
fisher.genes <- fisher.genes[padj < 0.05, ]
fwrite(fisher.genes, "../dat/fisher.genes.tsv", sep = "\t")


### question 2: are genes differentially expressed in disease X enriched in genes genetically associated with disease X, compared to other diseases?
## perform Fisher's test
# loop through diseases
fisher.genes.diseases <- foreach (i = seq(opentargets.list), .combine = rbind, .errorhandling = "remove") %do% {
    efo.id.opentargets <- names(opentargets.list)[i]
    foreach (j = seq(stopgap.list), .combine = rbind, .errorhandling = "remove") %dopar% {
        efo.id.stopgap <- names(stopgap.list)[j]
        tmp <- checkOverlapSignificance(opentargets.list[[efo.id.opentargets]]$ensembl.id, stopgap.list[[efo.id.stopgap]]$ensembl.id, ensembl.ids)
    # annotate
    tmp <- cbind(efo.id.opentargets, efo.id.stopgap, tmp)
    tmp <- merge(tmp, unique(stopgap[, .(efo.id, efo.term)]), by.x = "efo.id.stopgap", by.y = "efo.id", all.x = TRUE, all.y = FALSE)
    tmp <- merge(tmp, unique(opentargets[, .(efo.id, efo.term)]), by.x = "efo.id.opentargets", by.y = "efo.id", all.x = TRUE, all.y = FALSE, suffixes = c(".stopgap", ".opentargets"))
    setcolorder(tmp, c("efo.id.opentargets", "efo.term.opentargets", "efo.id.stopgap", "efo.term.stopgap", "set1", "set2", "overlap", "universe", "odds.ratio", "p.value"))
    }
}

# export
fwrite(fisher.genes.diseases, "../dat/fisher.genes.diseases.tsv", sep = "\t")

# remove zero and infinite values
fisher.genes.diseases[p.value == 0, p.value := 3e-324]

# correct p-values
fisher.genes.diseases <- fisher.genes.diseases[order(efo.id.opentargets, p.value), ]
fisher.genes.diseases[, padj := p.adjust(p.value, method = "fdr"), by = efo.id.opentargets]

# add dummy variable
fisher.genes.diseases[, same.efo := ifelse(efo.id.opentargets == efo.id.stopgap, TRUE, FALSE)]

# plot adjusted p-value distributions
ggplot(fisher.genes.diseases, aes(x = same.efo, y = -log10(padj))) +
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = quantile(fisher.genes.diseases[, -log10(padj)], c(0.025, 0.975))) +
    xlab("Disease") +
    ylab("-log10(adjusted p-value)") +
    theme_bw(14) +
    scale_x_discrete(breaks = c(FALSE, TRUE), labels = c("Different", "Same"))
ggsave("../dat/fisher.genes.diseases.png", width = 20, height = 10, unit = "cm", dpi = 150)

# test if adjusted p-value  distributions significantly different with Mann-Whitney test (Wilcoxon rank sum test)
res <- wilcox.test(x = fisher.genes.diseases[same.efo == FALSE, -log10(padj)], y = fisher.genes.diseases[same.efo == TRUE, -log10(padj)])
print(res)
print(res$p.value)

### question 3: are genes differentially expressed after treatment with drug for disease X enriched for genes genetically associated with disease X?
# TODO: need a way to convert BRD IDs to ChEMBL
