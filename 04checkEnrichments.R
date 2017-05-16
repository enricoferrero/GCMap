library(EnsDb.Hsapiens.v75)
library(data.table)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)
library(BiocParallel)
register(MulticoreParam(workers = 1))
library(fgsea)
library(ggplot2)
library(pROC)

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
opentargets.drugs <- fread("../dat/opentargets.drugs.tsv")
opentargets.degs <- fread("../dat/opentargets.degs.tsv")
harmonizome <- fread("../dat/harmonizome.tsv")

# split by disease
stopgap.list <- split(stopgap, stopgap$efo.id)
opentargets.degs.list <- split(opentargets.degs, opentargets.degs$efo.id)
#opentargets.drugs.list <- split(opentargets.drugs, opentargets.drugs$efo.id)

# create universes for Fisher's test
ensembl.ids <- keys(EnsDb.Hsapiens.v75)
efo.ids <- unique(c(stopgap[, efo.id], opentargets.degs[, efo.id], opentargets.drugs[, efo.id]))
chembl.ids <- unique(opentargets.drugs[, chembl.id])

## create lists for GSEA
#stopgap.list.gsea <- sapply(stopgap.list, "[[", "ensembl.id")
#opentargets.degs.list.gsea <- sapply(opentargets.degs.list, function(x) setNames(x$lfc, x$ensembl.id))

## number of GSEA permutations
#n.perm = 10000

### question 1: are genes differentially expressed in disease X enriched in genes genetically associated with disease X, compared to other diseases?
## perform Fisher's test
# loop through diseases
fisher.genes <- foreach (i = seq(opentargets.degs.list), .combine = rbind, .errorhandling = "remove") %do% {
    efo.id.opentargets.degs <- names(opentargets.degs.list)[i]
    foreach (j = seq(stopgap.list), .combine = rbind, .errorhandling = "remove") %dopar% {
        efo.id.stopgap <- names(stopgap.list)[j]
        # test significance
        tmp <- checkOverlapSignificance(opentargets.degs.list[[efo.id.opentargets.degs]]$ensembl.id, stopgap.list[[efo.id.stopgap]]$ensembl.id, ensembl.ids)
        # annotate
        tmp <- cbind(efo.id.opentargets.degs, efo.id.stopgap, tmp)
        tmp <- merge(tmp, unique(stopgap[, .(efo.id, efo.term)]), by.x = "efo.id.stopgap", by.y = "efo.id", all.x = TRUE, all.y = FALSE)
        tmp <- merge(tmp, unique(opentargets.degs[, .(efo.id, efo.term)]), by.x = "efo.id.opentargets.degs", by.y = "efo.id", all.x = TRUE, all.y = FALSE, suffixes = c(".stopgap", ".opentargets.degs"))
        setcolorder(tmp, c("efo.id.opentargets.degs", "efo.term.opentargets.degs", "efo.id.stopgap", "efo.term.stopgap", "set1", "set2", "overlap", "universe", "odds.ratio", "p.value"))
    }
}
fisher.genes <- unique(fisher.genes)

# correct p-values
fisher.genes[p.value == 0, p.value := 3e-324]
fisher.genes <- fisher.genes[order(p.value), ]
fisher.genes[, p.adjusted := p.adjust(p.value, method = "fdr")]

# add dummy variable
fisher.genes[, same.disease := ifelse(efo.id.opentargets.degs == efo.id.stopgap, TRUE, FALSE)]

# plot adjusted p-value distributions
png("../dat/fisher.genes.boxplots.png", width = 10 * 150, height = 5 * 150, res = 150)
print(ggplot(fisher.genes, aes(x = same.disease, y = -log10(p.adjusted))) +
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = quantile(fisher.genes[, -log10(p.adjusted)], c(0.020, 0.980))) +
    xlab("Same disease") +
    ylab("-log10(adjusted p-value)") +
    theme_bw(14) +
    scale_x_discrete(breaks = c(FALSE, TRUE), labels = c("No", "Yes")))
dev.off()

# test if adjusted p-value  distributions significantly different with Mann-Whitney test (Wilcoxon rank sum test)
mw.res <- wilcox.test(x = fisher.genes[same.disease == FALSE, -log10(p.adjusted)], y = fisher.genes[same.disease == TRUE, -log10(p.adjusted)])
print(mw.res)
print(mw.res$p.value)

# perform ROC analysis
resp <- as.numeric(fisher.genes[, same.disease])
pred <- ifelse(fisher.genes[, p.adjusted] < 0.05, 1, 0)
roc.res <- roc(response = resp, predictor = pred)
print(roc.res$auc)
png("../dat/fisher.genes.roc.png", width = 6 * 150, height = 6 * 150, res = 150)
plot(roc.res, main = paste("AUC:", round(roc.res$auc, 2)), xlab = "False positive rate", ylab = "True positive rate", col = "firebrick", lwd = 2, identity.lty = 2)
dev.off()

# export
fwrite(fisher.genes, "../dat/fisher.genes.tsv", sep = "\t")

### perform GSEA
## loop through diseases
#gsea.genes <- foreach (i = seq(opentargets.degs.list.gsea), .combine = rbind, .errorhandling = "remove") %do% {
#    efo.id.opentargets.degs <- names(opentargets.degs.list.gsea)[i]
#    foreach (j = seq(stopgap.list.gsea), .combine = rbind, .errorhandling = "remove") %dopar% {
#        efo.id.stopgap <- names(stopgap.list.gsea)[j]
#        # test significance
#        tmp <- fgsea(pathways = stopgap.list.gsea[efo.id.stopgap], stats = opentargets.degs.list.gsea[[efo.id.opentargets.degs]], minSize = 5, maxSize = 1000, nperm = n.perm)
#        # annotate
#        if (nrow(tmp) > 0) {
#            tmp <- tmp[, .(efo.id.stopgap = pathway, p.value = pval, p.adjusted = padj, es = ES, nes = NES, n.more.extreme = nMoreExtreme, size)]
#            tmp <- cbind(efo.id.opentargets.degs, tmp)
#            tmp <- merge(tmp, unique(stopgap[, .(efo.id, efo.term)]), by.x = "efo.id.stopgap", by.y = "efo.id", all.x = TRUE, all.y = FALSE)
#            tmp <- merge(tmp, unique(opentargets.degs[, .(efo.id, efo.term)]), by.x = "efo.id.opentargets.degs", by.y = "efo.id", all.x = TRUE, all.y = FALSE, suffixes = c(".stopgap", ".opentargets.degs"))
#            setcolorder(tmp, c("efo.id.opentargets.degs", "efo.term.opentargets.degs", "efo.id.stopgap", "efo.term.stopgap", "p.value", "p.adjusted", "es", "nes", "n.more.extreme", "size"))
#        }
#    }
#}
#gsea.genes <- unique(gsea.genes)
#
## correct p-values
#gsea.genes[p.value == 0, p.value := 3e-324]
#gsea.genes <- gsea.genes[order(p.value), ]
#gsea.genes[, p.adjusted := p.adjust(p.value, method = "fdr")]
#
## add dummy variable
#gsea.genes[, same.disease := ifelse(efo.id.opentargets.degs == efo.id.stopgap, TRUE, FALSE)]
#
## plot adjusted p-value distributions
#ggplot(gsea.genes, aes(x = same.disease, y = -log10(p.adjusted))) +
#    geom_boxplot(outlier.shape = NA) +
#    coord_cartesian(ylim = quantile(gsea.genes[, -log10(p.adjusted)], c(0.025, 0.975))) +
#    xlab("Disease") +
#    ylab("-log10(adjusted p-value)") +
#    theme_bw(14) +
#    scale_x_discrete(breaks = c(FALSE, TRUE), labels = c("Different", "Same"))
#ggsave("../dat/gsea.genes.png", width = 20, height = 10, unit = "cm", dpi = 150)
#
## test if adjusted p-value  distributions significantly different with Mann-Whitney test (Wilcoxon rank sum test)
#mw.res <- wilcox.test(x = gsea.genes[same.disease == FALSE, -log10(p.adjusted)], y = gsea.genes[same.disease == TRUE, -log10(p.adjusted)])
#print(mw.res)
#print(mw.res$p.value)
#
## perform ROC analysis
#resp <- as.numeric(gsea.genes[, same.disease])
#pred <- ifelse(gsea.genes[, p.adjusted] < 0.05, 1, 0)
#roc.res <- roc(response = resp, predictor = pred)
#print(roc.res$auc)
#png("../dat/gsea.genes.roc.png", width = 6 * 150, height = 6 * 150, res = 150)
#plot(roc.res, main = paste("AUC:", round(roc.res$auc, 2)), xlab = "False positive rate", ylab = "True positive rate", col = "firebrick", lwd = 2, identity.lty = 2)
#dev.off()
#
## export
#fwrite(gsea.genes, "../dat/gsea.genes.tsv", sep = "\t")

### question 2: are genes differentially expressed after treatment with drug for disease X enriched for genes genetically associated with disease X, compared to other diseases?

## perform Fisher's test
# loop through diseases
fisher.drugs <- foreach (i = seq(efo.ids), .combine = rbind, .errorhandling = "remove") %dopar% {
    this.efo.id <- efo.ids[i]
    # loop through drugs
    foreach (j = seq(chembl.ids), .combine = rbind, .errorhandling = "remove") %do% {
             this.chembl.id <- chembl.ids[j]
             degs <- unique(harmonizome[chembl.id == this.chembl.id, ensembl.id])
             gags <- stopgap.list[[this.efo.id]]$ensembl.id
             if (length(degs) > 0 && length(gags) > 0) {
                # test significance
                tmp <- checkOverlapSignificance(degs, gags, ensembl.ids)
                # annotate
                tmp <- cbind(chembl.id = this.chembl.id, efo.id = this.efo.id, tmp)
                tmp <- merge(unique(opentargets.drugs[, .(efo.id, efo.term)]), tmp, by = "efo.id")
                tmp <- merge(unique(opentargets.drugs[, .(chembl.id, chembl.name)]), tmp, by = "chembl.id")
                # add dummy variable
                tmp[, existing.indication := ifelse(nrow(opentargets.drugs[efo.id == this.efo.id & chembl.id == this.chembl.id]) > 0, TRUE, FALSE)]
                setcolorder(tmp, c("chembl.id", "chembl.name", "efo.id", "efo.term", "set1", "set2", "overlap", "universe", "odds.ratio", "p.value", "existing.indication"))
             }
    }
}
fisher.drugs <- unique(fisher.drugs)

# correct p-values
fisher.drugs[p.value == 0, p.value := 3e-324]
fisher.drugs <- fisher.drugs[order(p.value), ]
fisher.drugs[, p.adjusted := p.adjust(p.value, method = "fdr")]

# plot adjusted p-value distributions
png("../dat/fisher.drugs.boxplots.png", width = 10 * 150, height = 5 * 150, res = 150)
print(ggplot(fisher.drugs, aes(x = existing.indication, y = -log10(p.adjusted))) +
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = quantile(fisher.drugs[, -log10(p.adjusted)], c(0.020, 0.980))) +
    xlab("Existing indication") +
    ylab("-log10(adjusted p-value)") +
    theme_bw(14) +
    scale_x_discrete(breaks = c(FALSE, TRUE), labels = c("No", "Yes")))
dev.off()

# test if adjusted p-value  distributions significantly different with Mann-Whitney test (Wilcoxon rank sum test)
mw.res <- wilcox.test(x = fisher.drugs[existing.indication == FALSE, -log10(p.adjusted)], y = fisher.drugs[existing.indication == TRUE, -log10(p.adjusted)])
print(mw.res)
print(mw.res$p.value)

# perform ROC analysis
resp <- as.numeric(fisher.drugs[, existing.indication])
pred <- ifelse(fisher.drugs[, p.adjusted] < 0.05, 1, 0)
roc.res <- roc(response = resp, predictor = pred)
print(roc.res$auc)
png("../dat/fisher.drugs.roc.png", width = 6 * 150, height = 6 * 150, res = 150)
plot(roc.res, main = paste("AUC:", round(roc.res$auc, 2)), xlab = "False positive rate", ylab = "True positive rate", col = "firebrick", lwd = 2, identity.lty = 2)
dev.off()

# export
fwrite(fisher.drugs, "../dat/fisher.drugs.tsv", sep = "\t")

### perform GSEA
## loop through diseases
#gsea.drugs <- foreach (i = seq(opentargets.drugs.list), .combine = rbind, .errorhandling = "remove") %do% {
#    efo.id <- names(opentargets.drugs.list)[i]
#    # loop through drugs
#    foreach (j = 1:nrow(opentargets.drugs.list[[i]]), .combine = rbind, .errorhandling = "remove") %dopar% {
#        this.chembl.id <- opentargets.drugs.list[[i]][j, chembl.id]
#        degs <- unique(harmonizome[chembl.id == this.chembl.id, .(ensembl.id, direction)])
#        degs <- setNames(degs$direction, degs$ensembl.id)
#            if (length(degs) > 0) {
#                # test significance
#                tmp <- fgsea(pathways = stopgap.list.gsea[efo.id], stats = degs, minSize = 5, maxSize = 1000, nperm = n.perm)
#                if (nrow(tmp) > 0) {
#                    # annotate
#                    tmp <- tmp[, .(efo.id = pathway, p.value = pval, p.adjusted = padj, es = ES, nes = NES, n.more.extreme = nMoreExtreme, size)]
#                    tmp <- cbind(chembl.id = this.chembl.id, tmp)
#                    tmp <- merge(unique(opentargets.drugs[, .(efo.id, efo.term, chembl.id, chembl.name)]), tmp)
#                    setcolorder(tmp, c("chembl.id", "chembl.name", "efo.id", "efo.term", "p.value", "p.adjusted", "es", "nes", "n.more.extreme", "size"))
#                }
#            }
#    }
#}
#gsea.drugs <- unique(gsea.drugs)
#
## correct p-values
#gsea.drugs[p.value == 0, p.value := 3e-324]
#gsea.drugs <- gsea.drugs[order(p.value), ]
#gsea.drugs[, p.adjusted := p.adjust(p.value, method = "fdr")]
#
## export
#fwrite(gsea.drugs, "../dat/gsea.drugs.tsv", sep = "\t")
