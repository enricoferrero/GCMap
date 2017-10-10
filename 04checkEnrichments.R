library(EnsDb.Hsapiens.v75)
library(data.table)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)
library(BiocParallel)
register(MulticoreParam(workers = 1))
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
stopgap <- unique(fread("../dat/stopgap.tsv"))
opentargets.drugs <- unique(fread("../dat/opentargets.drugs.tsv"))
opentargets.degs <- unique(fread("../dat/opentargets.degs.tsv"))
harmonizome <- unique(fread("../dat/harmonizome.tsv"))

# split by disease
stopgap.list <- split(stopgap, stopgap$efo.id)
opentargets.degs.list <- split(opentargets.degs, opentargets.degs$efo.id)

# create universes for Fisher's test
ensembl.ids <- keys(EnsDb.Hsapiens.v75)
efo.ids <- unique(c(stopgap[, efo.id], opentargets.degs[, efo.id], opentargets.drugs[, efo.id]))
chembl.ids <- unique(opentargets.drugs[, chembl.id])

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

## test if distributions significantly different with Mann-Whitney test (Wilcoxon rank sum test)
# adjusted p-values
mw.res <- wilcox.test(x = fisher.genes[same.disease == FALSE, -log10(p.adjusted)], y = fisher.genes[same.disease == TRUE, -log10(p.adjusted)])
print(mw.res)
print(mw.res$p.value)
png("../dat/fisher.genes.pvalues.boxplots.png", width = 6 * 150, height = 6 * 150, res = 150)
print(ggplot(fisher.genes, aes(x = same.disease, y = -log10(p.adjusted))) +
    geom_boxplot(fill = "#0066ff", outlier.shape = NA) +
    coord_cartesian(ylim = quantile(fisher.genes[, -log10(p.adjusted)], c(0.03, 0.97))) +
    xlab("Same disease") +
    ylab("-log10(adjusted p-value)") +
    theme_bw(18) +
    scale_x_discrete(breaks = c(FALSE, TRUE), labels = c("No", "Yes")) +
    ggtitle(paste("p-value =", sprintf("%.2e", mw.res$p.value))))
dev.off()
# odds ratios
mw.res <- wilcox.test(x = fisher.genes[same.disease == FALSE, odds.ratio], y = fisher.genes[same.disease == TRUE, odds.ratio]) 
print(mw.res)
print(mw.res$p.value)
png("../dat/fisher.genes.oddsratios.boxplots.png", width = 6 * 150, height = 6 * 150, res = 150)
print(ggplot(fisher.genes, aes(x = same.disease, y = odds.ratio)) +
    geom_boxplot(fill = "#0066ff", outlier.shape = NA) +
    coord_cartesian(ylim = quantile(fisher.genes[, odds.ratio], c(0.06, 0.94))) +
    xlab("Same disease") +
    ylab("Odds ratio") +
    theme_bw(18) +
    scale_x_discrete(breaks = c(FALSE, TRUE), labels = c("No", "Yes")) +
    ggtitle(paste("p-value =", sprintf("%.2e", mw.res$p.value))))
dev.off()

# perform ROC analysis
preds <- fisher.genes[, -log10(p.adjusted)]
labls <- as.numeric(fisher.genes[, same.disease])
roc.res <- roc(response = labls, predictor = preds, algorithm = 2, ci = TRUE, ci.method = "bootstrap", boot.n = 1000, parallel = TRUE, progress = "none")
print(roc.res)
sp.ci <- ci.sp(roc.res, sensitivities = seq(0, 1, 0.05), boot.n = 1000, parallel = TRUE, progress = "none")
se.ci <- ci.se(roc.res, specifities = seq(0, 1, 0.05), boot.n = 1000, parallel = TRUE, progress = "none")
png("../dat/fisher.genes.roc.png", width = 6 * 150, height = 6 * 150, res = 150)
par(pty = "s")
plot(roc.res, main = paste("AUC =", round(roc.res$auc, 2)), xlab = "False positive rate", ylab = "True positive rate", identity.lty = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, cex = 1.5)
plot(se.ci, type = "shape", col = "lightgrey", border = NA, no.roc = TRUE)
plot(sp.ci, type = "shape", col = "lightgrey", border = NA, no.roc = TRUE)
plot(roc.res, add = TRUE, col = "#0066ff", lwd = 3)
dev.off()

## export
# master file
fwrite(fisher.genes, "../dat/fisher.genes.tsv", sep = "\t")
# table
fwrite(fisher.genes[same.disease == TRUE & p.adjusted < 0.05, ], "../doc/TableS1.csv", sep = ",")


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

## test if distributions significantly different with Mann-Whitney test (Wilcoxon rank sum test)
# adjusted p-values
mw.res <- wilcox.test(x = fisher.drugs[existing.indication == FALSE, -log10(p.adjusted)], y = fisher.drugs[existing.indication == TRUE, -log10(p.adjusted)])
print(mw.res)
print(mw.res$p.value)
png("../dat/fisher.drugs.pvalues.boxplots.png", width = 6 * 150, height = 6 * 150, res = 150)
print(ggplot(fisher.drugs, aes(x = existing.indication, y = -log10(p.adjusted))) +
    geom_boxplot(fill = "#ff6600", outlier.shape = NA) +
    coord_cartesian(ylim = quantile(fisher.drugs[, -log10(p.adjusted)], c(0.06, 0.94))) +
    xlab("Existing indication") +
    ylab("-log10(adjusted p-value)") +
    theme_bw(18) +
    scale_x_discrete(breaks = c(FALSE, TRUE), labels = c("No", "Yes")) +
    ggtitle(paste("p-value =", sprintf("%.2e", mw.res$p.value))))
dev.off()
# odds ratios
mw.res <- wilcox.test(x = fisher.drugs[existing.indication == FALSE, odds.ratio], y = fisher.drugs[existing.indication == TRUE, odds.ratio]) 
print(mw.res)
print(mw.res$p.value)
png("../dat/fisher.drugs.oddsratios.boxplots.png", width = 6 * 150, height = 6 * 150, res = 150)
print(ggplot(fisher.drugs, aes(x = existing.indication, y = odds.ratio)) +
    geom_boxplot(fill = "#ff6600", outlier.shape = NA) +
    coord_cartesian(ylim = quantile(fisher.drugs[, odds.ratio], c(0.04, 0.96))) +
    xlab("Existing indication") +
    ylab("Odds ratio") +
    theme_bw(18) +
    scale_x_discrete(breaks = c(FALSE, TRUE), labels = c("No", "Yes")) +
    ggtitle(paste("p-value =", sprintf("%.2e", mw.res$p.value))))
dev.off()

# perform ROC analysis
preds <- fisher.drugs[, -log10(p.adjusted)]
labls <- as.numeric(fisher.drugs[, existing.indication])
roc.res <- roc(response = labls, predictor = preds, algorithm = 2, ci = TRUE, ci.method = "bootstrap", boot.n = 1000, parallel = TRUE, progress = "none")
print(roc.res)
sp.ci <- ci.sp(roc.res, sensitivities = seq(0, 1, 0.05), boot.n = 1000, parallel = TRUE, progress = "none")
se.ci <- ci.se(roc.res, specifities = seq(0, 1, 0.05), boot.n = 1000, parallel = TRUE, progress = "none")
png("../dat/fisher.drugs.roc.png", width = 6 * 150, height = 6 * 150, res = 150)
par(pty = "s")
plot(roc.res, main = paste("AUC =", round(roc.res$auc, 2)), xlab = "False positive rate", ylab = "True positive rate", identity.lty = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, cex = 1.5)
plot(se.ci, type = "shape", col = "lightgrey", border = NA, no.roc = TRUE)
plot(sp.ci, type = "shape", col = "lightgrey", border = NA, no.roc = TRUE)
plot(roc.res, add = TRUE, col = "#ff6600", lwd = 3)
dev.off()

## export
# master file
fwrite(fisher.drugs, "../dat/fisher.drugs.tsv", sep = "\t")
# tables
fwrite(fisher.drugs[existing.indication == TRUE & p.adjusted < 0.05, ], "../doc/TableS2.csv", sep = ",")
fwrite(fisher.drugs[existing.indication == FALSE & p.adjusted < 1e-10, ], "../doc/TableS3.csv", sep = ",")
