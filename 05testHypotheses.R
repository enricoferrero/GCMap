set.seed(16, kind = "L'Ecuyer-CMRG")
library(data.table)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)
library(BiocParallel)
register(MulticoreParam(workers = 1))
library(fgsea)

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

# reverse and mimic
for (search.mode in c("reverse", "mimic")) {

    if (search.mode == "reverse") {
        aggravate.value <- FALSE
    } else if (search.mode == "mimic") {
        aggravate.value <- TRUE
    }

    # read datasets
    stopgap <- fread("../dat/stopgap.tsv")
    opentargets <- fread("../dat/opentargets.tsv")
    harmonizome <- fread("../dat/harmonizome.tsv")
    lincs <- fread(paste0("../dat/lincs.", search.mode, ".tsv"))

    ## hypothesis 1: for each disease, check if the overlap between drugs predicted to reverse the genetic signature and drugs currently used to treat the disease is significant

    # split by disease
    opentargets.list <- split(opentargets, opentargets$efo.id)
    lincs.list <- split(lincs, lincs$efo.id)

    # create universe of compounds
    chembl <- unique(c(opentargets$chembl.id, lincs$chembl.id))

    # perform Fisher's test
    fisher.drugs <- foreach (i = seq(opentargets.list), .combine = rbind) %dopar% {
        efo.id <- names(opentargets.list)[i]
        fisher.drugs <- checkOverlapSignificance(opentargets.list[[efo.id]]$chembl.id, lincs.list[[efo.id]]$chembl.id, chembl)
        cbind(efo.id, fisher.drugs)
    }

    # correct p-values
    fisher.drugs <- fisher.drugs[order(p.value), ]
    fisher.drugs[, padj := p.adjust(p.value, method = "fdr")]

    # check results 
    fisher.drugs[, sum(padj < 0.05) / .N]
    fisher.drugs[padj < 0.05, ]

    # export
    fisher.drugs <- fisher.drugs[padj < 0.05, ]
    fwrite(fisher.drugs, paste0("../dat/fisher.drugs.", search.mode, ".tsv"), sep = "\t")


    ## hypothesis 2: for each disease, check if the overlap between genes differentially expressed after drug treatment and genes genetically associated with the disease is significant

    # add perturbation string to LINCS data to match it with the (adapted) Harmonizome data
    lincs[, perturbation := tolower(paste(lincs.id, lincs.name, lincs.cell, lincs.time, lincs.time.unit, lincs.dose, lincs.dose.unit, sep = "_"))]
    harmonizome[, perturbation := tolower(gsub("_([1-9]*[0-9])\\.0_", "_\\1_", perturbation))]

    # split by disease
    lincs.list <- split(lincs, lincs$efo.id)
    stopgap.list <- split(stopgap$ensembl.id, stopgap$efo.id)

    # loop thorugh diseases
    gsea.genes <- foreach (i = seq(lincs.list), .combine = rbind) %dopar% {
        
        # pick a disease
        efo.id <- names(lincs.list)[i]
        stopgap.efo.id <- stopgap.list[efo.id]

        # loop thorugh compounds
        foreach (j = seq(nrow(lincs.list[[efo.id]])), .combine = rbind) %dopar% {

            # get gene list
            tmp <- harmonizome[grep(lincs.list[[efo.id]][j, perturbation], perturbation, fixed = TRUE), ]
            genes <- tmp[, direction]
            names(genes) <- tmp[, ensembl.id]

            # perform GSEA
            gsea.genes <- fgsea(pathways = stopgap.efo.id, stats = genes, minSize = 5, maxSize = 1000, nperm = 10000)
            # annotate
            gsea.genes <- gsea.genes[, .(efo.id = pathway, gsea.pval = pval, gsea.padj = padj, gsea.es = ES, gsea.nes = NES, gsea.n.more.extreme = nMoreExtreme, gsea.size = size)]
            if (nrow(gsea.genes) > 0) {
                gsea.genes <- cbind(gsea.genes, lincs.list[[efo.id]][j, .(chembl.id, lincs.score, lincs.id, lincs.name, lincs.cell, lincs.dose, lincs.dose.unit, lincs.time, lincs.time.unit)])
                gsea.genes <- merge(unique(stopgap[, .(efo.id, efo.term)]), gsea.genes, by = "efo.id")
            }

        }

    }

    # correct p-values
    gsea.genes <- gsea.genes[order(gsea.pval), ]
    gsea.genes[, gsea.padj := p.adjust(gsea.pval, method = "fdr")]

    # check results 
    gsea.genes[, sum(gsea.padj < 0.05) / .N]
    gsea.genes[gsea.padj < 0.05, ]

    # export
    gsea.genes <- gsea.genes[gsea.padj < 0.05, ]
    fwrite(gsea.genes, paste0("../dat/gsea.genes.", search.mode, ".tsv"), sep = "\t")

}
