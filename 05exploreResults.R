library(data.table)
library(ggplot2)
library(riverplot)
library(RColorBrewer)

# read datasets
opentargets.tas <- fread("../dat/opentargets.therapeuticareas.txt")
relevant.tas <- c("cardiovascular system", "digestive system", "endocrine system", "hematological system", "immune system", "infection", "metabolic system", "neoplasm", "nervous system", "reproductive system", "respiratory system", "skeletal system")
fisher.genes <- fread("../dat/fisher.genes.tsv")
fisher.drugs <- fread("../dat/fisher.drugs.tsv")

## hypothesis 1 (fisher.genes)
# add therapeutic area
fisher.genes <- merge(fisher.genes, opentargets.tas[, .(ID, therapeutic.area = TA_LABEL)], by.x = "efo.id.stopgap", by.y = "ID", all = FALSE)
# calculate number of disease per TA
genes.numbers <- fisher.genes[same.disease == TRUE & therapeutic.area %in% relevant.tas, .(number = sum(p.adjusted < 0.05)), by = therapeutic.area]
# plot
png("../dat/fisher.genes.ta.barplot.png", res = 150, width = 9 * 150, height = 6 * 150)
print(ggplot(genes.numbers, aes(x = reorder(therapeutic.area, -number), y = number)) +
      geom_bar(stat = "identity", colour = "black", fill = "#0066ff") +
      xlab("Therapeutic area") +
      ylab("Number of significant diseases") +
      theme_bw(18) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))
dev.off()

## hypothesis 2 (fisher.drugs)
# add therapeutic area
fisher.drugs <- merge(fisher.drugs, opentargets.tas[, .(ID, therapeutic.area = TA_LABEL)], by.x = "efo.id", by.y = "ID", all = FALSE, allow.cartesian = TRUE)
# calculate number of disease per TA
drugs.numbers <- fisher.drugs[existing.indication == TRUE & therapeutic.area %in% relevant.tas, .(number = sum(p.adjusted < 0.05)), by = therapeutic.area]
# plot
png("../dat/fisher.drugs.ta.barplot.png", res = 150, width = 9 * 150, height = 6 * 150)
print(ggplot(drugs.numbers, aes(x = reorder(therapeutic.area, -number), y = number)) +
      geom_bar(stat = "identity", colour = "black", fill = "#ff6600") +
      xlab("Therapeutic area") +
      ylab("Number of significant diseases") +
      theme_bw(18) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))
dev.off()

## repurposing visualitation with sankey plot
# grab sources (repositioning from)
sources <- unique(fisher.drugs[therapeutic.area %in% relevant.tas & existing.indication == TRUE, .(chembl.id, from.ta = therapeutic.area)])
# grab targets (repurposing to)
targets <- unique(fisher.drugs[therapeutic.area %in% relevant.tas & existing.indication == FALSE & p.adjusted < 1e-10, .(chembl.id, to.ta = therapeutic.area)])
# merge
drugs.repos <- merge(sources, targets, by = "chembl.id", all = FALSE, allow.cartesian = TRUE)
# edges
edges <- drugs.repos[, .(Value = .N), by = .(N1 = paste(from.ta, 1), N2 = paste(to.ta, 2))]
# nodes
nodes <- data.table(ID = c(paste(relevant.tas, 1), paste(relevant.tas, 2)), x = c(rep(1, length(relevant.tas)), rep(2, length(relevant.tas))), labels = rep(relevant.tas, 2), col = rep(brewer.pal(length(relevant.tas), "Set3"), 2))
# create riverplot object
sankey <- makeRiver(as.data.frame(nodes), as.data.frame(edges))
png("../dat/fisher.drugs.repositioning.sankeyplot.png", res = 300, width = 12 * 300, height = 25 * 300)
plot(sankey, srt = "0", textcex = 1.4)
dev.off()

## explore trajectories
# most common
edges[head(order(-Value), 10)]
# largest sources
edges[, .(Value = sum(Value)), by = N1][order(-Value),]
# largest targets
edges[, .(Value = sum(Value)), by = N2][order(-Value),]



