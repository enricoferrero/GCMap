# GCMap
### Genetics Connectivity Map: Couple drug gene expression data with disease genetics data to predict genetics-backed indications for drug repositioning.

1. Test if genes differentially expressed after drug treatment are enriched for genes genetically associated with disease (for main drug indication).
	a. Get drug gene expression signatures (CMap, LINCS, BioMAP)
	b. Get drug indications (PharmaProjects, OpenTargets)
	c. Get genes genetically associated with disease from (GWAS catalog, GRASP, STOPGAP)
	d. Perform enrichment test (GSEA, others)

2.  Build machine/deep learning model that predicts indications based on enrichment of drug gene expression signatures in GWAS results.
	a. Observations are drug - disease pairs
	b1. Features are protein-coding genes with binary values (1: gene is in drug signature and genetically associated with disease; 0: it's not)
	b2. Features are protein-coding genes with factor values (A: gene is in drug signature and genetically associated with disease; B: gene is in drug signature; C: gene is genetically associated with disease; D: it's neither)
	c. Positive labels are succesfull drug - indication pairs, the rest are negative/unlabelled.


Alternative: instead of drug profiles use CRISPR KO profiles (GenomeCRISPR)
