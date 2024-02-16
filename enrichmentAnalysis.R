
# Load necessary libraries for enrichment analysis and visualization
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(ggplot2)

# Set working directory to the location of the gene expression data files
setwd("~/Documents/Studies/LVH_DeepLearning/TopGenes/")

# Read DEG data for GSE36961
degs_gse36961_data <- read.table("data/DEGs_GSE36961.txt", header = TRUE, sep = "\t")
# Display the dimensions and the first few rows of GSE36961 DEG data
dim(degs_gse36961_data)
head(degs_gse36961_data)

##### Analysis for GSE36961 ######
# Create a list of gene symbols for enrichment analysis
deg_list <- degs_gse36961_data$ID

# Map gene symbols to Entrez IDs using the org.Hs.eg.db package
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = deg_list,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
entrez_ids = unname(entrez_ids) # Remove names for further analysis

# Perform KEGG pathway enrichment analysis
kegg_result <- enrichKEGG(gene = entrez_ids,
                          organism = 'hsa',
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
# Visualize KEGG pathway enrichment results
p1 = dotplot(kegg_result) 

# Perform GO enrichment analysis for Biological Process
go_bp <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
# Visualize Biological Process enrichment results
p2 = dotplot(go_bp)

# Perform GO enrichment analysis for Cellular Component
go_cc <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
# Visualize Cellular Component enrichment results
p3 = dotplot(go_cc)

# Perform GO enrichment analysis for Molecular Function
go_mf <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
# Visualize Molecular Function enrichment results
p4 = dotplot(go_mf)

# Combine and display all dot plots with annotations
p = (p2 + p3 + p4 + p1) + plot_annotation(tag_levels = 'A')
# Save the combined plot to a file
ggsave('~/Documents/Studies/LVH_DeepLearning/Figures/GSE36961_kegg.png', p, width = 12, height = 8)


