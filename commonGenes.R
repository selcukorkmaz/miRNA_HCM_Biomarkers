###### Part 1: Venn diagrams show common genes between DEGs in GSE36961 and and target genes of DE-miRNAs from GSE36946 #####

# Load necessary libraries for data manipulation, miRNA-target interaction analysis, and plotting
library(GEOquery)   # For downloading and accessing GEO datasets
library(multiMiR)   # For querying miRNA-target interactions databases
library(ggvenn)     # For creating Venn diagrams to visualize overlaps
library(patchwork)


# Read differential expression genes (DEGs) file for GSE36946 with headers, separated by tabs
demirna = read.table("TopGenes/DEGs_GSE36946.txt", header=TRUE, sep="\t")

# Define a vector of miRNAs involved in the study
mirnas = c("hsa-miR-373", "hsa-miR-373", "hsa-miR-514", "hsa-miR-514", "hsa-miR-514", "hsa-miR-514", "hsa-miR-10a", "hsa-miR-10a",
           "hsa-miR-144", "hsa-miR-144")

# Define a corresponding vector of precursors for the miRNAs
precursor = c("hsa-miR-373-3p", "hsa-miR-373-5p", "hsa-miR-514a-3p", "hsa-miR-514a-5p", "hsa-miR-514b-3p", "hsa-miR-514b-5p",
              "hsa-miR-10a-3p", "hsa-miR-10a-5p", "hsa-miR-144-3p", "hsa-miR-144-5p")

# Define a vector indicating whether the miRNAs are up or down-regulated
upDown = c("up", "up", "up", "up", "up", "up", "down", "down", "down", "down")

# Combine miRNA, precursor, and expression vectors into a data frame
demirnas = cbind.data.frame(miRNA = mirnas, precursor = precursor, expr = upDown)

# Retrieve target genes for the specified miRNAs from the multiMiR database
getTargets = get_multimir(mirna = demirnas$precursor, summary = TRUE)

# Read DEGs data for GSE36961 and store it in a variable
degs_gse36961_data = read.table("TopGenes/DEGs_GSE36961.txt", header = TRUE, sep = "\t")
# Display dimensions and the first few rows of the dataset to understand its structure
dim(degs_gse36961_data)
head(degs_gse36961_data)

# Find up-regulated and down-regulated common DEGs based on log fold change (logFC)
upCommonDEGs = degs_gse36961_data[degs_gse36961_data$logFC > 0, "ID"]
downCommonDEGs = degs_gse36961_data[degs_gse36961_data$logFC < 0, "ID"]

# Filter miRNA targets from the retrieved target list based on their expression direction
upMiRNAtargets = getTargets@data[getTargets@data$mature_mirna_id %in% demirnas[demirnas$expr == "up", "precursor"], "target_symbol"]
downMiRNAtargets = getTargets@data[getTargets@data$mature_mirna_id %in% demirnas[demirnas$expr == "down", "precursor"], "target_symbol"]

# Identify and save up-regulated DEGs that are also targeted by up-regulated miRNAs
updegs = degs_gse36961_data[degs_gse36961_data$logFC > 0,]
upDEGs = updegs[updegs$ID %in% upMiRNAtargets,]
write.table(upDEGs, "TopGenes/common_upDEGs.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# Read top genes data from GSE141910 for further analysis
# topgenes_gse141910 = read.table("TopGenes/TopGenes_GSE141910.txt", 
                               # header=TRUE, sep="\t")

# Retrieve platform annotation data for GPL6104 to map gene symbols
# gpl6104 = getGEO("GPL6104")
# annotation = gpl6104@dataTable@table %>% dplyr::select(ID, ILMN_Gene)
# 
# # Merge top genes data with annotation data based on ID
# topgenes_gse141910 = left_join(topgenes_gse141910, annotation, by="ID")
# 
# # Find up-regulated DEGs in top genes list and annotate them
# upDEGs_test = topgenes_gse141910[topgenes_gse141910$ID %in% annotation[annotation$ILMN_Gene %in% updegs[updegs$ID %in% upMiRNAtargets,"ID"],"ID"],]

# Generate a list for Venn diagram comparison of DEGs and DE-miRNA Targets
x = list(
  DEGs = upCommonDEGs,
  "DE-miRNA Targets" = upMiRNAtargets
)

# Create a Venn diagram for up-regulated DEGs and miRNA targets
p1 = ggvenn(
  x, 
  fill_color = c("#FF2C21", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 6, text_size = 6
)

# Exclude a specific gene from down-regulated DEGs for further analysis
downdegs = degs_gse36961_data[degs_gse36961_data$logFC < 0,]
# downdegs = downdegs[downdegs$ID != "SERPINE1",]

# Save down-regulated DEGs that are targeted by down-regulated miRNAs
write.table(downdegs[downdegs$ID %in% downMiRNAtargets,], "TopGenes/common_downDEGs.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# Find down-regulated DEGs in top genes list and annotate them
downDEGs_test = topgenes_gse141910[topgenes_gse141910$ID %in% annotation[annotation$ILMN_Gene %in% downdegs[downdegs$ID %in% downMiRNAtargets,"ID"],"ID"],]
downDEGs_test$ILMN_Gene

# Prepare a list for Venn diagram of down-regulated DEGs excluding a specific gene
x = list(
  DEGs = downCommonDEGs,
  "DE-miRNA Targets" = downMiRNAtargets
)

# Identify overlaps between down-regulated DEGs and miRNA targets, excluding a specific gene
downCommonDEGs[downCommonDEGs %in% downMiRNAtargets]

# Create a Venn diagram for down-regulated DEGs and miRNA targets
p2 = ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 6, text_size = 6
)

# Combine both Venn diagrams for up and down-regulated genes and targets
p = p1 + p2
p = p + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = '18'))

# Save the combined Venn diagram plot to a file
ggsave('~/Documents/Studies/LVH_DeepLearning/IJC/Figures/Figure4.png', p, width = 12, height = 8)


###### Part 2: Expression patterns in GSE141910 of common miRNAs with targets from the DEGs #####

# Load necessary libraries for statistical analysis and enhanced plotting
library(rstatix)    # Helper functions for statistical tests and summaries
library(ggpubr)     # For enhancing 'ggplot2' plots for publication quality

# Read expression data for the study GSE141910
exprs_gse141910 = read.table("TopGenes/exprs_GSE141910.txt", header = TRUE, sep = "\t")
# Display the dimensions of the expression data
dim(exprs_gse141910)
# Preview a specific subset of the data
exprs_gse141910[1:5,13131:13132]

# Prepare data for upregulated DEGs for analysis
# Select relevant columns based on upregulated DEGs and group them
upDEGs_test_data = exprs_gse141910 %>% dplyr::select(upCommonDEGs[upCommonDEGs %in% upMiRNAtargets], group)
# Rename the first two columns based on gene symbols for clarity
# colnames(upDEGs_test_data)[1:2] = upDEGs_test$ILMN_Gene

# Convert the wide data to a long format for easier analysis
upDEGs_test_data2 <- upDEGs_test_data %>%
  gather(key = "Expression", value = "Value", upCommonDEGs[upCommonDEGs %in% upMiRNAtargets]) %>%
  convert_as_factor(Expression) %>% na.omit()

# Perform t-tests on the long-format data, adjust p-values for multiple testing using Bonferroni method, and add significance labels
stat.test <- upDEGs_test_data2 %>%
  group_by(Expression) %>%
  t_test(Value ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

# Plot the results of the upregulated DEGs analysis
myplot <- upDEGs_test_data2 %>%
  drop_na() %>%
  ggplot(aes(x=Expression, y=Value, color=group)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(), alpha = 0.3) + 
  scale_color_manual(values=c("#1B9E77", "#7570B3")) +
  ylab("Target Gene Expression") +
  xlab("") +
  theme_bw() + 
  theme(text = element_text(size = 20)) + 
  theme(legend.position="top", legend.title = element_blank())

# Add statistical annotations to the plot
stat.test <- stat.test %>% add_xy_position(x = "Expression")
p1 = myplot + stat_pvalue_manual(stat.test, label = "p.adj.signif")

# Repeat the process for downregulated DEGs
downDEGs_test_data = exprs_gse141910 %>% dplyr::select(downCommonDEGs[downCommonDEGs %in% downMiRNAtargets], group)
# colnames(downDEGs_test_data)[1:6] = downDEGs_test$ILMN_Gene

downDEGs_test_data2 <- downDEGs_test_data %>%
  gather(key = "Expression", value = "Value", downCommonDEGs[downCommonDEGs %in% downMiRNAtargets]) %>%
  convert_as_factor(Expression) %>% na.omit()

stat.test <- downDEGs_test_data2 %>%
  group_by(Expression) %>%
  t_test(Value ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

# Plot the results of the downregulated DEGs analysis
myplot <- downDEGs_test_data2 %>%
  drop_na() %>%
  ggplot(aes(x=Expression, y=Value, color=group)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(), alpha = 0.3) + 
  scale_color_manual(values=c("#1B9E77", "#7570B3")) +
  ylab("Target Gene Expression") +
  xlab("") +
  theme_bw() + 
  theme(text = element_text(size = 20)) + 
  theme(legend.position="top", legend.title = element_blank())

stat.test <- stat.test %>% add_xy_position(x = "Expression")
p2 = myplot + stat_pvalue_manual(stat.test, label = "p.adj.signif")

# Combine plots for both upregulated and downregulated DEGs
p = p1 + p2+ plot_layout(widths = c(1, 3))
p = p + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = '18'))

# Save the combined plot to a file
ggsave('~/Documents/Studies/LVH_DeepLearning/Figures/clinical_validation.png', p, width = 18, height = 8)

