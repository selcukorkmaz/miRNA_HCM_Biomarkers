# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0

# Load required libraries for data analysis and visualization
library(GEOquery)
library(limma)
library(umap)
library(readxl)
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(EnhancedVolcano)
library(ggrepel)


# Load series and platform data from GEO
gset = read.table("data/GSE36961_non-normalized.txt", header = TRUE, sep = "\t")
gset = gset[,c(1, seq(2,290,2))]
gset = gset[complete.cases(gset),]

# Load group information
groups_df = read_excel("data/GSE36961_groups.xlsx")
ids = as.data.frame(colnames(gset)[-1])
colnames(ids) = "Title"

# Merge group information with sample IDs
groups_df2 <- left_join(ids, groups_df, by = "Title")
groups = as.factor(groups_df2$`Sample type`)
group_names <- make.names(c("HCM","Control"))

# Log2 transformation check and application
ex <- gset[,-1]
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[ex <= 0] <- NaN
gset_exp <- log2(ex) }


# Preprocessing of transformed data
gset_exp[gset_exp <= 0] <- NaN
gset_exp <- gset_exp[complete.cases(gset_exp),]
gset_t = gset_exp  %>% t() %>% as.data.frame()
colnames(gset_t) = gset$ID_REF
gset_t$group <- groups


# Box plots before normalizing
ex <- as.matrix(gset_t[,-ncol(gset_t)] %>% t())
ex_df <- as.data.frame(ex)
ex_df$SampleID <- rownames(ex_df)
long_df <- melt(ex_df, id.vars = "SampleID", variable.name = "miRNA", value.name = "Expression")

long_df <- long_df %>%
  mutate(Group = rep(groups, each = nrow(long_df) / length(groups)))

long_df$miRNA_Group <- interaction(long_df$miRNA, long_df$Group, sep = " - ")
head(long_df)

before=ggplot(long_df, aes(x = miRNA_Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + # Rotate x labels for readability
  labs(x = "", y = "") +
  ylim(c(5,8))+
  scale_fill_manual(values = c("Control" = "#1B9E77", "HCM" = "#7570B3"))+
  theme(legend.text=element_text(size=10), legend.title=element_blank())+
  scale_x_discrete(breaks=long_df$miRNA_Group,
                   labels=long_df$miRNA)


# Normalize expression data
gset_exp_norm <- normalizeBetweenArrays(gset_exp) 

# Prepare for linear modeling
gset_exp2 = gset_exp_norm  %>% t() %>% as.data.frame()
colnames(gset_exp2) = gset$ID_REF
gset_exp2$group <- groups
design <- model.matrix(~group + 0, gset_exp2)
colnames(design) <- levels(groups)

# Vooma transformation and linear modeling
transposed_gset_exp2 <- as.matrix(gset_exp2[,-ncol(gset_exp2)]) %>% t()
v <- voomaByGroup(transposed_gset_exp2, group=group_names, design, plot=T, cex=0.1, pch=".", col=1:nlevels(groups))
v$genes <- gset$ID_REF

# Fit linear model to the voom-transformed data
fit  <- lmFit(v)

# Set up contrasts for differential expression analysis
cts <- c(paste(group_names[1],"-",group_names[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# Compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","AveExpr"))
write.table(tT, file="data/TopGenes_GSE36961.txt", row.names=F, sep="\t", quote = F)


# Summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=1)

# Filter and save differentially expressed genes
degs = dT@.Data[dT@.Data[,1]!=0,]
degs_data = tT[abs(tT$logFC)>=1,]
degs_data = degs_data[degs_data$ID %in% names(degs),]
degs_data <- subset(degs_data, select=c("ID","adj.P.Val","P.Value","t","B","logFC"))
write.table(degs_data, file="data/DEGs_GSE36961.txt", row.names=F, sep="\t", quote = F)


# Creating plots using layout: voom:mean-variance trend plot, moderated t statistic plot, 
# p-adj value distribuution plot, expression value distribuution plot
ex <- as.matrix(gset_exp2[,-ncol(gset_exp2)] %>% t())
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates

t.good <- which(!is.na(fit2$F)) 

layout <- layout(matrix(c(1,2,3,4), 2, 2, byrow = F)) 
layout <- voomaByGroup(transposed_gset_exp2, group=group_names, design, plot=T, cex=0.1, pch=".", col=1:nlevels(groups))
layout <- hist(tT$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
               ylab = "Number of genes", main = "P-adj value distribution")
layout <- qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

layout <- plotDensities(ex, group=groups, main="GSE36961: Expression value distribution", legend ="topright")

# Box plot after notmalization
ex_df <- as.data.frame(ex)
ex_df$SampleID <- rownames(ex_df)

# Melt the data frame to long format
long_df <- melt(ex_df, id.vars = "SampleID", variable.name = "miRNA", value.name = "Expression")

# Add group information to the long data frame
long_df <- long_df %>%
  mutate(Group = rep(groups, each = nrow(long_df) / length(groups)))

# Create a new variable for interaction between miRNA and Group to plot side by side
long_df$miRNA_Group <- interaction(long_df$miRNA, long_df$Group, sep = " - ")
head(long_df)

# Now, use ggplot2 to create the box plot with Control and HCM side by side
after=ggplot(long_df, aes(x = miRNA_Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylim(c(5,8)) +
  labs(x = "", y = "") +
  scale_fill_manual(values = c("Control" = "#1B9E77", "HCM" = "#7570B3"))+
  theme(legend.text=element_text(size=10), legend.title=element_blank())+
  scale_x_discrete(breaks=long_df$miRNA_Group,
                   labels=long_df$miRNA)


p = before / after

p = p+ plot_annotation(tag_levels = 'A')

ggsave('~/Documents/Studies/LVH_DeepLearning/Figures/GSE36961_RLE.png', p, width = 16, height = 8)


# Volcano plot
adjPVal_cutoff = 0.05
logFC_cutoff = 1

down <- rownames(tT)[tT$adj.P.Val < adjPVal_cutoff & tT$logFC < -logFC_cutoff]
up <- rownames(tT)[tT$adj.P.Val < adjPVal_cutoff & tT$logFC > logFC_cutoff]

keyvals.colour1 <- ifelse(
  rownames(tT) %in% up, 'red',
  ifelse(rownames(tT) %in% down, 'blue',
         'black'))

keyvals.colour1[is.na(keyvals.colour1)] <- 'black'
names(keyvals.colour1)[keyvals.colour1 == 'black'] <- 'Not significant'
names(keyvals.colour1)[keyvals.colour1 == 'red'] <- 'Up-regulated'
names(keyvals.colour1)[keyvals.colour1 == 'blue'] <- 'Down-regulated'


p1 <- EnhancedVolcano(tT,
                      lab = NA,
                      x = 'logFC',
                      y = 'adj.P.Val',
                      xlim = c(-4,4),
                      title = NULL,  
                      subtitle = NULL, 
                      caption = NULL,
                      pCutoff = adjPVal_cutoff,
                      FCcutoff = logFC_cutoff,
                      pointSize = c(ifelse(tT$adj.P.Val < adjPVal_cutoff & abs(tT$logFC) > logFC_cutoff, 4, 2)),
                      colCustom = keyvals.colour1,
                      legendLabSize = 8,
                      legendIconSize = 4,
                      axisLabSize = 10,
                      legendPosition = 'bottom')

# Mean-Difference Plot of Expression Data
df <- data.frame(
  logFC = tT$logFC,
  AveExpr = tT$AveExpr,
  Significant = names(keyvals.colour1)   # Assuming 'dT' indicates significant genes as non-zero
)

# Create the Mean-Difference Plot using ggplot2
p2 <- ggplot(df, aes(x = AveExpr, y = logFC, color = Significant)) +
  geom_point(alpha = 0.6) +  # Adjust point transparency with 'alpha'
  scale_color_manual(values = c("Not significant" = "grey7", "Up-regulated" = "red", "Down-regulated" = "blue")) +  # Change colors as needed
  labs(x = "Average Expression", y = "Log Fold Change") +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend if not needed



# UAMP plot

df <- data.frame(UMAP1 = ump$layout[,1], UMAP2 = ump$layout[,2], Group = groups)

# Creating the UMAP plot with ggplot2 and ggrepel for label repulsion
p3 <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(aes(shape = Group), size = 3) +  # Adjust size as needed
  scale_shape_manual(values=c(19,19)) +  # Ensure unique shapes for groups if desired
  # geom_text_repel(aes(label = rownames(df)), size = 2, max.overlaps = Inf) +  # Adjust size for visibility
  scale_color_manual(values = c("#179f76", "#7570a7")) +  # Customize colors if needed
  theme_minimal() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "right") + 
  theme(legend.text=element_text(size=10), legend.title=element_blank())
# guides(color = guide_legend(title = "Group"), shape = guide_legend(title = "Group"))  # Customize legend



p = p1/p2|p3

p = p+ plot_annotation(tag_levels = 'A')

ggsave('~/Documents/Studies/LVH_DeepLearning/Figures/GSE36961_volcano_umap.png', p, width = 12, height = 8)



