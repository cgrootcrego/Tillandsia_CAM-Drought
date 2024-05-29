pacman::p_load("edgeR", "ggplot2", "ggpubr", "tidyverse", "gridExtra", "ggfortify", "pheatmap", "VennDiagram", "plotly", "ggrepel", "gprofiler2", "topGO", "dplyr", "xlsx", "tibble", "readxl")

Trimming_group <- function(dat,groups,minCount){
  #method 2
  keep <- rowSums(cpm(dat)>=minCount) >= 3 #number of biological replicates
  out2 <- dat[keep,]
  print(paste0("Number of genes after trimming using min count in min 3 samples:",nrow(out2)))
  out2 <- DGEList(counts=out2,group=groups)
  out2 <- calcNormFactors(out2, method ="TMM")
}

# 1. get info from previous file
info <- read.csv(file = "Sampling_info_for_DGE.txt", header=T)
info <- info %>% dplyr::select("ID","Treatment","Species","DayNight")

# 1.1. Count data input
feature_counts <- read.table("Tfas_gene22.counts.csv", sep = "\t", row.names = 1, header=T)
colnames(feature_counts) <- paste(info$ID, info$Treatment, info$Species, info$DayNight, sep = "_")
combined <- colnames(feature_counts)
rownames(info) <- colnames(feature_counts)

group <- paste(info$Species, info$DayNight, info$Treatment, sep = "_")
y <- Trimming_group(feature_counts, group, minCount = 1) # 19,213 genes kept
logCPM <- cpm(y, log=TRUE)  # Log-transformed CPM

#### MAKE PCA FIGURE ####
pca <- prcomp(t(feature_counts))
pca <- prcomp(t(logCPM))
pca_df <- as.data.frame(pca$x)
# Merge pca_df and info_df by row names
pca_df <- cbind(pca_df, info[rownames(pca_df), ])
pca_df$group <- factor(
  interaction(pca_df$Species, pca_df$Treatment),
  levels = unique(interaction(pca_df$Species, pca_df$Treatment))
)
# Convert 'group' to a factor and specify the level order
pca_df$group <- factor(pca_df$group, levels = c("Tio.control", "Tio.drought", "Tle.control", "Tle.drought"))

variance_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

p <-ggplot(pca_df, aes(x = PC1, y = PC2, fill = group, shape = DayNight)) +
  geom_point(size = 3, color = "black", stroke = 0.5) +
  theme_bw() +
  labs(fill = "Watering regime", shape = "Timepoint") +
  ggtitle("") +
  xlab(paste0("PC1 (", format(variance_explained[1], digits = 4), "% )")) +
  ylab(paste0("PC2 (", format(variance_explained[2], digits = 4), "% )")) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(labels = c("Tio.control" = "T. vanhyningii, control", "Tio.drought"="T. vanhyningii, drought", "Tle.control" = "T. leiboldiana, control", "Tle.drought"="T. leiboldiana, drought"),
                    values = c("Tio.control" = "#ffa400", "Tio.drought"="#FFE4B4", "Tle.control" = "#2a9d8f", "Tle.drought"="#BFEEE8")) +
  theme(axis.text = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = 21, fill = c("#ffa400", "#FFE4B4", "#2a9d8f", "#BFEEE8"),
                                                 color = "black")),
         shape = guide_legend(override.aes = list(fill = "white")))

# Save the plot as a PDF
ggsave("FigureS2_PCA_ReadCounts_VANLEI.Normalized.pdf", p, width = 6, height = 4)

# Save the plot as a PNG
ggsave("FigureS2_PCA_ReadCounts_VANLEI.Normalized.png", p, width = 6, height = 4, dpi = 300)
