# DGE script for calling DE genes in T. leiboldiana and T. vanhyningii between timepoints and watering regimes

# Originally written by S. Saadain, rewritten and reconstructed by C. Groot Crego in September 2023

# Load libraries
pacman::p_load("edgeR", "ggplot2", "ggpubr", "tidyverse", "gridExtra", "ggfortify", "pheatmap", "VennDiagram", "plotly", "ggrepel", "gprofiler2", "topGO", "dplyr", "xlsx", "tibble", "readxl", "statmod", "reshape2", "rtracklayer", "GenomicFeatures")

# functions

# Added on 14.11.2023, this is quite similar to my initial filtering strategy but a bit less strict.
Trimming_group <- function(dat,groups,minCount){
  #method 2
  keep <- rowSums(cpm(dat)>=minCount) >= 3 #number of biological replicates
  out2 <- dat[keep,]
  print(paste0("Number of genes after trimming using min count in min 4 samples:",nrow(out2)))
  out2 <- DGEList(counts=out2,group=groups)
  out2 <- calcNormFactors(out2, method ="TMM")
}

# Perform the QLFTest on fitted GLM's for each gene to obtain DGE results for a given contrast (e.g Day vs night or Drought vs. Control), will also adjust p-values for multiple testing
makeContrast <- function(matrix, my_contrast, contrast) {
  glmQLF_res <- glmQLFTest(matrix, contrast = my_contrast[, contrast])
  glmQLF_res$table$padj <- p.adjust(glmQLF_res$table$PValue, method = "BH")
  res_table <- glmQLF_res$table
  return(res_table)
}

# Performs filtering of DE genes based on adjusted p-value and logfoldchange. Additionally, it marks whether a gene is up- or downregulated
filtering <- function(DGE, res_table_name_list, nameslist, pvalue=0.05, minchange=2) {
  gene_list <- vector("list", length(res_table_name_list) * 2)
  name <- vector("list", length(res_table_name_list) * 2)
  i <- 1
  for (y in nameslist) {
    name[[i]] <- paste0(y, "_up")
    i <- i + 1
    name[[i]] <- paste0(y, "_dn")
    i <- i + 1
  }
  i <- 1
  for (res in res_table_name_list) {
    gene_list[[i]] <- rownames(DGE)[res$padj < pvalue & res$logFC > minchange] # use log2(minchange) to have a foldchange of 1.5!
    i <- i + 1
    gene_list[[i]] <- rownames(DGE)[res$padj < pvalue & res$logFC < -minchange]
    i <- i + 1
  }
  names(gene_list) <- name
  gene_list
}

# This will compile lists of up- and down-regulated DE genes and write out
create_df_and_write <- function(up_vector, dn_vector, tio, name, species) {
  # Create a dataframe for gene_id and direction
  df_direction <- data.frame(
    Gene_ID = c(up_vector, dn_vector),
    Direction = c(rep("up", length(up_vector)), rep("down", length(dn_vector)))
  )
  # Add a new column to tio that contains the row names
  tio$Gene_ID <- rownames(tio)
  # Extract the rows from tio that match the Gene_ID in df_direction
  df_info = tio[tio$Gene_ID %in% df_direction$Gene_ID, ]
  # Add the Direction column to df_info
  df_info = merge(df_info, df_direction, by = "Gene_ID")
  # Write out the info dataframe
  write.table(df_info, paste0(species, "_", name, "_DEGenes.Info.txt"), sep = "\t", row.names = FALSE, quote = F)
}

listgoanal <- function(list, geneID2GO) {
  j <- 1
  res_list <- list()
  for (i in list){
    allgenes <- factor(as.numeric(rownames(feature_counts)  %in% list[[j]]))
    names(allgenes) <- rownames(feature_counts)
    print(summary(allgenes))
    res_list[[j]] <- topGoAnalysis(allGenes = allgenes, gene2GO = geneID2GO)
    j <- j + 1
  }
  names(res_list) <- names(list)
  return(res_list)
}

to_ID <- function(query) {
  id <- strsplit(query, split = "_")
  return(unlist(id)[2 * (1:length(query))]) # da die Ursprungsdatei eine Liste ist wird der Output auch eine Liste, um auf das richtige Element zuzugreifen acessen wir jedes 2. element
}

# 1. Read sampling info
info <- read.csv(file = "Sampling_info_for_DGE.txt", header=T, sep = "\t") # File was created in the old EdgeR script from Sarah by combining Figure_FeatureCounts2022.csv and Tfas_gene22.counts.csv

# 1.1. Read count data
feature_counts <- read.table("counts_table.Mapped_to_Tfas.csv", sep = "\t", row.names = 1, header=T)
colnames(feature_counts) <- paste(info$ID, info$Treatment, info$Species, info$DayNight, sep = "_")
combined <- colnames(feature_counts)

####
# For PCA, check Make_PCA_Figure.R
####

# 2. Model Design
# Tio
counts_tio <- feature_counts[,info$Species == "Tio"]
tio_info <- info[info$Species == "Tio",]
model_tio <- data.frame(row.names = combined[info$Species == "Tio"],
                        Time = as.factor(tio_info$DayNight),
                        Treatment = as.factor(tio_info$Treatment))
group <- paste(tio_info$DayNight, tio_info$Treatment, sep = "_")
design_tio <- model.matrix(~0 + group)

# Tle
counts_tle <- feature_counts[,info$Species == "Tle"]
tle_info <- info[info$Species == "Tle",]
model_tle <- data.frame(row.names = combined[info$Species == "Tle"],
                        Time = as.factor(tle_info$DayNight),
                        Treatment = as.factor(tle_info$Treatment))
group <- paste(tle_info$DayNight, tle_info$Treatment, sep = "_")
design_tle <- model.matrix(~0 + group)
dge_tle <- DGEList(counts_tle)

# 3. DGE analysis

# Tio

# Filter lowly expressed genes

# Old method where the average cpm per row has to be > 1
#counts_tio_filtered <- counts_tio[apply(cpm(counts_tio),1,function(x){!(mean(x)<1)}),]
#dim(counts_tio_filtered) # 17,231 genes kept

# New method where the cpm has to be > 1 in at least 3 samples. This method automatically converts the counts to a EdgeR DGEList object and applies calcNormFactors normalization
group <- paste(tio_info$DayNight, tio_info$Treatment, sep = "_")
counts_tio_filtered <- Trimming_group(counts_tio, group, minCount = 1) # 17,899 genes kept

# Estimate dispersion
DGE_tio <- estimateDisp(counts_tio_filtered, design = design_tio, robust = T) # uses Cox-Reid profile-adjusted likelihood (CR) method to estimate dispersions. Given the DGEList object (y here as its already filered and normalized) and the design matrix, GLMs are fitted. This allows valid estimation of the dispersion, since all systematic sources of variation are accounted for.

# Fit GLM's
GLM_tio <- glmQLFit(DGE_tio, design = design_tio) # glmQLFit() fits neg. binomial GLM for each tag and produces on object of class DGEGLM with new components, which then can be passed to glmQLFTest() to carry out the QL F-test (see makeContrasts function above). Here the coefficients can be chosen from the full design matrix. This gives the null model against which the full model is compared.
# QL F-test (not LR test) is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It provides more robust and reliable error rate control when the number of replicates is small.

# Create contrasts between the groups
tiocontrasts <- makeContrasts("DayNight" = 0.5 * (groupDay_control + groupDay_drought) - 0.5 * (groupNight_control + groupNight_drought),
                              "DroughtControl" = 0.5 * (groupDay_drought + groupNight_drought) - 0.5 * (groupDay_control + groupNight_control),
                              "DayNight_atControl" = groupDay_control - groupNight_control,
                              "DayNight_atDrought" = groupDay_drought - groupNight_drought,
                              "DroughtControl_atDay" = groupDay_drought - groupDay_control,
                              "DroughtControl_atNight" = groupNight_drought - groupNight_control, levels = design_tio)
# DayNightcontrol contains everything differentially expressed between Day/Night at control
# droughtcontrolNight contains everything which is differentially expressed between drought/control at Night

# Create a list containing glmQLFTest() output for each contrast (see function makeContrast above)
# output is list of genes for each contrast
res_tio <- list(DayNight = makeContrast(GLM_tio, tiocontrasts, "DayNight"),
                DroughtControl = makeContrast(GLM_tio, tiocontrasts, "DroughtControl"),
                DayNight_atControl = makeContrast(GLM_tio, tiocontrasts, "DayNight_atControl"),
                DayNight_atDrought = makeContrast(GLM_tio, tiocontrasts, "DayNight_atDrought"),
                DroughtControl_atDay = makeContrast(GLM_tio, tiocontrasts, "DroughtControl_atDay"),
                DroughtControl_atNight = makeContrast(GLM_tio, tiocontrasts, "DroughtControl_atNight"))

# filters the output from before for p-value below 0.05, and a log2-fold-change below or above 1.5
res_tio_filt <- filtering(GLM_tio, res_tio, c("DayNight", "DroughtControl", "DayNight_atControl", "DayNight_atDrought", "DroughtControl_atDay", "DroughtControl_atNight"), minchange=1.5)

# Obtain the number of significant DE genes for each contrast
sumtio <- summary(res_tio_filt)[,1]

# Write the list of DE genes for each contrast, with a second column marking if they are down or upregulated
for (name in unique(names(res_tio_filt))) {
  # Replace "_up" or "_dn" with an empty string
  name = sub("_up", "", name)
  name = sub("_dn", "", name)
  print(name)
  up_vector = res_tio_filt[[paste0(name, "_up")]]
  dn_vector = res_tio_filt[[paste0(name, "_dn")]]
  tio = res_tio[[name]]
  create_df_and_write(up_vector, dn_vector, tio, name, "VAN")
}

# Tle

# Filter
group <- paste(tle_info$DayNight, tle_info$Treatment, sep = "_")
counts_tle_filtered <- Trimming_group(counts_tle, group, minCount = 1) # 17,676 genes kept
# Estimate dispersion
DGE_tle <- estimateDisp(counts_tle_filtered, design = design_tle, robust = T)
# Fit GLM's
GLM_tle <- glmQLFit(DGE_tle, design = design_tle)
# Make contrasts
tlecontrasts <- makeContrasts("DayNight" = 0.5 * (groupDay_control + groupDay_drought) - 0.5 * (groupNight_control + groupNight_drought),
                              "DroughtControl" = 0.5 * (groupDay_drought + groupNight_drought) - 0.5 * (groupDay_control + groupNight_control),
                              "DayNight_atControl" = groupDay_control - groupNight_control,
                              "DayNight_atDrought" = groupDay_drought - groupNight_drought,
                              "DroughtControl_atDay" = groupDay_drought - groupDay_control,
                              "DroughtControl_atNight" = groupNight_drought - groupNight_control, levels = design_tle)

# Compile GLMQLFtest results
res_tle <- list(DayNight = makeContrast(GLM_tle, tlecontrasts, "DayNight"),
            DroughtControl = makeContrast(GLM_tle, tlecontrasts, "DroughtControl"),
            DayNight_atControl = makeContrast(GLM_tle, tlecontrasts, "DayNight_atControl"),
            DayNight_atDrought = makeContrast(GLM_tle, tlecontrasts, "DayNight_atDrought"),
            DroughtControl_atDay = makeContrast(GLM_tle, tlecontrasts, "DroughtControl_atDay"),
            DroughtControl_atNight = makeContrast(GLM_tle, tlecontrasts, "DroughtControl_atNight"))

# filter
res_tle_filt <- filtering(GLM_tle, res_tle, c("DayNight", "DroughtControl", "DayNight_atControl", "DayNight_atDrought", "DroughtControl_atDay", "DroughtControl_atNight"), minchange = 1.5)
# Obtain counts of genes for each contrast
sumtle <- summary(res_tle_filt)[,1]

# Write the list of DE genes for each contrast, with a second column marking if they are down or upregulated
for (name in unique(names(res_tle_filt))) {
  # Replace "_up" or "_dn" with an empty string
  name = sub("_up", "", name)
  name = sub("_dn", "", name)
  print(name)
  up_vector = res_tle_filt[[paste0(name, "_up")]]
  dn_vector = res_tle_filt[[paste0(name, "_dn")]]
  tle = res_tle[[name]]
  create_df_and_write(up_vector, dn_vector, tle, name, "LEI")
}

# Overlap of trimmed reads for hypergeometric test
genesTio <- rownames(DGE_tio) # 17899
genesTle <- rownames(DGE_tle) # 17676
commonGenes <- length(intersect(genesTio, genesTle)) # 16462 (92 %, 93 %)

# 4. Compile sums of DE genes for all contrasts and write to output
sum_all <- as.data.frame(cbind(sumtio, sumtle))
colnames(sum_all) <- c("T.ion (CAM)", "T. lei (C3)")
# Turn down-regulated into negative values
rownames_neg <- grepl("_dn", rownames(sum_all))
sum_all$`T.ion (CAM)` <- as.numeric(sum_all$`T.ion (CAM)`)
sum_all$`T. lei (C3)` <- as.numeric(sum_all$`T. lei (C3)`)
sum_all[rownames_neg, ] <- -sum_all[rownames_neg, ]
# write to barplot directory
write.csv(sum_all, "II. Barplots UpDown/DE_Gene_counts.ForBarplots.csv")

# 5. Compile overlap for venn diagrams
# venn diagramms
venn_DayNight <- list(c(res_tio_filt$DayNightcontrol_up, res_tio_filt$DayNightcontrol_dn), c(res_tio_filt$DayNightdrought_up, res_tio_filt$DayNightdrought_dn), c(res_tle_filt$DayNightcontrol_up, res_tle_filt$DayNightcontrol_dn), c(res_tle_filt$DayNightdrought_up, res_tle_filt$DayNightdrought_dn))
saveRDS(venn_DayNight, file = "III. Venn Diagrams/VennDiagram_Input_Day_vs_Night.20230914.rds")

venn_DroughtControl <- list(c(res_tio_filt$droughtcontrolDay_up, res_tio_filt$droughtcontrolDay_dn), c(res_tle_filt$droughtcontrolDay_up, res_tle_filt$droughtcontrolDay_dn), c(res_tle_filt$droughtcontrolNight_up, res_tle_filt$droughtcontrolNight_dn), c(res_tio_filt$droughtcontrolNight_up, res_tio_filt$droughtcontrolNight_dn))
saveRDS(venn_DroughtControl, file = "III. Venn Diagrams/VennDiagram_Input_Drought_vs_Control.20230914.rds")

# 6. Create separate gene lists for each section of the Venn Diagram

# First for Day vs Night Venn diagram
droughttio <- c(res_tio_filt$DayNightdrought_up, res_tio_filt$DayNightdrought_dn)
droughttle <- c(res_tle_filt$DayNightdrought_up, res_tle_filt$DayNightdrought_dn)
controltle <- c(res_tle_filt$DayNightcontrol_up, res_tle_filt$DayNightcontrol_dn)
controltio <- c(res_tio_filt$DayNightcontrol_up, res_tio_filt$DayNightcontrol_dn)

# Unique to one group
a = setdiff(controltio, union(union(droughttle,droughttio), controltle))
b = setdiff(droughttio,union(union(droughttle,controltio), controltle))
c = setdiff(droughttle, union(union(controltio,droughttio), controltle))
d = setdiff(controltle, union(union(controltio,droughttio), droughttle))
# Shared between two groups
ab = setdiff(intersect(controltio, droughttio), union(controltle, droughttle))
ac = setdiff(intersect(controltio, droughttle), union(droughttio, controltle))
ad = setdiff(setdiff(intersect(controltio, controltle), droughttio), droughttle)
bc = setdiff(intersect(droughttle,droughttio), union(controltle, controltio))
bd = setdiff(setdiff(intersect(droughttio, controltle), droughttle), controltio)
cd = setdiff(intersect(controltle, droughttle), union(controltio, droughttio))
# Shared between three groups
abc = intersect(setdiff(intersect(droughttio, droughttle), controltle), controltio)
abd = intersect(setdiff(intersect(controltio, controltle), droughttle), droughttio)
acd = setdiff(intersect(intersect(controltio,controltle),droughttle), droughttio)
bcd = setdiff(intersect(intersect(droughttio,droughttle),controltle), controltio)
# Shared between all groups
abcd = intersect(intersect(intersect(droughttio,droughttle),controltio),controltle)

# Write lists

# Create a list where the names of the list elements are 'ab', 'ac', 'ad', 'bc'
venn_sections_DayNight <- list(a = a, b = b, c = c, d = d, ab = ab, ac = ac, ad = ad, bc = bc, bd = bd, cd = cd, abc = abc, abd = abd, acd = acd, bcd = bcd, abcd = abcd)

# Loop over the names of the list
for (name in names(venn_sections_DayNight)) {
  print(name)
  # Check if the name is 'c', 'd', or 'cd'
  if (name %in% c("c", "d", "cd")) {
    # If it is, use res_tle to create the subset
    subset <- res_tle$DayNight[venn_sections_DayNight[[name]],]
  } else {
    # If it's not, use res_tio
    subset <- res_tio$DayNight[venn_sections_DayNight[[name]],]
  }
  # Convert row names to a column named 'Gene_ID'
  subset <- rownames_to_column(subset, var = "Gene_ID")
  head(subset)
  # Write the subset to a file
  write.table(subset, file = paste0("Venn_GeneList_DayNight_Intersection_", name, ".txt"), sep = "\t", quote = F, row.names = F)
}

# Now for Drought vs control
daytio <- c(res_tio_filt$droughtcontrolDay_up, res_tio_filt$droughtcontrolDay_dn)
daytle <- c(res_tle_filt$droughtcontrolDay_up, res_tle_filt$droughtcontrolDay_dn)
nighttle <- c(res_tle_filt$droughtcontrolNight_up, res_tle_filt$droughtcontrolNight_dn)
nighttio <- c(res_tio_filt$droughtcontrolNight_up, res_tio_filt$droughtcontrolNight_dn)

# Unique to one group
a = setdiff(daytio,union(union(daytle,nighttio), nighttle))
b = setdiff(nighttio, union(union(daytle,daytio), nighttle))
c = setdiff(nighttle, union(union(nighttio,daytio), daytle))
d = setdiff(daytle, union(union(nighttio,daytio), nighttle))
# Shared between two groups
ab = setdiff(intersect(nighttio, daytio), union(nighttle, daytle))
ac = setdiff(setdiff(intersect(daytio, nighttle), daytle), nighttio)
ad = setdiff(intersect(daytle,daytio), union(nighttle, nighttio))
bd = setdiff(intersect(nighttio, daytle), union(daytio, nighttle))
bc = setdiff(setdiff(intersect(nighttio, nighttle), daytio), daytle)
cd = setdiff(intersect(nighttle, daytle), union(nighttio, daytio))
# Shared between three groups
abc = intersect(setdiff(intersect(nighttio, nighttle), daytle), daytio)
abd = intersect(setdiff(intersect(daytio, daytle), nighttle), nighttio)
acd = setdiff(intersect(intersect(daytio,daytle),nighttle), nighttio)
bcd = setdiff(intersect(intersect(nighttio,nighttle),daytle), daytio)
# Shared between all groups
abcd = intersect(intersect(intersect(daytio,daytle),nighttio),nighttle)

# Write lists

# Create a list where the names of the list elements are 'ab', 'ac', 'ad', 'bc'
venn_sections_DroughtControl <- list(a = a, b = b, c = c, d = d, ab = ab, ac = ac, ad = ad, bc = bc, bd = bd, cd = cd, abc = abc, abd = abd, acd = acd, bcd = bcd, abcd = abcd)

# Loop over the names of the list
for (name in names(venn_sections_DroughtControl)) {
  print(name)
  # Check if the name is 'c', 'd', or 'cd'
  if (name %in% c("c", "d", "cd")) {
    # If it is, use res_tle to create the subset
    subset <- res_tle$DayNight[venn_sections_DroughtControl[[name]],]
  } else {
    # If it's not, use res_tio
    subset <- res_tio$DayNight[venn_sections_DroughtControl[[name]],]
  }
  # Convert row names to a column named 'Gene_ID'
  subset <- rownames_to_column(subset, var = "Gene_ID")
  # Write the subset to a file
  write.table(subset, file = paste0("Venn_GeneList_DroughtControl_Intersection_", name, ".txt"), sep = "\t", quote = F, row.names = F)
}

### For gene-specific expression analyses, write the normalized counts
cpm_tio <- cpm(counts_tio)
cpm_tle <- cpm(counts_tle)
cpm_counts <- cbind(cpm_tio, cpm_tle)
write.table(cpm_counts, "Tfas_gene22.counts.NormalizedCPM.WithinSpecies.csv", row.names = T, quote = F, sep = "\t")

# For visualization purposes, perform TPM on all (to adjust for exonic length per gene)
cpm_counts_all <- cpm(feature_counts)

# Exon lengths obtained with following code:
# awk '$3 == "mRNA" || $3 == "exon" {print $0}' Tillandsia_leiboldiana_v1.2.edited_allfeatures.26chrom.gff | cut -f 3,4,5,9 | sed 's/;/\t/g' | cut -f 1,2,3,4 | sed 's/:/\t/g' | cut -f 1,2,3,4 | sed 's/ID=//g' > Tlei_intron_lengths.txt
# dat <- read.delim("Tfas_gff_digest.txt", sep = "\t", header = F)
# colnames(dat) <- c("feature", "start", "end", "name")
# # calculate full transcript length
# dat$length <- dat$end - dat$start
# # sum up all exon lengths
# dat2 <- aggregate(dat$length, by=list(dat$feature, dat$name), FUN=sum)
# colnames(dat2) <- c("feature", "name", "length")
# # split exon and mrna
# date <- dat2[dat2$feature == "exon",]
# colnames(date) <- c("feature", "name", "exonic_length")
# datg <- dat2[dat2$feature == "mRNA",]
# colnames(datg) <- c("feature", "name", "full_length")
# # merge again
# dat3 <- merge(date[,c(2,3)], datg[,c(2,3)], by = "name", all = T, )
# dat3[is.na(dat3$exonic_length),2] <- 0
# # calculate intron length
# dat3$intronic_length <- dat3$full_length - dat3$exonic_length
#write.table(dat3, file = "Tfas_gene_exon_intron_length.txt", quote = F, sep = "\t", row.names = F)
lengths <- read.table("Tfas_gene_exon_intron_length.txt", sep = "\t", header = T)

# Function to transform counts
tpm3 <- function(counts, len) {
  # Filter out genes with 0 length
  valid_indices <- which(len != 0)
  counts <- counts[valid_indices, , drop = FALSE]
  len <- len[valid_indices]

  rpk <- counts / len  # Normalizing by gene length to get RPK (Reads Per Kilobase)
  col_sums <- colSums(rpk)
  if (any(col_sums == 0, na.rm = TRUE)) {
    warning("Some columns have a sum of zero; these will result in NaN values")
  }
  tpm <- rpk / colSums(rpk) * 1e6  # Normalizing by total counts to get TPM
  return(tpm)
}

# Make a named vector of lengths, with names as gene IDs
lengths_vector <- setNames(lengths$exonic_length, lengths$name)
# Re-order lengths_vector to match the order of genes in counts
lengths_ordered <- lengths_vector[rownames(feature_counts)]
# Convert lengths from base pairs to kilobases
lengths_kb <- lengths_ordered / 1000

# Calculate TPM
tpm_counts <- tpm3(feature_counts, lengths_kb)
write.table(tpm_counts, file = "Tfas_gene22.counts.NormalizedTPM.All.csv", row.names = T, quote = F,
            sep = "\t")
