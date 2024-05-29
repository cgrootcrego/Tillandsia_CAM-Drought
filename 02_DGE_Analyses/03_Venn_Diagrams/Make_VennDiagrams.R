pacman::p_load("ggplot2","VennDiagram", "devtools", "ggvenn", "ggrepel", "SuperExactTest", "UpsetR")
#devtools::install_github("yanlinlin82/ggvenn")

# Coordinates for labels
coord <- read.table('coordinates_venn.txt', header = T)

venn_time <- readRDS("VennDiagram_Input_Day_vs_Night.20230914.rds")

# Make Venn Diagram for timewise comparisons
names(venn_time) <- c("T. vanhyningii,\n control",
                      "T. vanhyningii, drought",
                      "T. leiboldiana,\n control",
                      "T. leiboldiana, drought")

venn_time <- list(
  "T. vanhyningii,\n control" = venn_time[["T. vanhyningii,\n control"]],
  "T. vanhyningii, drought" = venn_time[["T. vanhyningii, drought"]],
  "T. leiboldiana, drought" = venn_time[["T. leiboldiana, drought"]],
  "T. leiboldiana,\n control" = venn_time[["T. leiboldiana,\n control"]]
)

# Hypergeometric test

# This was done followingthe tutorial / example at https://github.com/mw201608/SuperExactTest/blob/master/examples/set_html.Md

# The number of trimmed reads for the two species is 17,899 and 17,676 for Tio and Tlei respectively. These have a 92 % and 93 % of those respectively overlap with each other (16,462). Given the high similarity in starting gene set where the DE analysis is performed, I decided to take the average between the two gene sets, 17,788 as the gene population size for the hypergeometric test.

# Calculate the expected overlap between groups based on the total gene population they were sampled from
total=17788
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))
# The expected overlap between groups is 27 %

res=supertest(venn_time, n=total)
plot(res, Layout="landscape", degree=2:4, sort.by="size", margin=c(0.5,5,1,2))
write.table(summary(res)$Table, file="HyperGeometricTest_Venn_Time.table.csv", row.names=FALSE, sep = "\t")

# Number of unique elements in all groups
VAN_list <- c(venn_time[[1]], venn_time[[2]])
LEI_list <- c(venn_time[[3]], venn_time[[4]])
total_count <- length(Reduce(union, venn_time))

control_list <- c(venn_time[[1]], venn_time[[4]])
drought_list <- c(venn_time[[2]], venn_time[[3]])

VAN_unique_genes_time <- unique(VAN_list) # 1,952 genes
LEI_unique_genes_time <- unique(LEI_list) # 1,299 genes
control_unique_genes_time <- unique(control_list) # 1,972 genes
drought_unique_genes_time <- unique(drought_list) # 1,973 genes

# Find the intersection (shared genes) between the two species
shared_genes <- intersect(VAN_unique_genes_time, LEI_unique_genes_time)
length(shared_genes)

common_genes_VAN <- intersect(venn_time[[1]], venn_time[[2]]) # 868 genes
common_genes_LEI <- intersect(venn_time[[3]], venn_time[[4]]) # 455 genes
common_genes_control <- intersect(venn_time[[1]], venn_time[[4]]) # 284 genes
common_genes_drought <- intersect(venn_time[[2]], venn_time[[3]]) # 345 genes

p <- ggvenn(venn_time, set_name_size = 6, text_size = 8) +
  scale_fill_manual(values = c("#0C7C59", "#5B4B49", "#E1CA96", "#A7E8BD"))+  geom_label(data = coord, aes(x, y, label = s), size = 6)

# Save the plot as a PDF
ggsave("Figure4A_VennDiagram_Time.pdf", p, width = 14, height = 12)

# Save the plot as a PNG
ggsave("Figure4A_VennDiagram_Time.png", p, width = 14, height = 12, dpi = 500)


# Make Venn Diagram for watering comparisons
venn_water <- readRDS("VennDiagram_Input_Drought_vs_Control.20230914.rds")

names(venn_water) <- c("T. vanhyningii, day",
                      "T. leiboldiana, day",
                      "T. leiboldiana, night",
                      "T. vanhyningii, night")

venn_water <- list(
  "T. vanhyningii, day" = venn_water[["T. vanhyningii, day"]],
  "T. vanhyningii, night" = venn_water[["T. vanhyningii, night"]],
  "T. leiboldiana, night" = venn_water[["T. leiboldiana, night"]],
  "T. leiboldiana, day" = venn_water[["T. leiboldiana, day"]]
)

total=17788
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))
# The expected overlap between groups is 27 %

res=supertest(venn_water, n=total)
summary(res)
write.table(summary(res)$Table, file="HyperGeometricTest_Venn_WateringRegime.table.csv", row.names=FALSE, sep = "\t")
plot(res, Layout="landscape", degree=2:4, sort.by="size", margin=c(0.5,5,1,2))

# Total number of unique genes in the venn-diagram
total_count <- length(Reduce(union, venn_water)) # 1348 genes

# Number of unique elements in all groups
VAN_list <- c(venn_water[[1]], venn_water[[2]])
LEI_list <- c(venn_water[[3]], venn_water[[4]])

day_list <- c(venn_water[[1]], venn_water[[4]])
night_list <- c(venn_water[[2]], venn_water[[3]])

VAN_unique_genes_water <- unique(VAN_list) # 123 genes
LEI_unique_genes_water <- unique(LEI_list) # 1,257 genes
day_unique_genes_water <- unique(day_list) # 459 genes
night_unique_genes_water <- unique(night_list) # 1,185 genes

# Find the intersection (shared genes) between the two species
shared_genes <- intersect(VAN_unique_genes, LEI_unique_genes)
length(shared_genes)

common_genes_VAN <- intersect(venn_water$`T. vanhyningii,
day`, venn_water$`T. vanhyningii, night`) # 10 genes
common_genes_LEI <- intersect(venn_water$`T. leiboldiana, night`, venn_water$`T. leiboldiana,
day`) # 276 genes
common_genes_day <- intersect(venn_water$`T. vanhyningii,
day`, venn_water$`T. leiboldiana,
day`) # 11 genes
common_genes_night <- intersect(venn_water$`T. vanhyningii, night`, venn_water$`T. leiboldiana, night`) # 11 genes

p_w <- ggvenn(venn_water, set_name_size = 7, text_size = 8, show_elements = F) +
  scale_fill_manual(values = c("#F18F01", "#76415E", "#BAA0AE", "#EFEBCE"))+  geom_label(data = coord, aes(x, y, label = s), size = 6)
p_w

# Save the plot as a PDF
ggsave("Figure4B_VennDiagram_WateringRegime.pdf", p_w, width = 14, height = 12)

# Save the plot as a PNG
ggsave("Figure4B_VennDiagram_WateringRegime.png", p_w, width = 14, height = 12, dpi = 500)

# Find overlap between DE genes time-wise and condition-wise
# Combine all lists within venn_time to get a unique list of genes
combined_venn_time <- Reduce(union, venn_time)
length(combined_venn_time)
combined_all <- union(combined_venn_time, combined_venn_water)
length(combined_all)
# Combine all lists within venn_water to get a unique list of genes
combined_venn_water <- Reduce(union, venn_water)
# Find the intersection (overlap) between the two combined lists
overlap_genes <- intersect(combined_venn_time, combined_venn_water)
# Count the number of genes in the overlap
length(unique(overlap_genes)) # 582 genes are both DE in time and in condition across both species

# Initialize a list to keep track of the DE status for each gene
DE_status <- list()
# For each overlapping gene, check DE status in each dataset
for (gene in overlap_genes) {
  DE_status[[gene]] <- c(
    DE_time_VAN = gene %in% VAN_unique_genes_time,
    DE_time_LEI = gene %in% LEI_unique_genes_time,
    DE_water_VAN = gene %in% VAN_unique_genes_water,
    DE_water_LEI = gene %in% LEI_unique_genes_water
  )
}
# DE_status now contains the DE status for each overlapping gene
# You can convert it to a data frame for easier viewing and analysis
DE_status_df <- data.frame(t(sapply(DE_status, unlist)))
write.table(DE_status_df, "Genes_DE_in_time_AND_watering_regime.txt", quote = F, sep = "\t", row.names = T)
# Print out the DE status data frame
summary(DE_status_df)
list_data <- lapply(DE_status_df, function(x) row.names(DE_status_df)[x])
# Of the 582 genes that are DE in both comparisons, 225 are DE only in LEI, 41 are DE only in VAN. 151 genes are DE in VAN for time but not water and DE in LEI in water but not time. 137 genes are DE in time in both species and DE in water in LEI but not VAN
library(UpSetR)
pdf("Overlap_Time_Watering_Analysis.pdf", height = 6, width = 8)
upset(fromList(list_data), order.by = "freq")
dev.off()
### CAM-related only ###
cam_genes <- scan("../V. Search CAM genes/Genes_of_interest_all.in_Tfas.ID.txt", what = "", sep = "\n")
venn_time_cam <- lapply(venn_time, function(x) intersect(x, cam_genes))

# Number of unique elements in all groups
VAN_list <- c(venn_time_cam[[1]], venn_time_cam[[2]])
LEI_list <- c(venn_time_cam[[3]], venn_time_cam[[4]])
total_count <- length(Reduce(union, venn_time_cam))

control_list <- c(venn_time_cam[[1]], venn_time_cam[[4]])
drought_list <- c(venn_time_cam[[2]], venn_time_cam[[3]])

VAN_unique_genes <- unique(VAN_list) # 101 genes
LEI_unique_genes <- unique(LEI_list) # 86 genes
control_unique_genes <- unique(control_list) # 116 genes
drought_unique_genes <- unique(drought_list) # 102 genes

# Find the intersection (shared genes) between the two species
shared_genes <- intersect(VAN_unique_genes, LEI_unique_genes)
length(shared_genes)

common_genes_VAN <- intersect(venn_time_cam$`T. vanhyningii,
control`, venn_time_cam$`T. vanhyningii, drought`) # 46 genes
common_genes_LEI <- intersect(venn_time_cam$`T. leiboldiana, drought`, venn_time_cam$`T. leiboldiana,
control`) # 31 genes
common_genes_control <- intersect(venn_time_cam$`T. vanhyningii,
control`, venn_time_cam$`T. leiboldiana,
control`) # 24 genes
common_genes_drought <- intersect(venn_time_cam$`T. vanhyningii, drought`, venn_time_cam$`T. leiboldiana, drought`) # 22 genes

p <- ggvenn(venn_time_cam, set_name_size = 6, text_size = 8) +
  scale_fill_manual(values = c("#0C7C59", "#5B4B49", "#E1CA96", "#A7E8BD"))+  geom_label(data = coord, aes(x, y, label = s), size = 6)
# Save the plot as a PDF
ggsave("Figure4A_VennDiagram_Time_CAM.pdf", p, width = 14, height = 12)
# Save the plot as a PNG
ggsave("Figure4A_VennDiagram_Time_CAM.png", p, width = 14, height = 12, dpi = 500)


venn_water_cam <- lapply(venn_water, function(x) intersect(x, cam_genes))

# Number of unique elements in all groups
VAN_list <- c(venn_water_cam[[1]], venn_water_cam[[2]])
LEI_list <- c(venn_water_cam[[3]], venn_water_cam[[4]])

day_list <- c(venn_water_cam[[1]], venn_water_cam[[4]])
night_list <- c(venn_water_cam[[2]], venn_water_cam[[3]])
total_count <- length(Reduce(union, venn_water_cam))

VAN_unique_genes <- unique(VAN_list) # 7 genes
LEI_unique_genes <- unique(LEI_list) # 78 genes
day_unique_genes <- unique(day_list) # 4 genes
night_unique_genes <- unique(night_list) # 77 genes

common_genes_VAN <- intersect(venn_water_cam[[1]], venn_water_cam[[2]]) # 1 genes
common_genes_LEI <- intersect(venn_water_cam[[3]], venn_water_cam[[4]]) # 22 genes
common_genes_day <- intersect(venn_water_cam[[1]], venn_water_cam[[4]]) # 1 genes
common_genes_night <- intersect(venn_water_cam[[2]], venn_water_cam[[3]]) # 0 genes

shared_genes <- intersect(VAN_unique_genes, LEI_unique_genes)
length(shared_genes)

p_w <- ggvenn(venn_water_cam, set_name_size = 7, text_size = 8, show_elements = F) +
  scale_fill_manual(values = c("#F18F01", "#76415E", "#BAA0AE", "#EFEBCE"))+  geom_label(data = coord, aes(x, y, label = s), size = 6)
p_w
# Save the plot as a PDF
ggsave("Figure4B_VennDiagram_WateringRegime_CAM.pdf", p_w, width = 14, height = 12)
# Save the plot as a PNG
ggsave("Figure4B_VennDiagram_WateringRegime_CAM.png", p_w, width = 14, height = 12, dpi = 500)
