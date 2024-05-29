pacman::p_load("ggplot2", "tidyr", "splitstackshape")

setwd('/home/botanik/Documents/GitHub/Tillandsia-CAM-Drought/2. DGE/II. Barplots UpDown/')

tab <- read.table("DE_Gene_counts.ForBarplots.csv", sep = ",", header = T)[,1:3]

plot <- cbind(as.numeric(c(tab$T.ion..CAM., tab$T..lei..C3.)), c(rep("tio", 12), rep("tle", 12)))
plot <- cbind(rep(tab$X,2), plot)
#plot <- cbind(cSplit(plot, "V1", sep = "_")[,3], plot)
colnames(plot) <- c("compall", "comp", "value", "species")
plot <- plot[-c(1,2,3,4,13,14,15,16),]
write.table(plot, "GeneCounts_UpDown_AllComparisons_VANLEI.txt", sep = '\t', quote = FALSE, row.names = FALSE)

p <- ggplot(as.data.frame(plot), aes(compall, as.numeric(value), fill = species)) +
  geom_col(position = "dodge", color = "black") +
  geom_hline(yintercept = 0) +
  #geom_vline(xintercept = 2.5) +
  geom_text(aes(label = abs(as.numeric(value)), y = as.numeric(value)), 
            position = position_dodge(width = 0.9), 
            vjust = ifelse(as.numeric(plot$value) > 0, -0.25, 1.25)) + scale_y_continuous(breaks=c(-1500,-1250,-1000, -750, -500, -250, 0, 250, 500, 750, 1000, 1250, 1500)) + 
  scale_fill_manual(labels = c("tio" = "T. vanhyningii",
                               "tle" = "T. leiboldiana"),
                    values = c("tio" = "#ffa400", 
                               "tle" = "#2a9d8f")) +
  scale_x_discrete(labels = c("DayNightcontrol" = "Day vs. Night at control",
                              "DayNightdrought" = "Day vs. Night at drought",
                              "droughtcontrolDay" = "Drought vs. Control at Day",
                              "droughtcontrolNight" = "Drought vs. Control at Night")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + ylab("") +
  labs(fill = "Species", y="Number of upregulated (>0) and downregulated (<0) genes", x="") +
  theme(legend.text = element_text(face = "italic", size = 11),
        legend.title = element_text(size=12),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size=11))

p
# Save the plot as a PDF
ggsave("Figure3_Barplots_UpDown.pdf", p, width = 8, height = 6)

# Save the plot as a PNG
ggsave("Figure3_Barplots_UpDown.png", p, width = 9, height = 7, dpi = 500)
