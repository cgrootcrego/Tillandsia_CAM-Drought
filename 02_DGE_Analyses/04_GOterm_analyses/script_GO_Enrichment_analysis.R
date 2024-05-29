#!/usr/bin/Rscript --vanilla

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("BiocManager", "topGO", "stringr", "GOplot")

#setwd('/home/botanik/Documents/GitHub/Tillandsia-CAM-Drought/2. DGE/')

# Load functions
 change_names <- function(data, name_list){
   colnames(data) <- name_list
   return(data)
 }

 rename <- function(table, geneNames){
   names(table) <- geneNames
   return(table)
 }

 attach_enriched_go_genes <- function(enriched_go_with_my_genes){
   enriched_go_with_my_genes.list = c()
   for (i in 1:length(enriched_go_with_my_genes)){
     enriched_go_with_my_genes.list = c(enriched_go_with_my_genes.list, enriched_go_with_my_genes[[i]])
   }
   return(enriched_go_with_my_genes.list)
 }

 circle_dat <- function(terms, genes){

  colnames(terms) <- tolower(colnames(terms))
  terms$genes <- toupper(terms$genes)
  genes$ID <- toupper(genes$ID)
  tgenes <- strsplit(as.vector(terms$genes), ', ')
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  logFC <- sapply(unlist(tgenes), function(x) genes$logFC[match(x, genes$ID)])
  if(class(logFC) == 'factor'){
    logFC <- gsub(",", ".", gsub("\\.", "", logFC))
    logFC <- as.numeric(logFC)
  }
  s <- 1; zsc <- c()
  for (c in 1:length(count)){
    value <- 0
    e <- s + count[c] - 1
    value <- logFC[s:e] # This will perform a zscore on the raw logFC values, therefore the z-score reflects the combined differences in expression
    #value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, -1)) # this performs a z-score on the count of genes that are under and over expressed, so not taking into account the abundance of expression change
    zsc <- c(zsc, sum(value, na.rm = F) / sqrt(count[c])) #### HERE : na.rm = TRUE, takes all the genes into account !
    s <- e + 1
  }
  if (is.null(terms$id)){
    df <- data.frame(category = rep(as.character(terms$category), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }else{
    df <- data.frame(category = rep(as.character(terms$category), count), ID = rep(as.character(terms$id), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }
  return(df)
}

 GOBar <- function(data, display, order.by.zscore = T, title, zsc.col, contrast){
   id <- adj_pval <- zscore <- NULL
   if (contrast == "DayNight"){
     low <- "Night"
     high <- "Day"
     zsc.col <- c('#F18F01', '#EFEBCE', '#76415E')
   } else if (contrast == "DroughtControl") {
     low <- "Control"
     high <- "Drought"
     zsc.col <- c('#5B4B49', '#E1CA96', '#0C7C59')
   } else {
     stop("Contrast not correctly formulated. Must be either DayNight or DroughtControl")
   }
   colnames(data) <- tolower(colnames(data))
   data$adj_pval <- -log(data$adj_pval, 10)
   sub <- data[!duplicated(data$term), ]
   sub <- sub[order(sub$zscore, decreasing = T), ]
   print(min(sub$zscore))
   print(max(sub$zscore))
   g <-  ggplot(sub, aes(y = factor(term, levels = rev(unique(stats::reorder(term, adj_pval)))), x = adj_pval, fill = zscore)) +
     geom_bar(stat = 'identity', position=position_dodge()) +
     theme_bw() +
     scale_fill_gradient2('Log(FC)\n(z-score)\n', space = 'Lab', low = zsc.col[3], mid = zsc.col[2], high = zsc.col[1], guide = guide_colourbar(title.position = "top", title.hjust = 0),
                          breaks = c(min(sub$zscore), max(sub$zscore)), labels = c(paste0(low," > ", high), paste0(high," > ", low))) +
     scale_color_manual(values=c("black", "black", "black"), guide = "none") +
     labs(title = '', y = '', x = '-log (adj p-value)') +
     geom_text(aes(y=term, x=-0.25, label=count, color = "black"))
   g
 }

# Load arguments
# 1 is GotoGenes.map, 2 is the subset of genes to test, 3 is contrast, 4 is the species and 5 is the background condition
args <- commandArgs(trailingOnly = TRUE)

# Note: gotogenes.map for these analyses was created like so:
# awk '$3=="mRNA" {print $0}' Tillandsia_fasciculata_v1.gff | cut -f 9 | awk -F ';' -v OFS="\t" '{id="NA"; ont="NA"; for(i=1;i<=NF;i++) if($i ~ /^ID=/) {id=substr($i, 4)} else if ($i ~ /^Ontology_id=/) {ont=substr($i, 13)}} {print id, ont}' | sort -k 1

# This contains the full gene annotation (non-robust and non-orthologous genes included), i.e. 34,886 genes of which 10,567 have no GO terms

# Read input
geneID2GO <- readMappings(file = args[[1]])
DE_genes <- read.table(args[[2]], header = T)
contrast = args[[3]] # "DroughtControl" or "DayNight"
print(contrast)
species = args[[4]] # VAN or LEI
print(species)
background = args[[5]] # The condition in which the contrast was made, e.g. in Drought versus Control at Day, it is day
print(background)
outputname = paste0("GOterm_enrichment_", species, "_", contrast, "_at", background, ".txt")

# How to run consecutively in command line
# for i in ../*_DEGenes.Info.txt; do
#   echo $i
#   species=`echo $i | sed 's,../,,g' | cut -f 1 -d'_'`
#   contrast=`echo $i | sed 's,../,,g' | cut -f 2 -d'_'`
#   background=`echo $i | sed 's,../,,g' | cut -f 3 -d'_' | sed 's/at//g'`
#   echo $species $contrast $background
#   if [ "$background" == "DEGenes.Info.txt" ]; then
#     continue
#   fi
#   Rscript script_GO_Enrichment_analysis.R genes_to_GO_Tfas.exonic.FullAnnotation.map $i $contrast $species $background
# done

 # Run enrichment
GO2geneID <- inverseList(geneID2GO)
geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% DE_genes[,1]))
names(geneList) <- geneNames

name_list = c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")
table = as.factor(geneNames) %in% DE_genes[,1]
int_table = as.integer(table)
int_fac_table = factor(int_table)
fac_table = rename(table = int_fac_table, geneNames = geneNames)

GOdata.BP = new("topGOdata", ontology = "BP", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.MF = new("topGOdata", ontology = "MF", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.CC = new("topGOdata", ontology = "CC", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultWeight01.BP = runTest(GOdata.BP, statistic = "fisher")
resultWeight01.MF = runTest(GOdata.MF, statistic = "fisher")
resultWeight01.CC = runTest(GOdata.CC, statistic = "fisher")

allRes.BP1 = GenTable(GOdata.BP, weight01_pval=resultWeight01.BP, orderBy = "weight01", ranksOf = "weight01",topNodes = 100, numChar=1000)
allRes.BP2 = cbind(allRes.BP1,"BP")
allRes.BP = change_names(data = allRes.BP2, name_list = name_list)

allRes.MF1 = GenTable(GOdata.MF, weight01_pval=resultWeight01.MF, orderBy = "weight01", ranksOf = "weight01",topNodes = 100, numChar=1000)
allRes.MF2 = cbind(allRes.MF1,"MF")
allRes.MF = change_names(data = allRes.MF2, name_list = name_list)

allRes.CC1 = GenTable(GOdata.CC, weight01_pval=resultWeight01.CC, orderBy = "weight01", ranksOf = "weight01",topNodes = 100, numChar=1000)
allRes.CC2 = cbind(allRes.CC1,"CC")
allRes.CC = change_names(data = allRes.CC2, name_list = name_list)

allRes1 = rbind(allRes.BP,allRes.MF)
allRes = rbind(allRes1, allRes.CC)
allRes[startsWith(allRes$weight01_pval, "<"),6] <- "1e-30"

allGO.BP = genesInTerm(GOdata.BP)
allGO.MF = genesInTerm(GOdata.MF)
allGO.CC = genesInTerm(GOdata.CC)
allGO = c(allGO.BP, allGO.MF, allGO.CC)

# Create final table
SAM_ANOTATION = lapply(allGO,function(x) x[x %in%  DE_genes[,1]])
enriched_go_with_my_genes = lapply(SAM_ANOTATION[allRes[,1]], paste0, collapse = ", ")
enriched_go_with_my_genes.list = attach_enriched_go_genes(enriched_go_with_my_genes)
go.dataframe = data.frame("Category" = allRes$branch, "ID" = allRes$GO.ID, "Term" = allRes$Term,
                          "Genes" = as.vector(enriched_go_with_my_genes.list),
                          "adj_pval" = as.numeric(sub(",", ".", allRes$weight01_pval, fixed = TRUE)))

# Write table that reports all significant go terms and the underlying genes
go.dataframe.significant <- go.dataframe[go.dataframe$adj_pval <= 0.05,]
write.table(go.dataframe.significant, file = outputname, row.names = F, quote = F, sep = "\t")

# Make barplot of up and down-regulated functions
EC.genelist = data.frame("ID" = DE_genes$Gene_ID, "logFC" = DE_genes$logFC, "logCPM" = DE_genes$logCPM, "P.Value" = DE_genes$PValue, "adj.P.Val" = DE_genes$padj)

circ = circle_dat(go.dataframe, EC.genelist)

circ.significant <- circ[circ$adj_pval <= 0.05,]

reduced_circ <- reduce_overlap(circ.significant, overlap = 0.9)

#if (contrast == "DayNight"){
#  reduced_circ <- subset(reduced_circ, count >= 3)
#}

pdf(file = paste0("Barplot_TopGO_LogFC_Zscore.", species, "_", contrast, "_at", background, ".pdf"), width = 12, height = 6)
GOBar(subset(reduced_circ, category == 'BP'), contrast = contrast)
dev.off()
