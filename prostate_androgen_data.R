library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
library(curatedTCGAData)
library(pheatmap)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(DESeq2)
library(magrittr)
library(GEOquery)
library(tidyverse)

# extract and prepare the raw data
raw.data <- read.csv("/home/admin/Sung-prostate-cancer/GSE128749_raw_counts_GRCh38.p13_NCBI.tsv", sep = '\t')
gene_info <- read.csv("/home/admin/Sung-prostate-cancer/Human.GRCh38.p13.annot (2).tsv", sep = '\t')
rownames(raw.data) <- make.unique(coldata$Symbol)
raw.data <- raw.data[,-1]

# prepare the sample info 
gse <- getGEO("GSE128749", GSEMatrix = TRUE)
coldata <- gse$GSE128749_series_matrix.txt.gz@phenoData@data
coldata <- coldata[,c(2,11)]
names(coldata)[2] <- "Treatment"

coldata <- coldata %>%
  mutate(Treatment= sub('^treatment: ', '', Treatment))

# making the rownames and column names identical
all(rownames(coldata) %in% colnames(raw.data))
all(rownames(coldata) == colnames(raw.data))

count_matrix <- round(raw.data)
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ 1) # not spcifying model

## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 11,]
nrow(dds) #  11143 genes


# perform variance stabilization
dds_norm <- varianceStabilizingTransformation(dds)

# get normalized counts
normalized.count.matrix <- assay(dds_norm) %>% 
  t()

# norm.counts <- norm.counts[tumor$ids,]
normalized.count.matrix <- t(normalized.count.matrix) %>% as.data.frame()

# write.csv(normalized.count.matrix, "normalized.count.matrix.androgenPC.csv")

##================================================================
# Add  the metabolon number column for normalized.count.matrix 
##================================================================

list_D <- read.csv("/home/admin/Sunchul/List_D.csv")
# Get the intersected gene symbols
intersected.genes <- intersect(rownames(normalized.count.matrix), list_D$Gene_name)
length(intersected.genes)
# Get the indices of intersected genes in metabolin data frame
intersected.indices <- which(list_D$Gene_name %in% intersected.genes)
list_D <- list_D[intersected.indices,]
length(list_D$Gene_name)
library(dplyr)

normalized.count.matrix<- normalized.count.matrix[which(rownames(normalized.count.matrix) %in% intersected.genes),]
normalized.count.matrix$metabolon_number <- list_D$List.C.metabolon
normalized.count.matrix <- normalized.count.matrix %>% relocate(metabolon_number, .before = 1)

normalized.count.matrix$metabolon_number %>% table()
#  M1 M12 M15 M16 M19 M20 M35 M37 M38 M46 M47 M48 M49  M5 M50 M56 M57 M58 M59  M6 M60 M62 M69 M78 M80 
#  12   3  12   6  10  14   8   4   3  26  18   9  12   7  12  20  15  18  14   9  32   3  27   9  20 

#-------------------------------------------------------------------------------
# 01 subset the first androgen treated samples 
#-------------------------------------------------------------------------------

treated.androgen <- coldata[which(coldata$Treatment == "Synthetic Androgen (R1881)"),]

treated.counts <- normalized.count.matrix[,treated.androgen$geo_accession]
treated.counts$metabolon_number <- normalized.count.matrix$metabolon_number
treated.counts <- treated.counts  %>% relocate(metabolon_number, .before = 1)

untreated.androgen <- coldata[which(coldata$Treatment == "96 % Ethanol"),]

untreated.counts <- normalized.count.matrix[,untreated.androgen$geo_accession]
untreated.counts$metabolon_number <- normalized.count.matrix$metabolon_number
untreated.counts <- untreated.counts  %>% relocate(metabolon_number, .before = 1)
untreated.counts <- untreated.counts[,-7]





##==============================================================================
# --------------------- one function foe correlations --------------------------
#-------------------------------------------------------------------------------
##==============================================================================
library(ggplot2)
library(reshape2)

##==============================================================================
# --------------------- one function foe correlations --------------------------
#-------------------------------------------------------------------------------
##==============================================================================
library(ggplot2)
library(reshape2)

generate_heatmap_pdf <- function(count_matrix, pdf_filename) {
  # Load necessary libraries
  library(ggplot2)
  library(reshape2)
  
  # Split the data by metabolin
  metabolin_data <- split(count_matrix, count_matrix$metabolon_number)
  
  # Function to calculate R^2 matrix for a single metabolin
  calculate_r_squared <- function(metabolin_df) {
    metabolin_df <- metabolin_df[, -1]
    cor_matrix <- cor(metabolin_df, use = "pairwise.complete.obs")
    r_squared_matrix <- cor_matrix^2
    return(r_squared_matrix)
  }
  
  # Calculate R^2 matrices for each metabolin
  r_squared_list <- lapply(metabolin_data, calculate_r_squared)
  
  # Function to plot R^2 heatmap
  plot_r_squared_heatmap <- function(r_squared_matrix, metabolin) {
    melted_r_squared <- melt(r_squared_matrix)
    melted_r_squared$is_diagonal <- melted_r_squared$Var1 == melted_r_squared$Var2
    
    ggplot(melted_r_squared, aes(x = Var1, y = Var2)) +
      geom_tile(aes(fill = ifelse(is_diagonal, "Diagonal", 
                                  ifelse(value > 0.55, "Above 0.55", "Below 0.55"))), 
                color = "black", size = 0.5) +
      geom_text(aes(label = ifelse(value > 0.55, sprintf("%.2f", value), "")), 
                size = 4, color = "white") +
      scale_fill_manual(values = c("Above 0.55" = "#FF0000", 
                                   "Below 0.55" = "#FFFFFF",
                                   "Diagonal" = "#808080"), 
                        name = "R^2 Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12),
            legend.position = "bottom",
            legend.text = element_text(size = 10, face = "bold"),
            legend.title = element_text(size = 12),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            panel.grid.major = element_line(color = "grey80"),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1)) +
      labs(title = paste("Metabolon", metabolin),
           x = NULL, y = NULL) +
      coord_fixed(ratio = 1)
  }
  
  # Generate plots for each metabolin and save to PDF
  pdf(pdf_filename, width = 15, height = 15)
  
  for (metabolin in names(r_squared_list)) {
    print(plot_r_squared_heatmap(r_squared_list[[metabolin]], metabolin))
  }
  
  dev.off()
}


generate_heatmap_pdf(untreated.counts, "Heatmap_Corr_Prostate_untreated_5Samples.pdf")













