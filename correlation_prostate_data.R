# load the required libraries
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
library(curatedTCGAData)
library(pheatmap)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(DESeq2)
library(magrittr)  # or library(dplyr)

count_matrix <- read.csv("/home/admin/Sung-prostate-cancer/prostate_matrix.csv")
rownames(count_matrix) <- count_matrix$X
count_matrix <- count_matrix[,-1]

############### Deseq2 Normalizations ####################

# my phenodata 
phdata <- colnames(count_matrix)
phdata <- as.data.frame(phdata)

# Assuming your data frame is called 'df'
phdata$disease_state <- ifelse(grepl("01A", phdata$phdata), "tumor", 
                               ifelse(grepl("11A",  phdata$phdata), "normal", NA))

names(phdata)[1] <- "ids"
rownames(phdata) <- phdata$ids

# making the rownames and column names identical
all(rownames(phdata) %in% colnames(count_matrix))
all(rownames(phdata) == colnames(count_matrix))

count_matrix <- round(count_matrix)
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = phdata,
                              design = ~ 1) # not spcifying model

## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 534,]
nrow(dds) #  10763 genes


# perform variance stabilization
dds_norm <- varianceStabilizingTransformation(dds)

# get normalized counts
normalized.count.matrix <- assay(dds_norm) %>% 
  t()

# norm.counts <- norm.counts[tumor$ids,]
normalized.count.matrix <- t(normalized.count.matrix) %>% as.data.frame()
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
# 12  14  11  10  12  14  13   8   4  27  17   9  14   7  18  19  16  18  16  13  31   6  31  10  19 

#-------------------------------------------------------------------------------
# 01 subset the first 50 tumor samples 
#-------------------------------------------------------------------------------

tumor <- phdata[which(phdata$disease_state == "tumor"),]

tumor.counts <- normalized.count.matrix[,tumor$ids]
tumor.counts$metabolon_number <- normalized.count.matrix$metabolon_number
tumor.counts <- tumor.counts %>% relocate(metabolon_number, .before = 1)
tumor.count.sub5<- tumor.counts[,c(1:6)]  # subse1t the (21:30, 31:40, 41:50, 51:60) 
ncol(tumor.count.sub5)
names(tumor.count.sub5)

library(ggplot2)
library(reshape2)
library(gridExtra)
#-------------------------------------------------------------------------------
# 02 subset the first 50 normal samples 
#-------------------------------------------------------------------------------

normal <- phdata[which(phdata$disease_state == "normal"),]

normal.counts <- normalized.count.matrix[,normal$ids]
normal.counts$metabolon_number <- normalized.count.matrix$metabolon_number
normal.counts <- normal.counts %>% relocate(metabolon_number, .before = 1)
normal.count.sub5<- normal.counts[,c(1:6)]  # subse1t the (21:30, 31:40, 41:50, 51:60) 
ncol(normal.count.sub5)
names(normal.count.sub5)

#-------------------------------------------------------------------------------
##==============================================================================
# --------------------- one function foe correlations --------------------------
#-------------------------------------------------------------------------------
##==============================================================================
library(ggplot2)
library(reshape2)

generate_heatmap_pdf <- function(count_matrix, pdf_filename) {
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
                size = 2, color = "white") +
      scale_fill_manual(values = c("Above 0.55" = "#FF0000", 
                                   "Below 0.55" = "#FFFFFF",
                                   "Diagonal" = "#808080"), 
                        name = "R^2 Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
            axis.text.y = element_text(size =10),
            legend.position = "bottom",
            legend.text = element_text(size = 8, face = "bold"),
            legend.title = element_text(size = 10),
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            panel.grid.major = element_blank(),
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

generate_heatmap_pdf(tumor.count.sub5, "Heatmap_Corr_5Prostate_Tumor_Samples.pdf")
generate_heatmap_pdf(normal.count.sub5, "Heatmap_Corr_5Prostate_Normal_Samples.pdf")

# subset 1:10 samples noraml and tumor samples 

tumor.count.sub1_51<- tumor.counts[,c(1:51)]  # subse1t the (21:30, 31:40, 41:50, 51:60) 
ncol(tumor.count.sub1_51)
names(tumor.count.sub1_51)

generate_heatmap_pdf(tumor.count.sub1_51, "Heatmap_Corr_1-50_Prostate_Tumor_Samples.pdf")


normal.count.sub41_50 <- normal.counts[,c(1,41:50)]  # subse1t the (21:30, 31:40, 41:50, 51:60) 
ncol(normal.count.sub41_50)
names(normal.count.sub41_50)
generate_heatmap_pdf(normal.count.sub41_50, "Heatmap_Corr_41-50_Prostate_Normal_Samples.pdf")













