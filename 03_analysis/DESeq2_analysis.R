# import library
library(DESeq2)
library(tidyverse)
library(here)
library(rtracklayer)

# Function to rename columns using metadata mapping
rename_columns <- function(data, metadata) {
  name_map <- setNames(metadata$Sample.name, metadata$SRA_Label)
  colnames(data) <- name_map[colnames(data)]
  return(data)
}

# Function to create study design
create_study_design <- function(sample_names) {
  data.frame(row.names = sample_names, sample = sample_names) %>%
    mutate(condition = if_else(str_detect(sample, "Cholestasis"), "Cholestasis", "BA"))
}

# Read data and metadata
dt <- read.table(here("data/counts.txt"), header = TRUE, row.names = 1)
metadata <- read.table(here("data/samples.txt"), sep = "\t", header = TRUE)

# Apply column renaming function
dt <- rename_columns(dt, metadata)

# Generate study design
StudyDesign <- create_study_design(colnames(dt))

# create DESeqDataSetObject
dds <- DESeqDataSetFromMatrix(countData = dt, colData = StudyDesign, design = ~ condition)

# Prefiltering: remove low expression gene < 10 lead counts in total data
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Set reference level for condition
dds$condition <- relevel(dds$condition, ref = "Cholestasis")

# estimate size factor
dds <- estimateSizeFactors(dds)

# Run deseq2
dds <- DESeq(dds)
res <- results(dds)

# normalized data
ntd <- normTransform(dds)

# save normalized counts. 
write.csv(counts(dds, normalized = TRUE), file = here("data/normalizedCounts.csv"))

# gene annotation
gtf <- readGFF(here("reference/gencode.v46.chr_patch_hapl_scaff.annotation.gtf"))
gtf_gene <- gtf %>% subset(gtf$type == "gene") 
gtf_gene <- gtf_gene[,c("gene_id", "gene_type", "gene_name")]

# Define splitLeft function 
removeEnsemblVersion <- function(x) {
  return(strsplit(x, "[.]")[[1]][[1]])
}
gtf_gene[,"gene_id"] <- sapply(gtf_gene[,"gene_id"], removeEnsemblVersion)

# merge data
deseq_res <- data.frame(res) %>% 
  rownames_to_column(var = "ensembl") %>% 
  left_join(gtf_gene, c("ensembl" = "gene_id")) %>% 
  distinct(ensembl, .keep_all = T)

# save table
write.csv(deseq_res, file=here("results/deseq2_result.csv"))
write.csv(subset(deseq_res, padj < 0.05), file=here("results/filtered_deseq2_result.csv"))

# Function to calculate sample distances
calculate_sample_distances <- function(ntd) {
  # Log2 transformation to normalize large count differences
  log_data <- log2(assay(ntd) + 1)
  
  # Compute Pearson correlation and convert to distance matrix
  sample_cor <- cor(log_data, method = "pearson")
  sample_dist_matrix <- as.matrix(as.dist(1 - sample_cor))
  
  # Assign row and column names
  rownames(sample_dist_matrix) <- colnames(ntd)
  colnames(sample_dist_matrix) <- colnames(ntd)
  
  return(sample_dist_matrix)
}

# Compute sample distance matrix
sample_dist_matrix <- calculate_sample_distances(ntd)
