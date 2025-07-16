# load library
library(pathview)
library(here)
library(org.Hs.eg.db)

# Read DESeq2 results
deseq_res <- read.csv(here("results/deseq2_result.csv"), row.names = 1)

# Map Ensembl IDs to Entrez IDs
deseq_res$entrez <- mapIds(org.Hs.eg.db,
                          keys = deseq_res$ensembl,
                          column="ENTREZID",
                          keytype = "ENSEMBL",
                          multiVals = "first")

# Prepare log2 fold change data
logFC <- deseq_res$log2FoldChange
names(logFC) <- deseq_res$entrez

# Visualize KEGG pathway hsa00190
pathview(gene.data = logFC, pathway.id = "hsa00190")

# Visualize KEGG pathway hsa04146
pathview(gene.data = logFC, pathway.id = "hsa04146")
