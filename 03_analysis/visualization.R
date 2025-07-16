# load library
library(ggplot2)
library(pheatmap)
library(viridis)
library(here)
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(ggsignif)
library(RColorBrewer)
library(patchwork)

# load deseq2 analysis results
source(here("scripts/03_analysis/DESeq2_analysis.R"))
# load plot function
source(here("scripts/03_analysis/plot_function.R"))

# Function to save plots
save_plot <- function(plot, filename, width = 7, height = 7, dpi = 600) {
  ggsave(filename = here("results/figures", filename), plot = plot, width = width, height = height, dpi = dpi)
}

# Function to get ensembl gene id from gmt file
getHeatmapData <- function(tidy_data, deseq_res, gmt_file_path) {
  # read gmt file
  gmt <- clusterProfiler::read.gmt(gmt_file_path)
  
  # get ensembl
  ENSEMBL <- deseq_res %>% 
    filter(gene_name %in% gmt$gene) %>% 
    pull(ensembl)
  
  # make df
  data <- tidy_data %>% 
    filter(ensembl %in% ENSEMBL) %>% 
    dplyr::select(-diagnosis) %>%
    pivot_wider(names_from = sample, values_from = expression) %>% 
    column_to_rownames(var = "ensembl")
  
  return(data)
}

# Plot total read counts per sample
barplot(colSums(dt)/1e6, main="Total read counts (millions)")

# Heatmap: sample distance
p_heatmap <- pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = as.dist(sample_dist_matrix),
  clustering_distance_cols = as.dist(sample_dist_matrix),
  clustering_method = "complete",
  col=viridis(100),
  border_color = NA, 
  legend_breaks = c(0, 0.09, 0.18),
  legend_labels = c("Close", "Medium", "Distant"),
  angle_col = 45,
  fontsize = 12,
  legend = TRUE,  
  family = "Arial"
  )

save_plot(p_heatmap, "heatmap.png")

# MAplot
plotMA(res, ylim = c(-5, 5))

# Count outliers. 
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

# Dispersion plot. 
plotDispEsts(dds)

# ---------------------
# PCA analysis
# ---------------------
normalized_count <- read.csv(here("data/normalizedCounts.csv"), row.names = 1)
data <- t(normalized_count)

# 分散を求めて、上位2000の遺伝子を抽出
score <- apply(data, 2, var) 
id <- sort.list(as.vector(score), decreasing = T)[1:2000] 

# 主成分分析を行う
rpca <- prcomp(x=data[, id], scale. = TRUE)

# summary
summary(rpca)

## data.frameを作成
rpca_df <- data.frame(rpca$x) %>% 
  mutate(group = if_else(str_detect(rownames(rpca$x), pattern = "BA"), "BA", "Cholestasis"))

# グラフのために%を取得
rpca_var <- rpca$sdev^2
rpca_var_percent <- round(rpca_var/sum(rpca_var)*100, digits = 1) 

# PC1, PC2で描画
p_pca <- pca_plot(df = rpca_df, rpca_var_percent = rpca_var_percent) 
# save plot
save_plot(plot = p_pca, filename = "pcaPlot.png", width = 7, height = 5) 

# ---------------------
# Volcano plot.
# ---------------------
deseq_res <- read.csv(here("results/deseq2_result.csv"), header = TRUE, row.names = 1)

# upregulated gene, downregulated geneを設定
log2FC_thresh <- log(2, 1.5) 
padj_thresh <- 0.1 

# 描画のためにexpression列を追加
deseq_res <- deseq_res %>% 
  mutate(Expression = case_when(log2FoldChange >= log2FC_thresh & padj <= padj_thresh ~ "Up-regulated",
                                log2FoldChange <= -log2FC_thresh & padj <= padj_thresh ~ "Down-regulated",
                                TRUE ~ "Unchanged"))

# labelをつけたい遺伝子をpick up
genes_to_label <- c("CYBB","NOX4","DUOX1","DUOX2","CAT",
                    "GPX1", "GSS","SOD1","SOD2","TXN", "CYGB",
                    "MAT1A", "GLS2","PKLR","MAT2A","GLS","PKM", 
                    "ESRP2", "SRSF3", "SLU7","HNF4A","HNF1A","CEBPA","FOXA2",
                    "FOXO1", "PCK1", "G6PC1")

# plot
save_plot(plot = volcano_plot(df = deseq_res), 
          filename = "volcanoPlot.png", height = 7, width = 5.5)
save_plot(plot = volcano_plot(df = deseq_res, label_gene = genes_to_label), 
          filename = "volcanoPlotwithLabel.png", height = 7, width = 5.5)

# ---------------------
# DEG heatmap
# ---------------------
deg <- deseq_res %>% 
  filter(Expression != "Unchanged") %>% 
  pull(ensembl)

# heatmap用のdata作成
data_heatmap <- data.frame(t(data)) %>% 
  filter(rownames(.) %in% deg) 
# 順番を入れ替える
desired_order <- c("BA_1","BA_2","BA_3","BA_4","BA_5","BA_6","BA_7",
                   "Cholestasis_1","Cholestasis_2","Cholestasis_3","Cholestasis_4","Cholestasis_5")

# z-score normalization (mean = 0, dispersion = 1), i.e. standardization
data_heatmap_z <- genefilter::genescale(data_heatmap, axis = 1, method = "Z") 

mat_orderd <- data_heatmap_z[, desired_order]

# 色設定
annotation_col <- data.frame(
  Group = factor(c("BA","BA","BA","BA","BA","BA","BA","Cholestasis",
                   "Cholestasis","Cholestasis","Cholestasis","Cholestasis"),
                 levels = c("BA", "Cholestasis"))
)
rownames(annotation_col) <- colnames(mat_orderd)

# アノテーションカラー設定
ann_colors <- list(Group = c("BA" = "#E41A1C", "Cholestasis" = "#4DAF4A"))

# カラースケール
breaksList <- seq(-2, 2, by = 0.01)
my_color <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(breaksList))

heatmap_plot <- plot_heatmap(data = mat_orderd, color = my_color, annotation_col = annotation_col, annotation_colors = ann_colors)

# save plot
save_plot(plot = heatmap_plot, filename = "DEGheatmap.png", width = 6, height = 9)

# ---------------------
# Create boxplot
# ---------------------
OUTPUT_DIR <- here("results/figures/boxplot")

if(!dir.exists(OUTPUT_DIR)){
  dir.create(OUTPUT_DIR)
}

# make tidy data
tidy_data <- normalized_count %>% 
  # rownamesをensembl列にする
  mutate(ensembl = rownames(.)) %>% 
  # tidy dataへ
  pivot_longer(cols = -ensembl, names_to = "sample", values_to = "expression") %>%
  # diagnosis列を作成
  mutate(diagnosis = case_when(
    str_detect(sample, pattern = "BA") ~ "BA",
    TRUE ~ "Cholestasis"
  ))


# legendのみ取り出す
legend_p <- tidy_data %>% 
  filter(ensembl == "ENSG00000079999") %>% 
  dplyr::select(-sample) %>% 
  group_by(diagnosis, ensembl) %>% 
  summarise_all(list(mean = mean, sd = sd)) %>% 
  ggplot(aes(x=diagnosis, y=mean, fill=diagnosis))+
  geom_bar(width=0.5, color="gray30", stat = "identity")+
  scale_fill_manual(values = c("BA" = "gray70", "Cholestasis" = "gray30"))+
  theme(legend.title = element_blank(),
        legend.text = element_text(family = "Arial", size=12))

# legendを取り出して保存
legend <- cowplot::get_legend(legend_p)
save_plot(plot = legend, filename = "legend.png", width = 7, height = 7) 

# genelists
GENE_LISTS <- list(NOX = c("ENSG00000165168","ENSG00000086991","ENSG00000137857", "ENSG00000140279"),
                   ROS = c("ENSG00000121691","ENSG00000233276","ENSG00000100983","ENSG00000142168", "ENSG00000291237","ENSG00000136810", "ENSG00000161544"),
                   ADULT_ISOZYME = c("ENSG00000151224", "ENSG00000135423", "ENSG00000143627"),
                   FETAL_ISOZYME = c("ENSG00000168906", "ENSG00000115419", "ENSG00000067225"),
                   SPLICING = c("ENSG00000103067", "ENSG00000112081", "ENSG00000164609"),
                   TRANSCRIPT = c("ENSG00000101076","ENSG00000135100", "ENSG00000245848", "ENSG00000125798"), 
                   GLUCOSE = c("ENSG00000150907", "ENSG00000101076", "ENSG00000124253", "ENSG00000131482")
                   )

names(GENE_LISTS$NOX) <- c("CYBB", "NOX4", "DUOX1", "DUOX2")
names(GENE_LISTS$ROS) <-  c("CAT","GPX1", "GSS","SOD1", "SOD2", "TXN", "CYGB")
names(GENE_LISTS$ADULT_ISOZYME) <- c("MAT1A", "GLS2", "PKLR")
names(GENE_LISTS$FETAL_ISOZYME) <- c("MAT2A", "GLS", "PKM")
names(GENE_LISTS$SPLICING) <- c("ESRP2", "SRSF3", "SLU7")
names(GENE_LISTS$TRANSCRIPT) <- c("HNF4A", "HNF1A", "CEBPA", "FOXA2")
names(GENE_LISTS$GLUCOSE) <- c("FOXO1", "HNF4A","PCK1", "G6PC1")

# plot
boxplots <- lapply(GENE_LISTS, function(genes) plot_boxplot (deseq_res = deseq_res, tidy_data = tidy_data, target = genes, output_dir = OUTPUT_DIR))

# 横並びのグラフ
# プロット保存用関数
save_merge_plot <- function(plot_list, filename) {
  merged_plot <- wrap_plots(plot_list) + plot_layout(nrow = 1)
  ggsave(filename = filename, plot = merged_plot, dpi = 600, height = 6, width = 3 * length(plot_list))
}

# save plot
for(i in 1:length(GENE_LISTS)) {
  save_merge_plot(plot_list = boxplots[[i]], 
                  filename = paste0(here("results/figures/boxplot/"), "merged_", names(GENE_LISTS[i]), ".png"))
}
