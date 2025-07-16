library(fgsea)
library(here)
library(ggplot2)

# Read DESeq2 results
deseq_res <- read.csv(here("results/deseq2_result.csv"), row.names = 1)

# Prepare log2 fold change data
logFC <- deseq_res$log2FoldChange
names(logFC) <- deseq_res$gene_name

# logFC を降順にソート
ranked_genes <- sort(logFC, decreasing = TRUE)
# 重複する遺伝子名を削除（最初の出現を残す）
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]

# MSigDB の gene sets 
msigdb_files <- c(
  here("reference/BIOCARTA_FXR_PATHWAY.v2024.1.Hs.gmt"), 
  here("reference/PPARG_01.v2024.1.Hs.gmt"), 
  here("reference/WP_NRF2_PATHWAY.v2024.1.Hs.gmt")
  )  

# GMT ファイルをリストとして読み込む
msigdb_lists <- lapply(msigdb_files, read.gmt)

# パスウェイごとに遺伝子名だけを抽出し、fgseaに適した形式に変換
pathways <- lapply(msigdb_lists, function(x) x$gene)
names(pathways) <- c("BIOCARTA_FXR_PATHWAY", "PPARG_01", "WP_NRF2_PATHWAY")

# fgsea 解析の実行
fgsea_results <- fgsea(
  pathways = pathways,  # 各遺伝子セット
  stats = ranked_genes,  # 遺伝子のランキング（log2FoldChange など）
  minSize = 5,  # 最小遺伝子数
  maxSize = 500,  # 最大遺伝子数
  eps = 1e-20  # 計算精度
)

# plot関数
plot_gsea <- function(pathway_name) {
  # GSEA結果から NES と FDR を取得
  fgsea_res_filtered <- fgsea_results[fgsea_results$pathway == pathway_name, ]
  
  # NES と FDR を取得
  NES_value <- fgsea_res_filtered$NES
  FDR_value <- fgsea_res_filtered$padj
  
  # プロットを作成
  p <- plotEnrichment(pathways[[pathway_name]], ranked_genes) +
    ggtitle(sprintf("%s\nNES = %.2f, FDR = %.3f", pathway_name, NES_value, FDR_value)) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      panel.grid = element_blank(),
      axis.line = element_line(size = 1, color = "black")
      ) +
    scale_color_manual(values = c("red"))
  
  # プロット表示
  print(p)
  
  # save
  ggsave(plot = p, width = 8, height = 6, 
         filename = here(paste0("results/figures/", sprintf("Enrichment_%s.png", pathway_name)))
  )
}

# plot
plot_gsea(pathway_name = names(pathways)[1])
plot_gsea(pathway_name = names(pathways)[2])
plot_gsea(pathway_name = names(pathways)[3])
