# 読み込み
library(dorothea)
library(viper)
library(dplyr)
library(here)
library(decoupleR)
library(tidyverse)

# ヒトの場合
data(dorothea_hs, package = "dorothea")

dorothea_regulon <- get_dorothea(organism = "human") %>%
  filter(confidence %in% c("A", "B", "C"))

expr_mat <- read.csv(here("data/normalizedCounts.csv"), row.names = 1)
sample_order <- c("BA_1", "BA_2", "BA_3","BA_4", "BA_5", "BA_6", "BA_7",
                  "Cholestasis_1", "Cholestasis_2", "Cholestasis_3", "Cholestasis_4", "Cholestasis_5")
sample_info <- data.frame(
  sample = sample_order,
  group = c(rep("BA", 7), rep("Cholestasis", 5))
)

expr_mat <- expr_mat[sample_order]

deseq_dt <- read.csv(here("results/deseq2_result.csv"), row.names = 1)

expr_df <- data.frame(expr_mat) %>% 
  tibble::rownames_to_column("ensembl") %>% 
  left_join(deseq_dt[, c("ensembl", "gene_name")]) %>% 
  mutate(gene = gene_name)

expr_df1 <- expr_df %>%
  distinct(gene, .keep_all = TRUE) %>%
  column_to_rownames(var = "gene") %>% 
  dplyr::select(-c(gene_name, ensembl))

tf_activity <- run_viper(
  mat = expr_df1,
  network = dorothea_regulon,  
  .source = "source",
  .target = "target",
  .likelihood = "likelihood"
)

tf_matrix <- tf_activity %>%
  select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>% 
  column_to_rownames("source")

# transposeして行列に（サンプル×TF）
tf_t <- t(tf_matrix)

# 群情報を統合
tf_df <- as.data.frame(tf_t) %>%
  rownames_to_column("sample") %>%
  left_join(sample_info, by = "sample")

# 群間比較（例：NFE2L2のt検定）
t.test_result <- t.test(NFE2L2 ~ group, data = tf_df)
print(t.test_result)


tf_names <- colnames(tf_df)[!(colnames(tf_df) %in% c("sample", "group"))]

result_list <- lapply(tf_names, function(tf) {
  t_res <- t.test(tf_df[[tf]] ~ tf_df$group)
  data.frame(
    TF = tf,
    p_value = t_res$p.value,
    mean_group1 = t_res$estimate[1],
    mean_group2 = t_res$estimate[2],
    diff = diff(t_res$estimate)
  )
})

tf_stats <- bind_rows(result_list) %>%
  arrange(p_value)
tf_stats$padj <- p.adjust(tf_stats$p_value, method = "BH")

# p値でソートして上位TFを確認
tf_stats[tf_stats$padj < 0.05, ]


pval_matrix <- tf_activity %>%
  select(source, condition, p_value) %>%
  pivot_wider(names_from = condition, values_from = p_value)

tf_wide <- tf_activity %>%
  pivot_wider(names_from = condition, values_from = score)



library(tidyverse)
library(pheatmap)

# 例として特定のTFだけ抽出
tfs_of_interest <- c("NFE2L2", "YAP1", "STAT3")  # Nrf2はNFE2L2で表記されること多いです

tf_matrix <- tf_activity %>%
  filter(source %in% tfs_of_interest) %>%
  select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source")

# NAがあれば0や平均などに補完（例は0に置換）
tf_matrix[is.na(tf_matrix)] <- 0

# ヒートマップ作成
pheatmap(tf_matrix,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         main = "TF Activity Scores for Selected TFs",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))





nrf2_targets <- c("HMOX1", "NQO1", "GCLC", "GCLM", "TXNRD1", "SLC7A11")

nrf2_module_score <- colMeans(expr_mat[nrf2_targets, ])

glutathione_genes <- c("GCLC", "GCLM", "GSS", "GSR", "GPX1", "GPX2","GPX4","GGT5","GSTT2", 
                       "GSTP1", "GSTA1", "GSTM1", "GGT1", "MGST1", "MGST2")


deseq_dt %>% filter(gene_name %in% glutathione_genes)

deseq_dt %>% filter(gene_name %in% c("GPX4", "SLC7A11", "SLC3A2", "ACSL4", "FTH1", "DHCR7", "BACH1", 
                                     "AIFM2", "SLC40A1", "HMOX1", "FTH1", "FTL1", "FSP1", "CHMP4B", "TMEM164")) 
