
# PCA
pca_plot <- function(df, rpca_var_percent) {
  ggplot(data= df, aes(x=PC1, y=PC2, color = group, shape = group)) +
    geom_point(size = 4, alpha=0.9) +
    geom_text_repel(aes(label=rownames(df)), size = 4, family = "Arial", max.overlaps = Inf, color = "black") +
    labs(x=paste0("PC1: ", rpca_var_percent[1], "%"),
         y=paste0("PC2: ", rpca_var_percent[2], "%")) +
    scale_color_manual(values = c("BA" = "gray70", "Cholestasis" = "gray30")) +
    scale_shape_manual(values = c("BA" = 16, "Cholestasis" = 17)) +
    theme_bw(base_family = "Arial") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "right",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)),
           shape = guide_legend(override.aes = list(size = 4)))
}

# Volcano plot

volcano_plot <- function(df, label_gene = NULL) {
  plot <- ggplot(data=df, aes(log2FoldChange, -log(padj, 10))) +
    geom_point(aes(color=Expression), size=2/5, alpha = 0.8, na.rm = TRUE) +
    scale_color_manual(values=c("#377EB8", "gray80", "#E41A1C")) +
    guides(color = guide_legend(override.aes = list(size = 3))) +  
    xlab(expression("log"[2]*" Fold Change")) +
    ylab(expression("-log"[10]*" Adjusted p-value")) +
    theme_minimal(base_family = "Arial")+
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "top",
      legend.text = element_text(size = 11),
      panel.grid = element_blank(), 
      panel.border = element_rect(color = "black", fill = NA)
    ) +
    scale_x_continuous(limits = c(-8, 8))
  
  # if label present
  if(!is.null(label_gene)) {
    plot <- plot +geom_text_repel(
      data = subset(df, gene_name %in% label_gene),
                    aes(label = gene_name),
                    size=4,
                    box.padding = 0.5,
                    max.overlaps = Inf,
                    segment.color = "grey50",
                    family = "Arial",
                    fontface = "bold")
  }
  
  return(plot)
}

# heatmap

plot_heatmap <- function(data, color, annotation_col, annotation_colors) {
  pheatmap(as.matrix(data),
           color = my_color,
           breaks = breaksList,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           scale = "row",
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           show_rownames = FALSE, 
           show_colnames = TRUE,
           fontsize = 10,
           border_color = NA,
           legend = TRUE,
           legend_labels = "Expression",
           fontsize_row = 6,
           fontsize_col = 10,
           family="Arial")
}

# boxplotを描画する関数
plot_boxplot <- function(deseq_res, tidy_data, target, output_dir){
  # 有意差を示す関数作成
  sig <- function(a){
    return(
      ifelse(a>0.05, "",
             ifelse(a>0.01, "*",
                    ifelse(a>0.001, "**","***")))
    )
  }
  
  p <- list()
  
  for (i in 1:length(target)){
    gene_ensembl <- target[i]
    gene_name <- names(target[i])
    
    # 有意差判定のためにpadj
    padj <- deseq_res %>% 
      filter(ensembl == gene_ensembl) %>% 
      dplyr::select(padj)
    
    # 有意差barの高さのためにyRoof
    yRoof <- tidy_data %>% 
      filter(ensembl == gene_ensembl) %>% 
      dplyr::select(expression) %>% 
      max()*1.2
    
    # plot
    p[[i]] <- tidy_data %>% 
      filter(ensembl == gene_ensembl) %>% 
      dplyr::select(-sample) %>% 
      group_by(diagnosis, ensembl) %>% 
      summarise_all(list(mean = mean, sd = sd)) %>% 
      ggplot(aes(x=diagnosis, y=mean, fill=diagnosis))+
      geom_bar(width=0.5, color="gray30", stat = "identity") +
      geom_errorbar(aes(ymin=mean - sd, ymax = mean + sd), color="black", width=0.2)+
      scale_fill_manual(values = c("BA" = "gray70", "Cholestasis" = "gray30"))+
      geom_signif(y_position = yRoof * 0.9, xmin = 1, xmax = 2, annotation = sig(padj),
                  textsize = 5, tip_length = 0.01, color="black") +
      labs(x = "", y = "Expression level", title = gene_name)+
      theme_classic(base_size = 14, base_family = "Arial")+
      theme(plot.title = element_text(hjust=0.5, face="bold", size=15),
            axis.title.x =element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text= element_text(size=12),
            axis.title=element_text(size=12),
            legend.position="none")
    
    ggsave(paste0(output_dir, "/", gene_name,".png"), plot =p[[i]],dpi=600, height=6, width=3)
  }
  
  return(p)
}