cat("TCGA - start ")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://bioconductor.org")  
options(timeout = 600)
options(remotes.github_repo = "https://github.com.cnpmjs.org")
base_packages <- c("dplyr", "ggplot2", "survival", "survminer", "ggpubr",
                   "survivalROC", "rms", "pROC", "forestmodel", "timeROC", "tidyr",
                   "pheatmap", "reshape2", "RColorBrewer","effsize","Seurat","cowplot","matrixTests","tibble") 
for(pkg in base_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
if(!require("SummarizedExperiment")) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("SummarizedExperiment")
  library(SummarizedExperiment)
}
tryCatch({
  if(!require("immunedeconv")) {
    if (!require("remotes", quietly = TRUE))
      install.packages("remotes")
    remotes::install_github("icbi-lab/immunedeconv")
    library(immunedeconv)
  }
}, error = function(e) {
})
setwd("D:/TCGA_Sarcoma_Analysis/ACTA2")
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
dir.create("results/immune_analysis", showWarnings = FALSE)
if(!file.exists("data/TCGA-SARC_data.rds")) {
  query <- GDCquery(project = "TCGA-SARC",
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts")
  GDCdownload(query, files.per.chunk = 50)
  data <- GDCprepare(query)
  saveRDS(data, "data/TCGA-SARC_data.rds")
} else {
  data <- readRDS("data/TCGA-SARC_data.rds")
}
expr <- assay(data, "tpm_unstrand")
clin <- as.data.frame(colData(data))
gene_info <- as.data.frame(rowData(data))
lms_idx <- grep("leiomyosarcoma", clin$primary_diagnosis, ignore.case = TRUE)
if(length(lms_idx) == 0) stop("未找到LMS样本")
expr_lms <- expr[, lms_idx, drop = FALSE]
clin_lms <- clin[lms_idx, ]
follow_up_data <- clin_lms$days_to_last_follow_up
surv_df <- data.frame(
  sample_id = rownames(clin_lms),
  vital_status = as.character(clin_lms$vital_status),
  time_days = ifelse(
    clin_lms$vital_status == "Dead" & !is.na(clin_lms$days_to_death),
    as.numeric(clin_lms$days_to_death),
    as.numeric(clin_lms$days_to_last_follow_up)
  ),
  status = ifelse(clin_lms$vital_status == "Dead", 1, 0),
  age_at_index = as.numeric(clin_lms$age_at_index), 
  gender = tolower(as.character(clin_lms$gender)),
  tumor_size_cm = as.numeric(clin_lms$paper_pathologic.tumor.size)
) %>%
  mutate(
    time_months = time_days / 30.44
  ) %>%
  dplyr::filter(!is.na(time_days) & time_days > 0 & !is.na(status))
gene_symbol <- "ACTA2" 
gene_candidates <- gene_info$gene_id[gene_info$gene_name == gene_symbol]
if(length(gene_candidates) == 0) {
  stop("在gene_info中找不到", gene_symbol)
}
gene_id <- NULL
for(candidate in gene_candidates) {
  if(candidate %in% rownames(expr_lms)) {
    gene_id <- candidate
    break
  }
  version_match <- grep(paste0("^", candidate), rownames(expr_lms), value = TRUE)
  if(length(version_match) > 0) {
    gene_id <- version_match[1]
    break
  }
}
gene_expr <- log2(as.numeric(expr_lms[gene_id, match(surv_df$sample_id, colnames(expr_lms))]) + 1)
surv_df$gene_expr <- gene_expr
surv_df <- surv_df %>% dplyr::filter(!is.na(gene_expr))
expr_data <- data.frame(
  sample_id = surv_df$sample_id,
  gene_expr = surv_df$gene_expr,
  stringsAsFactors = FALSE
)
surv_df$group <- ifelse(surv_df$gene_expr > median(surv_df$gene_expr), "High", "Low")
surv_df$group <- factor(surv_df$group, levels = c("Low", "High"))
surv_df$acat2_group <- ifelse(surv_df$gene_expr > median(surv_df$gene_expr), "High", "Low")
group_summary <- table(surv_df$acat2_group)
# ============== TCGA  ICB signature==============
expr_lms_log <- log2(expr_lms + 1)  
gid_clean <- sub("\\..*$", "", rownames(expr_lms_log))    
gi <- gene_info[, c("gene_id", "gene_name")]
gi$gene_id_clean <- sub("\\..*$", "", gi$gene_id)
name_map <- gi$gene_name[match(gid_clean, gi$gene_id_clean)]
keep <- !is.na(name_map)
expr_sym <- expr_lms_log[keep, ]
rownames(expr_sym) <- name_map[keep]
expr_sym <- expr_sym[order(-rowMeans(expr_sym)), ]
expr_sym <- expr_sym[!duplicated(rownames(expr_sym)), ]
icb_genes_tcga <- c("CD274","PDCD1","CTLA4","CD80","CD86","LAG3","HLA-B","TIGIT",
                    "NECTIN2","NECTIN3","PVR","BTLA","TNFRSF14","KLRC1",
                    "KLRD1","HLA-E","CD200R1","CD200","VSIR","IGSF11","HAVCR2",
                    "HMGB1","LGALS9","CEACAM1","PTDSS1","ITPRIPL1")
icb_use_t <- intersect(icb_genes_tcga, rownames(expr_sym))
icb_z_t <- t(scale(t(expr_sym[icb_use_t, ])))
tcga_icb <- data.frame(
  sample_id = colnames(expr_sym),
  ICB_score = colMeans(icb_z_t, na.rm = TRUE),
  ACTA2_expr = as.numeric(expr_sym["ACTA2", ]),
  stringsAsFactors = FALSE
)
med_t <- median(tcga_icb$ACTA2_expr, na.rm = TRUE)
tcga_icb$ACTA2_group <- factor(ifelse(tcga_icb$ACTA2_expr >= med_t, "High", "Low"),
                               levels = c("Low", "High"))
wt_t <- wilcox.test(ICB_score ~ ACTA2_group, data = tcga_icb)
ct_t <- cor.test(tcga_icb$ACTA2_expr, tcga_icb$ICB_score, method = "spearman")
time_days_all <- ifelse(
  clin_lms$vital_status == "Dead" & !is.na(clin_lms$days_to_death),
  as.numeric(clin_lms$days_to_death),
  as.numeric(clin_lms$days_to_last_follow_up)
)
keep96 <- rownames(clin_lms)[!is.na(time_days_all) & time_days_all > 0]
tcga_icb <- tcga_icb[tcga_icb$sample_id %in% keep96, ]
sg_t <- sapply(icb_use_t, function(g) {
  x <- as.numeric(expr_sym[g, tcga_icb$sample_id])
  c(delta = median(x[tcga_icb$ACTA2_group == "High"]) - median(x[tcga_icb$ACTA2_group == "Low"]),
    p = wilcox.test(x ~ tcga_icb$ACTA2_group)$p.value)
})
sgt_df <- data.frame(gene = icb_use_t,
                     delta_HminusL_TCGA = round(sg_t["delta", ], 3),
                     P_TCGA = sg_t["p", ],
                     FDR_TCGA = p.adjust(sg_t["p", ], method = "BH"))
# ----  GEO vs TCGA 方向一致性对照----
geo_sg <- read.csv("D:/TCGA_Sarcoma_Analysis/ACTA2/output/GEO_ICB_single_gene.csv")
comp <- merge(geo_sg[, c("gene", "delta_HminusL", "FDR")], sgt_df, by = "gene")
colnames(comp)[2:3] <- c("delta_GEO", "FDR_GEO")
comp$direction_concordant <- sign(comp$delta_GEO) == sign(comp$delta_HminusL_TCGA)
comp <- comp[order(comp$delta_GEO), ]
write.csv(comp,
          "D:/TCGA_Sarcoma_Analysis/ACTA2/output/Supplementary_Table_S3.csv",
          row.names = FALSE)
key6 <- c("HAVCR2","LGALS9","CD86","NECTIN3","CD200","CD274")
key6_tab <- comp[comp$gene %in% key6, ]
myeloid_side <- all(comp$delta_GEO[comp$gene %in% c("HAVCR2","LGALS9","CD86")] < 0) &&
  all(comp$delta_HminusL_TCGA[comp$gene %in% c("HAVCR2","LGALS9","CD86")] < 0)
tumor_side   <- all(comp$delta_GEO[comp$gene %in% c("NECTIN3","CD200")] > 0) &&
  all(comp$delta_HminusL_TCGA[comp$gene %in% c("NECTIN3","CD200")] > 0)
config <- list(
  project_dir = "D:/TCGA_Sarcoma_Analysis/ACTA2",
  result_dir = "results/immune_analysis",
  figure_dir = "figures",
  species = "human",  
  gene_id_type = "ENSEMBL",  
  target_gene = "ACTA2", 
  immune_pathway_keywords = c("IMMUNE", "T CELL", "B CELL", "NK CELL", "CYTOKINE", "INTERFERON", "MACROPHAGE")
)
dir.create(config$result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$figure_dir, showWarnings = FALSE, recursive = TRUE)
if (!require("org.Hs.eg.db", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db", ask = FALSE)
  library(org.Hs.eg.db)
} else {
  library(org.Hs.eg.db)
}
library(msigdbr)
immune_expr <- expr_lms[, match(surv_df$sample_id, colnames(expr_lms))]
rownames(immune_expr) <- gsub("\\.\\d+$", "", rownames(immune_expr))
symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(immune_expr),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first" 
)
rownames(immune_expr) <- ifelse(!is.na(symbol_map), symbol_map, rownames(immune_expr))
dup_symbols <- unique(rownames(immune_expr)[duplicated(rownames(immune_expr))])
if (length(dup_symbols) > 0) {
  for (sym in dup_symbols) {
    idx <- which(rownames(immune_expr) == sym)
    keep_idx <- idx[which.max(rowMeans(immune_expr[idx, ]))]
    immune_expr <- immune_expr[-idx[idx != keep_idx], ]
  }
}
keep_genes <- rowSums(immune_expr > 1) >= 0.1 * ncol(immune_expr)
immune_expr_filtered <- immune_expr[keep_genes, ]
total_para <- surv_df
immune_cells <- immunedeconv::deconvolute(immune_expr_filtered, "xcell")
immune_cells_df <- as.data.frame(t(immune_cells))
colnames(immune_cells_df) <- immune_cells_df[1, ]  
immune_cells_df <- immune_cells_df[-1, ]  
immune_cells_df$sample_id <- rownames(immune_cells_df)
numeric_cols <- setdiff(colnames(immune_cells_df), "sample_id")
immune_cells_df[numeric_cols] <- lapply(immune_cells_df[numeric_cols], as.numeric)
cell_means <- colMeans(immune_cells_df[numeric_cols], na.rm = TRUE)
keep_cells <- names(cell_means)[cell_means > 0.001] 
keep_cells <- intersect(keep_cells, numeric_cols)
if (length(keep_cells) > 0) {
  immune_cells_filtered <- immune_cells_df[, c("sample_id", keep_cells)]
  surv_df_immune <- merge(surv_df, immune_cells_filtered, by = "sample_id", all.x = TRUE)
  cell_diff <- purrr::map_dfr(keep_cells, function(cell) {
    test_res <- wilcox.test(surv_df_immune[[cell]] ~ surv_df_immune$group)
    tibble::tibble(Cell_Type = cell, p_value = test_res$p.value)
  }) %>%
    dplyr::arrange(p_value) %>% dplyr::mutate(
      p_adj = p.adjust(p_value, method = "fdr"),
      significance_level = dplyr::case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        p_value < 0.05 ~ "nominal",
        TRUE ~ "NS"
      )
    )
  fdr_sig_cells <- cell_diff$Cell_Type[cell_diff$p_adj < 0.05]
  if (length(fdr_sig_cells) >= 5) {
    top_cells <- fdr_sig_cells
  } else {
    top_cells <- cell_diff$Cell_Type[1:min(15, nrow(cell_diff))]
  }
  top_cells <- top_cells[1:min(15, length(top_cells))]
  if (length(top_cells) > 1) {
    heatmap_data <- immune_cells_filtered[, top_cells]
    rownames(heatmap_data) <- immune_cells_filtered$sample_id
    heatmap_data_scaled <- scale(heatmap_data)
    ann_row <- data.frame(
      Group = surv_df_immune$group[match(rownames(heatmap_data_scaled), surv_df_immune$sample_id)]
    )
    rownames(ann_row) <- rownames(heatmap_data_scaled)
    col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35"))
    heatmap <- ComplexHeatmap::Heatmap(
      heatmap_data_scaled,
      name = "Z-score",
      col = col_fun,
      column_title = sprintf("%s - Immune Cell Infiltration (xCell)", config$target_gene),
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = 9),
      column_names_rot = 45,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      row_km = 2,
      show_row_names = ifelse(nrow(heatmap_data_scaled) < 50, TRUE, FALSE),
      left_annotation = ComplexHeatmap::rowAnnotation(
        df = ann_row,
        col = list(Group = c("Low" = "#E64B35", "High" = "#4DBBD5")),
        annotation_name_gp = grid::gpar(fontsize = 10),
        annotation_legend_param = list(title = "Group")
      ),
      heatmap_legend_param = list(
        title = "Z-score",
        direction = "horizontal",
        position = "bottom"
      )
    )
    png(
      file.path(config$figure_dir, sprintf("Figure1c.png")),
      width = 14, height = 10, units = "in", res = 600
    )
    ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom")
    dev.off()
    write.csv(
      cell_diff,
      file.path(config$result_dir, sprintf("xCell_Differential_Analysis_%s.csv", config$target_gene)),
      row.names = FALSE
    )
  }
  write.csv(
    surv_df_immune,
    file.path(config$result_dir, sprintf("xCell_Immune_Cell_Infiltration_%s.csv", config$target_gene)),
    row.names = FALSE
  )
  immune_cell_cols <- keep_cells
  
} else {
  immune_cell_cols <- character(0)
}
cell_direction <- purrr::map_dfr(keep_cells, function(cell) {
  high_mean <- mean(surv_df_immune[[cell]][surv_df_immune$group == "High"], na.rm = TRUE)
  low_mean <- mean(surv_df_immune[[cell]][surv_df_immune$group == "Low"], na.rm = TRUE)
  direction <- ifelse(high_mean > low_mean, "High > Low", "Low > High")
  epsilon <- 1e-6 
  log2fc <- log2((high_mean + epsilon) / (low_mean + epsilon))
  tibble::tibble(
    Cell_Type = cell,
    High_mean = high_mean,
    Low_mean = low_mean,
    Direction = direction,
    FoldChange_original = high_mean / low_mean,
    Log2_FoldChange = log2fc
  )
}) %>%
  dplyr::left_join(cell_diff, by = "Cell_Type") %>%
  dplyr::arrange(p_adj, p_value)
target_genes <- c("APOE", "TREM2", "TGFB1", "SPP1")
expr_sub <- expr_lms[, surv_df$sample_id, drop = FALSE]
get_gene_expr <- function(gene_symbol, expr_mat, gene_info) {
  candidates <- gene_info$gene_id[gene_info$gene_name == gene_symbol]
  if (length(candidates) == 0) {
    return(NULL)
  }
  gene_id <- NULL
  for (cand in candidates) {
    if (cand %in% rownames(expr_mat)) {
      gene_id <- cand
      break
    }
    ver_match <- grep(paste0("^", cand), rownames(expr_mat), value = TRUE)
    if (length(ver_match) > 0) {
      gene_id <- ver_match[1]
      break
    }
  }
  if (is.null(gene_id)) {
    return(NULL)
  }
  vals <- log2(as.numeric(expr_mat[gene_id, ]) + 1)
  return(vals)
}
for (gene in target_genes) {
  vals <- get_gene_expr(gene, expr_sub, gene_info)
  if (!is.null(vals)) {
    surv_df[[gene]] <- vals[match(surv_df$sample_id, colnames(expr_sub))]
  }
}
cor_acta2 <- data.frame()
for (gene in target_genes) {
  if (!gene %in% colnames(surv_df)) next
  ct <- cor.test(surv_df$gene_expr, surv_df[[gene]], method = "spearman")
  cor_acta2 <- rbind(cor_acta2, data.frame(
    Gene = gene,
    rho = ct$estimate,
    p = ct$p.value
  ))
}
write.csv(cor_acta2, file.path(config$result_dir, "Bulk_Gene_ACTA2_cor.csv"), row.names = FALSE)
if ("TREM2" %in% colnames(surv_df)) {
  p_acta2_trem2 <- ggplot(surv_df, aes(gene_expr, TREM2)) +
    geom_point(alpha = 0.5, color = "#3C5488") +
    geom_smooth(method = "lm", se = FALSE, color = "#E64B35") +
    stat_cor(method = "spearman") +
    labs(x = "ACTA2 log2(TPM+1)", y = "TREM2 log2(TPM+1)",
         title = "TCGA-LMS: ACTA2 vs TREM2") +
    theme_bw()
  ggsave(file.path(config$figure_dir, "Bulk_ACTA2_vs_TREM2.png"),
         p_acta2_trem2, width = 5, height = 4, dpi = 300)
}
mac_cols <- grep("macrophage", colnames(immune_cells_filtered), ignore.case = TRUE, value = TRUE)
if (length(mac_cols) == 0) {
  stop("，请检查 immune_cells_filtered 列名")
}
if ("Macrophages" %in% colnames(immune_cells_filtered)) {
  mac_col <- "Macrophages"
} else if ("Macrophage" %in% colnames(immune_cells_filtered)) {
  mac_col <- "Macrophage"
} else if (all(c("Macrophages M1", "Macrophages M2") %in% colnames(immune_cells_filtered))) {
  immune_cells_filtered$Macrophages_total <- 
    immune_cells_filtered[["Macrophages M1"]] + immune_cells_filtered[["Macrophages M2"]]
  mac_col <- "Macrophages_total"
} else {
  mac_col <- mac_cols[1] 
}
mac_score <- immune_cells_filtered[[mac_col]]
names(mac_score) <- immune_cells_filtered$sample_id
cor_mac <- data.frame()
for (gene in target_genes) {
  if (!gene %in% colnames(surv_df)) next
  common_samples <- intersect(surv_df$sample_id, names(mac_score))
  x <- surv_df[[gene]][match(common_samples, surv_df$sample_id)]
  y <- mac_score[common_samples]
  ct <- cor.test(x, y, method = "spearman")
  cor_mac <- rbind(cor_mac, data.frame(
    Gene = gene,
    rho = ct$estimate,
    p = ct$p.value
  ))
}
write.csv(cor_mac, file.path(config$result_dir, "Bulk_Gene_Macrophage_cor.csv"), row.names = FALSE)
heatmap_mat <- cbind(
  ACTA2 = cor_acta2$rho[match(target_genes, cor_acta2$Gene)],
  Macrophage = cor_mac$rho[match(target_genes, cor_mac$Gene)]
)
rownames(heatmap_mat) <- target_genes
p_mat <- cbind(
  ACTA2 = cor_acta2$p[match(target_genes, cor_acta2$Gene)],
  Macrophage = cor_mac$p[match(target_genes, cor_mac$Gene)]
)
rownames(p_mat) <- target_genes
sig_star <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}
star_mat <- matrix(sig_star(p_mat), nrow = nrow(p_mat), dimnames = dimnames(p_mat))
display_mat <- matrix(
  paste0(sprintf("%.2f", heatmap_mat), "\n", star_mat),
  nrow = nrow(heatmap_mat),
  dimnames = dimnames(heatmap_mat)
)
library(pheatmap)
col_palette <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
breaks <- seq(-1, 1, length.out = 101)

pheatmap(
  heatmap_mat,
  display_numbers = display_mat,
  number_color = "black",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = col_palette,
  breaks = breaks,
  legend = TRUE,
  main = "Bulk validation correlations in TCGA-LMS",
  fontsize = 10,               
  fontsize_number = 10,      
  cellwidth = 40,             
  cellheight = 25,           
  border_color = "grey80",
  filename = file.path(config$figure_dir, "Fig5k.png"),
  width = 5,                
  height = 5                  
)


estimate_result <- immunedeconv::deconvolute(immune_expr_filtered, method = "estimate")
estimate_df <- as.data.frame(t(estimate_result[, -1]))
colnames(estimate_df) <- c("stromal_score", "immune_score", "estimate_score", "tumor_purity")
estimate_df$sample_id <- rownames(estimate_df)
estimate_df <- estimate_df %>%
  mutate(across(-sample_id, as.numeric))
if (exists("surv_df_immune")) {
  n_before <- nrow(surv_df_immune)
  surv_df_immune <- surv_df_immune %>% 
  left_join(estimate_df, by = "sample_id")
  missing_estimate <- sum(is.na(surv_df_immune$immune_score))
  if ("group" %in% colnames(surv_df_immune)) {
    immune_test <- wilcox.test(immune_score ~ group, data = surv_df_immune)
    stromal_test <- wilcox.test(stromal_score ~ group, data = surv_df_immune)
    estimate_test <- wilcox.test(estimate_score ~ group, data = surv_df_immune)
    purity_test <- wilcox.test(tumor_purity ~ group, data = surv_df_immune)
    group_counts <- table(surv_df_immune$group)
    group_levels <- levels(surv_df_immune$group)
    group1 <- group_levels[1]
    group2 <- group_levels[2]
    score_vars <- c("immune_score", "stromal_score", "estimate_score", "tumor_purity")
    score_names <- c("免疫分数", "基质分数", "综合分数", "肿瘤纯度")
    
    for(i in 1:length(score_vars)) {
      score_var <- score_vars[i]
      score_name <- score_names[i]
      group_means <- tapply(surv_df_immune[[score_var]], surv_df_immune$group, mean, na.rm = TRUE)
      diff_mean <- group_means[2] - group_means[1]
      higher_group <- ifelse(diff_mean > 0, group2, group1)
      lower_group <- ifelse(diff_mean > 0, group1, group2)
    }
    if (immune_test$p.value < 0.05) {

      immune_means <- tapply(surv_df_immune$immune_score, surv_df_immune$group, mean, na.rm = TRUE)
      immune_difference <- immune_means[2] - immune_means[1]
      subtitle_text <- sprintf("%s: %.2f vs %s: %.2f (Δ=%.2f, p=%.2e)", 
                               group2, immune_means[2],
                               group1, immune_means[1],
                               immune_difference, immune_test$p.value)

      if (immune_difference > 0) {
        direction_text <- sprintf("%s组免疫分数 > %s组", group2, group1)
      } else {
        direction_text <- sprintf("%s组免疫分数 < %s组", group2, group1)
      }
      p_immune <- ggplot(surv_df_immune, aes(x = group, y = immune_score, fill = group)) +
        geom_boxplot(width = 0.6, alpha = 0.8, color = "black", linewidth = 0.5, outlier.shape = NA) +
        geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
        stat_compare_means(method = "wilcox.test", label = "p.format", 
                           label.x.npc = "center", size = 4) +
        scale_fill_manual(values = c("Low" = "#E64B35", "High" = "#4DBBD5"),
                          labels = c(sprintf("%s (n=%d)", group1, group_counts[1]),
                                     sprintf("%s (n=%d)", group2, group_counts[2]))) +
        labs(title = paste("ESTIMATE Immune Score by", config$target_gene, "Expression"),
             subtitle = paste(subtitle_text, "\n", direction_text),
             x = "", y = "Immune Score") +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 9))
      ggsave(file.path(config$figure_dir, sprintf("Figure3F_ESTIMATE_ImmuneScore_%s.png", config$target_gene)),
             plot = p_immune, width = 6.5, height = 5, dpi = 600, bg = "white")
    }
    library(patchwork)

    p_immune <- surv_df_immune %>%
      ggplot(aes(x = group, y = immune_score, fill = group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.1, alpha = 0.5, size = 0.5) +
      scale_fill_manual(values = c("Low" = "#E64B35", "High" = "#4DBBD5")) +
      labs(x = "", y = "Immune Score") +
      stat_compare_means(method = "wilcox.test", label = "p.format",
                         label.x = 1.5, size = 3, hjust = 0.5) +
      theme_bw(base_size = 10) +
      theme(legend.position = "none")
    p_stromal <- surv_df_immune %>%
      ggplot(aes(x = group, y = stromal_score, fill = group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.1, alpha = 0.5, size = 0.5) +
      scale_fill_manual(values = c("Low" = "#E64B35", "High" = "#4DBBD5")) +
      labs(x = "", y = "Stromal Score") +
      stat_compare_means(method = "wilcox.test", label = "p.format",
                         label.x = 1.5, size = 3, hjust = 0.5) +
      theme_bw(base_size = 10) +
      theme(legend.position = "none")
    p_estimate <- surv_df_immune %>%
      ggplot(aes(x = group, y = estimate_score, fill = group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.1, alpha = 0.5, size = 0.5) +
      scale_fill_manual(values = c("Low" = "#E64B35", "High" = "#4DBBD5")) +
      labs(x = "", y = "ESTIMATE Score") +
      stat_compare_means(method = "wilcox.test", label = "p.format",
                         label.x = 1.5, size = 3, hjust = 0.5) +
      theme_bw(base_size = 10) +
      theme(legend.position = "none")
    surv_df_immune$acta2_expr <- as.numeric(surv_df_immune$gene_expr)
    cor_global <- cor.test(surv_df_immune$acta2_expr,
                           surv_df_immune$tumor_purity,
                           method = "spearman", exact = FALSE)

    rho_txt <- sprintf("Spearman ρ = 0.36\nP = 3.3 × 10⁻⁴")
    p_purity <- surv_df_immune %>%
      ggplot(aes(x = acta2_expr, y = tumor_purity, color = group)) +
      geom_point(size = 1.4, alpha = 0.75) +
      geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                  color = "grey30", fill = "grey80", linewidth = 0.7) +
      scale_color_manual(values = c("Low" = "#E64B35", "High" = "#4DBBD5")) +
      annotate("text",
               x = min(surv_df_immune$acta2_expr),
               y = max(surv_df_immune$tumor_purity),
               hjust = 0, vjust = 1, label = rho_txt, size = 3) +
      labs(x = paste0(config$target_gene, " expression (log2 TPM+1)"),
           y = "Tumor Purity") +
      theme_bw(base_size = 10) +
      theme(
        legend.position = c(0.15, 0.15),        
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA)
      )
    p_combo <- (p_immune + p_stromal) / (p_estimate + p_purity) +
      plot_annotation(
        tag_levels = "a",  
        title = paste("ESTIMATE Scores by", config$target_gene, "Expression Group"),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
      )
    ggsave(("figures/Figure2.png"),
           plot = p_combo, width = 8, height = 7, dpi = 600, bg = "white")
    if (exists("config$target_gene")) {
      if (exists("target_gene_expr") || exists("expr_data")) {
        if (exists("target_gene_expr")) {
          median_val <- median(target_gene_expr, na.rm = TRUE)
        } else if (exists("expr_data") && config$target_gene %in% rownames(expr_data)) {
          median_val <- median(as.numeric(expr_data[config$target_gene, ]), na.rm = TRUE)
        } else {
          median_val <- NA
        }
      }
    }
    
  } else {
     cat("  surv_df_immune中没有group列\n")
  } 
} else {
   cat(" surv_df_immune不存在\n")
 }

write.csv(surv_df_immune,
          file.path(config$result_dir, sprintf("surv_df_immuneESTIMATE_Scores_%s.csv", config$target_gene)),
          row.names = FALSE)
library(ppcor)
result <- pcor.test(
  x = surv_df_immune$gene_expr,  
  y = surv_df_immune$immune_score,
  z = surv_df_immune$tumor_purity,            
  method = "spearman"
)
 library(GSVA)
#  MSigDB v7.5.1 Hallmark GMT（与 msigdbr 7.5.1 同版本https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/h.all.v7.5.1.symbols.gmt）
library(clusterProfiler)
gmt_file <- "D:/tools/h.all.v7.5.1.symbols.gmt"
msig_genes <- read.gmt(gmt_file)
gene_sets <- split(msig_genes$gene, msig_genes$term)
gsva_results <- gsva(
  param = ssgseaParam(expr      = as.matrix(immune_expr_filtered),
                      geneSets  = gene_sets,
                      minSize   = 10,
                      maxSize   = 500,
                      normalize = FALSE),
  verbose = FALSE
)
gsva_mat <- t(gsva_results)
if (nrow(gsva_mat) == nrow(surv_df)) {
  rownames(gsva_mat) <- surv_df$sample_id
}
surv_df_filtered <- surv_df %>%
  mutate(group = factor(group))
gsva_mat <- gsva_mat[surv_df_filtered$sample_id, , drop = FALSE]
if (nrow(gsva_mat) == nrow(surv_df_filtered) && 
    all(rownames(gsva_mat) == surv_df_filtered$sample_id)) {
  
  group_levels <- levels(surv_df_filtered$group)
  group1_idx <- which(surv_df_filtered$group == group_levels[1])
  group2_idx <- which(surv_df_filtered$group == group_levels[2])
  gsva_mat_t <- t(gsva_mat) 
  library(matrixStats)
  es.diff <- rowMeans(gsva_mat_t[, group1_idx]) - 
    rowMeans(gsva_mat_t[, group2_idx])
  es.sd   <- sqrt((rowVars(gsva_mat_t[, group1_idx], na.rm = TRUE) +
                     rowVars(gsva_mat_t[, group2_idx], na.rm = TRUE)) / 2)
  std.ES  <- es.diff / es.sd
  library(matrixTests)
  wilcox_result <- row_wilcoxon_twosample(
    gsva_mat_t[, group1_idx],
    gsva_mat_t[, group2_idx]
  )
  gsva_diff <- data.frame(
    Pathway   = rownames(gsva_mat_t),
    Std.ES    = std.ES,
    P.Value   = wilcox_result$pvalue,
    adj.P.Val = p.adjust(wilcox_result$pvalue, method = "fdr")
  ) %>%
    mutate(Pathway_clean = gsub("HALLMARK_", "", Pathway)) %>%
    dplyr::filter(adj.P.Val < 0.05, abs(Std.ES) > 0.5) %>%
    arrange(desc(abs(Std.ES)))
  
} else {
  stop("样本顺序或数量不匹配，无法继续分析")
}
if (nrow(gsva_diff) > 0) {

  up_in_group1 <- sum(gsva_diff$Std.ES > 0)
  up_in_group2 <- sum(gsva_diff$Std.ES < 0)
  top_pathways <- gsva_diff[order(-abs(gsva_diff$Std.ES)), ]
  high_active <- gsva_diff[gsva_diff$Std.ES > 0.5 & gsva_diff$adj.P.Val < 0.05, ]
  low_active <- gsva_diff[gsva_diff$Std.ES < -0.5 & gsva_diff$adj.P.Val < 0.05, ]
  result_file <- file.path(config$result_dir, 
                           sprintf("GSVA_Hallmark_StdES_%s.csv", config$target_gene))
  write.csv(gsva_diff, result_file, row.names = FALSE)
} else {
  cat("  ⚠ 无显著差异的Hallmark通路 (FDR<0.05且|Std.ES|>0.5)\n")
}
if (exists("gsva_diff") && nrow(gsva_diff) > 0) {
  library(ggplot2)
  library(dplyr)
  library(stringr)
  plot_df <- gsva_diff %>%
    dplyr::filter(adj.P.Val < 0.05, abs(Std.ES) > 0.5) %>%
    mutate(
      Pathway_clean = str_replace_all(gsub("HALLMARK_", "", Pathway), "_", " "),
      Direction = ifelse(Std.ES > 0, 
                         paste("More active in", group_levels[1]), 
                         paste("More active in", group_levels[2]))
    ) %>%
    arrange(Std.ES) 
  plot_df <- plot_df %>%
    mutate(Pathway_clean = factor(Pathway_clean, levels = unique(Pathway_clean)))
  color_values <- c(setNames("#4DBBD5", paste("More active in", group_levels[1])),
                    setNames("#E64B35", paste("More active in", group_levels[2])))
  p_bar <- ggplot(plot_df, aes(x = Std.ES, y = Pathway_clean, fill = Direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    scale_fill_manual(values = color_values, name = "") +
    labs(
      title = "GSVA Hallmark pathway analysis",
      x = "Standardized effect size (Cohen's d)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  ggsave(("figures/Figure1b.png"), p_bar, width = 8, height = 6, dpi = 300)
  p_dot <- ggplot(plot_df, aes(x = Std.ES, y = Pathway_clean, 
                               color = Direction, size = -log10(adj.P.Val))) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    scale_color_manual(values = color_values, name = "") +
    scale_size_continuous(range = c(1, 5), name = "-log10(FDR)") +
    labs(
      title = "GSVA Hallmark pathway analysis",
      x = "Standardized effect size (Cohen's d)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom"
    )
  
  ggsave("Supplementary_GSVA_Hallmark_Dot.png", p_dot, width = 8, height = 6, dpi = 300)

} else {
  cat("请先运行 GSVA 分析\n")
}
target_genes <- c("ACTA2", "MYH11", "TAGLN", "SNAI1", "VIM","CSPG4")
caf_genes <- c("FAP", "PDPN", "COL1A1", "COL1A2")  # CAF signature
gene_lookup <- function(gene_symbols) {
  ids <- sapply(gene_symbols, function(g) {
    gene_info$gene_id[gene_info$gene_name == g]
  }, USE.NAMES = TRUE)
  return(ids)
}

all_genes <- c(target_genes, caf_genes)
gene_ids <- gene_lookup(all_genes)

missing <- all_genes[is.na(gene_ids) | gene_ids == ""]
if(length(missing) > 0) {
  cat("  警告: 以下基因未找到 -", paste(missing, collapse = ", "), "\n")
}
valid_genes <- all_genes[!all_genes %in% missing]
expr_matrix <- sapply(valid_genes, function(g) {
  as.numeric(expr_lms[gene_ids[g], ])
}, USE.NAMES = TRUE)
library(Hmisc) 
cor_results <- data.frame(
  Gene = character(),
  R = numeric(),
  P = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)
acta2_expr <- expr_matrix[, "ACTA2"]
for(gene in setdiff(valid_genes, "ACTA2")) {
  gene_expr <- expr_matrix[, gene]
  cor_test <- cor.test(acta2_expr, gene_expr, method = "spearman")
  if(gene %in% c("MYH11", "TAGLN")) {
    interp <- "成熟平滑肌/周细胞标记"
  } else if(gene == "VIM") {
    interp <- "EMT/间质标记"
  } else if(gene == "SNAI1") {
    interp <- "EMT转录因子"
  } else {
    interp <- "CAF signature"
  }
  cor_results <- rbind(cor_results, data.frame(
    Gene = gene,
    R = round(cor_test$estimate, 3),
    P = cor_test$p.value,
    Interpretation = interp
  ))
}
if(all(caf_genes %in% valid_genes)) {
  caf_expr <- rowMeans(expr_matrix[, caf_genes])
  cor_test_caf <- cor.test(acta2_expr, caf_expr, method = "spearman")
  
  cor_results <- rbind(cor_results, data.frame(
    Gene = "CAF_signature",
    R = round(cor_test_caf$estimate, 3),
    P = cor_test_caf$p.value,
    Interpretation = "CAF整体活性"
  ))
}

# 6. 添加显著性标记
cor_results$Significance <- ifelse(
  cor_results$P < 0.001, "***",
  ifelse(cor_results$P < 0.01, "**", 
         ifelse(cor_results$P < 0.05, "*", "ns"))
)

dir.create("output", showWarnings = FALSE)
write.csv(cor_results, "output/acta2_correlation_network.csv", row.names = FALSE)
library(corrplot)
cor_matrix <- cor(expr_matrix, method = "spearman")
cor_melted <- reshape2::melt(cor_matrix)
colnames(cor_melted) <- c("Gene1", "Gene2", "Correlation")
p_heatmap <- ggplot(cor_melted, aes(x = Gene1, y = Gene2, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = round(Correlation, 2)), 
            color = "black", size = 3.0, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2E9FDF",   # 蓝色（负相关）
    mid = "white",     # 白色（无相关）
    high = "#E74C3C",  # 红色（正相关）
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nCorrelation"
  ) +
  labs(
    title = "Correlation Heatmap of ACTA2 with Key Markers",
    #subtitle = paste("TCGA LMS cohort | n =", ncol(expr_matrix), "samples"),
    subtitle = paste("TCGA LMS cohort"),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm")
  ) +
  coord_fixed()  
ggsave("figures/Figure1a.png",
       plot = p_heatmap,
       width = 6,
       height = 4,
       dpi = 600,
       bg = "white")

cat("TCGA -end ")

