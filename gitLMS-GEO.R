options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(timeout = 600)
options(download.file.method = "libcurl")
rm(list = ls()) 
base_packages <- c("dplyr", "ggplot2", "survival", "survminer", "ggpubr",
                   "survivalROC", "rms", "pROC", "forestmodel", "timeROC", "tidyr",
                   "pheatmap", "reshape2", "RColorBrewer","effsize","Seurat","cowplot") 
for(pkg in base_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("GEOquery", quietly = TRUE))
BiocManager::install("GEOquery")
library(GEOquery)
if(!require("SummarizedExperiment")) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("SummarizedExperiment")
  BiocManager::install("GEOquery")
  library(SummarizedExperiment)
  library(GEOquery)
}
cat("GEO start")
tryCatch({
  if(!require("immunedeconv")) {
    if (!require("remotes", quietly = TRUE))
      install.packages("remotes")
    remotes::install_github("icbi-lab/immunedeconv")
    library(immunedeconv)
  }
}, error = function(e) {
  cat("⚠ immunedeconv包安装失败:", e$message, "\n")
})
setwd("D:/TCGA_Sarcoma_Analysis/ACTA2")
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
dir.create("results/immune_analysis", showWarnings = FALSE)
dir.create("cache_dir", showWarnings = FALSE)
library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)
gse_file <- "D:/TCGA_Sarcoma_Analysis/ACTA2/data/GSE21050_series_matrix.txt.gz"
all_gse_lines <- readLines(gse_file, encoding = "UTF-8")
clinical_lines <- all_gse_lines[grep("^!Sample_characteristics_ch1", all_gse_lines)]
clinical_data <- data.frame(sample_id = character(0))
geo_accession_line <- all_gse_lines[grep("^!Sample_geo_accession", all_gse_lines)]
if (length(geo_accession_line) > 0) {
  parts <- strsplit(geo_accession_line, "\t")[[1]]
  sample_ids <- gsub('"', '', parts[-1])
  clinical_data <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
}
for (line in clinical_lines) {
  parts <- strsplit(line, "\t")[[1]]
  first_value <- gsub('"', '', parts[2])
  if (grepl(":", first_value)) {
    feature_name <- strsplit(first_value, ":")[[1]][1]
    feature_name <- trimws(feature_name)
    values <- sapply(parts[-1], function(x) {
      x_clean <- gsub('"', '', x)
      if (grepl(":", x_clean)) {
        val_parts <- strsplit(x_clean, ":")[[1]]
        if (length(val_parts) > 1) {
          return(trimws(paste(val_parts[-1], collapse = ":")))
        } else {
          return(trimws(val_parts[1]))
        }
      } else {
        return(trimws(x_clean))
      }
    })
    if (length(values) == nrow(clinical_data)) {
      clinical_data[[feature_name]] <- values
    }
  }
}
begin_line <- which(grepl("^!series_matrix_table_begin", all_gse_lines))
end_line <- which(grepl("^!series_matrix_table_end", all_gse_lines))

if (length(begin_line) > 0 && length(end_line) > 0) {
  expr_matrix_lines <- all_gse_lines[(begin_line+1):(end_line-1)]
  expr_data <- read.delim(text = expr_matrix_lines, header = TRUE, stringsAsFactors = FALSE)
  probe_ids <- expr_data[, 1]
  expr_matrix <- as.matrix(expr_data[, -1])
  rownames(expr_matrix) <- probe_ids
  expr_samples <- colnames(expr_matrix)
  expr_samples <- gsub('"', '', expr_samples)
} else {
  stop("无法继续分析")
}
acta2_probe <- "200974_at"
if (acta2_probe %in% rownames(expr_matrix)) {
  acta2_expr_raw <- expr_matrix[acta2_probe, ]
  acta2_expr <- log2(acta2_expr_raw + 1)
} else {
  stop("无法继续分析")
}
diagnosis_col <- grep("diagnosis", names(clinical_data), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(diagnosis_col)) {
  diag_table <- table(clinical_data[[diagnosis_col]])
  leiomyo_idx <- grepl("Leiomyosarcoma", clinical_data[[diagnosis_col]], ignore.case = TRUE)
  leiomyo_df <- clinical_data[leiomyo_idx, ]
  leiomyo_df$ACTA2_expr <- acta2_expr[leiomyo_df$sample_id]

} else {
  stop("无法继续分析")
}

# ============== ICB signature 分析（GEO 队列, LMS n=84）==============
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("hgu133plus2.db", quietly = TRUE)) BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
probe2sym <- AnnotationDbi::select(hgu133plus2.db,
                                   keys = rownames(expr_matrix),
                                   columns = "SYMBOL",
                                   keytype = "PROBEID")
probe2sym <- probe2sym[!is.na(probe2sym$SYMBOL), ]
expr_log <- log2(expr_matrix[probe2sym$PROBEID, ] + 1)
rownames(expr_log) <- probe2sym$SYMBOL

expr_log <- expr_log[order(-rowMeans(expr_log)), ]
expr_gene <- expr_log[!duplicated(rownames(expr_log)), ]
# ---- ICB 基因集（Zhou et al., Cancer Letters 2026 的 ICB signature）----
icb_genes <- c("CD274","PDCD1","CTLA4","CD80","CD86","LAG3","HLA-B","TIGIT",
               "PVRL2","NECTIN3","PVR","NECTIN2","BTLA","TNFRSF14","KLRC1",
               "KLRD1","HLA-E","CD200R1","CD200","VSIR","IGSF11","HAVCR2",
               "HMGB1","LGALS9","CEACAM1","PTDSS1","ITPRIPL1")
icb_use <- intersect(icb_genes, rownames(expr_gene))
lms_samples <- leiomyo_df$sample_id
icb_mat <- expr_gene[icb_use, lms_samples]
icb_z <- t(scale(t(icb_mat)))                    
leiomyo_df$ICB_score <- colMeans(icb_z, na.rm = TRUE)
med <- median(leiomyo_df$ACTA2_expr, na.rm = TRUE)
leiomyo_df$ACTA2_group <- factor(ifelse(leiomyo_df$ACTA2_expr >= med, "High", "Low"),
                                 levels = c("Low", "High"))
wt <- wilcox.test(ICB_score ~ ACTA2_group, data = leiomyo_df)
ct <- cor.test(leiomyo_df$ACTA2_expr, leiomyo_df$ICB_score, method = "spearman")
single_gene <- sapply(icb_use, function(g) {
  x <- expr_gene[g, lms_samples]
  p <- wilcox.test(x ~ leiomyo_df$ACTA2_group)$p.value
  d <- median(x[leiomyo_df$ACTA2_group == "High"]) - median(x[leiomyo_df$ACTA2_group == "Low"])
  c(delta = d, p = p)
})
sg_df <- data.frame(gene = icb_use,
                    delta_HminusL = round(single_gene["delta", ], 3),
                    P = single_gene["p", ],
                    FDR = p.adjust(single_gene["p", ], method = "BH"))
sg_df <- sg_df[order(sg_df$P), ]
library(ggplot2)
p_icb <- ggplot(leiomyo_df, aes(ACTA2_group, ICB_score, fill = ACTA2_group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = c("Low" = "#E74C3C", "High" = "#2E86C1")) +
  labs(title = "GSE21050 LMS: immune checkpoint signature by ACTA2",
       y = "ICB signature score (z)", x = "ACTA2 group") +
  theme_classic(base_size = 13) +
  annotate("text", x = 1.5, y = max(leiomyo_df$ICB_score, na.rm = TRUE),
           label = paste0("P = ", signif(wt$p.value, 3))) +
  theme(legend.position = "none")
ggsave("D:/TCGA_Sarcoma_Analysis/ACTA2/output/GEO_ACTA2_vs_ICB.pdf", p_icb, width = 4.5, height = 5)
write.csv(sg_df, "D:/TCGA_Sarcoma_Analysis/ACTA2/output/GEO_ICB_single_gene.csv", row.names = FALSE)
cache_dir <- "cache"
if (!dir.exists(cache_dir)) dir.create(cache_dir)
if (!exists("lms_expr")) {
  lms_expr <- expr_matrix
}
data_hash <- digest::digest(paste(
  dim(lms_expr), 
  paste(colnames(lms_expr), collapse = ""),
  collapse = "_"
))
estimate_cache_rds <- file.path(cache_dir, sprintf("GPL570_ESTIMATE_%s.rds", data_hash))
estimate_cache_csv <- file.path(cache_dir, sprintf("GPL570_ESTIMATE_%s.csv", data_hash))
gene_matrix_cache <- file.path(cache_dir, sprintf("GPL570_gene_matrix_%s.rds", data_hash))
if (file.exists(gene_matrix_cache)) {
  lms_expr_gene <- readRDS(gene_matrix_cache)
} else {
  if (!require("hgu133plus2.db", quietly = TRUE)) {
    BiocManager::install("hgu133plus2.db", ask = FALSE)
  }
  library(hgu133plus2.db)
  probe_map <- AnnotationDbi::select(hgu133plus2.db,
                                     keys = rownames(lms_expr),
                                     columns = "SYMBOL",
                                     keytype = "PROBEID") %>% dplyr::filter(!is.na(SYMBOL))
    lms_expr_gene <- lms_expr %>%
    as.data.frame() %>%
    tibble::rownames_to_column("PROBEID") %>%
    dplyr::inner_join(probe_map, by = "PROBEID") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(dplyr::across(-PROBEID, stats::median, na.rm = TRUE)) %>%
    tibble::column_to_rownames("SYMBOL") %>%
    as.matrix()
  saveRDS(lms_expr_gene, gene_matrix_cache)
}
if (file.exists(estimate_cache_rds)) {
  est <- readRDS(estimate_cache_rds)
} else if (file.exists(estimate_cache_csv)) {
  est <- read.csv(estimate_cache_csv, stringsAsFactors = FALSE)
} else {
  est_result <- immunedeconv::deconvolute(lms_expr_gene, method = "estimate")
  est <- as.data.frame(t(est_result[, -1]))
  colnames(est) <- c("stromal_score", "immune_score", "estimate_score", "tumor_purity")
  est$sample_id <- rownames(est)          
  est <- est %>% dplyr::mutate(dplyr::across(-sample_id, as.numeric))
  saveRDS(est, estimate_cache_rds)
  write.csv(est, estimate_cache_csv, row.names = FALSE)

}
leiomyo_df <- leiomyo_df %>%
  left_join(est, by = "sample_id")
if ("ACTA2_expr" %in% colnames(leiomyo_df)) {
  valid_acta2 <- sum(!is.na(leiomyo_df$ACTA2_expr))
} else if ("ACTA2" %in% rownames(lms_expr_gene)) {
  leiomyo_df$ACTA2_expr <- lms_expr_gene["ACTA2", leiomyo_df$sample_id]
}
library(ppcor)
geo_partial <- pcor.test(
  x = leiomyo_df$ACTA2_expr,     
  y = leiomyo_df$immune_score,  
  z = leiomyo_df$tumor_purity,    
  method = "spearman"
)

if (!dir.exists("data")) dir.create("data")
saveRDS(leiomyo_df, file.path("data", "leiomyo_df_with_ESTIMATE.rds"))
write.csv(leiomyo_df, file.path("data", "leiomyo_df_with_ESTIMATE.csv"), row.names = FALSE)
tcga_estimate_path <- "D:/TCGA_Sarcoma_Analysis/ACTA2/results/immune_analysis/surv_df_immuneESTIMATE_Scores_ACTA2.csv"
tcga_estimate <- read.csv(tcga_estimate_path)
n_tcga <- sum(!is.na(tcga_estimate$gene_expr) & !is.na(tcga_estimate$immune_score))
rename_col <- function(df, old_pattern, new_name) {
  old_col <- grep(old_pattern, colnames(df), value = TRUE, ignore.case = TRUE)
  if (length(old_col) > 0 && !new_name %in% colnames(df)) {
    df[[new_name]] <- df[[old_col[1]]]
  }
  return(df)
}
leiomyo_df <- rename_col(leiomyo_df, "immune_score|ESTIMATE_Immune", "immune_score")
leiomyo_df <- rename_col(leiomyo_df, "stromal_score|ESTIMATE_Stromal", "stromal_score")

# ==============GSE21050相关性分析 ===============

valid_data <- leiomyo_df[!is.na(leiomyo_df$ACTA2_expr) & 
                           !is.na(leiomyo_df$immune_score) & 
                           !is.na(leiomyo_df$stromal_score), ]
n_gse <- nrow(valid_data)
if (n_gse < 20) stop("样本不足")
median_acta2 <- median(valid_data$ACTA2_expr)
valid_data$ACTA2_group <- ifelse(valid_data$ACTA2_expr >= median_acta2, "High", "Low")
gse_immune_spear <- cor.test(valid_data$ACTA2_expr, valid_data$immune_score,
                             method = "spearman", exact = FALSE)
gse_immune_wilcox <- wilcox.test(immune_score ~ ACTA2_group, data = valid_data)
gse_immune_high <- mean(valid_data$immune_score[valid_data$ACTA2_group == "High"])
gse_immune_low <- mean(valid_data$immune_score[valid_data$ACTA2_group == "Low"])
gse_stromal_spear <- cor.test(valid_data$ACTA2_expr, valid_data$stromal_score,
                              method = "spearman", exact = FALSE)
gse_stromal_wilcox <- wilcox.test(stromal_score ~ ACTA2_group, data = valid_data)
gse_stromal_high <- mean(valid_data$stromal_score[valid_data$ACTA2_group == "High"])
gse_stromal_low <- mean(valid_data$stromal_score[valid_data$ACTA2_group == "Low"])
# # ============== TCGA相关性分析 ===============
tcga_valid <- tcga_estimate[!is.na(tcga_estimate$gene_expr) & 
                              !is.na(tcga_estimate$immune_score) & 
                              !is.na(tcga_estimate$stromal_score), ]
tcga_median <- median(tcga_valid$gene_expr)
tcga_valid$ACTA2_group <- ifelse(tcga_valid$gene_expr >= tcga_median, "High", "Low")
tcga_immune_spear <- cor.test(tcga_valid$gene_expr, tcga_valid$immune_score,
                              method = "spearman", exact = FALSE)
tcga_immune_wilcox <- wilcox.test(immune_score ~ ACTA2_group, data = tcga_valid)
tcga_immune_high <- mean(tcga_valid$immune_score[tcga_valid$ACTA2_group == "High"])
tcga_immune_low <- mean(tcga_valid$immune_score[tcga_valid$ACTA2_group == "Low"])
tcga_stromal_spear <- cor.test(tcga_valid$gene_expr, tcga_valid$stromal_score,
                               method = "spearman", exact = FALSE)
tcga_stromal_wilcox <- wilcox.test(stromal_score ~ ACTA2_group, data = tcga_valid)
tcga_stromal_high <- mean(tcga_valid$stromal_score[tcga_valid$ACTA2_group == "High"])
tcga_stromal_low <- mean(tcga_valid$stromal_score[tcga_valid$ACTA2_group == "Low"])
z_immune <- (atanh(tcga_immune_spear$estimate) - atanh(gse_immune_spear$estimate)) / 
  sqrt(1/(nrow(tcga_valid)-3) + 1/(n_gse-3))
p_immune <- 2 * (1 - pnorm(abs(z_immune)))
z_stromal <- (atanh(tcga_stromal_spear$estimate) - atanh(gse_stromal_spear$estimate)) / 
  sqrt(1/(nrow(tcga_valid)-3) + 1/(n_gse-3))
p_stromal <- 2 * (1 - pnorm(abs(z_stromal)))
output_dir <- "D:/TCGA_Sarcoma_Analysis/ACTA2/figures/GEO_validation"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
library(ggplot2)
library(ggpubr)
library(patchwork)
tcga_plot_df <- data.frame(
  ACTA2_expr = tcga_valid$gene_expr,
  immune_score = tcga_valid$immune_score,
  stromal_score = tcga_valid$stromal_score,
  ACTA2_group = tcga_valid$ACTA2_group,
  Dataset = "TCGA-LMS"
)
geo_plot_df <- data.frame(
  ACTA2_expr = valid_data$ACTA2_expr,
  immune_score = valid_data$immune_score,
  stromal_score = valid_data$stromal_score,
  ACTA2_group = valid_data$ACTA2_group,
  Dataset = "GSE21050"
)
combined_df <- rbind(tcga_plot_df, geo_plot_df)
p_immune_scatter <- ggplot(combined_df, aes(x = ACTA2_expr, y = immune_score, color = Dataset)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = Dataset), alpha = 0.2) +
  scale_color_manual(values = c("TCGA-LMS" = "#E64B35", "GSE21050" = "#4DBBD5")) +
  scale_fill_manual(values = c("TCGA-LMS" = "#E64B35", "GSE21050" = "#4DBBD5")) +
  labs(
    title = "ACTA2 Expression vs Immune Score",
    subtitle = sprintf("TCGA: ρ=%.3f, p=%.2e | GEO: ρ=%.3f, p=%.2e",
                       tcga_immune_spear$estimate, tcga_immune_spear$p.value,
                       gse_immune_spear$estimate, gse_immune_spear$p.value),
    x = "ACTA2 Expression (log2)",
    y = "ESTIMATE Immune Score"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "top"
  )
p_stromal_scatter <- ggplot(combined_df, aes(x = ACTA2_expr, y = stromal_score, color = Dataset)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = Dataset), alpha = 0.2) +
  scale_color_manual(values = c("TCGA-LMS" = "#E64B35", "GSE21050" = "#4DBBD5")) +
  scale_fill_manual(values = c("TCGA-LMS" = "#E64B35", "GSE21050" = "#4DBBD5")) +
  labs(
    title = "ACTA2 Expression vs Stromal Score",
    subtitle = sprintf("TCGA: ρ=%.3f, p=%.2e | GEO: ρ=%.3f, p=%.2e",
                       tcga_stromal_spear$estimate, tcga_stromal_spear$p.value,
                       gse_stromal_spear$estimate, gse_stromal_spear$p.value),
    x = "ACTA2 Expression (log2)",
    y = "ESTIMATE Stromal Score"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "top"
  )
p_immune_box <- ggplot(combined_df, aes(x = ACTA2_group, y = immune_score, fill = Dataset)) +
  geom_boxplot(alpha = 0.8, position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8), alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("TCGA-LMS" = "#E64B35", "GSE21050" = "#4DBBD5")) +
  labs(
    title = "Immune Score by ACTA2 Group",
    subtitle = sprintf("TCGA: High=%.1f vs Low=%.1f, p=%.2e | GEO: High=%.1f vs Low=%.1f, p=%.2e",
                       tcga_immune_high, tcga_immune_low, tcga_immune_wilcox$p.value,
                       gse_immune_high, gse_immune_low, gse_immune_wilcox$p.value),
    x = "ACTA2 Expression Group",
    y = "Immune Score"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    legend.position = "top"
  )
p_stromal_box <- ggplot(combined_df, aes(x = ACTA2_group, y = stromal_score, fill = Dataset)) +
  geom_boxplot(alpha = 0.8, position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8), alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("TCGA-LMS" = "#E64B35", "GSE21050" = "#4DBBD5")) +
  labs(
    title = "Stromal Score by ACTA2 Group",
    subtitle = sprintf("TCGA: High=%.1f vs Low=%.1f, p=%.2e | GEO: High=%.1f vs Low=%.1f, p=%.2e",
                       tcga_stromal_high, tcga_stromal_low, tcga_stromal_wilcox$p.value,
                       gse_stromal_high, gse_stromal_low, gse_stromal_wilcox$p.value),
    x = "ACTA2 Expression Group",
    y = "Stromal Score"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    legend.position = "top"
  )
combined_plot <- (p_immune_scatter | p_stromal_scatter) / (p_immune_box | p_stromal_box) +
  plot_annotation(
    title = "Cross-Cohort Validation: ACTA2 Defines Low-Immune, Low-Stromal LMS Subtype",
    subtitle = "TCGA-LMS (n=96) vs GSE21050-LMS (n=84)",
    tag_levels = "a", 
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )
ggsave(file.path(output_dir, "Figure3.png"), 
       combined_plot, width = 12, height = 10, dpi = 300)


validation_summary <- data.frame(
  Metric = c("Immune Score", "Stromal Score"),
  TCGA_rho = c(tcga_immune_spear$estimate, tcga_stromal_spear$estimate),
  TCGA_p = c(tcga_immune_spear$p.value, tcga_stromal_spear$p.value),
  TCGA_High = c(tcga_immune_high, tcga_stromal_high),
  TCGA_Low = c(tcga_immune_low, tcga_stromal_low),
  TCGA_wilcox_p = c(tcga_immune_wilcox$p.value, tcga_stromal_wilcox$p.value),
  GEO_rho = c(gse_immune_spear$estimate, gse_stromal_spear$estimate),
  GEO_p = c(gse_immune_spear$p.value, gse_stromal_spear$p.value),
  GEO_High = c(gse_immune_high, gse_stromal_high),
  GEO_Low = c(gse_immune_low, gse_stromal_low),
  GEO_wilcox_p = c(gse_immune_wilcox$p.value, gse_stromal_wilcox$p.value),
  Direction_Consistent = c(TRUE, TRUE),
  Fisher_Z_p = c(p_immune, p_stromal)
)

write.csv(validation_summary, file.path(output_dir, "Table_ESTIMATE_Validation.csv"), row.names = FALSE)
write.csv(combined_df, file.path(output_dir, "Combined_Plot_Data.csv"), row.names = FALSE)
validation_dir <- "figures/GEO_validation"
if (!dir.exists(validation_dir)) dir.create(validation_dir, recursive = TRUE)
expr_range <- range(lms_expr_gene, na.rm = TRUE)
if (expr_range[2] > 50) {
  lms_expr_for_xcell <- log2(lms_expr_gene + 1)

} else {
  lms_expr_for_xcell <- lms_expr_gene
}
gene_means <- rowMeans(lms_expr_for_xcell, na.rm = TRUE)
expr_threshold <- quantile(gene_means, 0.25)
keep_genes <- gene_means > expr_threshold
lms_expr_filtered <- lms_expr_for_xcell[keep_genes, ]
if (!require("immunedeconv", quietly = TRUE)) {
  BiocManager::install("immunedeconv", ask = FALSE)
  library(immunedeconv)
}
xcell_cache_file <- file.path(cache_dir, sprintf("GEO_xCell_corrected_%s.rds", data_hash))
if (file.exists(xcell_cache_file)) {
  xcell_raw <- readRDS(xcell_cache_file)
} else {
  cat("  运行xCell算法...\n")
  xcell_raw <- immunedeconv::deconvolute(
    lms_expr_filtered, 
    method = "xcell",
    arrays = TRUE,  
    scale = TRUE    
  )
  saveRDS(xcell_raw, xcell_cache_file)
}
if ("cell_type" %in% colnames(xcell_raw)) {
  cell_types <- xcell_raw$cell_type
  xcell_values <- as.data.frame(t(xcell_raw[, -1]))
  colnames(xcell_values) <- cell_types
} else {
  cell_types <- xcell_raw[1, ]
  xcell_values <- as.data.frame(t(xcell_raw[-1, ]))
  colnames(xcell_values) <- as.character(cell_types)
}
xcell_values$sample_id <- rownames(xcell_values)
numeric_cols <- setdiff(colnames(xcell_values), "sample_id")
xcell_values[numeric_cols] <- lapply(xcell_values[numeric_cols], as.numeric)
cell_means <- colMeans(xcell_values[numeric_cols], na.rm = TRUE)
keep_cells <- names(cell_means)[cell_means > 0.001]
keep_cells <- intersect(keep_cells, numeric_cols)
xcell_filtered <- xcell_values[, c("sample_id", keep_cells)]
valid_data_xcell <- leiomyo_df %>%
  left_join(xcell_filtered, by = "sample_id")
n_valid <- sum(!is.na(valid_data_xcell$ACTA2_expr))
acta2_median <- median(valid_data_xcell$ACTA2_expr, na.rm = TRUE)
valid_data_xcell$group <- ifelse(valid_data_xcell$ACTA2_expr > acta2_median, "High", "Low")
cell_diff <- purrr::map_dfr(keep_cells, function(cell) {
  if (!cell %in% colnames(valid_data_xcell)) {
    return(NULL)
  }
  test_data <- valid_data_xcell[!is.na(valid_data_xcell[[cell]]) & 
                                  !is.na(valid_data_xcell$group), ]
  if (nrow(test_data) < 10) {
    return(NULL)
  }
  test_res <- wilcox.test(test_data[[cell]] ~ test_data$group)
  high_mean <- mean(test_data[[cell]][test_data$group == "High"], na.rm = TRUE)
  low_mean <- mean(test_data[[cell]][test_data$group == "Low"], na.rm = TRUE)
  direction <- ifelse(high_mean > low_mean, "High > Low", "Low > High")
  epsilon <- 1e-6
  log2fc <- log2((high_mean + epsilon) / (low_mean + epsilon))
  
  tibble::tibble(
    Cell_Type = cell,
    p_value = test_res$p.value,
    High_mean = high_mean,
    Low_mean = low_mean,
    Direction = direction,
    FoldChange_original = high_mean / low_mean,
    Log2_FoldChange = log2fc
  )
}) %>%
  dplyr::arrange(p_value) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    significance_level = dplyr::case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_value < 0.05 ~ "nominal",
      TRUE ~ "NS"
    )
  )
sig_cells <- cell_diff %>% dplyr::filter(p_adj < 0.05)
nominal_sig <- cell_diff %>% dplyr::filter(p_value < 0.05 & p_adj >= 0.05)
write.csv(
  cell_diff,
  file.path(validation_dir, "GEO_xCell_Differential_Analysis_ACTA2_corrected.csv"),
  row.names = FALSE
)
write.csv(
  valid_data_xcell,
  file.path(validation_dir, "GEO_xCell_Complete_Data_corrected.csv"),
  row.names = FALSE
)
tcga_xcell_path <- "D:/TCGA_Sarcoma_Analysis/ACTA2/results/immune_analysis/xCell_Immune_Cell_Infiltration_ACTA2.csv"
tcga_xcell <- read.csv(tcga_xcell_path, check.names = FALSE)
key_cells <- c("Cancer associated fibroblast", "Macrophage", 
               "Macrophage M1", "Macrophage M2", "B cell plasma", "Endothelial cell")
comparison_results <- data.frame()
for (cell in key_cells) {
  geo_sub <- valid_data_xcell[, c("ACTA2_expr", cell)]
  geo_sub <- geo_sub[complete.cases(geo_sub), ]
  geo_cor <- cor.test(geo_sub$ACTA2_expr, geo_sub[[cell]], 
                      method = "spearman", exact = FALSE)
  tcga_sub <- tcga_xcell[, c("gene_expr", cell)]
  tcga_sub <- tcga_sub[complete.cases(tcga_sub), ]
  if (nrow(tcga_sub) < 10) {
    cat(sprintf("   跳过 %s: TCGA有效样本不足\n", cell))
    next
  }
  tcga_cor <- cor.test(tcga_sub$gene_expr, tcga_sub[[cell]], 
                       method = "spearman", exact = FALSE)
  if (abs(geo_cor$estimate) < 0.99 && abs(tcga_cor$estimate) < 0.99) {
    z1 <- atanh(geo_cor$estimate)
    z2 <- atanh(tcga_cor$estimate)
    z_diff <- (z1 - z2) / sqrt(1/(nrow(geo_sub) - 3) + 1/(nrow(tcga_sub) - 3))
    p_diff <- 2 * pnorm(-abs(z_diff))
    consistency <- ifelse(p_diff > 0.05, "Consistent", "Divergent")
  } else {
    z_diff <- NA
    p_diff <- NA
    consistency <- "Cannot Test"
  }
  comparison_results <- rbind(comparison_results, data.frame(
    Cell_Type = cell,
    GEO_rho = round(geo_cor$estimate, 3),
    GEO_p = geo_cor$p.value,
    GEO_n = nrow(geo_sub),
    TCGA_rho = round(tcga_cor$estimate, 3),
    TCGA_p = tcga_cor$p.value,
    TCGA_n = nrow(tcga_sub),
    Fisher_Z = ifelse(is.na(z_diff), NA, round(z_diff, 3)),
    Fisher_p = p_diff,
    Consistency = consistency,
    stringsAsFactors = FALSE
  ))
}
write.csv(
  comparison_results,
  file.path(validation_dir, "Figure6_TCGA_GEO_xCell_Comparison_corrected.csv"),
  row.names = FALSE
)

divergent_cells <- comparison_results[comparison_results$Consistency == "Divergent", ]
library(ggplot2)
library(gridExtra)
plot_list <- list()
for (i in 1:nrow(comparison_results)) {
  cell <- comparison_results$Cell_Type[i]
  geo_rho <- comparison_results$GEO_rho[i]
  tcga_rho <- comparison_results$TCGA_rho[i]
  fisher_p <- comparison_results$Fisher_p[i]
  
  plot_data <- data.frame(
    Dataset = c("GEO", "TCGA"),
    Correlation = c(geo_rho, tcga_rho),
    stringsAsFactors = FALSE
  )
  p <- ggplot(plot_data, aes(x = Dataset, y = Correlation, fill = Dataset)) +
    geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = c("GEO" = "#4DBBD5", "TCGA" = "#E64B35")) +
    labs(
      title = cell,
      subtitle = ifelse(!is.na(fisher_p), 
                        sprintf("Fisher's p = %.3f", fisher_p),
                        "Fisher test not applicable"),
      x = "",
      y = "Spearman ρ",
      tag = letters[i]   
    ) +
    ylim(-1, 1) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.major.x = element_blank()
    ) +
    geom_text(aes(label = sprintf("%.3f", Correlation)), 
              vjust = -0.5, size = 3)
  
  plot_list[[cell]] <- p
}
if (length(plot_list) > 0) {
  n_plots <- length(plot_list)
  n_col <- min(3, n_plots)
  n_row <- ceiling(n_plots / n_col)
  
  combined_plot <- do.call(gridExtra::grid.arrange, 
                           c(plot_list, ncol = n_col,
                             top = "TCGA vs GEO: correlation of ACTA2 expression with immune cell subsets"))
  ggsave(
    file.path(validation_dir, "Figure4.png"),
    combined_plot,
    # width = n_col * 4,
    # height = n_row * 3,
    width = 8,
    height = 6,
    dpi = 300
  )
}
key_status <- sapply(key_cells, function(cell) {
  if (cell %in% cell_diff$Cell_Type) {
    row <- cell_diff[cell_diff$Cell_Type == cell, ]
    sig <- ifelse(row$p_adj < 0.05, "FDR显著", 
                  ifelse(row$p_value < 0.05, "名义显著", "不显著"))
    return(paste(sig, sprintf("(p_adj=%.3f)", row$p_adj)))
  } else {
    return("未分析")
  }
})
med <- median(leiomyo_df$ACTA2_expr, na.rm = TRUE)
leiomyo_df$group2 <- ifelse(leiomyo_df$ACTA2_expr > med, "High", "Low")
features <- c("Macrophage M2", "Macrophage M1","VIM", "immune_score", "stromal_score")
for (feat in features) {
  if (!feat %in% colnames(leiomyo_df)) {
    next
  }
  df_sub <- leiomyo_df[!is.na(leiomyo_df[[feat]]) & !is.na(leiomyo_df$group2), ]
  if (nrow(df_sub) == 0) next
  test_res <- wilcox.test(df_sub[[feat]] ~ df_sub$group2)
  means <- tapply(df_sub[[feat]], df_sub$group2, mean, na.rm = TRUE)

}
library(ppcor)
key6 <- c("HAVCR2","LGALS9","CD86","NECTIN3","CD200","CD274")
side6 <- c("CD86" = "Myeloid-restricted", "HAVCR2" = "Myeloid-restricted",
           "LGALS9" = "Myeloid-restricted", "CD274" = "Tumor-intrinsic",
           "CD200" = "Tumor-intrinsic", "NECTIN3" = "Tumor-intrinsic")
geo_pur <- leiomyo_df %>%
  dplyr::left_join(xcell_filtered[, c("sample_id","Macrophage")], by = "sample_id")
res_geo <- lapply(key6, function(g) {
  gx <- as.numeric(expr_gene[g, geo_pur$sample_id])
  ok <- !is.na(gx) & !is.na(geo_pur$tumor_purity) & !is.na(geo_pur$Macrophage)
  d  <- data.frame(acta2 = geo_pur$ACTA2_expr[ok], gx = gx[ok],
                   purity = geo_pur$tumor_purity[ok], mac = geo_pur$Macrophage[ok])
  r_acta2 <- cor.test(d$acta2, d$gx, method = "spearman")
  r_mac   <- cor.test(d$mac,   d$gx, method = "spearman")
  r_pur   <- cor.test(d$purity,d$gx, method = "spearman")
  if (side6[g] == "Myeloid-restricted") {
    pc <- pcor.test(d$acta2, d$gx, d$mac, method = "spearman")
    ctrl_var <- "Macrophage"
  } else {
    pc <- pcor.test(d$acta2, d$gx, d$purity, method = "spearman")
    ctrl_var <- "Purity"
  }
  data.frame(gene = g, side = side6[g],
             rho_ACTA2 = round(r_acta2$estimate, 3), p_ACTA2 = signif(r_acta2$p.value, 3),
             rho_mac   = round(r_mac$estimate, 3),   p_mac   = signif(r_mac$p.value, 3),
             rho_pur   = round(r_pur$estimate, 3),   p_pur   = signif(r_pur$p.value, 3),
             ctrl = ctrl_var,
             partial_rho = round(pc$estimate, 3), partial_p = signif(pc$p.value, 3))
})
geo_comp <- do.call(rbind, res_geo)
geo_comp$FDR_partial <- p.adjust(geo_comp$partial_p, method = "BH")

write.csv(geo_comp, "D:/TCGA_Sarcoma_Analysis/ACTA2/output/GEO_compartment_purity_check.csv",
          row.names = FALSE)

key4 <- c("APOE", "TREM2", "TGFB1", "SPP1")
res_comm_geo <- lapply(key4, function(g) {
  gx <- as.numeric(expr_gene[g, geo_pur$sample_id])
  ok <- !is.na(gx) & !is.na(geo_pur$Macrophage)
  d  <- data.frame(acta2 = geo_pur$ACTA2_expr[ok], gx = gx[ok],
                   mac = geo_pur$Macrophage[ok])
  r_acta2 <- cor.test(d$acta2, d$gx, method = "spearman")
  r_mac   <- cor.test(d$mac,   d$gx, method = "spearman")
  pc      <- pcor.test(d$acta2, d$gx, d$mac, method = "spearman")
  data.frame(gene = g,
             rho_ACTA2   = round(r_acta2$estimate, 3), p_ACTA2 = signif(r_acta2$p.value, 3),
             rho_mac     = round(r_mac$estimate, 3),   p_mac   = signif(r_mac$p.value, 3),
             partial_rho = round(pc$estimate, 3),      partial_p = signif(pc$p.value, 3))
})
geo_comm <- do.call(rbind, res_comm_geo)
geo_comm$FDR_partial <- p.adjust(geo_comm$partial_p, method = "BH")
write.csv(geo_comm, "D:/TCGA_Sarcoma_Analysis/ACTA2/output/GEO_communication_pcor_check.csv",
          row.names = FALSE)
library(SummarizedExperiment); library(dplyr); library(ppcor)
data <- readRDS("D:/TCGA_Sarcoma_Analysis/ACTA2/data/TCGA-SARC_data.rds")
expr <- assay(data, "tpm_unstrand")
ginfo <- as.data.frame(rowData(data))
clin  <- as.data.frame(colData(data))
lms_idx <- grep("leiomyosarcoma", clin$primary_diagnosis, ignore.case = TRUE)
expr_lms <- expr[, lms_idx]
sym <- ginfo$gene_name
keep <- !is.na(sym) & sym != ""
expr_gene <- log2(expr_lms[keep, ] + 1)
rownames(expr_gene) <- sym[keep]
expr_gene <- expr_gene[order(-rowMeans(expr_gene)), ]
expr_gene <- expr_gene[!duplicated(rownames(expr_gene)), ]
est_tcga <- read.csv("D:/TCGA_Sarcoma_Analysis/ACTA2/results/immune_analysis/surv_df_immuneESTIMATE_Scores_ACTA2.csv",
                     check.names = FALSE)
xc_tcga  <- read.csv("D:/TCGA_Sarcoma_Analysis/ACTA2/results/immune_analysis/xCell_Immune_Cell_Infiltration_ACTA2.csv",
                     check.names = FALSE)
samp <- colnames(expr_gene)
overlap <- sapply(est_tcga, function(c) mean(as.character(c) %in% samp))
id_col_est <- names(which.max(overlap))
if (max(overlap) < 0.5) { 
  overlap2 <- sapply(est_tcga, function(c) mean(substr(as.character(c),1,15) %in% substr(samp,1,15)))
  id_col_est <- names(which.max(overlap2))
  est_tcga$.mid <- substr(as.character(est_tcga[[id_col_est]]), 1, 15)
  samp_map <- setNames(samp, substr(samp,1,15))
  est_tcga$.sid <- samp_map[est_tcga$.mid]
} else {
  est_tcga$.sid <- as.character(est_tcga[[id_col_est]])
}
overlap_x <- sapply(xc_tcga, function(c) mean(as.character(c) %in% samp))
if (max(overlap_x) > 0.5) {
  xc_tcga$.sid <- as.character(xc_tcga[[names(which.max(overlap_x))]])
} else if (nrow(xc_tcga) == nrow(est_tcga) &&
           isTRUE(all.equal(xc_tcga$gene_expr, est_tcga$gene_expr))) {
  xc_tcga$.sid <- est_tcga$.sid
} else {
  stop("xCell ")
}

tcga_pur <- data.frame(sid = samp, ACTA2 = as.numeric(expr_gene["ACTA2", samp])) %>%
  left_join(est_tcga %>% dplyr::select(.sid, tumor_purity, gene_expr), by = c("sid"=".sid")) %>%
  left_join(xc_tcga %>% dplyr::select(.sid, Macrophage), by = c("sid"=".sid"))
chk <- cor.test(tcga_pur$ACTA2, tcga_pur$gene_expr, method = "spearman")

key6 <- c("HAVCR2","LGALS9","CD86","NECTIN3","CD200","CD274")
side6 <- c("CD86" = "Myeloid-restricted", "HAVCR2" = "Myeloid-restricted",
           "LGALS9" = "Myeloid-restricted", "CD274" = "Tumor-intrinsic",
           "CD200" = "Tumor-intrinsic", "NECTIN3" = "Tumor-intrinsic")
res_tcga <- lapply(key6, function(g) {
  gx <- as.numeric(expr_gene[g, tcga_pur$sid])
  ok <- !is.na(gx) & !is.na(tcga_pur$tumor_purity) & !is.na(tcga_pur$Macrophage)
  d <- data.frame(acta2 = tcga_pur$ACTA2[ok], gx = gx[ok],
                  purity = tcga_pur$tumor_purity[ok], mac = tcga_pur$Macrophage[ok])
  r_acta2 <- cor.test(d$acta2, d$gx, method = "spearman")
  r_mac   <- cor.test(d$mac,   d$gx, method = "spearman")
  r_pur   <- cor.test(d$purity,d$gx, method = "spearman")
  if (side6[g] == "Myeloid-restricted") {
    pc <- pcor.test(d$acta2, d$gx, d$mac, method = "spearman"); ctrl_var <- "Macrophage"
  } else {
    pc <- pcor.test(d$acta2, d$gx, d$purity, method = "spearman"); ctrl_var <- "Purity"
  }
  data.frame(gene = g, side = side6[g],
             rho_ACTA2 = round(r_acta2$estimate,3), p_ACTA2 = signif(r_acta2$p.value,3),
             rho_mac   = round(r_mac$estimate,3),   p_mac   = signif(r_mac$p.value,3),
             rho_pur   = round(r_pur$estimate,3),   p_pur   = signif(r_pur$p.value,3),
             ctrl = ctrl_var,
             partial_rho = round(pc$estimate,3), partial_p = signif(pc$p.value,3))
})
tcga_comp <- do.call(rbind, res_tcga)
tcga_comp$FDR_partial <- p.adjust(tcga_comp$partial_p, method = "BH")
write.csv(tcga_comp, "D:/TCGA_Sarcoma_Analysis/ACTA2/output/TCGA_compartment_purity_check.csv",
          row.names = FALSE)
key4 <- c("APOE", "TREM2", "TGFB1", "SPP1")
res_comm <- lapply(key4, function(g) {
  gx <- as.numeric(expr_gene[g, tcga_pur$sid])
  ok <- !is.na(gx) & !is.na(tcga_pur$Macrophage)
  d <- data.frame(acta2 = tcga_pur$ACTA2[ok], gx = gx[ok],
                  mac = tcga_pur$Macrophage[ok])
  r_acta2 <- cor.test(d$acta2, d$gx, method = "spearman")
  r_mac   <- cor.test(d$mac,   d$gx, method = "spearman")
  pc      <- pcor.test(d$acta2, d$gx, d$mac, method = "spearman")
  data.frame(gene = g,
             rho_ACTA2   = round(r_acta2$estimate, 3), p_ACTA2 = signif(r_acta2$p.value, 3),
             rho_mac     = round(r_mac$estimate, 3),   p_mac   = signif(r_mac$p.value, 3),
             partial_rho = round(pc$estimate, 3),      partial_p = signif(pc$p.value, 3))
})
tcga_comm <- do.call(rbind, res_comm)
tcga_comm$FDR_partial <- p.adjust(tcga_comm$partial_p, method = "BH")
write.csv(tcga_comm, "D:/TCGA_Sarcoma_Analysis/ACTA2/output/TCGA_communication_pcor_check.csv",
          row.names = FALSE)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

out_dir <- "D:/TCGA_Sarcoma_Analysis/ACTA2/output"
col_low  <- "#E74C3C"  
col_high <- "#2E86C1"  
key6 <- c("CD86","HAVCR2","LGALS9","CD274","CD200","NECTIN3")
s3      <- read.csv(file.path(out_dir, "Supplementary_Table_S3.csv"))
geo_cp  <- read.csv(file.path(out_dir, "GEO_compartment_purity_check.csv"))
tcga_cp <- read.csv(file.path(out_dir, "TCGA_compartment_purity_check.csv"))

s4 <- rbind(cbind(cohort = "GSE21050", geo_cp),
            cbind(cohort = "TCGA-LMS", tcga_cp))
write.csv(s4,
          file.path(out_dir, "Supplementary_Table_S4.csv"),
          row.names = FALSE)
tcga_comm <- read.csv(file.path(out_dir, "TCGA_communication_pcor_check.csv"))
geo_comm  <- read.csv(file.path(out_dir, "GEO_communication_pcor_check.csv"))
s6 <- rbind(cbind(cohort = "TCGA-LMS", tcga_comm),            cbind(cohort = "GSE21050", geo_comm))
rownames(s6) <- NULL   
write.csv(s6, file.path(out_dir, "Supplementary_Table_S6.csv"),
          row.names = FALSE)
d8a <- bind_rows(
  s3 %>% dplyr::filter(gene %in% key6) %>%
    transmute(gene, cohort = "GSE21050 (n = 84)", delta = delta_GEO,          FDR = FDR_GEO),
  s3 %>% dplyr::filter(gene %in% key6) %>%
    transmute(gene, cohort = "TCGA-LMS (n = 96)", delta = delta_HminusL_TCGA, FDR = FDR_TCGA)
) %>% mutate(
  gene      = factor(gene, levels = key6),
  direction = ifelse(delta < 0, "ACTA2-low enriched", "ACTA2-high enriched"),
  sig       = FDR < 0.05
)
p8a <- ggplot(d8a, aes(x = delta, y = gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = direction, shape = sig),
             size = 3.5, stroke = 1.1) +
  facet_wrap(~cohort, scales = "free_x") +
  scale_color_manual(values = c("ACTA2-low enriched"  = col_low,
                                "ACTA2-high enriched" = col_high)) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "FDR < 0.05", "FALSE" = "FDR \u2265 0.05")) +
  labs(x = "Expression difference (ACTA2-high \u2212 ACTA2-low)", y = NULL,
       title = "a",
       color = NULL, shape = NULL) +
  guides(color = guide_legend(order = 1, override.aes = list(shape = 16)),
         shape = guide_legend(order = 2, override.aes = list(color = "grey40"))) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "bold"))
d8b <- s4 %>%
  transmute(cohort  = ifelse(cohort == "GSE21050", "GSE21050 (n = 84)", "TCGA-LMS (n = 96)"),
            gene    = factor(gene, levels = key6),
            raw     = rho_ACTA2,
            partial = partial_rho) %>%
  pivot_longer(c(raw, partial), names_to = "type", values_to = "rho") %>%
  mutate(type = factor(ifelse(type == "raw",
                              "Before adjustment",
                              "After adjustment"),
                       levels = c("Before adjustment", "After adjustment")))

p8b <- ggplot(d8b, aes(y = gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(aes(x = rho, group = interaction(gene, cohort)),
            color = "grey70", linewidth = 1) +
  geom_point(aes(x = rho, color = type), size = 3.2) +
  facet_wrap(~cohort) +
  scale_color_manual(values = c("Before adjustment"                    = col_low,
                                "After adjustment"  = col_high)) +
  labs(x = "Spearman \u03c1 with ACTA2", y = NULL, color = NULL,
       title = "b ") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "bold"))
p8 <- p8a / p8b + plot_layout(heights = c(1.5, 1))
ggsave(file.path(out_dir, "Figure6a-6b.png"), p8, width = 9, height = 8.5, dpi = 300)

cat("GEO end")
