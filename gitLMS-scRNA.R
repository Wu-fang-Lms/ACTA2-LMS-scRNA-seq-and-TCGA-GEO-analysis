cat("scRNA start")
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
  cat("⚠ immunedeconv:", e$message, "\n")
})
while (sink.number() > 0) {
  sink()
}
setwd("D:/TCGA_Sarcoma_Analysis/ACTA2")
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
dir.create("results/immune_analysis", showWarnings = FALSE)
  set.seed(42)
  options(future.globals.maxSize = 14 * 1024^3, timeout = 600)
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
  library(CellChat)
  if(!requireNamespace("infercnv", quietly = TRUE)) {
    BiocManager::install("infercnv", ask = FALSE, update = FALSE)
  }
  must_have <- c("Seurat","harmony","infercnv")
  check <- sapply(must_have, requireNamespace, quietly = TRUE)
  if(any(!check)) {
    stop("以下包仍未安装成功，请检查:\n",
         paste(names(check)[!check], collapse = ", "))
  }
  HAS_DF <- requireNamespace("DoubletFinder", quietly = TRUE)
  suppressPackageStartupMessages({
    library(Matrix); library(Seurat); library(ggplot2); library(dplyr)
    library(patchwork); library(ggrepel); library(data.table)
    library(harmony); library(infercnv)
  })
  dirs <- c("Figures/SingleCell","Figures/SingleCell/PDF",
            "Results/SingleCell","Results/SingleCell/inferCNV","Results/SingleCell/QC")
  for(d in dirs) dir.create(d, recursive=TRUE, showWarnings=FALSE)
  samples <- list(
    SRC147            = "data/GSM6534012_SRC147/GSM6534012_SRC147",
    SRC241            = "data/GSM6534015_SRC241/GSM6534015_SRC241",
    SRC242_biceps     = "data/GSM6534016_SRC242_biceps/GSM6534016_SRC242_biceps",
    SRC242_pectoralis = "data/GSM6534017_SRC242_pectoralis/GSM6534017_SRC242_pectoralis"
  )
  PARAM <- list(
    min_umi = 500, min_features = 300, max_mt = 15,
    hvg = 2000, harmony_dims = 20, cluster_res = 0.6,
    infercnv_cutoff = 0.1, cnv_ref_q = 0.95, mal_cluster_frac = 0.5,
    infercnv_mesenchymal_cap = 18000,
    infercnv_ref_cap = 1500
  )
  qc_stats <- data.frame(); qc_plots <- list(); knee_plots <- list()
  open_maybe_gz <- function(p) if(file.exists(paste0(p,".gz"))) gzfile(paste0(p,".gz")) else p
  analyze_scRNA <- function(sample_name, base_path) {
    f_mtx <- paste0(base_path,"_matrix.mtx"); f_feat <- paste0(base_path,"_features.tsv"); f_bc <- paste0(base_path,"_barcodes.tsv")
    if(!all(sapply(c(f_mtx,f_feat,f_bc), function(f) file.exists(f) || file.exists(paste0(f,".gz"))))) {
      cat("缺失:", sample_name, "\n"); return(NULL)
    }
    gene_symbols <- read.delim(open_maybe_gz(f_feat), header=FALSE)[,2]
    if(anyDuplicated(gene_symbols)) gene_symbols <- make.unique(gene_symbols, sep="_")
    counts <- readMM(open_maybe_gz(f_mtx))
    rownames(counts) <- gene_symbols
    colnames(counts) <- read.delim(open_maybe_gz(f_bc), header=FALSE)[,1]
    counts <- as(counts, "CsparseMatrix")
    if(min(counts@x) < 0) counts@x[counts@x < 0] <- 0
    n_raw <- ncol(counts)
    totals <- colSums(counts)
    ts <- sort(totals, decreasing=TRUE); ts <- ts[ts>0]
    knee_plots[[sample_name]] <<- ggplot(data.frame(rank=seq_along(ts), total=ts), aes(rank,total)) +
      geom_point(size=0.2) + scale_x_log10() + scale_y_log10() +
      geom_hline(yintercept=PARAM$min_umi, linetype=2, color="red") +
      theme_bw() + ggtitle(sample_name) + labs(x="Barcode rank", y="Total UMI")
    counts_real <- counts[, totals >= PARAM$min_umi]
    n_called <- ncol(counts_real)
    rm(counts, totals, ts); gc()
    seu <- CreateSeuratObject(counts_real, project=sample_name, min.cells=5)
    rm(counts_real); gc()
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")
    n0 <- ncol(seu)
    upper_q <- if(HAS_DF) 0.99 else 0.975
    q_up <- quantile(seu$nFeature_RNA, upper_q)
    seu <- subset(seu, subset = nFeature_RNA >= PARAM$min_features &
                    nFeature_RNA < q_up & percent.mt < PARAM$max_mt)
    if(ncol(seu) < 100) { cat("❌ 细胞不足100,跳过\n"); return(NULL) }
    n_dbl <- 0
    if(HAS_DF) {
      seu <- NormalizeData(seu, verbose=FALSE)
      seu <- FindVariableFeatures(seu, nfeatures=PARAM$hvg, verbose=FALSE)
      seu <- ScaleData(seu, verbose=FALSE)
      seu <- RunPCA(seu, npcs=20, verbose=FALSE)
      seu <- FindNeighbors(seu, dims=1:20, verbose=FALSE)
      seu <- FindClusters(seu, resolution=0.5, verbose=FALSE)
      nExp <- round(0.075 * ncol(seu))
      seu <- doubletFinder(seu, PCs=1:20, pN=0.25, pK=0.09, nExp=nExp)
      dbl_col <- grep("^DF.classifications", colnames(seu@meta.data), value=TRUE)
      seu$doublet <- seu@meta.data[[dbl_col[1]]]
      n_dbl <- sum(seu$doublet == "Doublet")
      seu <- subset(seu, subset = doublet == "Singlet")
    } else {
      cat(" HAS_DF\n")
    }
    seu <- NormalizeData(seu, verbose=FALSE)
    seu <- FindVariableFeatures(seu, nfeatures=PARAM$hvg, verbose=FALSE)
    
    qc_stats <<- rbind(qc_stats, data.frame(
      sample=sample_name, raw_barcodes=n_raw, called_cells=n_called,
      doublets=n_dbl, final_cells=ncol(seu), median_umi=median(seu$nCount_RNA),
      median_genes=median(seu$nFeature_RNA), mean_mt=round(mean(seu$percent.mt),2)))
    qc_plots[[sample_name]] <<- VlnPlot(seu, c("nFeature_RNA","nCount_RNA","percent.mt"),
                                        ncol=3, pt.size=0, log=TRUE) +
      plot_annotation(title=paste0(sample_name," (post-QC, n=",ncol(seu),")"))
    saveRDS(seu, paste0("Results/SingleCell/", sample_name, "_seurat_v3.rds"))
    return(seu)
  }
  
  run_with_cache <- function(sample_name, base_path) {
    rds_path <- paste0("Results/SingleCell/", sample_name, "_seurat_v3.rds")
    if(file.exists(rds_path)) { return(readRDS(rds_path)) }
  }
  seurat_objects <- list()
  for(s in names(samples)) {
    obj <- run_with_cache(s, samples[[s]])
    if(!is.null(obj)) seurat_objects[[s]] <- obj
    rm(obj); gc()
  }
  seurat_objects <- seurat_objects[!sapply(seurat_objects, is.null)]
  write.csv(qc_stats, "Results/SingleCell/QC/QC_summary.csv", row.names=FALSE)
  ok <- intersect(names(seurat_objects), names(qc_plots))
  if(length(ok) > 0) {
    p_s1 <- wrap_plots(c(qc_plots[ok], knee_plots[ok]), ncol=1) + plot_annotation(tag_levels="a")
    ggsave("Figures/SingleCell/FigureS1.png", p_s1, width=11, height=3*length(ok)*2,  dpi=600, bg="white", limitsize=FALSE)
  }
  # ====================  Harmony====================
  integrated_rds <- "Results/SingleCell/Integrated_v3_Harmony.rds"
  if(file.exists(integrated_rds)) {
   integrated <- readRDS(integrated_rds)
  } else {
    integrated <- merge(seurat_objects[[1]], y=seurat_objects[-1], add.cell.ids=names(seurat_objects))
    rm(seurat_objects); gc()
    integrated <- tryCatch(JoinLayers(integrated), error=function(e) integrated)
    integrated <- FindVariableFeatures(integrated, selection.method="vst", nfeatures=PARAM$hvg, verbose=FALSE)
    integrated <- ScaleData(integrated, features=VariableFeatures(integrated), verbose=FALSE)
    integrated <- RunPCA(integrated, features=VariableFeatures(integrated), npcs=30, verbose=FALSE)
    integrated <- RunHarmony(integrated, group.by.vars="orig.ident", dims.use=1:PARAM$harmony_dims, verbose=FALSE)
    integrated <- RunUMAP(integrated, reduction="harmony", dims=1:PARAM$harmony_dims,
                          reduction.name="umap.harmony", reduction.key="UMAPHARMONY_", verbose=FALSE)
    integrated <- FindNeighbors(integrated, reduction="harmony", dims=1:PARAM$harmony_dims, verbose=FALSE)
    integrated <- FindClusters(integrated, resolution=PARAM$cluster_res, verbose=FALSE)
    integrated$patient <- sapply(strsplit(as.character(integrated$orig.ident), "_"), `[`, 1)
    saveRDS(integrated, integrated_rds)
  }
  cell_markers <- list(
    SmoothMuscle = c("ACTA2","MYH11","TAGLN","CNN1"),
    Fibro_CAF    = c("COL1A1","DCN","LUM","FAP"),
    Endothelial  = c("PECAM1","CDH5","VWF","CLDN5"),
    Macrophages  = c("CD68","CD163","CSF1R","C1QA"),
    T_cells      = c("CD3D","CD3E","CD8A","CD4"),
    NK_cells     = c("NCAM1","NCR1","KLRD1"),
    B_cells      = c("CD79A","MS4A1","CD19"),
    Mast_cells   = c("TPSAB1","CPA3","KIT"),
    Pericytes    = c("RGS5","PDGFRB","CSPG4")
  )
  DefaultAssay(integrated) <- "RNA"
  if(!"data" %in% Layers(integrated, assay = "RNA"))
    integrated <- JoinLayers(integrated, assay = "RNA")
  stopifnot("layers " = "data" %in% Layers(integrated, assay = "RNA"))
  avg <- AverageExpression(integrated, assays = "RNA", layer = "data",
                           group.by = "seurat_clusters", verbose = FALSE)$RNA
  stopifnot("AverageExpression " = !is.null(avg) && nrow(avg) > 0)
  cl_mean <- sapply(cell_markers, function(g) {
    g <- intersect(g, rownames(avg))
    if(length(g) == 0) return(rep(0, ncol(avg)))
    colMeans(avg[g, , drop = FALSE])
  })
    cluster_anno <- colnames(cl_mean)[apply(cl_mean, 1, which.max)]
  names(cluster_anno) <- sub("^g", "", rownames(cl_mean))  
   integrated$celltype_coarse <- unname(cluster_anno[as.character(integrated$seurat_clusters)])
  cnv_rds <- "Results/SingleCell/cnv_score_v3.rds"
  if(file.exists(cnv_rds)) {
    cnv_score <- readRDS(cnv_rds)
  } else {
    go_path <- "Results/SingleCell/gene_order.txt"
    if(!file.exists(go_path)) {
      gtf_files <- list.files("data", pattern="\\.gtf(\\.gz)?$", full.names=TRUE)
      con <- if(grepl("\\.gz$", gtf_files[1])) gzfile(gtf_files[1]) else file(gtf_files[1])
      gtf_lines <- readLines(con); close(con)
      gtf_lines <- gtf_lines[!startsWith(gtf_lines, "#")]
      gtf_lines <- gtf_lines[grepl("\tgene\t", gtf_lines)]
      gtf <- fread(text=gtf_lines, sep="\t", header=FALSE,
                   col.names=c("chr","source","type","start","end","score","strand","frame","attr"))
      gtf <- gtf[chr %in% paste0("chr", c(1:22,"X","Y"))]
      gtf$gene <- sub('.*gene_name "([^"]+)".*', "\\1", gtf$attr)
      gtf <- gtf[!duplicated(gene)]                   
      go <- data.frame(chr=gtf$chr, start=gtf$start, end=gtf$end, row.names=gtf$gene)
      write.table(go, go_path, sep="\t", quote=FALSE, col.names=FALSE)
      rm(gtf_lines, gtf); gc()
    }
    fin_rd <- "Results/SingleCell/inferCNV/run.final.infercnv_obj"
    if(!file.exists(fin_rd)) {
      cat("\n==================== inferCNV ====================\n")
      ct <- integrated$celltype_coarse; names(ct) <- colnames(integrated)
      ref_types <- intersect(c("Macrophages","T_cells","NK_cells","B_cells","Endothelial","Mast_cells"),
                             unique(ct))
      mes_cells <- colnames(integrated)[ct %in% c("SmoothMuscle","Fibro_CAF")]
      if(length(mes_cells) > PARAM$infercnv_mesenchymal_cap) mes_cells <- sample(mes_cells, PARAM$infercnv_mesenchymal_cap)
      ref_cells <- unlist(lapply(ref_types, function(rt) {
        x <- colnames(integrated)[ct == rt]
        if(length(x) > PARAM$infercnv_ref_cap) sample(x, PARAM$infercnv_ref_cap) else x
      }))
      cells_cnv <- c(mes_cells, ref_cells)
      cat("  inferCNV输入:", length(cells_cnv), "细胞(间质", length(mes_cells),
          "+ 参考", length(ref_cells), ")\n")
      stopifnot("间质/参考细胞数异常,检查上面的注释表" = length(mes_cells) > 100 & length(ref_cells) > 100)
      anno <- data.frame(group = ifelse(ct[cells_cnv] %in% ref_types, ct[cells_cnv], "query"),
                         row.names = cells_cnv)
      cnv_obj <- CreateInfercnvObject(
        raw_counts_matrix = GetAssayData(integrated, layer="counts")[, cells_cnv],
        annotations_file  = anno, gene_order_file = go_path,
        ref_group_names   = ref_types)
      cnv_obj <- infercnv::run(cnv_obj, cutoff=PARAM$infercnv_cutoff,
                               out_dir="Results/SingleCell/inferCNV",
                               cluster_by_groups=TRUE, denoise=TRUE, HMM=FALSE,
                               num_threads=2)
      cat(" Figure5b cnv is saved after infercnv::run \n")
      rm(cnv_obj); gc()
    }
    txt_d <- "Results/SingleCell/inferCNV/infercnv.observations_denoised.txt"
    txt_p <- "Results/SingleCell/inferCNV/infercnv.observations.txt"
    if(file.exists(txt_d) || file.exists(txt_p)) {
      f <- if(file.exists(txt_d)) txt_d else txt_p
      tmp <- fread(f, data.table=FALSE)
      obs <- as.matrix(tmp[,-1]); rownames(obs) <- tmp[[1]]; rm(tmp)
    } else {
      stopifnot(file.exists(fin_rd))
      cobj <- readRDS(fin_rd)
      obs <- cobj@expr.data       
      rm(cobj)
    }
    gc()
    base_chk <- mean(obs[, seq_len(min(200, ncol(obs)))])
    stopifnot(abs(base_chk - 1) < 0.5)
    cnv_score <- colMeans((obs-1)^2)
    saveRDS(cnv_score, cnv_rds); rm(obs); gc()
    }
  ct <- integrated$celltype_coarse; names(ct) <- colnames(integrated)
  ref_types <- intersect(c("Macrophages","T_cells","NK_cells","B_cells","Endothelial","Mast_cells"),
                         unique(ct))
  sc <- rep(NA_real_, ncol(integrated)); names(sc) <- colnames(integrated)
  sc[names(cnv_score)] <- cnv_score
  integrated$cnv_score <- unname(sc)
  cnv_thr <- quantile(cnv_score[ct[names(cnv_score)] %in% ref_types],
                      PARAM$cnv_ref_q, na.rm=TRUE)
  integrated$cnv_high <- !is.na(integrated$cnv_score) & integrated$cnv_score > cnv_thr
  mesenchymal <- integrated$celltype_coarse %in% c("SmoothMuscle","Fibro_CAF")
  cl_frac <- tapply(integrated$cnv_high[mesenchymal], integrated$seurat_clusters[mesenchymal], mean)
  mal_clusters <- names(cl_frac)[which(cl_frac > PARAM$mal_cluster_frac)]
  integrated$celltype_final <- integrated$celltype_coarse
  integrated$celltype_final[mesenchymal &  integrated$seurat_clusters %in% mal_clusters] <- "Malignant"
  integrated$celltype_final[mesenchymal & !integrated$seurat_clusters %in% mal_clusters] <- "CAF"
  saveRDS(integrated, "Results/SingleCell/Integrated_v3_annotated.rds")
  acta2_expr <- GetAssayData(integrated, layer="data")["ACTA2",]
  src_tab <- table(integrated$celltype_final[acta2_expr > 0])
  p_a1 <- DimPlot(integrated, reduction="umap.harmony", group.by="celltype_final",
               label=TRUE, repel=TRUE, pt.size=0.3, raster=TRUE) +
  ggtitle("Cell types (CNV-adjudicated)") + theme(plot.title=element_text(hjust=0.5, face="bold"))
  ggsave("Figures/SingleCell/Figure5a.png", p_a1, width = 7, height = 6, dpi = 600, bg = "white")
  p_a2 <- FeaturePlot(integrated, features="ACTA2", reduction="umap.harmony",
                      cols=c("lightgrey","red"), pt.size=0.3, raster=TRUE) +
    ggtitle("ACTA2 expression") + theme(plot.title=element_text(hjust=0.5, face="bold"))
  ggsave("Figures/SingleCell/Figure5c.png", p_a2, width = 7, height = 6, dpi = 600, bg = "white")
  p_a3 <- VlnPlot(integrated, features="ACTA2", group.by="celltype_final", pt.size=0) +
    NoLegend() + ggtitle("ACTA2 by cell type") +
    theme(axis.text.x=element_text(angle=45,hjust=1), plot.title=element_text(hjust=0.5, face="bold"))
  ggsave("Figures/SingleCell/Figure5d.png", p_a3, width = 8, height = 5, dpi = 600, bg = "white")
  df_v <- integrated@meta.data[!is.na(integrated$cnv_score), , drop = FALSE]
  p_a4 <- ggplot(df_v, aes(celltype_final, cnv_score, fill = celltype_final)) +
    geom_violin(scale = "width") +
    geom_hline(yintercept = cnv_thr, linetype = 2, color = "red") +
    theme_classic() + NoLegend() +
    labs(x = NULL, y = "CNV score") +
    ggtitle("CNV burden (inferCNV)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave("Figures/SingleCell/Figure5e.png", p_a4, width = 8, height = 5, dpi = 600, bg = "white")
  p_s2 <- DotPlot(integrated, features=unique(unlist(cell_markers)), group.by="celltype_final") +
    RotatedAxis() + ggtitle("Canonical markers")
  ggsave("Figures/SingleCell/FigureS2.png", p_s2, width=12, height=5, dpi=600, bg="white")
  mal_cells <- colnames(integrated)[integrated$celltype_final=="Malignant"]
  mal <- subset(integrated, cells=mal_cells)
  mal <- AddModuleScore(mal, features=list(c("ACTA2","MYH11","TAGLN","CNN1")), name="DIFF", ctrl=50, verbose=FALSE)
  ds <- rep(NA_real_, ncol(integrated)); names(ds) <- colnames(integrated)
  tmp_ds <- mal$DIFF1; names(tmp_ds) <- colnames(mal)
  ds[names(tmp_ds)] <- tmp_ds
  integrated$diff_score <- unname(ds)
  diff_med <- tapply(integrated$diff_score[!is.na(integrated$diff_score)],
                     integrated$orig.ident[!is.na(integrated$diff_score)], median)
  mac_frac <- prop.table(table(integrated$orig.ident, integrated$celltype_final),1)[,"Macrophages"]
  df_b <- data.frame(sample=names(diff_med), diff_score=as.numeric(diff_med),
                     mac_frac=as.numeric(mac_frac[names(diff_med)]))
  rho <- cor.test(df_b$diff_score, df_b$mac_frac, method="spearman")
  p_b1 <- ggplot(df_b, aes(diff_score, mac_frac, label=sample)) +
    geom_point(size=3, color="#C0392B") + geom_text_repel(size=3) +
    geom_smooth(method="lm", se=FALSE, linetype=2, color="grey40") +
    labs(x="Malignant-cell differentiation score", y="Macrophage fraction") + theme_bw()
  ggsave("Figures/SingleCell/Figure5f.png", p_b1, width=6, height=4, dpi=600, bg="white")
  df_apoe <- FetchData(mal, vars = "APOE")
  df_apoe$diff_score <- mal$DIFF1
  df_apoe$sample     <- mal$orig.ident
  colnames(df_apoe)[1] <- "APOE"
  ct <- suppressWarnings(cor.test(df_apoe$APOE, df_apoe$diff_score, method = "spearman"))
  df_s <- do.call(rbind, lapply(split(df_apoe, df_apoe$sample), function(d)
    data.frame(apoe_mean = mean(d$APOE),
               apoe_pos  = 100 * mean(d$APOE > 0),
               diff_med  = median(d$diff_score))))
  df_s$sample <- rownames(df_s)
  mac_cells <- colnames(integrated)[integrated$celltype_final == "Macrophages"]
  df_mac <- FetchData(integrated, vars = "APOE", cells = mac_cells)
  df_mac$sample <- integrated$orig.ident[mac_cells]
  p1 <- ggplot(df_apoe, aes(diff_score, APOE, color = sample)) +
    geom_point(size = 0.3, alpha = 0.3) +
    geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.8) +
    labs(x = "Differentiation score", y = "APOE (log-normalized)",
         title = paste0("Malignant cells, rho = ", round(unname(ct$estimate), 3))) +
    theme_bw()
  p2 <- ggplot(df_s, aes(diff_med, apoe_pos, label = sample)) +
    geom_point(size = 3, color = "#C0392B") + ggrepel::geom_text_repel(size = 3) +
    labs(x = "Sample median differentiation score", y = "APOE+ malignant cells (%)") +
    theme_bw()
  ggsave("Figures/SingleCell/FigB2_Malignant_APOE_vs_Diff.png", p1 | p2,
         width = 12, height = 5, dpi = 600, bg = "white")
  write.csv(df_s, "Results/SingleCell/APOE_vs_diffscore_by_sample.csv", row.names = FALSE)
  cellchat_rds <- "Results/SingleCell/CellChat_v3.rds"
  mal_idx <- integrated$celltype_final == "Malignant"
  cut_med <- median(integrated$diff_score[mal_idx], na.rm = TRUE)
  integrated$cc_group <- integrated$celltype_final
  integrated$cc_group[mal_idx & integrated$diff_score >= cut_med] <- "Malignant_ACTA2hi"
  integrated$cc_group[mal_idx & integrated$diff_score <  cut_med] <- "Malignant_ACTA2lo"
  if(file.exists(cellchat_rds)) {
   cellchat <- readRDS(cellchat_rds)
  } else {
    suppressPackageStartupMessages(library(CellChat))
    data.input <- GetAssayData(integrated, layer = "data")  
    meta_cc <- data.frame(group = factor(integrated$cc_group),
                          row.names = colnames(integrated))
    cellchat <- createCellChat(object = data.input, meta = meta_cc, group.by = "group")
    cellchat@DB <- CellChatDB.human                         
    cellchat <- subsetData(cellchat)
    future::plan("multisession", workers = 2)               
    cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat)                  
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    saveRDS(cellchat, cellchat_rds); gc()
  }
  
  df_cc <- subsetCommunication(cellchat,
                               sources.use = c("Malignant_ACTA2hi","Malignant_ACTA2lo"),
                               targets.use = "Macrophages")
  df_cc <- df_cc[order(-df_cc$prob), ]
  write.csv(df_cc, "Results/SingleCell/Supplementary_Table_S1.csv", row.names = FALSE)
  p_c1 <- netVisual_bubble(cellchat,
                           sources.use = c("Malignant_ACTA2hi","Malignant_ACTA2lo"),
                           targets.use = "Macrophages",
                           remove.isolate = FALSE,
                           angle.x = 45,
                           return.data = FALSE) +
    ggtitle("Malignant → Macrophage signaling") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave("Figures/SingleCell/FigC_a_bubble.png", p_c1, width = 6, height = 4, dpi = 600, bg = "white")
  library(ggplot2)
  library(dplyr)
  df_forward <- subsetCommunication(cellchat,
                                    sources.use = c("Malignant_ACTA2hi","Malignant_ACTA2lo"),
                                    targets.use = "Macrophages",
                                    slot.name = "net")
  df_forward <- df_forward %>%
    mutate(pathway_name = factor(pathway_name),
           pair = paste0(ligand, " – ", receptor))
  df_f <- df_forward %>% filter(prob > 0.05) %>%
    arrange(desc(prob))
  p_forward <- ggplot(df_f, aes(x = source, y = pair, size = prob, color = pathway_name)) +
    geom_point(alpha = 0.4) +
    scale_size_continuous(range = c(2, 6), name = "Communication\nprobability") +
    labs(x = "", y = "Ligand – Receptor pair",
         title = "Malignant → Macrophage signaling") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(face = "italic"),
          panel.grid.major.y = element_line(linewidth = 0.1, color = "grey90"),
          legend.position = "right")
  ggsave("Figures/SingleCell/Figure5g.png", p_forward, width = 6, height = 6, dpi = 600)
  png("Figures/SingleCell/Figure5i.png", width = 2200, height = 2200, res = 300)
  netVisual_aggregate(cellchat, signaling = "MIF", layout = "circle")
  dev.off()
  df_rev <- subsetCommunication(cellchat, sources.use = "Macrophages",
                          targets.use = c("Malignant_ACTA2hi","Malignant_ACTA2lo"))
  df_rev <- df_rev[order(-df_rev$prob), ]
  write.csv(df_rev, "Results/SingleCell/Supplementary_Table_S2.csv", row.names = FALSE)
  p_c2 <- netVisual_bubble(cellchat,
                           sources.use = "Macrophages",
                           targets.use = c("Malignant_ACTA2hi","Malignant_ACTA2lo"),
                           remove.isolate = FALSE,
                           angle.x = 45,
                           return.data = FALSE) +
    ggtitle("Macrophage → Malignant signaling") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave("Figures/SingleCell/Figue5h.png", p_c2, width = 6, height = 6, dpi = 600, bg = "white")
  df_reverse <- subsetCommunication(cellchat,
                                    sources.use = "Macrophages",
                                    targets.use = c("Malignant_ACTA2hi","Malignant_ACTA2lo"),
                                    slot.name = "net")
  df_reverse <- df_reverse %>%
    mutate(pathway_name = factor(pathway_name),
           pair = paste0(ligand, " – ", receptor))
  df_r <- df_reverse %>% filter(prob > 0.005) %>%
    arrange(desc(prob))
  
  p_reverse <- ggplot(df_r, aes(x = target, y = pair, size = prob, color = pathway_name)) +
    geom_point(alpha = 0.4) +
    scale_size_continuous(range = c(2, 6), name = "Communication\nprobability") +
    labs(x = "", y = "Ligand – Receptor pair",
         title = "Macrophage → Malignant signaling") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(face = "italic"),
          panel.grid.major.y = element_line(linewidth = 0.1, color = "grey90"),
          legend.position = "right")

  png("Figures/SingleCell/Figure5j.png", width = 2200, height = 2200, res = 300)
  netVisual_aggregate(cellchat, signaling = "SPP1", layout = "circle")
  dev.off()
  
  if(!exists("integrated")) {
    integrated <- readRDS("Results/SingleCell/Integrated_v3_annotated.rds")
  } else {
    # cat("使用当前环境中的 integrated 对象\n")
  }
  checkpoint_genes <- c("HAVCR2","LGALS9","CD86","NECTIN3","CD200","CD274","TIGIT")
  present_genes <- checkpoint_genes[checkpoint_genes %in% rownames(integrated)]
  missing_genes <- setdiff(checkpoint_genes, present_genes)
  expr_data <- GetAssayData(integrated, assay = "RNA", layer = "data")
  cell_types <- integrated$celltype_final
  checkpoint_stats <- data.frame(
    Gene = character(),
    Cell_type = character(),
    Mean_expression = numeric(),
    Percent_positive = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(g in present_genes) {
    g_expr <- expr_data[g, ]
    for(ct in sort(unique(cell_types))) {
      idx <- which(cell_types == ct)
      avg_exp <- mean(g_expr[idx])
      pct_pos <- mean(g_expr[idx] > 0) * 100
      checkpoint_stats <- rbind(checkpoint_stats, data.frame(
        Gene = g,
        Cell_type = ct,
        Mean_expression = round(avg_exp, 4),
        Percent_positive = round(pct_pos, 1),
        stringsAsFactors = FALSE
      ))
    }
  }
  write.csv(checkpoint_stats, 
            "Results/SingleCell/Supplementary_Table_S5.csv", 
            row.names = FALSE)
  for(g in present_genes) {
    g_expr <- expr_data[g, ]
    for(ct in sort(unique(cell_types))) {
      idx <- which(cell_types == ct)
      avg_exp <- mean(g_expr[idx])
      pct_pos <- mean(g_expr[idx] > 0) * 100
     
    }
  }
  library(Seurat)
  library(ggplot2)
  p_dot_main <- DotPlot(integrated, 
                        features = present_genes, 
                        group.by = "celltype_final", 
                        scale = FALSE) +
    RotatedAxis() +
    scale_colour_gradientn(
      colours = c("#EEEEEE", "#4A90D9", "#D62728")
      # name = "Mean expression"      
    ) +
    scale_size_continuous(
      range = c(1, 6)
      # name = "Percent positive"     
    ) +
    
    guides(
      colour = guide_colourbar(title = "Mean expression"),
      size   = guide_legend(title = "Percent positive")
    ) +
    ggtitle("Checkpoint gene expression by cell type") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )
  
  ggsave("Figures/SingleCell/Figure6c.png", p_dot_main, 
         width = 8, height = 5, dpi = 600, bg = "white")
  
  cat("scRNA  end")

  
  
  
  
  
  
 