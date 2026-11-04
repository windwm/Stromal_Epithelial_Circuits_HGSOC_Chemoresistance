samName       <- "Ovarian_CR_all"
gtf_path      <- "gencode.v32.primary_assembly.annotation.gtf.gz"
out_dir       <- "ATAC"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cd55_chr      <- "chr1"
cd55_start    <- 207320678
cd55_end      <- 207400966

ido1_chr      <- "chr8"
ido1_start    <- 399000
ido1_end      <- 399300

subtype_list <- c(
  "EC-Mix1","EC-Mix2","EC-Mix3","EC-Mix4",
  "EC-PS1","EC-PS2","EC-PS3","EC-PS4","EC-PS5",
  "EC-PS6","EC-PS7","EC-PS8","EC-PS9","EC-PS10"
)

cicero_sample_num <- 100
coaccess_cutoff_plot <- 0.25

library(Seurat)
library(Signac)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(cicero)
library(rtracklayer)
library(SeuratWrappers)
library(EnsDb.Hsapiens.v86)


topTFs <- function(markers_rna, markers_motifs, celltype, padj.cutoff = 1e-2) {
  ctmarkers_rna <- markers_rna %>%
    dplyr::filter(RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>%
    arrange(-RNA.auc)

  ctmarkers_motif <- markers_motifs %>%
    dplyr::filter(motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>%
    arrange(-motif.auc)

  top_tfs <- inner_join(
    x = ctmarkers_rna[, c("RNA.group", "gene", "RNA.padj", "RNA.logFC")],
    y = ctmarkers_motifs[, c("motif.group", "gene", "motif.padj", "motif.logFC", "motif.auc")],
    by = "gene"
  )
  top_tfs$avg_auc <- (ctmarkers_rna$RNA.auc + ctmarkers_motifs$motif.auc) / 2
  top_tfs <- arrange(top_tfs, -avg_auc)
  return(top_tfs)
}

message("==== 加载hg38注释信息 ====")
annotations_ensdb <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations_ensdb) <- "hg38"
genome_sizes <- seqlengths(annotations_ensdb)
genome.df <- data.frame(chr = names(genome_sizes), length = genome_sizes)
genome.df <- genome.df[1:24, ] # 保留标准染色体1-22,X,Y

gene_anno <- rtracklayer::readGFF(gtf_path)
gene_anno$chromosome <- gene_anno$seqid
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name
gene_anno1 <- gene_anno[, c("chromosome","start","end","strand","type","gene","transcript","symbol")]

# -------------------------- 5. 【可选】整体Cicero + GeneActivity（所有细胞一起跑） --------------------------
# ds_use <- readRDS("All-ATAC-motif.rds")
#
# message("==== 全局Cicero：构建cds ====")
# cds <- as.cell_data_set(x = ds_use, assay = "ATAC", default.reduction = "wnn.umap.harmony")
# cicero_obj <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)$WNN.UMAP.HARMONY)
# saveRDS(cicero_obj, file.path(out_dir, paste0(samName,"_cicero.rds")))
#
# message("==== run_cicero 共开放连接 ====")
# conns_all <- run_cicero(cicero_obj, genomic_coords = genome.df, sample_num = cicero_sample_num)
# saveRDS(conns_all, file.path(out_dir, paste0(samName,"_conns.rds")))
#
# ccans_all <- generate_ccans(conns_all)
# saveRDS(ccans_all, file.path(out_dir, paste0(samName,"_ccans.rds")))
#
# # 构建promoter注释（TSS上下游1k）
# annotations_gtf <- import(gtf_path)
# geneAnno <- unique(annotations_gtf[
#   annotations_gtf$type=="gene" & seqnames(annotations_gtf) %in% standardChromosomes(annotations_gtf),
#   c("gene_name")
# ])
# promoter <- promoters(geneAnno, upstream=1000, downstream=1000)
# promoterDF <- data.frame(
#   chromosome = seqnames(promoter),
#   start = start(promoter),
#   end = end(promoter),
#   gene = promoter$gene_name
# ) %>% unique()
#
# cicero_obj <- annotate_cds_by_site(cicero_obj, promoterDF)
#
# # Cicero gene activity
# unnorm_ga <- build_gene_activity_matrix(cicero_obj, conns_all)
# unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga)==0, !Matrix::colSums(unnorm_ga)==0]
# num_genes <- pData(cicero_obj)$num_genes_expressed
# names(num_genes) <- rownames(pData(cicero_obj))
# cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
# saveRDS(cicero_gene_activities, file.path(out_dir, paste0(samName,"_cicero_gene_activities_normalized.rds")))
#
# # 新建Seurat对象存放Cicero GA
# ds_cicero <- CreateSeuratObject(
#   counts = cicero_gene_activities,
#   assay = "GA_cicero",
#   project = samName,
#   meta.data = ds_use@meta.data
# )
# ds_cicero[["GA_cicero_unnorm"]] <- CreateAssayObject(counts = unnorm_ga)
# ds_cicero <- NormalizeData(
#   object = ds_cicero,
#   assay = "GA_cicero_unnorm",
#   normalization.method = "LogNormalize",
#   scale.factor = median(ds_use$nCount_ATAC)
# )

# -------------------------- 6. 【重要】批量分subtype跑Cicero（EC各个亚群分开） --------------------------
# ds <- readRDS("你的主seurat对象.rds") # !!! 这里你需要自己读入ds
# ds@meta.data$maintype <- meta$maintype[match(ds@meta.data$cellname, meta$X)]
# ds@meta.data$subtype  <- meta$subtype[match(ds@meta.data$cellname, meta$X)]

for (ct in subtype_list) {
  message("\n>>>> Processing subtype: ", ct)
  out_cicero   <- file.path(out_dir, paste0(ct,"_cicero.rds"))
  out_conns    <- file.path(out_dir, paste0(ct,"_conns.rds"))
  out_ccans    <- file.path(out_dir, paste0(ct,"_ccans.rds"))
  out_links    <- file.path(out_dir, paste0(ct,"_links.rds"))
  out_pdf_cd55 <- file.path(out_dir, paste0(ct,"_CD55.pdf"))

  # 跳过：结果文件全部存在就不重复计算
  if(all(file.exists(out_cicero, out_conns, out_ccans, out_links, out_pdf_cd55))){
    message("文件已存在，跳过：", ct)
    next
  }

  ds_use <- subset(ds, subtype == ct)
  if(ncol(ds_use) < 20){
    message("细胞数太少，跳过：", ct)
    next
  }

  # cicero cds
  cds <- as.cell_data_set(x = ds_use, assay="ATAC", default.reduction = "wnn.umap.harmony")
  cicero_obj <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)$WNN.UMAP.HARMONY)
  saveRDS(cicero_obj, out_cicero)

  # run cicero
  conns <- run_cicero(cicero_obj, genomic_coords = genome.df, sample_num = cicero_sample_num)
  saveRDS(conns, out_conns)

  # plot CD55
  pdf(out_pdf_cd55, width = 8, height = 4)
  p <- plot_connections(
    conns,
    chr = cd55_chr,
    start = cd55_start,
    end = cd55_end,
    gene_model = gene_anno1,
    coaccess_cutoff = 0,
    connection_width = 1,
    collapseTranscripts = "longest"
  )
  print(p)
  dev.off()

  # CCAN and links
  ccans <- generate_ccans(conns)
  saveRDS(ccans, out_ccans)
  links <- ConnectionsToLinks(conns = conns, ccans = ccans)
  saveRDS(links, out_links)
}

# conns_Mix2 <- readRDS(file.path(out_dir, "EC-Mix2_conns.rds"))
# M1 <- readRDS(file.path(out_dir, "EC-Mix1_conns.rds"))
# M2 <- conns_Mix2
#
# pdf(file.path(out_dir, "EC-Mix2_IDO1.pdf"), width=8, height=4)
# plot_connections(
#   conns_Mix2,
#   chr = ido1_chr, start = ido1_start, end = ido1_end,
#   gene_model = gene_anno1,
#   coaccess_cutoff = coaccess_cutoff_plot,
#   connection_width = 1,
#   collapseTranscripts = "longest"
# )
# dev.off()
#
# pdf(file.path(out_dir, "EC-Mix1_2_CD55.pdf"), width=8, height=4)
# plot_connections(
#   M1,
#   chr = cd55_chr, start = cd55_start, end = cd55_end,
#   gene_model = gene_anno1,
#   coaccess_cutoff = 0.1,
#   connection_width = 0.5,
#   comparison_track = M2,
#   include_axis_track = FALSE,
#   collapseTranscripts = "longest"
# )
# dev.off()

# --------------------------  Seurat-GeneActivity（GA_seurat） --------------------------
# gene.activities <- GeneActivity(ds, extend.upstream = 500, biotypes = NULL)
# ds[["GA_seurat"]] <- CreateAssayObject(counts = gene.activities)
# ds <- NormalizeData(
#   object = ds,
#   assay = "GA_seurat",
#   normalization.method = "LogNormalize"
# )

# --------------------------  presto chromvar motif marker --------------------------
message("\n==== presto 差异检验：RNA + chromvar motif ====")
markers_rna <- presto::wilcoxauc.Seurat(
  X = ds_ps, group_by = 'subtype', assay = 'data', seurat_assay = 'RNA'
)
markers_motifs <- presto::wilcoxauc.Seurat(
  X = ds_ps, group_by = 'subtype', assay = 'data', seurat_assay = 'chromvar'
)

motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(ds_ps, id = motif.names)

write.table(
  markers_rna,
  file.path(out_dir, paste0(samName,"_presto-markerrna.txt")),
  sep = "\t", quote = F, row.names = F
)
write.table(
  markers_motifs,
  file.path(out_dir, paste0(samName,"_presto-markermotif.txt")),
  sep = "\t", quote = F, row.names = F
)

