
library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)

#future::plan("multicore",workers=6)
#options(future.globals.maxSize = 4000 * 1024^2)

ds_use <- readRDS("CR_All_harmony_cluster.rds")
samName="Immune_all"
ds_use <- subset(ds_use, subset = (seurat_clusters %in% c("4","21")))
  
DefaultAssay(ds_use) <- "RNA"

ds_use<-NormalizeData(ds_use, normalization.method = "LogNormalize", scale.factor = 10000)
ds_use<-FindVariableFeatures(ds_use,selection.method = "vst", nfeatures =5000)

all.genes <- rownames(ds_use)
ds_use<-ScaleData(ds_use,features = all.genes)
ds_use<-RunPCA(ds_use)
ds_use <- RunHarmony(ds_use, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony.rna",assay.use = "RNA")

ds_use<- RunUMAP(ds_use, dims=1:50, reduction = "pca",reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
ds_use <- RunTSNE(ds_use, dims=1:50,reduction = "pca", reduction.name = 'tsne.rna', reduction.key = 'rnaTSNE_')
ds_use<- RunUMAP(ds_use, dims=1:50, reduction = "harmony.rna",reduction.name = 'umap.rna.harmony', reduction.key = 'rnaUMAPHarmony_')
ds_use <- RunTSNE(ds_use, dims=1:50,reduction = "harmony.rna", reduction.name = 'tsne.rna.harmony', reduction.key = 'rnaTSNEHarmony_')

DefaultAssay(ds_use) <- "ATAC"
ds_use<- RunTFIDF(ds_use)
ds_use <- FindTopFeatures(ds_use, min.cutoff='q0')
ds_use <- RunSVD(ds_use)
ds_use <- RunHarmony(ds_use, group.by.vars = "orig.ident", reduction = "lsi", reduction.save = "harmony.atac",assay.use = "ATAC", project.dim = F)

ds_use <- RunUMAP(ds_use, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
ds_use <- RunTSNE(ds_use, reduction = 'lsi', dims = 2:50, reduction.name = 'tsne.atac', reduction.key = 'atacTSNE_')
ds_use <- RunUMAP(ds_use, reduction = 'harmony.atac', dims = 2:50, reduction.name = "umap.atac.harmony", reduction.key = "atacUMAPHarmony_")
ds_use <- RunTSNE(ds_use, reduction = 'harmony.atac', dims = 2:50, reduction.name = 'tsne.atac.harmony', reduction.key = 'atacTSNEHarmony_')

ds_use <- FindMultiModalNeighbors(ds_use, reduction.list = list("harmony.rna", "harmony.atac"), dims.list = list(1:50, 2:50))
ds_use <- RunTSNE(ds_use, nn.name = "weighted.nn", reduction.name = "wnn.tsne.harmony", reduction.key = "wnnTSNEHarmony_")
ds_use <- RunUMAP(ds_use, nn.name = "weighted.nn", reduction.name = "wnn.umap.harmony", reduction.key = "wnnUMAPHarmony_")
ds_use <- FindClusters(ds_use, graph.name = "wsnn", algorithm = 3, verbose = FALSE,resolution=0.5)

print("Finish cluster!") 
saveRDS(ds_use, file=paste0(samName,"_harmony_cluster.rds"))

df2 <- as.data.frame(ds_use@meta.data)
write.csv(df2, file=paste0(samName,"_meta.csv"))

 
markers <- FindAllMarkers(ds_use, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2)

write.table(markers,file=paste0(samName,"_allmarker_gene0.1.txt"),sep="\t")

top30 <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.table(top30,file=paste0(samName,"_top30markergene0.1.txt"),sep="\t")

top30 = CaseMatch(search = as.vector(top30$gene), match = rownames(ds_use))

pdf("top30markersgene0.1.pdf", height=50, width=20)  
DoHeatmap(ds_use, features = top30, group.by = "seurat_clusters")
dev.off()

##annotation

Tmarkers <- c("CD4","CD8A","CD8B","CD3D","CD3E","CD3G","CCR7","TCF7","LEF1","GZMK","GZMA","PRF1","GNLY","SLC4A10","TRGV9","KLRF1")

Tmarkers2 <- c("CD4","CD8A","CD8B","CD69","SELL","CCR7","LEF1","FOXP3","ICOS","NKG7","GZMA","GNLY","KLRD1","PDCD1","CXCL13","CCNB2","DLGAP5")

Bmarkers <- c("CD19","CD24","CD27","CD38","SDC1","MS4A1","CD74","HLA-DRA","HMGN2","H2AFZ","FTH1","TYROBP","FCER1G","AGR2","MT1G","TSPAN8","PAX5")

MHC2 <- c("HLA-DRA","HLA-DRB1","HLA-DRB3","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1")

NKmarkers <- c("NKG7","KLRC1","KLRC2","KLRD1","KLRK1")

TAMfeatures <- c("HLA-DR","CD14","CX3CR1")
M1features <- c("CD68","CD80","CD86","IL6","IL1B","STAT1","IFN","CXCL10","IRF1","IRF5","TNFA","NO","ROS","CXCL9","MHC11","TLR")
M2features <- c("CD163","MSR1","MRC1","STAT3","MMP2","MMP9","IL10","IL1RA","CSF1R","TNFB","CCL2","CCL17","CCL18","CCL24","TGFB1","TGFB2","CCL22")


ds <- subset(ds,seurat_clusters %in% c(0:19))

ds@meta.data$Anno <- factor(ds@meta.data$Anno, levels=c("TC_CD4_FOXP3","TC_CD4_CXCL13","TC_CD4_S100A6","TC_CD8_IL7R","TC_CD8_IL4R",
"BC_CD38","BC_MS4A1","M_MRC1","M_SPP1","M_HLA-DRB5","M_MKI67","M_IGKC","NK_GNLY","DC","unknow"))


markers <- read.table("AnnoTcell_allmarkergene0.25.txt")
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
exprMatrix <- tcell@assays$RNA@data[which(rownames(tcell) %in% top15$gene),]
df <- as.data.frame(exprMatrix)
df$gene <- rownames(df)
df_melt <- reshape2::melt(df, id.vars="gene",
                        variable.name = "cell",
                        value.name = "value")
info <- data.frame(cell = rownames(tcell@meta.data), cluster=tcell@meta.data$Anno)

df_merge <- dplyr::left_join(x=df_melt,y=info,by="cell")
df_use <- df_merge %>%
dplyr::group_by(cluster, gene) %>%
dplyr::summarise(avg=mean(value))

df_use2 <- reshape2::dcast(df_use,gene~cluster,
                        value.var = "avg")
rownames(df_use2) <- df_use2$gene
df_use2 <- df_use2[,-1]

pm <- as.matrix(df_use2)
pickann_col <- data.frame(cluster=colnames(pm), row.names=colnames(pm))

pickstrain_color <-c("#DFF6FF","#90E0EF","#00B4D8","#344CB7","#000957","#7E370C")
names(pickstrain_color) <- c("TC_CD4_FOXP3","TC_CD4_CXCL13","TC_CD4_S100A6","TC_CD8_IL7R","TC_CD8_IL4R","NK_GNLY")
pickann_colors <- list(cluster = pickstrain_color)
ordergene <-  intersect(top15$gene,rownames(pm))
pm2 <- pm[c(ordergene),]

tiff("Tdiff_heatmap.tiff", units="in", width=4,height=13,res=300)
fig <- pheatmap(pm2, color = colorRampPalette(c("#6496A4", "white", "#DA3D4E"))(10),
            fontsize_col = 11,
            annotation_col = pickann_col,
            annotation_colors = pickann_colors,
            scale = "row",
            border=FALSE,
            treeheight_row = 20,
            fontsize = 6,
            #angle_col = "45",
            cellwidth = 20, cellheight = 8,
            # gaps_row = c(10, 14),
            # cutree_rows = 6,
            main = "ANNO_heatmap_cluster",
            cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T)
dev.off()

