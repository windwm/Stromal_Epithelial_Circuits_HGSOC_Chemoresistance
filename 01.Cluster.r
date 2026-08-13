library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)

set.seed(123)
samName = "CR_All";
print(samName)
setwd("workpath")  
ds_use <- readRDS(seurat_obj)

DefaultAssay(ds_use) <- "RNA"

ds_use<-NormalizeData(ds_use, normalization.method = "LogNormalize", scale.factor = 10000)
ds_use<-FindVariableFeatures(ds_use,selection.method = "vst", nfeatures =3000)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.genes <- rownames(ds_use)
ds_use<- CellCycleScoring(object = ds_use, g2m.features = g2m.genes, s.features = s.genes)
ds_use<-ScaleData(ds_use,features = all.genes, vars.to.regress =  c("S.Score", "G2M.Score","percent.mt"))
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
ds_use <- FindClusters(ds_use, graph.name = "wsnn", algorithm = 3, verbose = FALSE,resolution=0.2)

print("Finish cluster!") 
saveRDS(ds_use, file=paste0(samName,"_harmony_cluster.rds"))

df <- ds_use@meta.data

write.table(df,"CR_All_res0.2_hvg5000_meta.txt")
p1 <- DimPlot(ds_use, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE,pt.size = 1) + ggtitle("RNA")+labs(tag = "A")
ggsave(paste0(samName,"_rna_umap.pdf"), plot = p1, width = 15, height = 12)
p2 <- DimPlot(ds_use, reduction = "umap.rna", group.by="Patient" ,label = TRUE, label.size = 2.5, repel = TRUE,pt.size = 1) + ggtitle("Patient")+labs(tag = "B")
ggsave(paste0(samName,"_rna_patient_umap.pdf"), plot = p2, width = 15, height = 12)    
p3 <- DimPlot(ds_use, reduction = "umap.rna.harmony", label = TRUE, label.size = 2.5, repel = TRUE,pt.size = 1) + ggtitle("RNA_Harmony")+labs(tag = "C")
ggsave(paste0(samName,"_rna_harmony_umap.pdf"), plot = p3, width = 15, height = 12)
p4 <- DimPlot(ds_use, reduction = "umap.rna.harmony", group.by="Patient" ,label = TRUE, label.size = 2.5, repel = TRUE,pt.size = 1) + ggtitle("Patient")+labs(tag = "D")
ggsave(paste0(samName,"_rna_harmony_patient_umap.pdf"), plot = p4, width = 15, height = 12)

p <- (p1|p2)/(p3|p4)
ggsave(paste0(samName,"_rna_tsne_all.pdf"), plot = p, width = 20, height = 20)

p1 <- DimPlot(ds_use, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE,pt.size = 1) + ggtitle("ATAC")+labs(tag = "A")
ggsave(paste0(samName,"_atac_umap.pdf"), plot = p1, width = 15, height = 12)
p2 <- DimPlot(ds_use, reduction = "umap.atac", group.by="Patient" ,label = TRUE, label.size = 2.5, repel = TRUE,pt.size = 1) + ggtitle("Patient")+labs(tag = "B")
ggsave(paste0(samName,"_atac_patient_umap.pdf"), plot = p2, width = 15, height = 12)    
p3 <- DimPlot(ds_use, reduction = "umap.atac.harmony", label = TRUE, label.size = 2.5, repel = TRUE,pt.size = 1) + ggtitle("ATAC_Harmony")+labs(tag = "C")
ggsave(paste0(samName,"_atac_harmony_umap.pdf"), plot = p3, width = 15, height = 12)
p4 <- DimPlot(ds_use, reduction = "umap.atac.harmony", group.by="Patient" ,label = TRUE, label.size = 2.5, repel = TRUE,pt.size = 1) + ggtitle("Patient")+labs(tag = "D")
ggsave(paste0(samName,"_atac_harmony_patient_umap.pdf"), plot = p4, width = 15, height = 12)

p <- (p1|p2)/(p3|p4)
ggsave(paste0(samName,"_atac_umap_all.pdf"), plot = p, width = 20, height = 20)

