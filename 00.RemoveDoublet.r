library(dplyr)
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(scDblFinder)
library(VennDiagram)

set.seed(123)
fileinfo <- read.csv("Rinfo.csv", header=T)

for(i in c(1:nrow(fileinfo))){

  samName = fileinfo$samName[i];
  fileName = fileinfo$fileName[i];
  print(samName)
  print(fileName)
  addr = paste0(samName,"_filted.RData")
  load(addr)
  archdoubpath = paste0(fileName,"_doubletcell_barcode.txt")
  archrdoublet <- read.table(archdoubpath)
  archrsinglpath = paste0(fileName,"_sinletcell_barcode.txt")
  archrsinglet <- read.table(archrsinglpath)
  
  ds_use <- ds.filted

  DefaultAssay(ds_use) <- "RNA"
  ds_use<-NormalizeData(ds_use, normalization.method = "LogNormalize", scale.factor = 10000)
  ds_use<-FindVariableFeatures(ds_use,selection.method = "vst", nfeatures =3000)
  
  all.genes <- rownames(ds_use)
  ds_use<-ScaleData(ds_use,features = all.genes)
  ds_use<-RunPCA(ds_use)
  
  ds_use<-RunUMAP(ds_use, dims=1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  ds_use <- RunTSNE(ds_use, dims=1:50, reduction.name = 'tsne.rna', reduction.key = 'rnaTSNE_')
  
  
  DefaultAssay(ds_use) <- "ATAC"
  ds_use <- RunTFIDF(ds_use)
  ds_use <- FindTopFeatures(ds_use, min.cutoff='q0')
  ds_use <- RunSVD(ds_use)
  ds_use <- RunUMAP(ds_use, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  ds_use <- RunTSNE(ds_use, dims=2:50, reduction.name = 'tsne.atac', reduction.key = 'atacTSNE_')
  
  ds_use <- FindMultiModalNeighbors(ds_use, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
  ds_use <- RunTSNE(ds_use, nn.name = "weighted.nn", reduction.name = "wnn.tsne", reduction.key = "wnnTSNE_")
  ds_use <- RunUMAP(ds_use, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  ds_use <- FindClusters(ds_use, graph.name = "wsnn", algorithm = 3, verbose = FALSE,resolution=0.5)
  
  print("Finish cluster!")
  
  
  DefaultAssay(ds_use) <- "RNA"
  
  sweep.res.list_ds <- paramSweep_v3(ds_use, PCs = 1:50, sct = F)
  sweep.stats_ds <- summarizeSweep(sweep.res.list_ds, GT = FALSE)
  bcmvn_ds <- find.pK(sweep.stats_ds)
  opt_pK <- as.numeric(as.vector(bcmvn_ds$pK[which.max(bcmvn_ds$BCmetric)]))
  print("opt_pK:")
  print(opt_pK)
  DoubletRate = ncol(ds_use)*8*1e-6  #DoubletRate = ncol(ds_use)*7.6*1e-6
  print("DoubletRate:")
  print(DoubletRate)
  annotations <- ds_use@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           
  nExp_poi <- round(DoubletRate*nrow(ds_use@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  pN_use = 0.25
  ds_use <- doubletFinder_v3(ds_use, PCs = 1:30, pN=pN_use, pK = opt_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  
  doubletname <- paste("DF.classifications",pN_use,opt_pK,nExp_poi,sep="_")
  
  ds_use <- doubletFinder_v3(ds_use, PCs = 1:30, pN = pN_use, pK = opt_pK, nExp = nExp_poi.adj, reuse.pANN = doubletname, sct = F)
  
  colnames(ds_use@meta.data)[ncol(ds_use@meta.data)] <- "DoubletFinder"
  print(samName)
  print("DoubletFinder:")
  print(table(ds_use@meta.data$DoubletFinder))
  
  
  print("Finish DoubletFinder!")
  
  sce <- scDblFinder(ds_use@assays$RNA@data,clusters=ds_use@meta.data$seurat_clusters, nfeatures=2000, dims=50, includePCs=50)
  #For 10x data, it is usually safe to leave the dbr empty, and it will be automatically estimated.
  ds_use@meta.data[['scDblFinder_RNA']] <- sce$scDblFinder.class
  print(samName)
  print("scDblFinder_RNA:")
  print(table(ds_use@meta.data$scDblFinder_RNA))
  
  p1 <- DimPlot(ds_use, reduction = "wnn.tsne", group.by = "scDblFinder_RNA")
  p2 <- DimPlot(ds_use, reduction = "wnn.umap", group.by = "scDblFinder_RNA")
  p <- p1+p2  
  ggsave(paste0(samName,"_scDblFinderRNA_precluster.pdf"), plot = p, width = 15, height = 12)
  
  print("Finish scDblFinder_RNA!")
  
  satac <- scDblFinder(ds_use@assays$ATAC@data, aggregateFeatures=TRUE, nfeatures=200, clusters=ds_use@meta.data$seurat_clusters, processing="normFeatures")
  ds_use@meta.data[['scDblFinder_ATAC']] <- satac$scDblFinder.class
  print(samName)
  print("scDblFinder_ATAC:")
  print(table(ds_use@meta.data$scDblFinder_ATAC))
  
  print("Finish scDblFinder_ATAC!")
  
  Cell <- colnames(ds_use)
  DoubletFinder_RNA <- ds_use@meta.data$DoubletFinder
  scDblFinder_RNA <- ds_use@meta.data$scDblFinder_RNA
  scDblFinder_ATAC <- ds_use@meta.data$scDblFinder_ATAC
  ArchR_ATAC <- c(1:length(Cell))
  
  allsoft <- data.frame(Cell, DoubletFinder_RNA, scDblFinder_RNA, scDblFinder_ATAC, ArchR_ATAC)
  allsoft$ArchR_ATAC[which(paste0(fileName,"#", allsoft$Cell) %in% archrdoublet[,1])] = "Doublet"
  allsoft$ArchR_ATAC[which(paste0(fileName,"#", allsoft$Cell) %in% archrsinglet[,1])] = "Singlet"
  allsoft$ArchR_ATAC[which(!(allsoft$ArchR_ATAC %in% c("Doublet","Singlet")) )] = "NotFind"
  
  write.csv(allsoft, file=paste0(fileName,"_douletinfo.csv"))
  
  ds_use@meta.data[['ArchR_ATAC']] <- allsoft$ArchR_ATAC
  p1 <- DimPlot(ds_use, reduction = "wnn.tsne", group.by = "ArchR_ATAC")
  p2 <- DimPlot(ds_use, reduction = "wnn.umap", group.by = "ArchR_ATAC")
  p <- p1+p2  
  ggsave(paste0(samName,"_ArchRATAC_precluster.pdf"), plot = p, width = 15, height = 12)
  
  saveRDS(ds_use,file = paste0(samName,"_all_precluster.Rds"))
  
  draw1 <- allsoft$Cell[which(allsoft$DoubletFinder_RNA == "Singlet")]
  draw2 <- allsoft$Cell[which(allsoft$scDblFinder_RNA == "singlet")]
  draw3 <- allsoft$Cell[which(allsoft$scDblFinder_ATAC == "singlet")]
  draw4 <- allsoft$Cell[which(allsoft$ArchR_ATAC == "Singlet")]
  
  venn_list<- list(DoubletFinder_RNA=draw1,scDblFinder_RNA=draw2, scDblFinder_ATAC=draw3, ArchR_ATAC=draw4)
  
  venn.diagram(venn_list, filename = paste0(samName,"_",fileName,"_singlet_overlap.png"),
             fill = c('#a6611a', '#dfc27d','#80cdc1','#018571'), alpha = 0.50, cat.col = rep('black', 4), 
             col = 'black', cex = 1.5, fontfamily = 'serif', 
             cat.cex = 0.7, cat.fontfamily = 'serif')
             
  draw1 <- allsoft$Cell[which(allsoft$DoubletFinder_RNA == "Doublet")]
  draw2 <- allsoft$Cell[which(allsoft$scDblFinder_RNA == "doublet")]
  draw3 <- allsoft$Cell[which(allsoft$scDblFinder_ATAC == "doublet")]
  draw4 <- allsoft$Cell[which(allsoft$ArchR_ATAC == "Doublet")]
  
  venn_list<- list(DoubletFinder_RNA=draw1,scDblFinder_RNA=draw2, scDblFinder_ATAC=draw3, ArchR_ATAC=draw4)
  
  venn.diagram(venn_list, filename = paste0(samName,"_",fileName,"_doublet_overlap.png"),
             fill = c('#a6611a', '#dfc27d','#80cdc1','#018571'), alpha = 0.50, cat.col = rep('black', 4), 
             col = 'black', cex = 1.5, fontfamily = 'serif', 
             cat.cex = 0.7, cat.fontfamily = 'serif')

}