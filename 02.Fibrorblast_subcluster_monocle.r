library(Seurat)
library(dplyr)
library(Matrix)
library(monocle)

ds <- readRDS("Fibro_harmony_cluster.rds")
ds <- subset(ds,seurat_clusters %in% c(0:18))

ds@meta.data$Anno <- factor( ds@meta.data$Anno, levels=c("Theca cell","CAF-1","CAF-2","CAF-3","Mesenchymal-fib","Mesothelial cell") )

df <- as.data.frame(ds@meta.data)
write.csv(df, file="Fibro_meta.csv")
dim <- ds@reductions$wnn.umap.harmony@cell.embeddings
write.csv(dim,file="Fibro_dim.csv")

fib <- c("COL1A1","COL4A1","COL6A1")
myCAF <- c("ACTA2","RGS5","MYH11","MCAM")
iCAF <- c("IL6","PDGFRA", "IL11","IL8","CXCL1", "CXCL2","CXCL6","CXCL8", "CXCL12")
apCAF <- c("HLA-DRB1","HLA-DRA","CD74")
meso <- c("WT1","KRT7","MSLN")


pdf("alldotplot.pdf", width=12)
DotPlot(ds, features = c(fib,"FAP",myCAF,iCAF,apCAF,meso), cols=c("blue","red"),scale=F) + RotatedAxis()
dev.off()

CAFfeatures <- c("COL1A1","COL6A1","ACTA2","FAP","PDPN","PDGFRA","PDGFRB","MCAM","CAV1","WT1","KRT7","MSLN")

pdf("anno-dotplot.pdf", width=7,height=4)
DotPlot(ds, features = CAFfeatures,group.by="Anno", cols=c("blue","red")) + RotatedAxis()
dev.off()

markers <- FindAllMarkers(ds, min.pct = 0.25, logfc.threshold = 0.25)
samName <- "Fibroblasts"
write.table(markers,file=paste0(samName,"_allmarkergene0.25.txt"),sep="\t")

top50 <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.table(top50,file=paste0(samName,"_top50markergene0.25.txt"),sep="\t")

####MARKER HEATMAP####
top50 = CaseMatch(search = as.vector(top50$gene), match = rownames(ds))

pdf("top50markersgene0.25.pdf", height=70, width=25)  
DoHeatmap(ds,features = top50)
dev.off()

exprMatrix <- ds@assays$RNA@data[which(rownames(ds) %in% top50),]
df <- as.data.frame(exprMatrix)
df$gene <- rownames(df)
df_melt <- reshape2::melt(df, id.vars="gene",
                          variable.name = "cell",
                          value.name = "value")
info <- data.frame(cell = rownames(ds@meta.data), cluster=ds@meta.data$Anno)

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
pickstrain_color <-  c("#F4EEA9","#519259","#3E7C17","#064635","#F0BB62","#A3423C")
names(pickstrain_color) <- c("Theca cell","CAF-1","CAF-2","CAF-3","Mesenchymal-fib","Mesothelial cell")
pickann_colors <- list(cluster = pickstrain_color)

ordergene <-  intersect(top50,rownames(pm))
pm2 <- pm[c(ordergene),]
write.csv(pm2,file=paste0(samName,"_expr.csv"))

tiff("heatmap_annocluster.tiff", units="in", width=6, height=8,res=300)
fig <- pheatmap(pm2, color = colorRampPalette(c("#6496A4", "white", "#DA3D4E"))(10),
                fontsize_col = 11,
                annotation_col = pickann_col,
                annotation_colors = pickann_colors,
                scale = "row",
                border=FALSE,
                treeheight_row = 20,
                fontsize = 6,
                #angle_col = "45",
                cellwidth = 10, cellheight = 1.5,
                # gaps_row = c(10, 14),
                # cutree_rows = 6,
                main = "ANNO_heatmap_cluster",
                cluster_rows = F, cluster_cols = F,show_rownames = F, show_colnames = T)
dev.off()


####MONOCLE####

df <- read.table("Fibroblasts_allmarkergene0.25.txt")
top30 <- df %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

data <- as(as.matrix(ds@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = ds@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#
cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size()); 
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
saveRDS(cds,file="fiboMonocle.Rds")

ordering_genes_used<-as.matrix(top30)
cds <- setOrderingFilter(cds, ordering_genes = ordering_genes_used)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree',auto_param_selection = F)##��ά
cds <- orderCells(cds)
saveRDS(cds,file="fiboMonocle-orderCells.rds")
write.table(pData(cds),file="fiboMonocle-my_pseudotime.txt")

my_pseudotime_de <- differentialGeneTest(cds,fullModelFormulaStr = "~sm.ns(Pseudotime)") 
saveRDS(my_pseudotime_de,file="FIBRO_pseudotime_de.Robj")

#coregene <-subset(my_pseudotime_de, qval < 10^-5)

write.table(my_pseudotime_de[order(my_pseudotime_de[,4]),],file="my_pseudotime_de1.txt",sep="\t")

####draw####
cds <- readRDS("fiboMonocle-orderCells.rds")
de <- read.table("my_pseudotime_de1.txt")
de <- de[which(de$qval <= 0.01),]
c <- de$gene_short_name[1:150]

plot_genes_branched_heatmap(cds[c,], num_clusters = 4,use_gene_short_name = TRUE,show_rownames = TRUE,cores=5)
plot_cell_trajectory(cds, color_by = "Pseudotime",show_branch_points=F)
plot_cell_trajectory(cds, color_by = "Anno",show_branch_points=F)+scale_color_manual(values=c("#519259","#3E7C17","#064635","#F0BB62","#A3423C"))
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, markers=c("ACTA2","FAP","PDGFRA","PDGFRB","MCAM","CAV1","WT1"), cell_size=0.8, use_color_gradient=T,show_branch_points=F,show_tree = F) + scale_color_gradient2(low="gray",mid="yellow",high="red")
plot_genes_in_pseudotime(cds[c("FAP","PDGFRA","PDGFRB")],color_by="Anno")
