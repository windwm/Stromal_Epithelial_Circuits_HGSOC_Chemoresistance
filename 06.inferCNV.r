library(dplyr)
library(Seurat)
library(ggplot2)
library(infercnv)

sampleinfo <- read.csv("Rinfo.csv")
geneFile = "/ref/gencode.v32.primary_assembly.annotation_gene_pos.txt"

for (i in 1:4){
  samName = sampleinfo$samName[i]
  fileName = sampleinfo$fileName[i]
  print("Now we begin:")
  print(samName)
  expFile=paste0(samName,'_expFile.txt')
  groupFile=paste0(samName,'_groupFile.txt')
  fpath = paste0(samName,"/inferCNV/")
  setwd(fpath)
  
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFile,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=c("Non-malign"))

  infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, 
                             out_dir="./",
                             cluster_by_groups=FALSE,
                             denoise=T,
                             HMM=T)
}


##########Draw Result##########3
  
for (i in 1:nrow(sampleinfo)){
  samName = sampleinfo$samName[i]
  fileName = sampleinfo$fileName[i]
  print("Now we begin:")
  print(samName)
  fpath=paste0(samName,"/inferCNV_result/")
  setwd(fpath)
  
  fileaddr = paste0(samName,"/inferCNV/run.final.infercnv_obj")
  fileanno = paste0(samName,"_groupFile.txt")
  
  infercnv_obj = readRDS(fileaddr)
  expr <- infercnv_obj@expr.data
  
  anno.df= read.table(fileanno)
  anno.df$CB <- anno.df$V1
  anno.df$class <- anno.df$V2
  anno.df <- anno.df[,c("CB", "class")]

  gn <- rownames(expr)
  rownames(geneFile)=geneFile$V1
  geneFile_draw <- geneFile[which(geneFile$V2  %in% c(paste("chr",1:22,sep = ""))),]
  sub_geneFile <-  geneFile[intersect(gn,geneFile_draw$V1),]
  expr=expr[intersect(gn,geneFile$V1),]
  expr_draw =expr[intersect(gn,geneFile_draw$V1),]
  
  kmeans.result <- kmeans(t(expr), 7)
  kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
  kmeans_df$CB=rownames(kmeans_df)
  kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") 
  kmeans_df_s=arrange(kmeans_df,kmeans_class) 
  rownames(kmeans_df_s)=kmeans_df_s$CB
  kmeans_df_s$CB=NULL
  kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class)
  kmeans_df_s$class <- as.factor(kmeans_df_s$class)
  
  top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
  color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:7] #类别数
  names(color_v)=as.character(1:7)
  
  color_v2= colorlist[1:length(levels(kmeans_df_s$class))] 
  names(color_v2)=levels(kmeans_df_s$class) 
  left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=color_v2,kmeans_class=color_v))
  
  filename = paste0(samName,"_",fileName,"_heatmap.pdf")
  pdf(filename,width = 15,height = 10)
  ht = Heatmap(t(expr_draw)[rownames(kmeans_df_s),], #绘图数据的CB顺序和注释CB顺序保持一致
                  col = colorRamp2(c(0.5,1,1.5), c("#377EB8","#F0F0F0","#E41A1C")), #如果是10x的数据，这里的刻度会有所变化
                  cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
                  column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), #这一步可以控制染色体顺序，即使你的基因排序文件顺序是错的
                  column_gap = unit(2, "mm"),
                  heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
                  top_annotation = top_anno,left_annotation = left_anno, #添加注释
                  row_title = NULL,column_title = NULL)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  write.table(kmeans_df_s, file = paste0(samName,"_",fileName,"_kmeans_df_s.txt"), quote = FALSE, sep = '\t', row.names = T, col.names = T)
  df <- table(kmeans_df_s$kmeans_class,kmeans_df_s$class)
  write.csv(df, file = paste0(samName,"_",fileName,"_test.csv"))

  
  
  expr_2=expr-1
  expr_2=expr_2 ^ 2
  CNV_score=as.data.frame(colMeans(expr_2))
  colnames(CNV_score)="CNV_score"
  CNV_score$CB=rownames(CNV_score)
  kmeans_df_s$CB=rownames(kmeans_df_s)
  CNV_score=CNV_score%>%inner_join(kmeans_df_s,by="CB")
  write.csv(CNV_score, file= paste0(samName,"_score.csv"))
  
  p <- ggplot(CNV_score, aes(kmeans_class,CNV_score))+
  geom_violin(aes(fill=kmeans_class),color="NA")+
    scale_fill_manual(values = color_v)+
    theme_bw()
    
  ggsave(file=paste0(samName,"_",fileName,"_CNV_level.pdf"),plot=p,width = 10,height = 6,units = "cm")

}