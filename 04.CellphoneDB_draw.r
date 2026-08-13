library(RColorBrewer)
library(scales)
library(igraph)
library(stringr)
library(dplyr)
library(ggplot2)

setwd("CellphoneDB/out/")
pvalues=read.table("pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05))
colnames(statdf)=c("number")

rownames(statdf)=gsub("\\.\\.","_",rownames(statdf))
rownames(statdf)=gsub("F\\.","F-",rownames(statdf))
rownames(statdf)=gsub("EC\\.","EC-",rownames(statdf))
rownames(statdf)=gsub("A\\.","A-",rownames(statdf))
rownames(statdf)=gsub("al\\.","al-",rownames(statdf))

statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")
rankname=sort(unique(statdf$indexa))
len=length(rankname)
A=c()
B=c()
C=c()
remaining=rankname
for (i in rankname[-len]) {
  remaining=setdiff(remaining,i)
  for (j in remaining) {
    count=statdf[statdf$indexa == i & statdf$indexb == j,"number"]+
      statdf[statdf$indexb == i & statdf$indexa == j,"number"]
    A=append(A,i)
    B=append(B,j)
    C=append(C,count)
  }
}

statdf2=data.frame(indexa=A,indexb=B,number=C)
statdf2=statdf2 %>% rbind(statdf[statdf$indexa==statdf$indexb,c("indexa","indexb","number")])
statdf2=statdf2[statdf2$number > 0,] 

statdraw <- statdf2[which(statdf2$indexa == "CAF-2" & statdf2$indexb %in% c("EC-Mix1","EC-Mix2","EC-Mix3","EC-Mix4","EC-Mix5","EC-PS1","EC-PS2","EC-PS3","EC-PS4","EC-PS5","EC-PS6","EC-PS7","EC-PS8","EC-PS9","EC-PS10")),]


color1=c("#651fff","#b39ddb","#9575cd","#7e57c2","#673ab7","#f8bbd0","#f48fb1","#f06292","#ec407a","#880e4f","#ce93d8","#ba68c8","#ab47bc","#e040fb","#493240",
"#8d6e63",
"#F4EEA9","#519259","#3E7C17","#064635","#F0BB62","#A3423C",
"#DFF6FF","#90E0EF","#00B4D8","#344CB7","#000957",
               "#39A388","#1C7947",
               "#F9E0AE","#FFDF91","#FFCE45","#D4AC2B","#B05E27","#7E370C",
               "#C38B8B"
)

names(color1)=c("EC-Mix1","EC-Mix2","EC-Mix3","EC-Mix4","EC-Mix5","EC-PS1","EC-PS2","EC-PS3","EC-PS4","EC-PS5","EC-PS6","EC-PS7","EC-PS8","EC-PS9","EC-PS10",
"Endothelial cell",
"Fibroblast","CAF-1","CAF-2","CAF-3","Mesenchymal-fib","Mesothelial-cell",
"TC_CD4_FOXP3","TC_CD4_CXCL13","TC_CD4_S100A6","TC_CD8_IL7R","TC_CD8_IL4R",
                      "BC_CD38","BC_MS4A1",
                      "M_MRC1","M_SPP1","M_HLA-DRB5","M_MKI67","M_IGKC","NK_GNLY",
                      "DC")

color2=colorRampPalette(brewer.pal(9, "Reds")[3:7])(20) 
names(color2)=1:20 

net <- graph_from_data_frame(statdraw[,c("indexa","indexb","number")])
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
group <-  igraph::cluster_optimal(net)
coords <- igraph::layout_in_circle(net, order = order(membership(group)))

E(net)$width <- E(net)$number / 16 #将数值映射到连线的宽度，有时还需要微调，这里除以2就是这个目的
E(net)$color <- color2[as.character(ifelse(E(net)$number > 20,20,E(net)$number))] #用前面设置好的颜色赋给连线，颜色深浅对应数值大小
E(net)$label = E(net)$number #连线的标注
E(net)$label.color <- "black" #连线标注的颜色
V(net)$label.color <- "black" #节点标注的颜色
V(net)$color <- color1[names(V(net))] #节点的填充颜色，前面已经设置了；V(net)返回节点信息

pdf("interaction-caf2.pdf",height=10,width=10)
plot(net,
     edge.arrow.size = 0, #连线不带箭头
     edge.curved = 0, #连线不弯曲
     vertex.frame.color = "black", #节点外框颜色
     layout = coords,
     vertex.label.cex = .7, #节点标注字体大小
     vertex.size = 20) #节点大小
dev.off()


