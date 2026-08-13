library(Seurat)
library(dplyr)
library(SCENIC)
library(AUCell)
library(RcisTarget)


mydbs <- c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
           "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")


scenicOptions <- initializeScenic(org = "hgnc", 
                                  nCores = 4,
                                  dbDir = mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "HNC")
                                  
#saveRDS(scenicOptions, "hg38-scenicOptions.rds")

ds <- readRDS("CR_All_harmony_cluster.rds")

exprMat  <-  as.matrix(ds@assays$RNA@data)
cellInfo <- ds@meta.data[,c("nCount_RNA","nFeature_RNA","seurat_clusters")]

genesKept <- geneFiltering(exprMat, scenicOptions,minCountsPerGene = 0.015 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

saveRDS(scenicOptions, file="ScenicOptions.Rds") 