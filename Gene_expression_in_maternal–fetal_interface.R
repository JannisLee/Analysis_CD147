rm(list = ls())

#please download data from https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6701/E-MTAB-6701.processed.1.zip and https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6701/E-MTAB-6701.processed.1.zip
#then uncompress to get the files mentioned below

setwd('~/project/Single/')
exprSet = read.table('~/project/Single/raw_data_10x.txt',stringsAsFactors = F,header = T,
                    row.names = 1)
meta <- read.table('meta_10x.txt',header = T,sep = '\t',fill = T ,quote = "",
                   na.strings = "NA", comment.char = "",stringsAsFactors = F,check.names = F)
library(dplyr)
library(Seurat)
library(patchwork)
exprSet = as.data.frame(na.omit(exprSet[which(rowSums(exprSet) > 0),]))
fminter <- CreateSeuratObject(counts = exprSet, project = "fminter", 
                              min.cells = 3, min.features = 200)

#reorder the sample
fminter@meta.data$id <- rownames(fminter@meta.data)
meta$id <- rownames(meta)
meta <- meta[order(meta$id,decreasing = F),]
fminter@meta.data <- fminter@meta.data[order(fminter@meta.data$id,decreasing = F),]
fminter@meta.data$cluster <- meta$final_cluster
fminter@meta.data$annotation <- meta$annotation

#Normalization
fminter <- NormalizeData(fminter)

#Calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others)
fminter <- FindVariableFeatures(fminter, selection.method = "vst", nfeatures = 3000)

#Scaling the data
all.genes <- rownames(fminter)
fminter <- ScaleData(fminter, features = all.genes)

#PCA
fminter <- RunPCA(fminter)

#UMAP
fminter <- RunUMAP(fminter, dims = 1:40, reduction = "pca")

library(ggplot2)
png('Celltypes.png',width = 1600, height = 900)
DimPlot(fminter, reduction = "umap", label = T,group.by = 'annotation')
dev.off()

library(stringr)
oldname <- unlist(fminter@assays$RNA@counts@Dimnames[1])
newname <- str_sub(oldname,end = nchar(oldname)-16)
fminter@assays$RNA@counts@Dimnames[[1]] <- newname

png('BSG.png',width = 1600, height = 900)
FeaturePlot(fminter, features = c("BSG"),min.cutoff = 0,max.cutoff = 150,
            slot = 'counts',reduction = 'umap') 
dev.off()

png('ITGB1.png',width = 1600, height = 900)
FeaturePlot(fminter, features = c("ITGB1"),min.cutoff = 0,max.cutoff = 1000,
            slot = 'counts',reduction = 'umap') 
dev.off()

