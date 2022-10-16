## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(SPECK)
library(Seurat)
library(ggplot2)
library(gridExtra)

## ----message=FALSE, warning=FALSE---------------------------------------------
data("pbmc.rna.mat")
dim(pbmc.rna.mat)

## ----message=FALSE, warning=FALSE---------------------------------------------
speck.full <- speck(counts.matrix = pbmc.rna.mat, rank.range.end = 100, 
                    min.consec.diff = 0.01, rep.consec.diff = 2, 
                    manual.rank = NULL, max.num.clusters = 4, 
                    seed.rsvd = 1, seed.ckmeans = 2)

speck.rank <- speck.full$rrr.rank
paste("Rank: ", speck.rank, sep = "")

plot(speck.full$component.stdev, ylab = "Stdev. of non-centered sample PCs", 
xlab = "Rank range", main = paste("Selected rank (k=", speck.rank, ")", sep=""))
abline(v = speck.rank, lty = 2, col = "red")

head(speck.full$clust.num); table(speck.full$clust.num)
head(speck.full$clust.max.prop)

speck.output <- speck.full$thresholded.mat
paste("# of samples in RRR object:", dim(speck.output)[1])
paste("# of genes in RRR object:", dim(speck.output)[2])
SPECK_assay <- CreateAssayObject(counts = t(speck.output))
pbmc.rna.seurat <- CreateSeuratObject(counts = t(as.matrix(pbmc.rna.mat)))
pbmc.rna.seurat[["SPECK"]] <- SPECK_assay

## ---- message=FALSE, warning=FALSE--------------------------------------------
DefaultAssay(pbmc.rna.seurat) <- "RNA"
pbmc.rna.seurat <- NormalizeData(pbmc.rna.seurat)
pbmc.rna.seurat <- FindVariableFeatures(pbmc.rna.seurat, 
                                        selection.method = "vst", 
                                        nfeatures = 2000)
all.genes <- rownames(pbmc.rna.seurat)
pbmc.rna.seurat <- ScaleData(pbmc.rna.seurat, features = all.genes)
pbmc.rna.seurat <- RunPCA(pbmc.rna.seurat, 
                          features = VariableFeatures(object = pbmc.rna.seurat))
pbmc.rna.seurat <- FindNeighbors(pbmc.rna.seurat, dims = 1:10)
pbmc.rna.seurat <- FindClusters(pbmc.rna.seurat, resolution = 0.5)
pbmc.rna.seurat <- RunUMAP(pbmc.rna.seurat, dims = 1:10)

## ---- fig.width=7, fig.height=9, message=FALSE, warning=FALSE-----------------
DefaultAssay(pbmc.rna.seurat) <- "RNA"
p1 <- FeaturePlot(pbmc.rna.seurat, "CD14", cols = c("lightgrey", "#007ece")) + 
  ggtitle("CD14 RNA")
DefaultAssay(pbmc.rna.seurat) <- "SPECK"
p2 <- FeaturePlot(pbmc.rna.seurat, "CD14", cols=c("lightgrey", "#E64B35CC")) +
  ggtitle("CD14 SPECK")

DefaultAssay(pbmc.rna.seurat) <- "RNA"
p3 <- FeaturePlot(pbmc.rna.seurat, "CD79B", cols = c("lightgrey", "#007ece")) + 
  ggtitle("CD79B RNA")
DefaultAssay(pbmc.rna.seurat) <- "SPECK"
p4 <- FeaturePlot(pbmc.rna.seurat, "CD79B", cols=c("lightgrey", "#E64B35CC")) + 
  ggtitle("CD79B SPECK")

DefaultAssay(pbmc.rna.seurat) <- "RNA"
p5 <- FeaturePlot(pbmc.rna.seurat, "CD19", cols = c("lightgrey", "#007ece")) + 
  ggtitle("CD19 RNA")
DefaultAssay(pbmc.rna.seurat) <- "SPECK"
p6 <- FeaturePlot(pbmc.rna.seurat, "CD19", cols=c("lightgrey", "#E64B35CC")) + 
  ggtitle("CD19 SPECK")

grid.arrange(p2, p1,
             p4, p3, 
             p6, p5,
             nrow = 3)

