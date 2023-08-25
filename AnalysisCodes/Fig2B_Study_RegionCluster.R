# Identify spatial regions on simulated data
#rm(list=ls())
library(tidyverse)
library(data.table)
library(raster)
library(spatstat)
library(rlist)
library(parallel)

# library(STsimulator)

setwd("/Users/songxiaoyu152/Dropbox/SpatialTranscriptomics/Paper_Simulator/Github")
source("RPackage/R/PointSimulator_NoData.R")
source("RPackage/R/PointSimulator_STData.R")
source("RPackage/R/ExprSimulator.R")
source("RPackage/R/scDesign2_fit_revised.R")
source("RPackage/R/scDesign2_simulate_revised.R")
source("RPackage/R/ParameterDigest.R")
source("RPackage/R/MultiCell.R")
# detach(para)

input="ParameterFile/fig2b.tsv"
ParaSimulation(input=input)


# Get true region
para=ParaDigest(input)
win=RandomRegionWindow(nRegion=para$num_regions, seed=para$parent_simulation_seed)
plot(win$window[[1]], col="pink")
plot(win$window[[2]], col="blue", add=T)




FileName="fig2b"
i=1
expr=fread(paste0("OutputData/", FileName, "_count_", i, ".tsv")) %>% as.data.frame %>%
  column_to_rownames("GeneName")
spot_feature=as.data.frame(fread(paste0("OutputData/", FileName, "_meta_",i,".tsv"))) 


dp2=data.frame(spot_feature, t(expr))
ggplot(dp2, aes(x=x, y=y, col=log(NOC2L+1))) + geom_point()













# BayesSpace ----------- 

## ----setup--------------------------------------------------------------------
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)

# NA remove
n_gene_NA <- apply(expr, 1, function(x) sum(is.na(x)))
expr <- expr[n_gene_NA == 0, ]
dim(expr)


spot_feature$Spot = factor(colnames(expr), levels = colnames(expr)) 

SCE <- SingleCellExperiment(
  assays = list(counts = as(as.matrix(expr), "dgCMatrix")), # wrap in a list
  rowData = rownames(expr),
  colData = spot_feature,
  metadata = list(name = "SCE")
)


# colData(SCE)$col <- spot_feature %>% 
#   arrange(Spot) %>% 
#   pull(x)
# # y coordinate
# colData(SCE)$row <- spot_feature %>% 
#   arrange(Spot) %>% 
#   pull(y)

SCE2 <- reducedDim(SCE, "PCA")

SCE <- spatialPreprocess(sce=SCE, platform="ST", n.PCs = 15,
                         n.HVGs=1000, skip.PCA=F,
                         log.normalize=T, assay.type="counts")
# Tuning the choice of q (number of clusters) before running spatialCluster
# SCE <- qTune(SCE, qs=seq(2, 3), platform="ST")


SCE <- spatialCluster(SCE, 
                      q=3, # no of clusters
                      d=15, # no of pc
                      platform="ST",
                      init.method="mclust", 
                      model="normal", 
                      gamma=2, # ST uses 2; Visium uses 3
                      nrep=2000, burn.in=1000,
                      save.chain=F)

dat=as.matrix(sparseMatrix(i = spot_feature$col, 
                           j= spot_feature$row, x= colData(SCE)$spatial.cluster))
dat2=as.matrix(sparseMatrix(i = spot_feature$col, 
                           j= spot_feature$row, x= spot_feature$region))
heatmap(dat2, Rowv = NA, Colv =NA)
heatmap(dat, Rowv = NA, Colv =NA)
View(dat)

#  plot ARI by effect size

# ----- SpaceFlow -----------  python package
#  Identifying multicellular spatiotemporal organization of cells with SpaceFlow






## ----Seurat-----------------------------------------------------
# Seurat
library(Seurat)
library(patchwork)


obj <- CreateSeuratObject(counts = dat$count, project = "fig2b", 
                           min.cells = 0, min.features = 0)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

obj <- FindVariableFeatures(obj, selection.method = "dispersion", nfeatures = 10)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs=5)
obj <- FindNeighbors(obj, dims = 1:2)
obj <- FindClusters(obj, resolution = 0.001)
a=obj$seurat_clusters

dp=data.frame(dat$spot_loc, a)
ggplot(dp, aes(x=x, y=y, col=a)) + geom_point()



