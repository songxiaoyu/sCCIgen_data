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
ParaSimulation(input=input, parallel=F)


# Get true region
para=ParaDigest(input)
win=RandomRegionWindow(nRegion=para$num_regions, seed=567)
plot(win$window[[1]], col="pink")
plot(win$window[[2]], col="blue", add=T)




FileName="fig2b"
i=1
expr=fread(paste0("OutputData/", FileName, "_count_", i, ".tsv")) %>% as.data.frame %>%
  column_to_rownames("GeneName")
cell_feature=as.data.frame(fread(paste0("OutputData/", FileName, "_meta_",i,".tsv"))) 
dat=multicell(expr, cell_feature, NoSpot=500)

spot_expr=dat$count
spot_feature=dat$spot_feature



dp2=data.frame(dat$spot_feature, t(dat$count))

ggplot(dp2, aes(x=x, y=y, col=log(Gene1+1))) + geom_point()


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






# BayesSpace
## ----setup--------------------------------------------------------------------
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)


a=rbind(dat$count, dat$count, dat$count, dat$count)


SCE <- SingleCellExperiment(
  assays = list(counts = as(a, "dgCMatrix")), # wrap in a list
  rowData = rownames(a),
  colData = dat$spot_loc,
  metadata = list(name = "SCE")
)

SCE <- spatialPreprocess(SCE, platform="ST", n.PCs = 3,
                          n.HVGs=10, 
                         log.normalize=T, assay.type="counts")
SCE <- qTune(SCE, qs=seq(2, 3))


SCE <- spatialCluster(SCE, q=2, platform="ST", d=3,
                           init.method="mclust", model="t", gamma=2,
                           nrep=1000, burn.in=100,
                           save.chain=TRUE)
clusterPlot(SCE)


## ----readVisium, eval=FALSE---------------------------------------------------
#  sce <- readVisium("path/to/spaceranger/outs/")

## ----download-----------------------------------------------------------------
melanoma <- getRDS(dataset="2018_thrane_melanoma", sample="ST_mel1_rep2")


## ----preprocess---------------------------------------------------------------
set.seed(102)
melanoma <- spatialPreprocess(melanoma, platform="ST", 
                              n.PCs=7, n.HVGs=2000, log.normalize=FALSE)

## ----tuning_q-----------------------------------------------------------------
melanoma <- qTune(melanoma, qs=seq(2, 10), platform="ST", d=7)
qPlot(melanoma)

## ----cluster------------------------------------------------------------------
set.seed(149)
melanoma <- spatialCluster(melanoma, q=4, platform="ST", d=7,
                           init.method="mclust", model="t", gamma=2,
                           nrep=1000, burn.in=100,
                           save.chain=TRUE)

## ----cluster.results----------------------------------------------------------
head(colData(melanoma))

## ----cluster.plot, fig.width=7, fig.height=5----------------------------------
clusterPlot(melanoma)

## ----cluster.plot.customize, fig.width=7, fig.height=5------------------------
clusterPlot(melanoma, palette=c("purple", "red", "blue", "yellow"), color="black") +
  theme_bw() +
  xlab("Column") +
  ylab("Row") +
  labs(fill="BayesSpace\ncluster", title="Spatial clustering of ST_mel1_rep2")

## ----enhance, eval=TRUE-------------------------------------------------------
melanoma.enhanced <- spatialEnhance(melanoma, q=4, platform="ST", d=7,
                                    model="t", gamma=2,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=1000, burn.in=100,
                                    save.chain=TRUE)

## ----enhance.results----------------------------------------------------------
head(colData(melanoma.enhanced))

## ----enhance.plot, eval=TRUE, fig.width=7, fig.height=5-----------------------
clusterPlot(melanoma.enhanced)

## ----enhanceFeatures----------------------------------------------------------
markers <- c("PMEL", "CD2", "CD19", "COL1A1")
melanoma.enhanced <- enhanceFeatures(melanoma.enhanced, melanoma,
                                     feature_names=markers,
                                     nrounds=0)

## ----enhanced.logcount--------------------------------------------------------
logcounts(melanoma.enhanced)[markers, 1:5]

## ----enhanced.rmse------------------------------------------------------------
rowData(melanoma.enhanced)[markers, ]

## ----enhanced.featurePlot-----------------------------------------------------
featurePlot(melanoma.enhanced, "PMEL")

## ----enhanced.markers, fig.width=12, fig.height=8-----------------------------
enhanced.plots <- purrr::map(markers, function(x) featurePlot(melanoma.enhanced, x))
patchwork::wrap_plots(enhanced.plots, ncol=2)

## ----compare.resolution, fig.width=16, fig.height=8---------------------------
spot.plots <- purrr::map(markers, function(x) featurePlot(melanoma, x))
patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol=4)

## ----mcmcChain, eval=TRUE-----------------------------------------------------
chain <- mcmcChain(melanoma)
chain[1:5, 1:5]




