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

input="ParameterFile/fig2d_seqfish.tsv"
ParaSimulation(input=input, parallel=F)




FileName="fig2d_seqfish"
i=1
expr=fread(paste0("OutputData/", FileName, "_count_", i, ".tsv")) %>% as.data.frame %>%
  column_to_rownames("GeneName")
cell_feature=as.data.frame(fread(paste0("OutputData/", FileName, "_meta_",i,".tsv"))) 
dat=multicell(expr, cell_feature, NoSpot=500)

spot_expr=dat$count
spot_feature=dat$spot_feature
spot_feature[1:3,]
spot_sum=apply(spot_feature[, -c(1:2)], 1, sum)

spot_prop=spot_feature[, -c(1:2)]/spot_sum

# decomposition method -- CARD
library(CARD)

spot_expr=as(spot_expr, "dgCMatrix")
spot_loc=spot_feature[, 1:2]
expr=as(as.matrix(expr), "dgCMatrix")
cell_meta=cell_feature %>% 
  dplyr::select(c(1,2)) %>% 
  mutate(sampleInfo="sample1")
rownames(cell_meta)=cell_meta$Cell
CARD_obj = createCARDObject(
  sc_count = expr,
  sc_meta = cell_meta,
  spatial_count = spot_expr,
  spatial_location = spot_loc,
  ct.varname = "annotation",
  ct.select = unique(cell_meta$annotation),
  sample.varname = "sampleInfo",
  minCountGene = 0,
  minCountSpot = 0) 

CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

est_prop=CARD_obj@Proportion_CARD

plot(spot_prop, est_prop)


