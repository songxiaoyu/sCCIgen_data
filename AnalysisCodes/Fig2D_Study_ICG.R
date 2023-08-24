# Identify ICG on simulated data
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
input="ParameterFile/fig2a1.tsv"
input="ParameterFile/fig2a2.tsv"
input="ParameterFile/fig2a3.tsv"
ParaSimulation(input=input, parallel=F)


# Giotto


FileName="fig2a1"
i=1
expr=fread(paste0("OutputData/", FileName, "_count_", i, ".tsv")) %>% as.data.frame %>%
  column_to_rownames("GeneName")
cell_feature=as.data.frame(fread(paste0("OutputData/", FileName, "_meta_",i,".tsv"))) 


# cell-cell distance cutoff 0.02
dist=pairdist(cell_feature[, 3:4])
dist2=dist[which(cell_feature$annotation=="CellType1"), 
     which(cell_feature$annotation=="CellType2")]
rownames(dist2)=cell_feature$Cell[which(cell_feature$annotation=="CellType1")]
colnames(dist2)=cell_feature$Cell[which(cell_feature$annotation=="CellType2")]

dlong=dist2 %>% as.data.frame() %>% rownames_to_column("CellType1") %>%
   pivot_longer(-CellType1, names_to="CellType2", values_to="Dist") %>%
  mutate(Group=ifelse(Dist<0.02, 1, 0)) %>%
  dplyr::select(CellType1, Group) %>%
  group_by(CellType1) %>%
  summarise(max=max(Group))
idx1=dlong$CellType1[which(dlong$max==1)]
idx0=dlong$CellType1[which(dlong$max==0)]
p1=NULL
for (i in 1:10) {
   tt=t.test(expr[i, idx1], expr[i, idx0])
   p1=c(p1,tt$p.value)
 }



