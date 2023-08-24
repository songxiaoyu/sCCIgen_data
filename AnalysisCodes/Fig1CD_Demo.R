# Identify spatial variable genes on simulated data
#rm(list=ls())
library(tidyverse)
library(data.table)
library(raster)
library(spatstat)
library(rlist)
library(parallel)
library(doParallel)
library(doRNG);
# library(STsimulator)

setwd("/Users/songxiaoyu152/Dropbox/SpatialTranscriptomics/Paper_Simulator/Github")
source("RPackage/R/PointSimulator_NoData.R")
source("RPackage/R/PointSimulator_STData.R")
source("RPackage/R/ExprSimulator.R")
source("RPackage/R/scDesign2_fit_revised.R")
source("RPackage/R/scDesign2_simulate_revised.R")
source("RPackage/R/ParameterDigest.R")
source("RPackage/R/MultiCell.R")

detach(para)
input="ParameterFile/fig1cd.tsv"
ParaSimulation(input=input)

# 

# Get true region
para=ParaDigest(input)
win=RandomRegionWindow(nRegion=para$num_regions, seed=123)

plot(win$window[[1]], col="pink")
plot(win$window[[2]], col="blue", add=T)
plot(win$window[[3]], col="orange", add=T)

  # region in ggplot
dat1=data.frame(x=win$window[[1]]$bdry[[1]]$x, y=win$window[[1]]$bdry[[1]]$y, r=1)
dat1=rbind(dat1, dat1[1,])
dat2=data.frame(x=win$window[[2]]$bdry[[1]]$x, y=win$window[[2]]$bdry[[1]]$y, r=2)
dat2=rbind(dat2, dat2[1,])
dat3=data.frame(x=win$window[[3]]$bdry[[1]]$x, y=win$window[[3]]$bdry[[1]]$y, r=3)
dat3=rbind(dat3, dat3[1,])
dat=rbind(dat1, dat2, dat3)
map <- ggplot() + 
  geom_path(aes(x = x, y = y), data = dat1)+ 
  geom_path(aes(x = x, y = y), data = dat2)+ 
  geom_path(aes(x = x, y = y), data = dat3) + 
  coord_equal()+theme_classic() +
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

map


# plot cell distribution by region 
FileName="fig1cd"
i=1
expr=fread(paste0("OutputData/", FileName, "_count_", i, ".tsv")) %>% as.data.frame %>%
  column_to_rownames("GeneName")
cell_feature=as.data.frame(fread(paste0("OutputData/", FileName, "_meta_",i,".tsv"))) 
gene_feature=as.data.frame(fread(paste0("OutputData/", FileName, "_expr_pattern_",i,".tsv")))


cell=map+geom_point(data=cell_feature, aes(x=x.loc, y=y.loc, color=annotation, shape=annotation))+
  scale_fill_viridis(discrete = TRUE)
cell
# ------------ plot spatial pattern ------------

gene_feature[1:3,]
# Type Region        CellType GeneID
# 1 SpatialChange      3 Epithelial cell   DOK3
# 2 SpatialChange      3 Epithelial cell   BIVM
# 3 SpatialChange      3 Epithelial cell  CCDC9
GeneID=gene_feature[which(gene_feature$Type=="SpatialChange"), "GeneID"] 
expd=data.frame(cell_feature, t(expr[rownames(expr) %in% GeneID,]))
expp=map+geom_point(data=expd[which(expd$annotation=="Epithelial cell"),], 
                    aes(x=x.loc, y=y.loc, color=BIVM))
expp
# ------------ plot expr-distance interaction ------------------
gene_feature[which(gene_feature$Type=="DistanceAssoGenes"), ] [1:3,]

GeneID=gene_feature[which(gene_feature$Type=="DistanceAssoGenes"), "GeneID"] 
expd=data.frame(cell_feature, t(expr[rownames(expr) %in% GeneID,]))

expp=map+
  geom_point(data=expd[which(expd$annotation=="Fibroblast"),], 
             aes(x=x.loc, y=y.loc, shape=annotation))+
  scale_shape_manual(values=c(2))+
  geom_point(data=expd[which(expd$annotation=="Epithelial cell"),], 
             aes(x=x.loc, y=y.loc, color=KLHL17, size=KLHL17*0.2))
expp
# ------------ plot expr-expr interaction ------------
gene_feature[which(gene_feature$Type=="DistanceAssoGenes"), ] [1:3,]

GeneID=gene_feature[which(gene_feature$Type=="DistanceAssoGenes"), "GeneID"] 
expd=data.frame(cell_feature, t(expr[rownames(expr) %in% GeneID,]))

expp=map+
  geom_point(data=expd[which(expd$annotation=="Fibroblast"),], 
             aes(x=x.loc, y=y.loc, shape=annotation))+
  scale_shape_manual(values=c(2))+
  geom_point(data=expd[which(expd$annotation=="Epithelial cell"),], 
             aes(x=x.loc, y=y.loc, color=KLHL17, size=1+(KLHL17>2)))
expp

