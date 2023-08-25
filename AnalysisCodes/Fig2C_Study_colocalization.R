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

input="ParameterFile/fig2c.tsv"
ParaSimulation(input=input)


FileName="fig2c"
i=1

spatial_feature=as.data.frame(fread(paste0("OutputData/", FileName, "_meta_",i,".tsv"))) 

# -- How to quantify cell-cell attraction and inhibition -- 
  # Kest, Fest, Gest, Jest
  # Gcross, Gdoc, Gmulti, Kcross Kdot, Kmulti, Jross, Jdot, Jmulti

# 1 Kest: Ripley's K-function - a  statistic summarising aspects of inter-point 
 # “dependence” and “clustering”. 
# 2 Fest: Estimate the Empty Space Function or its Hazard Rate - a  statistic 
 # summarising the sizes of gaps in the pattern.
# 3 Gest: Nearest Neighbour Distance Function G.

# 4 Jest: J = (1-Gest)/(1-Fest). J = 1 under randome point process. Deviations J(r) < 1J(r)<1 
 # or J(r) > 1J(r)>1 typically indicate spatial clustering or spatial regularity, respectively. 
 # Deviations between the empirical and theoretical KK curves may suggest spatial clustering 
 # or spatial regularity.

# Use `allstats` to estimates of all four functions F, G, J, K



# ---- ppp
W=simu.window(PointLoc=spatial_feature[,c("x.loc", "y.loc")], method="rectangle")
p=as.ppp(spatial_feature[,c("x.loc", "y.loc")], W=W)
marks(p)=as.factor(spatial_feature$annotation)

# Kest
K1=Kest(p[which(marks(p)=="Endothelial cell")], correction="isotropic")
plot(K1)
plot(K1, cbind(r, sqrt(iso/pi)) ~ r)

K1=Kest(p[which(marks(p)=="Immune (myeloid)")], correction="isotropic")
plot(K1)
plot(K1, cbind(r, sqrt(iso/pi)) ~ r)

K1=Kest(p[which(marks(p)=="Fibroblast")], correction="isotropic")
plot(K1)
plot(K1, cbind(r, sqrt(iso/pi)) ~ r)

# Fest

F1=Fest(p[which(marks(p)=="Endothelial cell")])
plot(F1)
plot(F1, cbind(km, trans, border) ~ theo)

F1=Fest(p[which(marks(p)=="Fibroblast")], correction="rs")
plot(F1)
# Gest
G1=Gest(p[which(marks(p)=="Endothelial cell")])
plot(G1)

#Jest
J1=Jest(p[which(marks(p)=="Endothelial cell")])
plot(J1)

J1=Jest(p[which(marks(p)=="Muscle")])
plot(J1)

# Kmulti

K <- Kmulti(p, marks(p) == "Epithelial cell", 
            marks(p) =="Adipocyte", correction="isotropic")
plot(K)
K <- Kmulti(p, marks(p) == "Muscle", marks(p) =="Fibroblast", 
            correction="isotropic")
plot(K)

# Jmulti
J=Jmulti(p, marks(p) == "Epithelial cell", 
         marks(p) =="Adipocyte", correction="rs")
plot(J)

J=Jmulti(p, marks(p) == "Muscle", marks(p) =="Fibroblast", correction="rs")
plot(J)




# Moran i 
## <https://stats.oarc.ucla.edu/r/faq/how-can-i-calculate-morans-i-in-r/>
## <https://mgimond.github.io/simple_moransI_example/>
library(ape)
# calculate inverse distance weights
dists <- as.matrix(dist(cbind(spatial_feature$x, spatial_feature$y)))
dists.inv <- 1/dists
diag(dists.inv) <- 0
Moran.I(, dists.inv)
