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

# ---- data generation 

input="ParameterFile/fig2c1.tsv"; para=ParaDigest(input);  
attach(para); p1=ParaCellsNoST(para=para, all_seeds=
                                 simulation_seed_for_each_dataset)[[1]][[1]]; detach(para)

input="ParameterFile/fig2c2.tsv"; para=ParaDigest(input);  
attach(para); p2=ParaCellsNoST(para=para, all_seeds=
                                 simulation_seed_for_each_dataset)[[1]][[1]]; detach(para)

input="ParameterFile/fig2c3.tsv"; para=ParaDigest(input);  
attach(para); p3=ParaCellsNoST(para=para, all_seeds=
                                 simulation_seed_for_each_dataset)[[1]][[1]]; detach(para)

input="ParameterFile/fig2c4.tsv"; para=ParaDigest(input);  
attach(para); p4=ParaCellsNoST(para=para, all_seeds=
                                 simulation_seed_for_each_dataset)[[1]][[1]]; detach(para)

input="ParameterFile/fig2c5.tsv"; para=ParaDigest(input);  
attach(para); p5=ParaCellsNoST(para=para, all_seeds=
                                 simulation_seed_for_each_dataset)[[1]][[1]]; detach(para)

# ParaSimulation(input=input)


# ---- summary plot 
F1=Fest(p1[which(marks(p1)=="Endothelial cell")], correction="km")
plot(F1)
F2=Fest(p2[which(marks(p2)=="Endothelial cell")], correction="km")
plot(F2)
F3=Fest(p3[which(marks(p3)=="Endothelial cell")], correction="km")
plot(F3)
F4=Fest(p4[which(marks(p4)=="Endothelial cell")], correction="km")
plot(F4)
F5=Fest(p5[which(marks(p5)=="Endothelial cell")], correction="km")
plot(F5)

dat=data.frame(rbind(cbind("Large Inhibition Effect", F1$r, F1$km),
                     cbind("Small Inhibition Effect", F2$r, F2$km),
                     cbind("No Effect", F3$r, F3$km),
                     cbind("Small Attraction Effect", F4$r, F4$km),
                     cbind("Large Attraction Effect", F5$r, F5$km)
                     ))

class(dat[,2])="numeric";class(dat[,3])="numeric"
colnames(dat)=c("Group", "Radius", "Fest" )
dat$Group=factor(dat$Group, levels=c("Large Inhibition Effect", 
                                        "Small Inhibition Effect",
                                        "No Effect",
                                        "Small Attraction Effect",
                                        "Large Attraction Effect"))
ggplot(dat, aes(x=Radius, y=Fest, color=Group)) + 
  geom_point() + geom_line() + xlim(0, 0.1) + xlab("Distance from Point") +
  ylab("Estiamte of Empty Space Function (Fest)")+ theme_classic()+ 
  theme(legend.position = c(0.8, 0.25), legend.title=element_blank())

a1=split(p1)[[2]]; dat1=data.frame(x=a1$x, y=a1$y)
ggplot(dat1, aes(x=x, y=y))+geom_point()
par(mfrow=c(1,3))
plot(split(p1)[[2]])
plot(split(p3)[[2]])
plot(split(p5)[[2]])


# FileName="fig2c"
# i=1
# spatial_feature=as.data.frame(fread(paste0("OutputData/", FileName, "_meta_",i,".tsv"))) 

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
# W=simu.window(PointLoc=spatial_feature[,c("x.loc", "y.loc")], method="rectangle")
# p=as.ppp(spatial_feature[,c("x.loc", "y.loc")], W=W)
# marks(p)=as.factor(spatial_feature$annotation)
# 
# Kest
K1=Kest(p[which(marks(p)=="Fibroblast")], correction="isotropic")
plot(K1)
plot(K1, cbind(r, sqrt(iso/pi)) ~ r)

K1=Kest(p1[which(marks(p1)=="Endothelial cell")], correction="isotropic")
plot(K1)
plot(K1, cbind(r, sqrt(iso/pi)) ~ r)

K1=Kest(p1[which(marks(p1)=="Immune (myeloid)")], correction="isotropic")
plot(K1)
plot(K1, cbind(r, sqrt(iso/pi)) ~ r)



# Fest
F1=Fest(p[which(marks(p)=="Fibroblast")], correction="rs")
plot(F1)
F1=Fest(p1[which(marks(p1)=="Endothelial cell")], correction="km")
plot(F1)
F2=Fest(p2[which(marks(p2)=="Endothelial cell")], correction="km")
plot(F2)
F3=Fest(p3[which(marks(p3)=="Endothelial cell")], correction="km")
plot(F3)
dat=data.frame(rbind(cbind("Attraction (effect = 0)", F1$r, F1$km),
                     cbind("Attraction (effect = 2) ", F2$r, F2$km),
                     cbind("Attraction (effect = 4) ", F3$r, F3$km)))
class(dat[,2])="numeric";class(dat[,3])="numeric"
colnames(dat)=c("Group", "Radius", "Fest" )
ggplot(dat, aes(x=Radius, y=Fest, color=Group)) + geom_point() + geom_line() + xlim(0, 0.1)



# plot(F1, cbind(km, trans, border) ~ theo)
F1=Fest(p1[which(marks(p1)=="Immune (myeloid)")], correction="km")
plot(F1)
F2=Fest(p2[which(marks(p2)=="Immune (myeloid)")], correction="km")
plot(F2)
F3=Fest(p3[which(marks(p3)=="Immune (myeloid)")], correction="km")
plot(F3)

dat=data.frame(rbind(cbind("Inhibition (effect = 0)", F1$r, F1$km),
 cbind("Inhibition (effect = 2) ", F2$r, F2$km),
 cbind("Inhibition (effect = 4) ", F3$r, F3$km)))
class(dat[,2])="numeric";class(dat[,3])="numeric"
colnames(dat)=c("Group", "Radius", "Fest" )
ggplot(dat, aes(x=Radius, y=Fest, color=Group)) + geom_point() + geom_line() + xlim(0, 0.1)


# Gest

G1=Gest(p[which(marks(p)=="Fibroblast")])
plot(G1) # null
G1=Gest(p[which(marks(p)=="Endothelial cell")])
plot(G1) # attraction



G1=Gest(p1[which(marks(p1)=="Immune (myeloid)")])
plot(G1) # inhibition
G2=Gest(p2[which(marks(p2)=="Immune (myeloid)")])
plot(G2) # inhibition
G3=Gest(p3[which(marks(p3)=="Immune (myeloid)")])
plot(G3) # inhibition

#Jest
J1=Jest(p[which(marks(p)=="Fibroblast")])
plot(J1) # null
J1=Jest(p1[which(marks(p1)=="Endothelial cell")])
plot(J1) # attraction
J1=Jest(p2[which(marks(p2)=="Endothelial cell")])
plot(J1)
J1=Jest(p3[which(marks(p3)=="Endothelial cell")])
plot(J1)
J1=Jest(p4[which(marks(p4)=="Endothelial cell")])
plot(J1)

J1=Jest(p5[which(marks(p5)=="Endothelial cell")])
plot(J1)


J1=Jest(p1[which(marks(p1)=="Immune (myeloid)")])
plot(J1) # inhibition



# Kmulti

K <- Kmulti(p, marks(p) == "Epithelial cell", 
            marks(p) =="Adipocyte", correction="isotropic")
plot(K)
K <- Kmulti(p, marks(p) == "Adipocyte", marks(p) =="Fibroblast", 
            correction="isotropic")
plot(K)
