# Identify spatial variable genes on simulated data
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

# detach(para)
input="ParameterFile/fig2a1.tsv"
input="ParameterFile/fig2a2.tsv"
input="ParameterFile/fig2a3.tsv"
ParaSimulation(input=input, parallel=F)

# Identify spatial variable genes on simulated data
# binSpect


library(boost)

Parameter_binSpect=function(FileName="fig2a1", NoSim=20) {

  pmatrix=NULL
  for (i in 1:NoSim) {
    # load data
    expr=fread(paste0("OutputData/", FileName, "_count_", i, ".tsv"))
    cell_feature=as.data.frame(fread(paste0("OutputData/", FileName, "_meta_",i,".tsv"))) 
    # clean
    expr2=expr%>% column_to_rownames("GeneName") %>% t() 
    PointLoc=cell_feature[, c("x.loc", "y.loc")] %>% as.matrix()
    # binSpect
    A <- get.neighbors(PointLoc, 4, method = "KNN")
    
    GeneID=expr$GeneName
    pvector=NULL
    for (j in 1:length(GeneID)) {
      g <- binarize.st(expr2, GeneID[j], cluster.method = "GMC")
      res <- binSpect(g, A, do.fisher.test = FALSE, gene.name =GeneID[j])
      p=res$measures$p.val
      pvector=cbind(pvector, p)
    }
    pmatrix=rbind(pmatrix, pvector)
  }
  colnames(pmatrix)=GeneID
  return(pmatrix=pmatrix)
}


pmatrix0=Parameter_binSpect(FileName="fig2a1", NoSim=20)
mean(pmatrix0[,1:2]<0.05)
mean(pmatrix0[,-c(1:2)]<0.05)


pmatrix1=Parameter_binSpect(FileName="fig2a2", NoSim=20)
mean(pmatrix1[,1:2]<0.05)
mean(pmatrix1[,-c(1:2)]<0.05)

pmatrix2=Parameter_binSpect(FileName="fig2a3", NoSim=20)
mean(pmatrix2[,1:2]<0.05)
mean(pmatrix2[,-c(1:2)]<0.05)

library(ggplot2)

a=rbind(cbind(0, "Power", mean(pmatrix0[,1:2]<0.05)),
        cbind(0, "Type I Error", mean(pmatrix0[,-c(1:2)]<0.05)),
        cbind(0.5, "Power", mean(pmatrix1[,1:2]<0.05)),
        cbind(0.5, "Type I Error", mean(pmatrix1[,-c(1:2)]<0.05)),
        cbind( (-0.5), "Power", mean(pmatrix2[,1:2]<0.05)),
        cbind((-0.5), "Type I Error", mean(pmatrix2[,-c(1:2)]<0.05))) %>% 
  as.data.frame() %>%
  type.convert(as.is =TRUE) 

ggplot(a, aes(x=V1, y=V3, col=V2)) + 
  geom_point()+ geom_line()+ 
  facet_grid(~V2) + 
  theme_classic() +
  labs(x="Effect Size", y = " ")+
  theme(legend.position = "none")
  












# slow ----- trendsceek
library(trendsceek)

# load data
expr=fread("OutputData/fig2a1_count_1.tsv")
cell_feature=as.data.frame(fread("OutputData/fig2a1_meta_1.tsv"))


PointLoc=cell_feature[, c("x.loc", "y.loc")]
cell_win=simu.window(PointLoc=PointLoc, 
                     method="convex5")
p=as.ppp(PointLoc, W=cell_win)
marks(p)=data.frame(CellType1=ifelse(cell_feature$annotation=="CellType1", 1, 0), 
                    Region=cell_feature$region)

##run trendsceek
set.seed(153)
idx=sample(1:10, p$n, replace=T)
p1=split(p, f=as.factor(idx))
trendstat_list = trendsceek_test(p1[[1]], nrand = 10, ncores = 1)
head(trendstat_list[['supstats_wide']])


##show significant genes
sig_list = extract_sig_genes(trendstat_list, alpha = 0.1)
sig_genes = sig_list[['markcorr']][, 'gene']
print(sig_genes)
plot_trendstats(trendstat_list, sig_genes)
pp_sig = pp_select(pp, sig_genes)
plot_pp_scatter(pp_sig, log_marks = FALSE, scale_marks = TRUE, pal.direction = -1)

##cells located in high-expressing regions of the significant genes
cellpeaks_siggenes = cellsceek_test(pp_sig)
sig_cells = get_sigcells(cellpeaks_siggenes)
plot_pp_density(pp_sig, log_marks = FALSE, cells2highlight = sig_cells)	 

