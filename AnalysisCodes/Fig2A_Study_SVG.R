# Identify spatial variable genes on simulated data
#rm(list=ls())
library(tidyverse)
library(data.table)
library(raster)
library(spatstat)
library(rlist)
library(parallel)
library(doParallel)

# library(STsimulator)

setwd("/Users/songxiaoyu152/Dropbox/SpatialTranscriptomics/Paper_Simulator/Github")
source("RPackage/R/PointSimulator_NoData.R")
source("RPackage/R/PointSimulator_STData.R")
source("RPackage/R/ExprSimulator.R")
source("RPackage/R/scDesign2_fit_revised.R")
source("RPackage/R/scDesign2_simulate_revised.R")
source("RPackage/R/ParameterDigest.R")
source("RPackage/R/MultiCell.R")

# use normal breast snRNAseq with two regions and different gene means in two regions. 
input="ParameterFile/fig2a0.tsv" # ormal breast snRNAseq with two regions with random 10% gene with high means in region 1. 
# input="ParameterFile/fig2a2.tsv"
# input="ParameterFile/fig2a3.tsv"
ParaSimulation(input=input)

# Identify spatial variable genes on simulated data
## two sample t-test and FDR control

Parameter_SVG=function(FileName="fig2a0"){
  count=fread("OutputData/fig2a0_count_1.tsv") 
  expr=count %>% column_to_rownames("GeneName")
  meta=fread("OutputData/fig2a0_meta_1.tsv")
  pattern=fread("OutputData/fig2a0_expr_pattern_1.tsv")
  
  gidx=which(count$GeneName %in% pattern$GeneID)
  gidx2=setdiff(1:nrow(expr), gidx)
  cidx=which(meta$region ==1)
  p_in=unlist(mclapply(gidx[1:10], function(f) t.test(expr[f,cidx], expr[f,-cidx])$p.value))
  p_out=unlist(mclapply(gidx2[1:10], function(f) t.test(expr[f,cidx], expr[f,-cidx])$p.value) ) 
  power=mean(p_in<0.05)
  type1=mean(p_out<0.05)
  summary=data.frame(type1, power)
  return(list(summary=summary, p_in=p_in,p_out= p_out))
}


pmatrix0=Parameter_SVG(FileName="fig2a0")$summary


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

