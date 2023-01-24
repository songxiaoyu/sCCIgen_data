rm(list=ls())
library(raster)
library(spatstat)
library(tidyverse)
library(rlist)
#library(scDesign2)
library(parallel, MASS, pscl)

# (1.3.1) --> No spatial data
# We need expression data

# Rpackage: /Users/anna/Projects/STsimulator/RPackage/R
source("R/PointSimulator_NoData.R") # This one to Anna
source("R/scDesign2_revised.R")
source("R/ExprSimulator.R")


# PointSimulator Functions for 2 regions---------
  # input scRNAseq  data:
  load("../DataToShare/scRNAseq_GTEx_breast.RData")

  # expr --> gene by cell  matrix of counts --> from datasets
  # anno --> colnames of expr
  anno=colnames(expr)
  # CopulaEst --> optional; estimated Gaussian copulas from expr data 
  
  # put in parameters --
  seed=124 # "Seed (to ensure reproducible simulation):" 
  nRegion=2  # No. of regions (suggested: 1-10)
  N=1000    # No. of cells 
  
  # more spatial parameters --
  # "No. of cells in each cell type in each region" 
    # Give default as below; allow users to chagne"
  cell.prop=vector("list", nRegion); 
  for (r in 1:nRegion) {cell.prop[[r]]=table(anno)/length(anno)} 
 
  
  # "Cell inhibition/attraction of cells from appearing its neighborhood"
    # create a L rows by 3 columns by nRegion array 
    # default NULL; Allow users to choose cell type 1, cell type 2, strength, and add rows. 
    # below is an example of 2 by 3 by 2 array
  cell.inh.attr.input=vector("list", nRegion); 
  cell.inh.attr.input[[1]]=data.frame(C1=c("Adipocyte", "Epithelial cell"), 
                              C2=c("Endothelial cell","Epithelial cell"), 
                               S=c(2, -1)) # interactive input by cell type 
  # 
  same.dis.cutoff =0.05
  even.distribution.coef=0.1
  
  # more expression parameters
  Copula=CopulaEst # Users do not need to choose; make Copula= NULL if an external data is used
  depth_simu_ref_ratio=1 # "Ratio of simulated/reference sequencing depth"; default=1
  
  # more parameters every time for adding a spatial pattern
  
  # EX1: Add.Spatial.Expr.Pattern
  # r=1
  # CellType="Adipocyte"
  # GeneID=NULL
  # PropOfGenes=0.1
  # delta.mean=1 
  # delta.sd=0.001
  
  # EX1: Add.Spatial.Expr.Pattern
  # r=1
  # CellType="Adipocyte"
  # GeneID="VHL"
  # PropOfGenes=NULL
  # delta.mean=1 
  # delta.sd=0.001  
  
  #  EX2: Add.Distance.Asso.Pattern
  # r=2,
  # perturbed.cell.type="Epithelial cell",
  # adjacent.cell.type="Immune (myeloid)",
  # int.dist.threshold=0.1,
  # delta.mean=1,
  # delta.sd=0.001,
  # GeneID=NULL, 
  # PropOfGenes=0.1
  
  # EX3: Add.Expr.Asso.Pattern
  # r=1,
  # perturbed.cell.type="Fibroblast",
  # adjacent.cell.type="Endothelial cell",
  # int.dist.threshold=0.1,
  # delta.mean=1,
  # delta.sd=0.001,
  # GenePairIDMatrix=NULL,
  # PropOfGenes=0.1,
  # Bidirectional=T
  


  # Step 1: Generate regions ----
    # Input: seed, nRegion
    # Output: window, area
  win=RandomRegionWindow(nRegion=nRegion, seed=seed)
  
  # Step 2: Generate cells in each region (repeat for `nRegion` times) ----
  # Input: N, win, cell.prop, cell.inh.attr.input, same.dis.cutoff, even.distribution.coef, seed
  # Output: cell.loc

  # R1 ---> make it inside of a function

  cell.loc=cell.loc.fc(N=N, win=win, cell.prop=cell.prop, 
                       cell.inh.attr.input=cell.inh.attr.input,
                         same.dis.cutoff =same.dis.cutoff,
                         even.distribution.coef=even.distribution.coef,
                         seed=seed)
 
  
  
  # # Step 3: Generate expression profiles for cells in each region (repeat it for `nRegion` times) ----
  # # Input: expr, anno, Copula (optional; if use our data, Copula is pre-estimated), 
  # #        cell.loc, depth_simu_ref_ratio
  # # Output: a list of length `nRegion`; each element is a G by N_r count matrix 
  # 
  # 
  # 
  # sim.count=Use_scDesign2(ppp.obj=cell.loc, 
  #                       expr=expr, anno=anno, 
  #                       Copula=CopulaEst, 
  #                       depth_simu_ref_ratio=depth_simu_ref_ratio,
  #                       seed=seed)
  # 
  # 
  # 
  # # Step 4: Add spatial patterns (repeatedly)
  # # Input: expr, r, CellType (interactive select), GeneID (upload) or PropOfGenes, 
  # #               delta.mean (averaged effect size), delta.sd (sd of effect size)
  # # Output: G by N.R1 beta.matrix
  # 
  # pattern1=Add.Spatial.Expr.Pattern(sim.count, r=1,
  #                          CellType="Adipocyte",
  #                          GeneID=NULL,
  #                          PropOfGenes=0.1,
  #                          delta.mean=1,
  #                          delta.sd=0.001, seed=seed)
  # 
  # 
  # # Step 5: Add distance to gene interactions
  # 
  # pattern2= Add.Distance.Asso.Pattern(ppp.obj=cell.loc,
  #                                   sim.count=sim.count, r=2,
  #                                   perturbed.cell.type="Epithelial cell",
  #                                   adjacent.cell.type="Immune (myeloid)",
  #                                   int.dist.threshold=0.1,
  #                                   delta.mean=1,
  #                                   delta.sd=0.001,
  #                                   GeneID=NULL, 
  #                                   PropOfGenes=0.1,
  #                                   seed=seed)
  # # Step 6: Add expr to gene interactions
  # pattern3= Add.Expr.Asso.Pattern(ppp.obj=cell.loc,
  #                           sim.count=sim.count,  r=1,
  #                            perturbed.cell.type="Fibroblast",
  #                            adjacent.cell.type="Endothelial cell",
  #                            int.dist.threshold=0.1,
  #                            delta.mean=1,
  #                            delta.sd=0.001,
  #                            GenePairIDMatrix=NULL,
  #                            PropOfGenes=0.1,
  #                            Bidirectional=T,
  #                            seed=seed)
  # 
  # # Step 7: Adjust total counts for Steps 4-6
  # sim.count.update=Pattern.Adj(sim.count, pattern.list=list(pattern1, pattern2, pattern3),
  #           bond.extreme=T, keep.total.count=T,
  #           integer=T)
  # 
  # # Step 5: combine points in different regions
  # output=MergeRegion(points.list=cell.loc, expr.list=sim.count.update)
  # 
