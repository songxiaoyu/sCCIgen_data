# Create Example Datasets for Users to Download
#rm(list=ls())
# library(tidyverse)
# library(data.table)
# library(raster)
# library(spatstat)
# library(rlist)
# library(parallel)
# library(doParallel)
# library(proxy)
library(STsimulator)

setwd("/Users/songxiaoyu152/Dropbox/SpatialTranscriptomics/Paper_Simulator/Github")
# source("RPackage/R/PointSimulator_NoData.R")
# source("RPackage/R/PointSimulator_STData.R")
# source("RPackage/R/ExprSimulator.R")
# source("RPackage/R/scDesign2_fit_revised.R")
# source("RPackage/R/scDesign2_simulate_revised.R")
# source("RPackage/R/ParameterDigest.R")
# source("RPackage/R/MultiCell.R")

# detach(para)
input="ParameterFile/example1.tsv"
ParaSimulation(input=input)

input="ParameterFile/example2.tsv"
ParaSimulation(input=input)

input="ParameterFile/example3.tsv"
ParaSimulation(input=input)









