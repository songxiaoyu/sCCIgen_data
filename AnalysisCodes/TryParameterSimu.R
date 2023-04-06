# Take input
rm(list=ls())
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

# No ST simple
# input="ParameterFile/fake1_input_simple.tsv"
# input="ParameterFile/fake1_input_expand.tsv"
# input="ParameterFile/fake1_input_1region.tsv"
# input="ParameterFile/fake1_input_1region_expand.tsv"

# ST simple existing cells
# input="ParameterFile/fake2_input_existing_cells_simple.tsv"
# input="ParameterFile/fake2_input_new_cells_simple.tsv"
# input="ParameterFile/fake3_input_existing_cells_simple.tsv"
# input="ParameterFile/fake3_input_new_cells_simple.tsv"

ParaSimulation(input=input, parallel=F)
