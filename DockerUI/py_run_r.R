# Create Example Datasets for Users to Download
rm(list=ls())
library(data.table)
library(raster)
library(spatstat)
library(rlist)
library(parallel)
library(doParallel)
library(proxy)
library(dplyr)
library(tibble)
library(readr)
library(STsimulator)

args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
ParaSimulation(input = input)