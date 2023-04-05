# Take input
#install.packages("rgeos")
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
input="ParameterFile/expr_input_simple2.tsv"
# ST simple existing cells
input="ParameterFile/MERFISH_input_existing_cells_simple.tsv"
# ST simple new cells
input="ParameterFile/ST_input_new_cells_simple2.tsv"
# No ST  - expanding expr patterns
input="ParameterFile/expr_input_expand_expr2.tsv"


# ----------------- ParaDigest ---------------
ParaDigest=function(input) {
  # digest parameters
  para1=para=read_tsv(input) %>% column_to_rownames("parameters") %>%
    t() %>% as.data.frame()
  suppressWarnings(class(para1) <-"numeric")
  para[,is.na(para1)==F]=para1[is.na(para1)==F]
  attach(para)

  # Path1: No ST data;
  # Path2: ST data - new cells
  # Path3: ST data - existing cells
  feature=CellFeatureLoad(para)
  para$path=ifelse(
    ncol(feature)==1, 1,
    ifelse(simulate_spatial_data=="TRUE", 2, 3)
  )

  # estimate copula
  if (gene_cor==T & copula_input=="NULL") {
    #  copula_name
    if (save_copula) {
      out_path_name=paste0(path_to_output_dir, output_name, "_Copula.RData")
      ParameterCopula(para=para, feature=feature, out_path_name=out_path_name)
    }
  }

  detach(para)
  return(para)
}

ParameterCopula=function(para, feature, out_path_name=NULL){

  # read  data
  expr=ExprLoad(para)
  # anno
  anno=feature[,1]
  colnames(expr)=anno

  # estimate copula
  expr=as.matrix(expr)
  Copula=Est_GeneCopula(expr=expr, anno=anno, zp_cutoff=0.8, ncores=10)
  save(Copula, file=out_path_name)
  return(Copula)
}

# ----------- expr load ---------------
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
ExprLoad=function(para){
  if(expression_data_file_type=="Rdata") {
    expr=as.data.frame(loadRData(paste0(path_to_input_dir, expression_data_file)))
  }
  if (expression_data_file_type=="tsv") {
    expr=as.data.frame(fread(paste0(path_to_input_dir, expression_data_file)))
  }
  return(expr)
}
CellFeatureLoad=function(para){
  type=tail(unlist(strsplit(cell_feature_data_file, "[.]")), 1)
  if(type=="Rdata") {
    feature=as.data.frame(loadRData(paste0(path_to_input_dir, cell_feature_data_file)))
  }
  if (type=="tsv") {
    feature=as.data.frame(fread(paste0(path_to_input_dir, cell_feature_data_file)))
  }
  return(feature)
}
CopulaLoad=function(para){
  # Copula
  if (copula_input!="NULL") {
    CopulaEst=loadRData(copula_input)
  }
  if (copula_input=="NULL" & gene_cor=='TRUE'){
    out_path_name=paste0(path_to_output_dir, output_name, "_Copula.RData")
    CopulaEst=loadRData(out_path_name)
  }
  return(CopulaEst)
}


# ----------------- ParaNoSTCells ---------------
ParaCellsNoST=function(para, all_seeds, parallel=F){

  # determine cell type proportion in each region
  cell_type_proportion=vector("list", num_regions);
  tmp=as.matrix(para[,grep("cell_type_proportion_", colnames(para))])
  tmp2=matrix(unlist(strsplit(tmp, ",")), ncol=3, byrow = T)
  for (i in 1:num_regions){
    cell_type_proportion[[i]]= tmp2%>%
      data.frame() %>% filter(X1==i, X3>0) %>%
      column_to_rownames("X2") %>% dplyr::select(X3) %>%
      transform(X3=as.numeric(X3)) %>% t()
  }

  # determine cell cell location interactions in each region
   cell_location_interactions=vector("list", num_regions);

  # parallel starts here:
  cell_loc=vector("list", num_simulated_datasets)
  for (m in 1: num_simulated_datasets) {
    seed=all_seeds[m]
    win=RandomRegionWindow(nRegion=num_regions, seed=seed)
    cell_loc[[m]]=cell.loc.fc(N=num_simulated_cells,
                           win=win,
                           cell.prop=cell_type_proportion,
                           cell.inh.attr.input=cell_location_interactions,
                           same.dis.cutoff =cell_overlap_cutoff,
                           even.distribution.coef=cell_even_distribution,
                           seed=seed)
  }
   return(cell_loc)
}



# ----------------- ParaSTNewCells ---------------
ParaCellsST=function(para, feature, all_seeds, parallel=F) {
  cell_loc=vector("list", num_simulated_datasets)
  if (ncol(feature)==4) {R=feature[,4]} else {R=rep(1, nrow(feature))}

  for (i in 1:num_simulated_datasets) {
    cell_loc[[i]]=cell.region.loc.model.fc(n=num_simulated_cells,
                                    PointLoc=feature[,c(2:3)],
                                    PointAnno=feature[,1],
                                    PointRegion=R,
                                    window_method=window_method,
                                    seed=all_seeds[[i]])
  }
  return(cell_loc)
}
# ----------------- ParaSTExistingCells ---------------
ParaExistingCellsST=function(m, feature) {
  if (ncol(feature)==4) {R=feature[,4]} else {R=rep(1, nrow(feature))}
  cell_loc1=cell.region.loc.existing.fc(PointLoc=feature[,c(2:3)],
                                       PointAnno=feature[,1],
                                       PointRegion=R,
                                       window_method="rectangle")
  cell_loc=rep(list(cell_loc1), times=num_simulated_datasets)
  return(cell_loc)
}

# ----------------- ParaExpr ---------------
ParaExpr=function(para, cell_loc_list, expr, feature,
                  CopulaEst, all_seeds, parallel=F){
  sim_method=ifelse(gene_cor=="TRUE", "copula", "ind")
  anno=feature[,1]
  colnames(expr)=anno
  # check No. of adding spatial patterns
  t1=length(grep("spatial_pattern_", colnames(para)))/6
  t2=length(grep("spatial_int_dist_", colnames(para)))/8
  t3=length(grep("spatial_int_expr_", colnames(para)))/10

  for (i in 1:num_simulated_datasets) {
    sim_count=Use_scDesign2(ppp.obj=cell_loc_list[[i]],
                            expr=expr, feature=feature,
                            Copula=CopulaEst,
                            depth_simu_ref_ratio=expr_depth_ratio,
                            seed=all_seeds[[i]], sim_method=sim_method)
    # add spatial
    if (t1>0) {
      for (tt1 in 1:t1) {
        temp=Add.Spatial.Expr.Pattern(sim.count = sim_count,
                                 r,
                                 CellType,
                                 GeneID=NULL,
                                 PropOfGenes=0.1,
                                 delta.mean=1,
                                 delta.sd=0.01,
                                 seed)
      }
    }





    sim_count_update=Pattern.Adj(sim.count=sim_count, pattern.list=NULL,
                                 bond.extreme=T, keep.total.count=T,
                                 integer=T)
    output=MergeRegion(points.list=cell_loc_list[[i]], expr.list=sim_count_update)
    print(paste("Finish simulating expression data", i))


    write_tsv(output$meta, file=paste0(path_to_output_dir, output_name, "_meta_", i, ".tsv"))
    write_tsv(as.data.frame(output$count)%>% rownames_to_column("GeneName"),
              file=paste0(path_to_output_dir, output_name, "_count_", i, ".tsv"))
    return(output)
  }
}


# ----------------- ParaSimulation ---------------
ParaSimulation=function(input, parallel=F) {
  # Digest parameters
  para=ParaDigest(input)
  attach(para)

  # Load  data
  expr=ExprLoad(para)
  feature=CellFeatureLoad(para)
  if (gene_cor=="TRUE") {
    CopulaEst=CopulaLoad(para)
  } else {CopulaEst=NULL}



  # parallel parameters
  if (num_simulated_datasets>1) {
    all_seeds=as.numeric(unlist(strsplit(simulation_seed_for_each_dataset, ",")))
  } else{all_seeds=simulation_seed_for_each_dataset}


  # Simulate cells
  if (path==1) {
    cell_loc_list=ParaCellsNoST(para=para,
                                all_seeds=all_seeds, parallel=parallel)
  }
  if (path==2) {
    cell_loc_list=ParaCellsST(para=para, feature=feature,
                              all_seeds=all_seeds, parallel=parallel)
  }
  if (path==3) {
    cell_loc_list=ParaExistingCellsST(m=num_simulated_datasets,
                                      feature=feature)
  }

  # Simulate Expr for these cells
  cell_expr=ParaExpr(para=para, cell_loc_list=cell_loc_list,
           expr=expr, feature=feature, CopulaEst=CopulaEst, all_seeds=all_seeds,
           parallel=parallel)

  detach(para)

}


# aa -----------------------------------

# sub functions





