# Take input

# ----------------- ParaDigest ---------------
ParaDigest=function(input) {
  # digest parameters
  para1=para=read_tsv(input) %>% column_to_rownames("parameters") %>%
    t() %>% as.data.frame()
  suppressWarnings(class(para1) <-"numeric")
  para[,is.na(para1)==F]=para1[is.na(para1)==F]

  attach(para)


  # clean seeds
  if (num_simulated_datasets>1) {
    para$all_seeds=as.numeric(unlist(strsplit(simulation_seed_for_each_dataset, ",")))
  } else{para$all_seeds=simulation_seed_for_each_dataset}


  # Path1: No ST data;
  # Path2: ST data - new cells
  # Path3: ST data - existing cells
  feature=CellFeatureLoad(para)
  para$path=ifelse(
    ncol(feature)==1, 1,
    ifelse(simulate_spatial_data=="TRUE", 2, 3)
  )

  detach(para)
  return(para)
}

# ----------- expr load ---------------
#' Load RData regardless of the name
#'
#' @param fileName  File name
#' @return Data
#' @export
#'
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
ExprLoad=function(para){
  if(expression_data_file_type=="Rdata" | expression_data_file_type=="RData") {
    expr=as.data.frame(loadRData(paste0(path_to_input_dir, expression_data_file)))
  }
  if (expression_data_file_type=="tsv") {
    expr=as.data.frame(fread(paste0(path_to_input_dir, expression_data_file)))
  }
  expr=as.matrix(expr)
  return(expr)
}
CellFeatureLoad=function(para){
  type=tail(unlist(strsplit(cell_feature_data_file, "[.]")), 1)
  if(type=="Rdata" | type=="RData" ) {
    feature=as.data.frame(loadRData(paste0(path_to_input_dir, cell_feature_data_file)))
  }
  if (type=="tsv") {
    feature=as.data.frame(fread(paste0(path_to_input_dir, cell_feature_data_file)))
  }
  return(feature)
}

# ----------- copula load/create ---------------
ParameterCopula=function(para, expr, feature, ncores=1){
  # Copula -- add region info
  if (gene_cor=="FALSE") {copula_input="NULL"; CopulaEst=NULL}
  if (gene_cor=="TRUE" & copula_input!="NULL") {CopulaEst=loadRData(copula_input)}
  if (gene_cor=="TRUE" & copula_input=="NULL") {

    expr=as.matrix(expr)
    anno=feature[,1]
    colnames(expr)=anno
    #L=length(unique(anno))
    if (region_specific_model=="FALSE") {
      CopulaEst=list(Est_GeneCopula(expr=expr,
                                    anno=anno, zp_cutoff=0.8, ncores=ncores))
    }
    if (region_specific_model=="TRUE") {
      Ridx=unique(feature[,4])
      R=length(Ridx)
      for (r in 1:R) {
        CopulaEst[[r]]=Est_GeneCopula(expr=expr,
                                      anno=anno, zp_cutoff=0.8, ncores=ncores)
        l=length( CopulaEst[[r]])

      }
      names(CopulaEst)=Ridx
    }
   # save
   out_path_name=paste0(path_to_output_dir, "_Copula.RData")
   save(CopulaEst, file=out_path_name)
   print("Finish estimating gene-gene correlation in expression data")
  }
  return(CopulaEst)
}

# ----------------- ParaNoSTCells ---------------
ParaCellsNoST=function(para, all_seeds){

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
  cell_type_proportion2=lapply(cell_type_proportion, function(f) f/sum(f))

  # determine cell cell location interactions in each region
   cell_location_interactions=vector("list", num_regions);

   if (custom_cell_location_interactions!="FALSE") {
     tmp1=unlist(strsplit(custom_cell_location_interactions, split=";"))
     tmp2=strsplit(tmp1, split = ",")
     tmp3=as.data.frame(matrix(unlist(tmp2),ncol=3,byrow=T))
     class(tmp3$V3)="numeric"
     for ( r in 1:num_regions) {cell_location_interactions[[r]]=tmp3}
   }



  # parallel starts here:
   cell_loc=foreach (m = 1: num_simulated_datasets) %dopar% {

     seed=all_seeds[m]
     win=RandomRegionWindow(nRegion=num_regions, seed=seed)
     cell.loc.fc(N=num_simulated_cells,
                               win=win,
                               cell.prop=cell_type_proportion2,
                               cell.inh.attr.input=cell_location_interactions,
                               same.dis.cutoff =cell_overlap_cutoff,
                               even.distribution.coef=cell_even_distribution,
                               seed=seed)

   }
   return(cell_loc)
}
# ----------------- ParaSTNewCells ---------------
ParaCellsST=function(para, feature, all_seeds) {

  if (ncol(feature)==4) {R=feature[,4]} else {R=rep(1, nrow(feature))}

  cell_loc=foreach (i = 1:num_simulated_datasets) %dopar% {
    cell.region.loc.model.fc(n=num_simulated_cells,
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

# ----------------- ParaPattern ---------------
ParaPattern=function(para, sim_count, cell_loc_list_i,
                     seed=NULL){

  # check No. of adding spatial patterns
  t1=length(grep("spatial_pattern_", colnames(para)))/6
  t2=length(grep("spatial_int_dist_", colnames(para)))/8
  t3=length(grep("spatial_int_expr_", colnames(para)))/10
  t0=sum(t1, t2, t3)
  beta.all=vector("list", num_simulated_datasets)
  for (i in 1:num_simulated_datasets) {

    # add spatial
    if (t0>0) {
      beta.all[[i]]=vector("list", t0)
    }
    for (tt1 in t1:1) {
      if (tt1==0) {break}
      #para[grep("spatial_pattern_", colnames(para))]
      r=eval(parse(text=paste0("spatial_pattern_",
                               tt1, "_region")))
      CellType=eval(parse(text=paste0("spatial_pattern_",
                                      tt1, "_cell_type")))
      GeneID1=eval(parse(text=paste0("spatial_pattern_",
                                     tt1, "_gene_id")))
      if (GeneID1=="NULL") {
        GeneID=eval(parse(text=GeneID1))
      } else {GeneID=unlist(strsplit(GeneID1, ","))}

      PropOfGenes=eval(parse(text=paste0("spatial_pattern_",
                                         tt1, "_gene_prop")))
      if (PropOfGenes=="NULL") {PropOfGenes=eval(parse(text=PropOfGenes))}

      delta.mean=eval(parse(text=paste0("spatial_pattern_",
                                        tt1, "_mean")))
      delta.sd=eval(parse(text=paste0("spatial_pattern_",
                                      tt1, "_sd")))
      beta.all[[i]][[tt1]]=Add.Spatial.Expr.Pattern(sim.count = sim_count,
                                    r=r,
                                    CellType=CellType,
                                    GeneID=GeneID,
                                    PropOfGenes=PropOfGenes,
                                    delta.mean=delta.mean,
                                    delta.sd=delta.sd,
                                    seed=seed)
    }

    # add expr-distance interaction
    for (tt1 in t2:1) {
      if (tt1==0) {break}
      #para[grep("spatial_int_dist_", colnames(para))]
      r=eval(parse(text=paste0("spatial_int_dist_",
                               tt1, "_region")))
      if (r=="NULL") {r=1}
      perturbed.cell.type=eval(parse(text=paste0("spatial_int_dist_",
                                      tt1, "_cell_type_perturbed")))
      adjacent.cell.type=eval(parse(text=paste0("spatial_int_dist_",
                                                 tt1, "_cell_type_adj")))
      int.dist.threshold=eval(parse(text=paste0("spatial_int_dist_",
                                                tt1, "_dist_cutoff")))
      GeneID1=eval(parse(text=paste0("spatial_int_dist_",
                                     tt1, "_gene_id1")))
      if (GeneID1=="NULL") {
        GeneID=eval(parse(text=GeneID1))
      } else {GeneID=unlist(strsplit(GeneID1, ","))}

      PropOfGenes=eval(parse(text=paste0("spatial_int_dist_",
                                         tt1, "_gene_prop")))
      if (PropOfGenes=="NULL") {PropOfGenes=eval(parse(text=PropOfGenes))}

      delta.mean=eval(parse(text=paste0("spatial_int_dist_",
                                        tt1, "_mean")))
      delta.sd=eval(parse(text=paste0("spatial_int_dist_",
                                      tt1, "_sd")))

      beta.all[[i]][[(t1+tt1)]]=Add.Distance.Asso.Pattern(ppp.obj=cell_loc_list_i,
                                       sim.count=sim_count,
                                       r=r,
                                       perturbed.cell.type=perturbed.cell.type,
                                       adjacent.cell.type=adjacent.cell.type,
                                       int.dist.threshold=int.dist.threshold,
                                       delta.mean=delta.mean,
                                       delta.sd=delta.sd,
                                       GeneID=GeneID, # Cell A Gene 1--> Cell B
                                       PropOfGenes=PropOfGenes,
                                       seed=seed)
    }
    # add expr-distance interaction
    for (tt1 in t3:1) {
      if (tt1==0) {break}
      #para[grep("spatial_int_expr_", colnames(para))]
      r=eval(parse(text=paste0("spatial_int_expr_",
                               tt1, "_region")))
      perturbed.cell.type=eval(parse(text=paste0("spatial_int_expr_",
                                                 tt1, "_cell_type_perturbed")))
      adjacent.cell.type=eval(parse(text=paste0("spatial_int_expr_",
                                                tt1, "_cell_type_adj")))
      int.dist.threshold=eval(parse(text=paste0("spatial_int_expr_",
                                                tt1, "_dist_cutoff")))

      GeneID1=eval(parse(text=paste0("spatial_int_expr_",
                                     tt1, "_gene_id1")))
      if (GeneID1=="NULL") {
        GeneID=eval(parse(text=GeneID1))
      } else {GeneID=unlist(strsplit(GeneID1, ","))}
      GeneID2=eval(parse(text=paste0("spatial_int_expr_",
                                     tt1, "_gene_id2")))
      if (GeneID2=="NULL") {
        GeneIDp=eval(parse(text=GeneID2))
      } else {GeneIDp=unlist(strsplit(GeneID2, ","))}

      if (is.null(GeneID)) {
        GenePairIDMatrix=NULL
      } else {GenePairIDMatrix=cbind(GeneID, GeneIDp)}

      PropOfGenes=eval(parse(text=paste0("spatial_int_expr_",
                                         tt1, "_gene_prop")))
      if (PropOfGenes=="NULL") {PropOfGenes=eval(parse(text=PropOfGenes))}

      Bidirectional1=eval(parse(text=paste0("spatial_int_expr_",
                                            tt1, "_bidirectional")))
      Bidirectional=eval(parse(text=Bidirectional1))

      delta.mean=eval(parse(text=paste0("spatial_int_expr_",
                                        tt1, "_mean")))
      delta.sd=eval(parse(text=paste0("spatial_int_expr_",
                                      tt1, "_sd")))

      beta.all[[i]][[(t1+t2+tt1)]]=Add.Expr.Asso.Pattern(ppp.obj=cell_loc_list_i,
                                                    sim.count=sim_count,
                                                    r=r,
                                                    perturbed.cell.type=perturbed.cell.type,
                                                    adjacent.cell.type=adjacent.cell.type,
                                                    int.dist.threshold=int.dist.threshold,
                                                    delta.mean=delta.mean,
                                                    delta.sd=delta.sd,
                                                    GenePairIDMatrix=GenePairIDMatrix,
                                                    PropOfGenes=PropOfGenes,
                                                    Bidirectional=Bidirectional,
                                                    seed=seed)

    }
  }
    return(beta.all)

}


# ----------------- ParaExpr ---------------
ParaFitExpr=function(para, expr, feature, CopulaEst, ncores){
  sim_method=ifelse(gene_cor=="TRUE", "copula", "ind")
  # fit by input data

   model_params=Use_scDesign2_model_params(expr=expr,
                                      feature=feature,
                                      Copula=CopulaEst,
                                      sim_method = sim_method,
                                      region_specific_model=region_specific_model,
                                      ncores=ncores)
    return(model_params)
}




ParaExpr=function(para, cell_loc_list, expr, feature,
                  CopulaEst, all_seeds, ncores=1){

  sim_method=ifelse(gene_cor=="TRUE", "copula", "ind")
  # fit by input data
  model_params=ParaFitExpr(para, expr, feature, CopulaEst, ncores=ncores)

  # simulate

  for (i in 1:num_simulated_datasets) {

    sim_count=Use_scDesign2(ppp.obj=cell_loc_list[[i]],
                            model_params=model_params,
                            expr=expr,
                            feature=feature,
                            depth_simu_ref_ratio=expr_depth_ratio,
                            sim_method=sim_method,
                            region_specific_model=region_specific_model,
                            seed=all_seeds[[i]])


    pattern_list=ParaPattern(para=para, sim_count=sim_count,
                             cell_loc_list_i=cell_loc_list[[i]],
                         seed=all_seeds[[i]])


    sim_count_update=Pattern.Adj(sim.count=sim_count,
                                 pattern.list=pattern_list[[i]],
                                 bond.extreme=T, keep.total.count=T,
                                 integer=T)

    output=MergeRegion(points.list=cell_loc_list[[i]],
                       expr.list=sim_count_update)

    expr_pattern=ExprPattern(pattern.list.i=pattern_list[[i]]) %>% as.data.frame()

    print(paste("Finished simulating data", i))
    # save

    save_name=paste0(path_to_output_dir, output_name)
    if (is.null(expr_pattern)==F) {
      write_tsv(expr_pattern, file=paste0(save_name, "_expr_pattern_", i, ".tsv"))
    }

    # multicell?
    if (num_spots=="NULL") {

      write_tsv(output$meta,
                file=paste0(save_name, "_meta_", i, ".tsv"))
      write_tsv(as.data.frame(output$count)%>% rownames_to_column("GeneName"),
                file=paste0(save_name, "_count_", i, ".tsv"))
    } else{
      output2=multicell(expr=output$count, cell_feature=output$meta, NoSpot=num_spots)

      write_tsv(output2$spot_feature,
                file=paste0(save_name, "_meta_", i, ".tsv"))
      write_tsv(as.data.frame(output2$count)%>% rownames_to_column("GeneName"),
                file=paste0(save_name, "_count_", i, ".tsv"))
    }


    print(paste("Finished saving data", i))

  }
}


# ----------------- ParaSimulation ---------------
ParaSimulation=function(input) {

  # parallel
  ncores=detectCores()-2; registerDoParallel(ncores)
  print(paste("No. of Cores in Use:", ncores))
  print("Start the simulation")

  # Digest parameters
  para=ParaDigest(input)
  attach(para)

  # Load  data
  expr=ExprLoad(para)
  feature=CellFeatureLoad(para)
  colnames(expr)=feature[,1]
  print("Finished loading data")

  # Simulate cells
  if (path==1) {
    cell_loc_list=ParaCellsNoST(para=para, all_seeds=all_seeds)
  }
  if (path==2) {
    cell_loc_list=ParaCellsST(para=para, feature=feature, all_seeds=all_seeds)
  }
  if (path==3) {
    cell_loc_list=ParaExistingCellsST(m=num_simulated_datasets, feature=feature)
  }
  print("Finished simulating the cell spatial maps")


  # Trim  data (no. of cells) in expr to make the expr estimation faster
  ctype=table(feature[,1])
  idxc=foreach (i = 1:length(ctype), .combine = "c") %dopar% {
    if (ctype[i]>2500) {
      sample(which(feature[,1]==names(ctype)[i]), 2500)
    } else {which(feature[,1]==names(ctype)[i])
    }
  }
  expr2=expr[,idxc]

  # Copula
  CopulaEst=ParameterCopula(para=para, expr=expr2, feature=feature, ncores=ncores)

  # Simulate Expr for these cells
  ParaExpr(para=para,
           cell_loc_list=cell_loc_list,
           expr=expr2, feature=feature,
           CopulaEst=CopulaEst, all_seeds=all_seeds,
           ncores=ncores)
  print("Finished simulating the expression of cells")
  detach(para)
  print("Finished the simulation")

}











