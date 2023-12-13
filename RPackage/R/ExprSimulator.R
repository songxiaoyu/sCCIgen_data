

Use_scDesign2_1region=function(ppp.obj1, Genes, model_params,
                       depth_simu_ref_ratio=1, cell_type_sel, seed,
                       sim_method = c('copula', 'ind')) {
  # cell types in simulated and reference data
  n.ordered=table(ppp.obj1$marks)
  exist.cell.type=names(n.ordered)
  cell_type_prop=n.ordered/ppp.obj1$n
  model_params_exist=model_params[exist.cell.type]
  set.seed(seed)
  sim_count <- scDesign2.revised(model_params=model_params_exist,
                                 n_cell_new=ppp.obj1$n,
                                 cell_type_prop = cell_type_prop,
                                 depth_simu_ref_ratio=depth_simu_ref_ratio,
                                 sim_method =sim_method)

  # Update the order of sim_count to match the cell type of ppp.obj1
  sim_count2=matrix(NA, ncol=ncol(sim_count), nrow=nrow(sim_count))
  colnames(sim_count2)=ppp.obj1$marks
  rownames(sim_count2) = Genes
  for (f in levels(ppp.obj1$marks)) {
    sim_count2[, which(colnames(sim_count2)==f)]=
      sim_count[,which(colnames(sim_count)==f)]
  }

  return(sim.count=sim_count2)
}




Use_scDesign2_model_params=function(expr,
                       feature,
                       Copula=NULL,
                       sim_method = c('copula', 'ind'),
                       region_specific_model,
                       ncores=1) {

  expr=as.matrix(expr)
  cell_type_sel=names(table(colnames(expr)))
  Genes=rownames(expr)

  if (region_specific_model!="TRUE") { # not region specific
    model_params=  fit_model_scDesign2(data_mat=expr,
                                       cell_type_sel=cell_type_sel,
                                       sim_method = 'ind',
                                       marginal='zinb',
                                       ncores = ncores)
    if (sim_method=="copula") {
      CellType=names(model_params)
      for (i in 1:length(CellType)){model_params[[CellType[i]]]$cov_mat=Copula[[1]][[i]]}
    }
  }

  #  region specific model
  if (region_specific_model=="TRUE") { #  region specific

    Region=feature[,4]
    Runiq=unique(Region)
    R=length(Runiq)

    model_params= foreach (r = 1:R) %dopar% {
      idx=which(Region==Runiq[r])

      model_params1=fit_model_scDesign2(data_mat=expr[,idx],
                                         cell_type_sel=cell_type_sel,
                                         sim_method = 'ind',
                                         marginal='zinb',
                                         ncores = ncores)
      if (sim_method=="copula") {
        CellType=names(model_params1)
        for (i in 1:length(CellType)){model_params1[[CellType[i]]]$cov_mat=Copula[[r]][[i]]}
      }
      model_params1
    }

  }

  return(model_params)
}


#' Simulate expression under no spatial patterns using revised scDesign2 functions
#'
#' This function uses input parameters to simulate expression for all regions.
#' @param ppp.obj Cells as points on a spatial map for all regions.
#' @param model_params Provide model parameters including marginal distributions and copula (if not NULL).
#' @param expr Gene expression (count) in reference data.
#' @param feature Cell features (e.g. cell type, spatial coordinates, regions) of reference data.
#' @param depth_simu_ref_ratio Relative sequencing depth in comparison of reference data.
#' @param sim_method Simulate independent genes using'ind' or correlated genes using 'copula'.
#' @param region_specific_model Whether estimation model differ in different regions.
#' @param seed Seed
#' @return Simulated expression count data for all regions.
#' @export

Use_scDesign2=function(ppp.obj,
                       model_params,
                       expr,
                       feature,
                       depth_simu_ref_ratio=1,
                       sim_method = c('copula', 'ind'),
                       region_specific_model,
                       seed) {

  expr=as.matrix(expr)
  R=length(ppp.obj)
  cell_type_sel=names(table(colnames(expr)))
  Genes=rownames(expr)

  if (region_specific_model!="TRUE") { # not region specific

    sim.count= foreach (r = 1:R) %dopar%{
      Use_scDesign2_1region(ppp.obj1=ppp.obj[[r]],
                                           Genes=Genes,
                                           model_params=model_params,
                                           depth_simu_ref_ratio=depth_simu_ref_ratio,
                                           cell_type_sel=cell_type_sel,
                                           seed=seed*31+r*931,
                                           sim_method = sim_method)
    }
  }

  #  region specific model
  if (region_specific_model=="TRUE") { #  region specific

    Region=feature[,4]
    Runiq=unique(Region)

    sim.count= foreach (r = 1:R) %in% {

      Use_scDesign2_1region(ppp.obj1=ppp.obj[[r]],
                                           Genes=Genes,
                                           model_params=model_params[[r]],
                                           depth_simu_ref_ratio=depth_simu_ref_ratio,
                                           cell_type_sel=cell_type_sel,
                                           seed=seed*31+r*931,
                                           sim_method = sim_method)
    }

  }

  return(sim.count)
}

# Spatial patterns -----------------------------------------------


#' adds spatial differential expressed pattern for one cell type
#'
#' This function adds spatial differential expressed patterns for one cell type.
#' @param sim.count Cells as points on a spatial map.
#' @param r Region index.
#' @param CellType Cell type index.
#' @param GeneID Gene(s) index.
#' @param PropOfGenes Proportion of genes with this pattern if GeneID is not provided.
#' @param delta.mean Mean effect
#' @param delta.sd SD of effect.
#' @param seed Seed
#' @return Adding spatial differential expressed patterns for one cell type in one region.
#' @export
#'
Add.Spatial.Expr.Pattern= function(sim.count,
                                   r,
                                   CellType,
                                   GeneID=NULL,
                                   PropOfGenes=0.1,
                                   delta.mean=1,
                                   delta.sd=0.01, seed) {
  set.seed(938*seed-142)
  R=length(sim.count)
  sim.count1=sim.count[[r]]
  G=nrow(sim.count1)
  N=ncol(sim.count1)
  GeneAll=rownames(sim.count1)

  # key matrix
  beta.matrix=vector("list", length=R)
  for (i in 1:R) {beta.matrix[[i]]=matrix(0, nrow=G, ncol=ncol(sim.count[[i]]))}
  colnames(beta.matrix[[r]])=colnames(sim.count1)
  rownames(beta.matrix[[r]])=GeneAll

  # GeneID
  if (is.null(GeneID)) {GeneID=sample(GeneAll, round(PropOfGenes * G))}

  CellID=which(colnames(beta.matrix[[r]])== CellType )
  beta=rnorm(length(GeneID), delta.mean, delta.sd)
  SignalSummary=data.frame(Type="SpatialChange", Region=r, CellType, GeneID,
                           AdjCellType="NA",
                           AdjGene="NA", beta)
  beta.matrix[[r]][GeneID,CellID] = beta +  beta.matrix[[r]][GeneID,CellID]

  return(list(SignalSummary=SignalSummary, beta.matrix=beta.matrix))
}

# Find.Neighbor.Pairs

Find.Neighbor.Pairs=function(ppp.obj,
                             interacting.cell.type.pair,
                             int.dist.threshold) {
  cell.loc=cbind(ppp.obj$x, ppp.obj$y)
  cell1.idx=which(ppp.obj$marks==interacting.cell.type.pair[1])
  cell2.idx=which(ppp.obj$marks==interacting.cell.type.pair[2])
  m=spatstat.geom::crossdist(cell.loc[cell1.idx,1], cell.loc[cell1.idx,2],
            cell.loc[cell2.idx,1],cell.loc[cell2.idx,2])
  # in neighbor or not?
  dmax=max( max(ppp.obj$x)-min(ppp.obj$x), max(ppp.obj$y)-min(ppp.obj$y))

  neighbo.loc.idx=which(m< (int.dist.threshold*dmax), arr.ind = TRUE)

  # index in original data
  nbr.idx=cbind(cell1.idx[neighbo.loc.idx[,1]],   cell2.idx[neighbo.loc.idx[,2]])

  colnames(nbr.idx)=interacting.cell.type.pair
  return(neighbor.idx=nbr.idx)
}

#' Add a type of cell-cell interaction pattern (expr-distance) to a pair of cell types
#'
#' This function add a type of cell-cell interactions to a pair of cell types. The expression in a cell type associated with the proximity of
#' the other cell type. One can repeat this function for multiple times to
#' add cell-cell interactions for many cell type pairs and regions.
#' @param ppp.obj Spatial info for cell type 1 (e.g. neuron)
#' @param sim.count Spatial info for cell type 2 (e.g. microglia)
#' @param r Expression info for cell type 1 (e.g. neuron)
#' @param perturbed.cell.type Expression info for cell type 1 (e.g. microglia)
#' @param adjacent.cell.type Int.Cell.Pair.Idx estimated from function `Find.Neighbor.Pairs`
#' @param int.dist.threshold effect.size
#' @param delta.mean effect.size
#' @param delta.sd effect.size
#' @param GeneID effect.size
#' @param PropOfGenes effect.size
#' @param seed effect.size
#' @return
#' \item{SignalSummary:}{SignalSummary}
#' \item{beta.matrix:}{beta.matrix}
#' @export

Add.Distance.Asso.Pattern = function(ppp.obj,
                                   sim.count, r,
                                   perturbed.cell.type,
                                   adjacent.cell.type,
                                   int.dist.threshold=0.1,
                                   delta.mean=1,
                                   delta.sd=0.001,
                                   GeneID=NULL, # Cell A Gene 1--> Cell B
                                   PropOfGenes=NULL,
                                   seed=NULL) {


  set.seed(seed*478-50194)


  R=length(sim.count)
  sim.count1=sim.count[[r]]
  N=ncol(sim.count1)
  G=nrow(sim.count1)
  GeneAll=rownames(sim.count1)

  # key matrix
  beta.matrix=vector("list", length=R)
  for (i in 1:R) {beta.matrix[[i]]=matrix(0, nrow=G,  ncol=ncol(sim.count[[i]]))}
  colnames(beta.matrix[[r]])=colnames(sim.count1)
  rownames(beta.matrix[[r]])=GeneAll

  # spatial info
  nbr.idx=Find.Neighbor.Pairs(ppp.obj=ppp.obj[[r]],
                              interacting.cell.type.pair=c(perturbed.cell.type,
                                                           adjacent.cell.type),
                              int.dist.threshold=int.dist.threshold)

  # GeneID
  if (is.null(GeneID)) {GeneID=sample(GeneAll, round(PropOfGenes * G))}

  beta=rnorm(length(GeneID), delta.mean, delta.sd)
  beta.matrix[[r]][GeneID, nbr.idx[,1]] = beta +  beta.matrix[[r]][GeneID, nbr.idx[,1]]

  SignalSummary=data.frame(Type="DistanceAssoGenes", Region=r, CellType=perturbed.cell.type,
                           GeneID, AdjCellType=adjacent.cell.type,
                           AdjGene="NA",beta)


  return(list(SignalSummary=SignalSummary, beta.matrix=beta.matrix))
}

#' Add a type of cell-cell interaction pattern (expr-expr) to a pair of cell types
#'
#' This function add cell-cell interactions to a pair of cell types (e.g.
#' neuron-microglia) for expression in a cell type associated with expression of
#' the neighboring other cell type.
#' One can repeat this function for multiple times to add cell-cell interactions for many cell types.
#' @param ppp.obj Spatial info for cell type 1 (e.g. neuron)
#' @param sim.count Spatial info for cell type 2 (e.g. microglia)
#' @param r Expression info for cell type 1 (e.g. neuron)
#' @param perturbed.cell.type Expression info for cell type 1 (e.g. microglia)
#' @param adjacent.cell.type Int.Cell.Pair.Idx estimated from function `Find.Neighbor.Pairs`
#' @param int.dist.threshold effect.size
#' @param delta.mean effect.size
#' @param delta.sd effect.size
#' @param GenePairIDMatrix effect.size
#' @param PropOfGenes effect.size
#' @param Bidirectional effect.size
#' @param seed effect.size
#' @return
#' \item{SignalSummary:}{SignalSummary}
#' \item{beta.matrix:}{beta.matrix}

Add.Expr.Asso.Pattern = function(ppp.obj, sim.count, r,
                           perturbed.cell.type,
                           adjacent.cell.type,
                           int.dist.threshold=0.1,
                           delta.mean=1,
                           delta.sd=0.001,
                           GenePairIDMatrix=NULL,
                           PropOfGenes=NULL,
                           Bidirectional=T,
                           seed=NULL) {


  set.seed(seed*3+194)


  R=length(sim.count)
  sim.count1=sim.count[[r]]
  N=ncol(sim.count1)
  G=nrow(sim.count1)
  GeneAll=rownames(sim.count1)

  # key matrix
  beta.matrix=vector("list", length=R)
  for (i in 1:R) {beta.matrix[[i]]=matrix(0, nrow=G,  ncol=ncol(sim.count[[i]]))}
  colnames(beta.matrix[[r]])=colnames(sim.count1)
  rownames(beta.matrix[[r]])=GeneAll

  # spatial info
  nbr.idx=Find.Neighbor.Pairs(ppp.obj=ppp.obj[[r]],
                              interacting.cell.type.pair=c(perturbed.cell.type, adjacent.cell.type),
                              int.dist.threshold=int.dist.threshold)

  # GeneID
  if (is.null(GenePairIDMatrix)) {GenePairIDMatrix=matrix(sample(GeneAll, 2*round(PropOfGenes * G)),
                                                          ncol=2)}

  beta=rnorm(nrow(GenePairIDMatrix), delta.mean, delta.sd)

  # 1 --> 2
  count2=sim.count1[GenePairIDMatrix[,2], nbr.idx[,2]]

  beta.matrix[[r]][GenePairIDMatrix[,1], nbr.idx[,1]]=
    beta.matrix[[r]][GenePairIDMatrix[,1], nbr.idx[,1]]+
    beta*log(count2+1)

  SignalSummary=data.frame(Type="ExprAssoGenes", Region=r, CellType=perturbed.cell.type,
                           GeneID=GenePairIDMatrix[,1], AdjCellType=adjacent.cell.type,
                           AdjGene=GenePairIDMatrix[,2], beta)


  # 2 --> 1
  if (Bidirectional==T) {
    count1=sim.count1[GenePairIDMatrix[,1], nbr.idx[,1]]

    beta.matrix[[r]][GenePairIDMatrix[,2], nbr.idx[,2]]=
      beta.matrix[[r]][GenePairIDMatrix[,2], nbr.idx[,2]]+
      beta*log(count1+1)

    SignalSummary=rbind(SignalSummary,
                        data.frame(Type="ExprAssoGenes", Region=r, CellType=adjacent.cell.type ,
                             GeneID=GenePairIDMatrix[,2], AdjCellType=perturbed.cell.type,
                             AdjGene=GenePairIDMatrix[,1], beta))

  }



  return(list(SignalSummary=SignalSummary, beta.matrix=beta.matrix))
}


# Pattern.adj.1region

Pattern.adj.1region= function(sim.count1, combined.beta.matrix,
                    bond.extreme=T, keep.total.count=T,
                    integer=T) {
  if (is.null(combined.beta.matrix)) {
    # integer
    if (integer==T) {
      sim.count1.update=round(sim.count1)
    } else {
      sim.count1.update=sim.count1
    }
  } else {
    sim.count1.update= exp(log(sim.count1+1) +combined.beta.matrix)-1

    # bond extreme values bonds at 10% reads
    if (bond.extreme==T) {
      TenPerReads=apply(sim.count1.update, 2, sum, na.rm=T)/10
      for (i in 1:ncol(sim.count1.update)) {
        sim.count1.update[,i][which(sim.count1.update[,i]>TenPerReads[i])]=TenPerReads[i]
      }
    }
    # scale by total count
    if (keep.total.count==T) {
      ratio=sum(sim.count1, na.rm=T)/sum(sim.count1.update, na.rm=T)
      sim.count1.update=sim.count1.update*ratio
    }
    # integer
    if (integer==T) {
      sim.count1.update=round(sim.count1.update)
    }

  }

  return(sim.count1.update)
}


#' Adjust the count data for all regions based on the input spatial patterns
#'
#' Adjust the count data for all regions based on the input spatial patterns
#' @param sim.count Spatial info for cell type 1 (e.g. neuron)
#' @param pattern.list Spatial info for cell type 2 (e.g. microglia)
#' @param bond.extreme Expression info for cell type 1 (e.g. neuron)
#' @param keep.total.count Expression info for cell type 1 (e.g. microglia)
#' @param integer Int.Cell.Pair.Idx estimated from function `Find.Neighbor.Pairs`
#' @return sim.count.update
#' @export

Pattern.Adj= function(sim.count, pattern.list=NULL,
                            bond.extreme=T, keep.total.count=T,
                            integer=T) {
  R=length(sim.count)

  sim.count.update=vector("list", length=R)
  for (i in 1:R) {
     if (is.null(pattern.list)) {
       combined.beta.matrix=NULL
     } else {
       K=length(pattern.list)
       beta.matrix.list=lapply(1:K, function(k) pattern.list[[k]]$beta.matrix[[i]])
       combined.beta.matrix=Reduce("+", beta.matrix.list)
     }

     sim.count.update[[i]]=Pattern.adj.1region(sim.count1=sim.count[[i]],
                                                combined.beta.matrix=combined.beta.matrix,
                                                bond.extreme=bond.extreme, keep.total.count=keep.total.count,
                                                integer=integer)

  }

  return(sim.count.update=sim.count.update)
}



ExprPattern=function(pattern.list.i){
  L=length(pattern.list.i)
  res=NULL
  for (l in 1:L) {
    res=rbind(res, pattern.list.i[[l]]$SignalSummary)
  }
  return(res)
}



#' Merge spatial and expression data from multiple regions
#'
#' Merge spatial and expression data from multiple regions
#' @import spatstat
#' @param points.list points.list is a list of points from multiple regions
#' @param expr.list expr.list is a list of expressions from multiple regions
#' @export
#' @return
#' \item{meta:}{meta}
#' \item{count:}{count}

MergeRegion=function(points.list, expr.list) {
  K=length(points.list)
  # points
  x.combine=unlist(lapply(1:K, function(f) points.list[[f]]$x))
  x.combine2=round(x.combine, digits=4)
  y.combine=unlist(lapply(1:K, function(f) points.list[[f]]$y))
  y.combine2=round(y.combine, digits=4)
  annotation=unlist(sapply(1:K, function(f) points.list[[f]]$marks))
  n=sapply(1:K, function(f) points.list[[f]]$n)
  Cell=paste0("Cell", seq(1, length(x.combine)))

  if (K>1) {
    region=rep(1:K, times=n)
    meta=data.frame(Cell=Cell, annotation=annotation,  x.loc=x.combine2,
                    y.loc=y.combine2,
                    region=region)
  } else {
    meta=data.frame(Cell=Cell, annotation=annotation,
                    x.loc=x.combine2, y.loc=y.combine2)
  }

  # expr
  expr.combine=Reduce(cbind, expr.list)
  colnames(expr.combine)=Cell

  return(list(meta=meta,count=as.data.frame(expr.combine)))
}


