

#' Simulate the window of spatial data
#'
#' This function esimates window using spatial location data of existing cells.
#' @param PointLoc PointLoc
#' @param method method=c("convex", "rectangle", "network")
#' @import spatstat
#' @export
# window=simu.window(PointLoc=NULL)
# window=simu.window(PointLoc=PointLoc, method="network")
simu.window=function(PointLoc=NULL, method="network") {

  if (is.null(PointLoc)==F & (method=="convex" | method=="rectangle")) {
    x0=PointLoc[,1]
    y0=PointLoc[,2]
    res=ripras(x0, y0, shape=method)
  }
  if (is.null(PointLoc)==F & method=="network") {

    # generate random perturbation
    n=nrow(PointLoc)

    dx=max(PointLoc[,1])-min(PointLoc[,1])
    dy=max(PointLoc[,2])-min(PointLoc[,2])
    dmax=max(dx, dy)
    x=PointLoc[,1]+0.005*dmax*cos(seq(0, 2*pi, length=n*4))
    y=PointLoc[,2]+0.005*dmax*sin(seq(0, 2*pi, length=n*4))
    PointLoc_noise=cbind(x,y)

    delaunay_triangle = geometry::delaunayn(PointLoc_noise)
    max_dist=sapply(1:nrow(delaunay_triangle), function(f)
      max(dist(PointLoc_noise[delaunay_triangle[f,],])))
    #trimesh(delaunay_triangle,cbind(x,y))

    ## get rid of the edges that connect distant points
    max_dist_standize=max_dist/max(max_dist)
    log_dist=log(max_dist_standize+1)
    idx=which(log_dist> mean(log_dist) + 3*sd(log_dist))
    delaunay_triangle2=delaunay_triangle[-idx,]
    #trimesh(delaunay_triangle2,PointLoc_noise)

    ## get the outer window
    delaunay_edges <- rbind(delaunay_triangle2[ ,c(1,2)],
                            delaunay_triangle2[ ,c(1,3)],
                            delaunay_triangle2[ ,c(2,3)])
    delaunay_edges_ordered=cbind(apply(delaunay_edges, 1, min),
                                 apply(delaunay_edges, 1, max))
    dup1= duplicated(delaunay_edges_ordered)
    dup2=duplicated(delaunay_edges_ordered, fromLast = TRUE)
    outer_edges=delaunay_edges[which(dup1==F & dup2 ==F),]

    points_ordered=outer_edges[1,]
    edges_remained=outer_edges[-1,]
    n=nrow(edges_remained)
    for (i in 1:(n-1)) {
      #print(i)
      idx=which(apply(edges_remained, 1, function(r)
        any(r==points_ordered[i+1])))
      points_ordered=c(points_ordered,
                       setdiff(edges_remained[idx,], points_ordered[i+1]))
      edges_remained=edges_remained[-idx,]
      #print(points_ordered)
    }

    a <- try(owin(poly=list(x=PointLoc_noise[points_ordered,1],
                            y=PointLoc_noise[points_ordered,2])), silent =T)
    if (inherits(a, "try-error")) {
      points_ordered2=rev(points_ordered)
      res=owin(poly=list(x=PointLoc_noise[points_ordered2,1],
                         y=PointLoc_noise[points_ordered2,2]))
    } else {res=a}

  }
  return(res)
}



#' simulate ST data location based on parametric model
#' @import spatstat
#' @export
cell.loc.model.fc=function(n,
                           PointLoc,
                           PointAnno,
                           window_method,
                           seed=NULL) {

  #
  if(is.null(seed)==F) {set.seed(seed)}
  cell_win=simu.window(PointLoc=PointLoc, method=window_method)
  p=as.ppp(PointLoc, W=cell_win)
  marks(p)=as.factor(PointAnno)
  # if too many cells
  if (p$n>10000) {
    idx=rbinom(p$n, 1, prob=10000/p$n)
    p2=subset(p, idx==1)
    p=p2
  }

  x=p$x
  y=p$y
  fit=ppm(p, ~marks+ marks:polynom(x,y,2),Poisson())
  nsim=ceiling(n/p$n)
  if (nsim>1) {
    a=rmh(model=fit, nsim=nsim)
    b=superimpose(a)
    marks(b)=marks(b)[,1]
  } else{b=rmh(model=fit, nsim=nsim)}

  # get rid of cells on the same location
  while(b$n-n>10) {
    delete.n=(b$n-n)
    dis=pairdist(b)
    dis[lower.tri(dis, diag=T)]=NA
    r= quantile(dis, probs=delete.n*2/b$n/(b$n-1), na.rm=T)
    dis2=dis< r
    same.loc.idx=which(dis2 == T, arr.ind = TRUE)
    b=b[setdiff(1:b$n, same.loc.idx[,1]), ]
  }

  return(b)
}




#' simulate ST data location based on parametric model
#' @import spatstat
#' @export
cell.region.loc.model.fc=function(n,
                           PointLoc,
                           PointAnno,
                           PointRegion,
                           window_method,
                           seed=NULL) {
  Rcat=unique(PointRegion)
  bb=vector("list", length(Rcat))
  for ( i in 1:length(Rcat)) {
    idx= which(PointRegion %in% Rcat[i])
    bb[[i]]=cell.loc.model.fc(n=length(idx),
                              PointLoc=PointLoc[idx,],
                              PointAnno=PointAnno[idx],
                             window_method=window_method,
                             seed=seed)
  }
  return(bb)
}



#' simulate ST data location based on parametric model
#' @import spatstat
#' @export
cell.loc.existing.fc=function(PointLoc,
                           PointAnno,
                           window_method="rectangle") {


  cell_win=simu.window(PointLoc=PointLoc, method=window_method)
  p=as.ppp(PointLoc, W=cell_win)
  marks(p)=as.factor(PointAnno)

  return(p)
}

#' simulate ST data location based on parametric model
#' @import spatstat
#' @export
cell.region.loc.existing.fc=function(PointLoc,
                              PointAnno,
                              PointRegion,
                              window_method="rectangle") {
  Rcat=unique(PointRegion)
  pp=vector("list", length=length(Rcat))
  for (i in 1:length(Rcat)) {
    idx= which(PointRegion %in% Rcat[i])
    pp[[i]]=cell.loc.existing.fc(PointLoc=PointLoc[idx,],
                         PointAnno=PointAnno[idx],
                         window_method=window_method)
  }

  return(pp)
}
