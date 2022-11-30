library(raster)
library(spatstat)
library(tidyverse)
library(rlist)
library(scDesign2)
library(ggplot2)
library(gridExtra)

#' connectUp
#'
#' This function assigns random connected regions on a square. Used within function `RandomRegionWindow`.
#' @param nPoly: nPoly is the No. of regions (e.g. nPoly=3)
#' @param poly: poly is a SpatialPolygonsDataFrame (e.g. 20 by 20 square).
#' @return
#' \item{selected.unique:}{A list of the selected polygons for each region.}

connectUp <- function(poly, nPoly, seed=NULL){
  if (is.null(seed)==F) {set.seed(seed)}
  nb = spdep::poly2nb(poly, queen=F) # find neighbors
  nb2=nb[which(sapply(nb, function(f) min(f))!=0)] # list neighbors for each spot
  g = igraph::graph.adjlist(nb2) # change to network
  nb.length=sapply(nb2, length)
  n=nrow(poly)
  selected = as.list(sample(which(nb.length<4),nPoly))
  selected.sum=unlist(selected)
  n0=length(selected.sum)

  while(n0 < n){
    for (p in 1:nPoly) {
      nbrs = unlist(igraph::ego(g, 1, selected[[p]], mindist=1)) # find neighbors
      selected.sum=unlist(selected)
      newnbrs = nbrs[!nbrs %in% selected.sum]
      if (length(newnbrs)!=0) {
        tb=table(newnbrs)
        cell.id=as.numeric(names(tb))
        prob=as.numeric(tb)/nb.length[cell.id]
        newnbr.freq=round((10^(4*prob-1)))
        newbr.seq=rep(cell.id, newnbr.freq)
        if (length(newbr.seq)>1) {
          selected[[p]] = c(selected[[p]], sample(newbr.seq, 1))
        } else { selected[[p]] = c(selected[[p]], newbr.seq)}
        n0=sum(sapply(selected, length))
      }
    }
    # print(selected)
  }
  length(unlist(selected))
  length(unique(unlist(selected)))
  return(selected)
}


#' RandomRegionWindow
#'
#' This function generates random regions on a unit square.
#' @export
#' @import rgeos
#' @param nRegion: nRegion is the No. of regions (e.g. nRegion=3)
#' @param nGrid: nGrid is the No. of spots on x and y.
#' @return
#' \item{selected.unique:}{A list of the selected polygons for each region.}

RandomRegionWindow <- function(nRegion=3, nGrid=20, seed=123){
  if (nRegion==1) {
    win=vector(mode = "list", length = 1); win[[1]] = unit.square()
  } else {
    # generate polygon
    r <- raster(ncols=nGrid, nrows=nGrid,xmn=0, xmx=1, ymn=0, ymx=1)
    poly <- rasterToPolygons(r)
    # connect
    connected = connectUp(poly=poly, nPoly=nRegion, seed=seed)

    region=lapply(1:nRegion, function(f) rgeos::gUnaryUnion(poly[connected[[f]],]))
    coords=lapply(region, function(f) data.frame(f@polygons[[1]]@Polygons[[1]]@coords))
    win=lapply(1:length(coords), function(f)
      owin(poly=list(x=rev(coords[[f]])[,1], y=rev(coords[[f]])[,2])))
  }

  # plot(win[[1]], xlim=c(0, 1), ylim=c(0,1), col=1)
  # plot(win[[2]], add=T, col=2)
  # plot(win[[3]], add=T, col=3)
  area=sapply(win, area.owin)
  return(list(window=win, area=area))
}

