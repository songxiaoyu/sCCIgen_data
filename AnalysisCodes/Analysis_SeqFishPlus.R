rm(list=ls())
setwd("/Users/songxiaoyu152/Library/CloudStorage/OneDrive-TheMountSinaiHospital/SpatialTranscriptomics/Paper_Simulator")
# source("Paper_Simulator/Github/R/PointSimulator_STData.R")
# source("Paper_Simulator/Github/R/ExprSimulator.R")
# library(splines)
# library(geometry)




## ----- MERFISH s260 ------
# load("/Users/songxiaoyu152/Dropbox/SpatialTranscriptomics/Paper_Simulator/DataToShare/MERFISH_s260.RData")
# loc=MERFISH_s260_loc
# anno=MERFISH_s260_anno
# expr=MERFISH_s260_exprs 
# expr=round(expr)
## ----- SeqFISH ------
load("Github/Data/SeqFishPlusCortexFilter_expr.Rdata")
dim(expr)
load("Github/Data/SeqFishPlusCortexFilter_loc.Rdata")
dim(loc)

# spatial
ft=cell.loc.model.fc(n=1000, PointLoc=loc, PointAnno=anno)
# expr
expr=Use_scDesign2(ppp.obj=ft,  expr=expr, anno=anno)















# ----- old ------
simulator.with.ST=function(n, loc, anno, exprs, min.no.cell=5) {
  # filter out input cell type that have too few cells
  anno_no=table(anno)
  anno_exclude=names(anno_no)[which(as.numeric(anno_no)<min.no.cell)]
  loc_filter=loc[which(!anno %in% anno_exclude),]
  anno_filter=anno[which(!anno %in% anno_exclude)]
  expr_filter=expr[,which(!anno %in% anno_exclude)]
  
  work.loc=scale.loc(PointLoc=loc_filter)
  # generate window
  window=simu.window(PointLoc=work.loc, method="network")
  ppp.obj=ppp.obj.fc(PointLoc=work.loc,PointAnno=anno_filter)
  plot(ppp.obj, pch=16, show.window=F, cols=1:11)
  # simulation location of cells from real data
  pts=cell.loc.model.fc(n=n,ppp.obj=ppp.obj) 
  plot.by.cell(ppp.obj, nrows=2)
  plot.by.cell(pts, nrows=2) 
  
  
  # estimate expression levels of the cells 
  ppp.obj=pts, cell1="Astrocyte", cell2="ExcitatoryNeuron"

  
  # simulation expression levels of the cells 

  
}
plot.by.cell(ppp.obj, nrows=1)
plot.by.cell(pts, nrows=2)

# ------ scRNAseq data ------- 
## ----- SeqFISH ------
load("DataToShare/SeqFishPlus.RData")
expr=SeqFishPlus_exprs
expr=expr[1:20,]
anno=SeqFishPlus_anno
simulator.with.scRNAseq=function(n=1000, G=30, K=3, cell.prop=c(0.5, 0.3, 0.2), expr=expr, anno=anno, seed=1453) {
  # n=1000;K=3;cell.prop=c(0.5, 0.3, 0.2);  G=30; seed=1453
  window=simu.window()
  n.cell=n.cell.type.fc(n=n, cell.prop=cell.prop, seed=seed)
  
  # simulation location of cells from real data
  da=common.density.fc(opt.x="random", opt.y="random")
  db=common.density.fc(opt.x="polynomial", opt.y="random", poly.coef.x=c(1, 4, 8))
  dc=common.density.fc(opt.x="Gaussian", opt.y="Gaussian", G.mu.x=0.2, G.sd.x=0.1,
                       G.mu.y=0.5, G.sd.y=0.2)
  den=list(da, db, dc)
  pts=cell.loc.random.fc(n.vec=n.cell, density.list=den, window=window,
                         same.loc.dist=0.005, 
                         inhibition.celltype.pair=c(1,2),
                         attraction.celltype.pair=c(1,3),
                         strength=0.5)
  
  # estimate expression levels of the cells 
  est.data=est.copula.fc(expr=expr, anno=anno, min.cell.no=50, zp_cutoff = 0.8)
  # simulation expression levels of the cells 
  dat=simu.ST.expr.fc(ppp.obj=pts$final.ppp, G=G, est.data=est.data, No.spatial.Gene=c(1,2,1), 
                  No.Neighbor.Cell.Gene=c(2,2,1), No.Neighbor.GenePair=c(2,1,3), 
                  n.read=50000)
  
  
}


library(ggplot2)
plot.spatial.gene.fc=function(ppp.obj=pts, expr=dat$expr, k=1, g=dat$Spatial.gene.idx[[1]]) {
  count=expr[g,]
  count[pts$final.ppp$marks!=k]=NA
  dat1=data.frame(loc.x=pts$final.ppp$x, loc.y=pts$final.ppp$y, 
                  mark=pts$final.ppp$marks, count=count)
  
  p=ggplot(dat1, aes(x=loc.x, y=loc.y, col=count, shape=factor(mark))) + geom_point() + 
    scale_color_gradient2(low="white", high="red")
  p
  return(p)
}

p1=plot.spatial.gene.fc(ppp.obj=pts, expr=dat$expr, k=1, g=dat$Spatial.gene.idx[[1]]);p1
p1=plot.spatial.gene.fc(ppp.obj=pts, expr=dat$expr, k=2, g=2);p1

library(ggplot2)
plot.neighbor.cell.gene.fc=function(ppp.obj=pts, expr=dat$expr, k=1, g=dat$Cell.gene.idx[[1]][1]) {
  count=expr[g,]
  count[pts$final.ppp$marks!=k]=NA
  dat1=data.frame(loc.x=pts$final.ppp$x, loc.y=pts$final.ppp$y, 
                  mark=pts$final.ppp$marks, count=count)

  p=ggplot(dat1, aes(x=loc.x, y=loc.y, col=count, shape=factor(mark))) + geom_point(size=2) + 
    scale_color_gradient2(low="white", high="red")
  return(p)
}
p2=plot.neighbor.cell.gene.fc();p2
p2=plot.neighbor.cell.gene.fc(g=3);p2

library(ggplot2)
plot.gene.pair.fc=function(ppp.obj=pts, counts=dat$expr, k1=1, k2=2, g=dat$gene.pair.idx[[1]][1,]) {

  count1=counts[g[1],]
  count1[pts$final.ppp$marks!=k1]=NA
  count2=counts[g[2],]
  count2[pts$final.ppp$marks!=k2]=NA
  count=ifelse(is.na(count1)==F, count1, -count2)
  dat1=data.frame(loc.x=pts$final.ppp$x, loc.y=pts$final.ppp$y, 
                  mark=pts$final.ppp$marks, count=count)
  
  p=ggplot(dat1) + geom_point(aes(x=loc.x, y=loc.y, col=count, shape=factor(mark)),size=2) + 
    scale_color_gradient2(low="blue",mid="white", high="red") 
  p
  return(p)
}

p3=plot.gene.pair.fc();p3

p3=plot.gene.pair.fc(g=c(5,7));p3

