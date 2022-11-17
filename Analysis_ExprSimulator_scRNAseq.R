rm(list=ls())

source("Paper_Simulator/codes/Function_ExprSimulator_scRNAseq.R")
## ----- SeqFISH ------
load("/Users/songxiaoyu152/Dropbox/SpatialTranscriptomics/Paper_Simulator/DataToShare/SeqFishPlus.RData")
loc=SeqFishPlus_loc
anno=SeqFishPlus_anno
ct=apply(SeqFishPlus_exprs, 1, sum)
expr=SeqFishPlus_exprs[which(ct>1000)[1:20],] #  10 variable genes

# first run analysis PointSimulator --------
seed=1523
work.loc=window()
n.cell=n.cell.type.fc(n=500, cell.prop=c(0.2, 0.2,0.2, 0.2, 0.2))
n.cell
# Option 1: random generate 
# den=random.common.density.fc(K=4, seed=round(seed/100))$density
# Option 2: commonly used density function
da=common.density.fc(opt.x="random", opt.y="random")
db=common.density.fc(opt.x="polynomial", opt.y="random", poly.coef.x=c(1, 4, 8))
dc=common.density.fc(opt.x="Gaussian", opt.y="Gaussian", G.mu.x=0.2, G.sd.x=0.2,
                     G.mu.y=0.5, G.sd.y=0.4)
dd=common.density.fc(opt.x="random", opt.y="polynomial", G.mu.x=0.6, G.sd.x=0.2,
                     poly.coef.y=c(0, 0.1, 0.5))  
de=common.density.fc(opt.x="random", opt.y="random")
den=list(da, db, dc,dd, de)
# Option 3: user specified functions

# simulation location of cells from real data
pts=cell.loc.function.fc(n.vec=n.cell, density.list=den, window=work.loc, 
                         same.loc.force.ratio =0.1, 
                         inhibition.attraction.celltype.pair=rbind(c(1,2), c(3,3),c(4,5)),
                         inhibition.attraction.coef=c(2, 0,0),
                         density.dist.flexibility.coef=c(1,1,1,1,1),
                         even.distribution.coef=0) 
pts1=pts$inhi.att.exclude.ppp

# use scDesign2
expr1=Use_scDesign2(ppp.obj=pts1, expr=expr, anno=anno)

expr1$cell.type.match
expr3 = Use_scDesign2(ppp.obj=pts1, expr, anno,
                               min.cell.type.num=10,
                               input.cell.type=c("Microglia", "Endothelial", 
                                                 "Neuroblast", "Oligodendrocyte", 
                                                 "Interneuron"),
                               total_count_new=50000) 
  
  

# add cell-cell interaction using SeqFISH+ data --------

expr2=Add.Cell.Cell.Int(ppp.obj=expr1$ppp.obj.update, sim.count=expr1$sim.count, 
                            dist.threshold=0.05,
                            Gaussian.effect.size.parameter=c(1,0.1),
                           interacting.cell.types=list(c("ExcitatoryNeuron", "Astrocyte")),
                          mech1.prop=0.1, mech2.prop=0.2,
                            mech3.prop=0.1, seed=1531) 

  
contrast.mech1.expr(ppp.obj=expr2$ppp.obj, gene1before=expr2$sim.count.before[1,], 
                    gene1after=expr2$sim.count.after[1,],
                    neighbor=expr2$neighbor.all.pairs[[1]]) 


contrast.mech2.expr.neighbor(ppp.obj=expr2$ppp.obj, gene1before=expr2$sim.count.before[2,], 
                             gene1after=expr2$sim.count.after[2,],
                             gene2before=expr2$sim.count.before[13,], 
                             gene2after=expr2$sim.count.after[13,],
                             cell.type=c("ExcitatoryNeuron", "Astrocyte"),
                             neighbor=expr2$neighbor.all.pairs[[1]]) 


# add Spatial pattern using SeqFISH+ data --------

  
  
  
  


