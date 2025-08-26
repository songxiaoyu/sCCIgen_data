# Identify spatial variable genes on simulated data
#rm(list=ls())
library(tidyverse)
library(data.table)
library(raster)
library(spatstat)
library(rlist)
library(parallel)
library(doParallel)
library(sCCIgen)
library(cowplot)
library(Seurat)
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_sCCIgen")


########################
# Simulate data: Region-specific genes  ----------------- 
######################## 
# run_interactive_sCCIgen()
# With CCI2
input="Github/sCCIgen_data/sample_parameter_file/snRNAseq_Region3_RVG_DisExpr_Fig7_param.yml"
ParaSimulation(input=input)
# Without CCI2
input="Github/sCCIgen_data/sample_parameter_file/snRNAseq_Region3_RVG_Fig7_param.yml"
ParaSimulation(input=input)

########################
# Load reference and simulated data  ----------------- 
########################
load("Github/sCCIgen_data/input_data/snRNAseq_breast_2025_expr.Rdata")
anno=colnames(expr)

expr2=fread("R1_outputs/snRNAseq_Region3_RVG_DisExpr_Fig7_count_1.tsv")  %>% column_to_rownames("GeneName")
pattern2=fread("R1_outputs/snRNAseq_Region3_RVG_DisExpr_Fig7_expr_pattern_1.tsv")
meta2=fread("R1_outputs/snRNAseq_Region3_RVG_DisExpr_Fig7_meta_1.tsv")


#####################################
############## plots #############
#####################################

# template ------------
input="Github/sCCIgen_data/sample_parameter_file/snRNAseq_Region3_RVG_DisExpr_Fig7_param.yml"
para=ParaDigest(input)
win=RandomRegionWindow(nRegion=para$num_regions, seed=para$all_seed[[1]])
# win=RandomRegionWindow(nRegion=3, seed=1111)
# plot(win$window[[1]], col="pink")
# plot(win$window[[2]], col="blue", add=T)
# plot(win$window[[3]], col="orange", add=T)

line1=data.frame(x=win$window[[1]]$bdry[[1]]$x, y=win$window[[1]]$bdry[[1]]$y, r=1)
line1=rbind(line1, line1[1,])
line2=data.frame(x=win$window[[2]]$bdry[[1]]$x, y=win$window[[2]]$bdry[[1]]$y, r=2)
line2=rbind(line2, line2[1,])
line3=data.frame(x=win$window[[3]]$bdry[[1]]$x, y=win$window[[3]]$bdry[[1]]$y, r=3)
line3=rbind(line3, line3[1,])
template=ggplot() +coord_fixed() +
  coord_equal()+theme_classic() +
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.title=element_blank(),
        plot.title=element_text(hjust = 0.5))  + 
  geom_path(aes(x = x, y = y), data = line1, color="gray30", linewidth=0.5)+ 
  geom_path(aes(x = x, y = y), data = line2, color="gray30", linewidth=0.5)+ 
  geom_path(aes(x = x, y = y), data = line3, color="gray30", linewidth=0.5) 

# Figure 7a: Plot cell type by region --------------------
library(RColorBrewer)
p1=template + 
  geom_point(data=meta2, aes(x=x.loc, y=y.loc, color=annotation, shape=annotation))+
  scale_color_manual(values = brewer.pal(6, "Set2")) +
  theme(legend.position="none")

p1
p1_legend=template + 
  geom_point(data=meta2, aes(x=x.loc, y=y.loc, color=annotation, shape=annotation))+
  scale_color_manual(values = brewer.pal(6, "Set2")) +
  theme(legend.position="bottom")


# Figure 7b: plot interaction pairs -----------------
d1=meta2 %>%filter(annotation=="Epithelial")
d2=meta2 %>% filter(annotation=="Immune")
cutoff=para$spatial_int_dist_1_dist_cutoff
distance=crossdist(d1$x.loc, d1$y.loc, d2$x.loc, d2$y.loc)
pair=which(distance<cutoff, arr.ind=T)
pair2=data.frame(rbind(cbind(d1[pair[,1],], idx=1:nrow(pair)),
                       cbind(d2[pair[,2],], idx=1:nrow(pair))))
meta2_sub=meta2[which(meta2$annotation=="Immune" | meta2$annotation=="Epithelial"),]

p3=template+
  geom_point(data=meta2_sub, aes(x=x.loc, y=y.loc, color=annotation, shape=annotation))+ 
  geom_line(data=pair2, colour="blue", aes(x=x.loc, y=y.loc, group=idx), linewidth=0.3)+
  theme(legend.position="none")
p3


# Figure 6c: Plot signature genes + ICG  expression by block ---------------
# - create four gene categories 
idx1=pattern2[which(Type=="DistanceAssoGenes"), ]$GeneID
idx2=pattern2[which(Type=="SpatialChange"), ]$GeneID
# four gene categories
x11=intersect(idx1, idx2) # SVG + SG
x10=setdiff(idx1, idx2) # ICG only
x01=setdiff(idx2, idx1) # SVG
x00=setdiff(rownames(expr2), union(idx1, idx2)) # non-ICG non-SVG

naf=apply(expr2, 1, function(f) mean(is.na(f)))
expr2r=expr2[which(naf==0),]
meta2$x11=apply(expr2r[rownames(expr2r) %in% x11,], 2, mean) %>% scale()%>% pmin(., 3) %>% pmax(., -3)
meta2$x10=apply(expr2r[rownames(expr2r) %in% x10,], 2, mean) %>% scale()%>% pmin(., 3) %>% pmax(., -3)
meta2$x01=apply(expr2r[rownames(expr2r) %in% x01,], 2, mean) %>% scale()%>% pmin(., 3) %>% pmax(., -3)
meta2$x00=apply(expr2r[rownames(expr2r) %in% x00,], 2, mean) %>% scale()%>% pmin(., 3) %>% pmax(., -3)

meta2_epi=meta2[which(meta2$annotation=="Epithelial"),]
meta2_nonepi=meta2[which(meta2$annotation!="Epithelial"),]
# plot 

p21=template + 
  geom_point(data=meta2_epi, aes(x=x.loc, y=y.loc, color=x11))+ 
  scale_colour_gradientn(colors=c("blue", "gray90","red"), limits=c(-3,3), breaks=c(-2, 0, 2)) 
#p21
p22=template + 
  geom_point(data=meta2_epi, aes(x=x.loc, y=y.loc, color=x10))+ 
  scale_colour_gradientn(colors=c("blue", "gray90","red"), limits=c(-3,3), breaks=c(-2, 0, 2)) 
#p22
p23=template + 
  geom_point(data=meta2_epi, aes(x=x.loc, y=y.loc, color=x01))+ 
  scale_colour_gradientn(colors=c("blue", "gray90","red"), limits=c(-3,3), breaks=c(-2, 0, 2)) 
p23
p24=template + 
  geom_point(data=meta2_epi, aes(x=x.loc, y=y.loc, color=x00))+ 
  scale_colour_gradientn(colors=c("blue", "gray90","red"), limits=c(-3,3), breaks=c(-2, 0, 2)) 
#p24
p2=plot_grid(p21+ theme(legend.position="none"), 
             p22+ theme(legend.position="none"), 
             p23+ theme(legend.position="none"),
             p24+ theme(legend.position="none"), nrow=1)
p2
p25=template + 
  geom_point(data=meta2_nonepi, aes(x=x.loc, y=y.loc, color=x11))+ 
  scale_colour_gradientn(colors=c("blue", "gray90","red"), limits=c(-3,3), breaks=c(-2, 0, 2)) 
#p25
p26=template + 
  geom_point(data=meta2_nonepi, aes(x=x.loc, y=y.loc, color=x10))+ 
  scale_colour_gradientn(colors=c("blue", "gray90","red"), limits=c(-3,3), breaks=c(-2, 0, 2)) 
#p26
p27=template + 
  geom_point(data=meta2_nonepi, aes(x=x.loc, y=y.loc, color=x01))+ 
  scale_colour_gradientn(colors=c("blue", "gray90","red"), limits=c(-3,3), breaks=c(-2, 0, 2)) 
#p27
p28=template + 
  geom_point(data=meta2_nonepi, aes(x=x.loc, y=y.loc, color=x00))+ 
  scale_colour_gradientn(colors=c("blue", "gray90","red"), limits=c(-3,3), breaks=c(-2, 0, 2)) 
#p28
p2non=plot_grid(p25+ theme(legend.position="none"), 
                p26+ theme(legend.position="none"), 
                p27+ theme(legend.position="none"),
                p28+ theme(legend.position="none"), nrow=1)
p2non
# Figure 6d: make it multi-cell resolution ---------------

multi=multicell(expr=expr2, spatial=meta2, NoSpot=400, cl=10)
expr3=multi$count
meta3=multi$spot_feature
save(expr3, meta3, file="R1_outputs/snRNAseq_Region3_SVG_DisExpr_FigS1_spot_1.Rdata")
load("R1_outputs/snRNAseq_Region3_SVG_DisExpr_FigS1_spot_1.Rdata")

naf=which(apply(expr3, 2, function(f) all(f==0)))
meta3r=meta3[-naf,]
expr3r=expr3[, -naf]

exprr3 =  NormalizeData(expr3r, normalization.method = "LogNormalize", scale.factor = 10000) 
temp=apply(exprr3, 2, mean, na.rm=T)
tm=mean(temp) 
ts=sd(temp)*2

meta3r$x11=apply(exprr3[rownames(exprr3) %in% x11,], 2, mean, na.rm=T) %>% scale(., scale=ts) %>% pmin(., 3) %>% pmax(., (-3))
meta3r$x10=apply(exprr3[rownames(exprr3) %in% x10,], 2, mean, na.rm=T) %>% scale(.,scale=ts) %>% pmin(., 3) %>% pmax(., (-3))
meta3r$x01=apply(exprr3[rownames(exprr3) %in% x01,], 2, mean, na.rm=T) %>% scale(., scale=ts) %>% pmin(., 3) %>% pmax(., (-3))
meta3r$x00=apply(exprr3[rownames(exprr3) %in% x00,], 2, mean, na.rm=T) %>% scale(., scale=ts) 


p41=template +  
  geom_point(data=meta3r, aes(x=x, y=y, color=x11), shape=15, size=3)+
  scale_colour_gradientn(colors=c("blue", "gray90","red"), limits=c(-3,3), breaks=c(-2, 0, 2), name="Set 1")
p41

p42=template +  
  geom_point(data=meta3r, aes(x=x, y=y, color=x10), shape=15, size=3)+
  scale_colour_gradientn(colors=c("blue", "gray90","red"), limits=c(-3,3), breaks=c(-2, 0, 2), name="Set 1")
p42
p43=template +  
  geom_point(data=meta3r, aes(x=x, y=y, color=x01), shape=15, size=3)+
  scale_colour_gradientn(colors=c("blue", "gray90","red"), limits=c(-3,3), breaks=c(-2, 0, 2), name="Set 1")
p43
p44=template +  
  geom_point(data=meta3r, aes(x=x, y=y, color=x00), shape=15, size=3)+
  scale_colour_gradientn(colors=c("blue", "gray90","red"), limits=c(-3,3), breaks=c(-2, 0, 2), name="Set 1")
p44
p4=plot_grid(p41+ theme(legend.position="none"),
             p42+ theme(legend.position="none"),
             p43+ theme(legend.position="none"),
             p44+ theme(legend.position="none"),nrow=1)
p4

#####################################
############## SVG #############
#####################################
# Figure 7e: SVG Identify spatial variable genes on simulated data using anova - compare between no CCI and CCI



Parameter_SVG=function(expr2, meta2){
  # NA remove first 
  naf=which(apply(expr2, 1, function(f) any(is.na(f))))
  expr2r=expr2[-naf, ]
  exprr2 =  NormalizeData(expr2r, normalization.method = "LogNormalize", scale.factor = 10000)
  p=rep(NA, nrow(expr2))
  for (i in 1: nrow(expr2)) {
    if (!i %in% naf) {
      y=exprr2[which(rownames(exprr2) == rownames(expr2)[i]),]
      p[i]=anova(lm(y~factor(meta2$region)))$`Pr(>F)`[1]
    }
  }
  return(p)
}


expr2=fread("R1_outputs/snRNAseq_Region3_RVG_DisExpr_FigS1_count_1.tsv")  %>% column_to_rownames("GeneName")
meta2=fread("R1_outputs/snRNAseq_Region3_RVG_DisExpr_FigS1_meta_1.tsv")
P_CCI=Parameter_SVG(expr2, meta2)

expr2=fread("R1_outputs/snRNAseq_Region3_RVG_FigS1_count_1.tsv")  %>% column_to_rownames("GeneName")
meta2=fread("R1_outputs/snRNAseq_Region3_RVG_FigS1_meta_1.tsv")
P_noCCI =Parameter_SVG(expr2, meta2)
hist(P_noCCI)

pdat=cbind(P_CCI, P_noCCI)

p5=ggplot() + geom_point(data=pdat, aes(x= -log10(P_noCCI), y= -log10(P_CCI)))+coord_fixed() +
  coord_equal()+theme_classic()+
  geom_abline(intercept=0, slope=1, col="red")+
  labs(x="-log10(P) w/o CCI", 
       y="-log10(P) with CCI")+
  xlim(0, -log10(min(pdat, na.rm=T)))+
  ylim(0, -log10(min(pdat, na.rm=T)))
p5



#####################################
############## SVG #############
#####################################


ggdraw() +
  draw_plot(p1, x = 0, y = .7, width = .3, height = .3) +
  draw_plot(p3, x = .35, y = .7, width = .3, height = .3) +
  draw_plot(p5, x = .7, y = .7, width = .3, height = .3) +
  draw_plot(p2, x = 0, y = 0.45, width = 1, height = 0.2) +
  draw_plot(p2non, x = 0, y = 0.25, width = 1, height = 0.2) +
  draw_plot(p4, x = 0, y = 0, width = 1, height = 0.25) + 
  draw_plot_label(label = c("a", "b", "e", "c", "d"))

ggsave("R2/fig7.pdf", width=9, height=12)
