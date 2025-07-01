# Identify ICG on simulated data
#rm(list=ls())
library(tidyverse)
library(data.table)
library(raster)
library(spatstat)
library(rlist)
library(parallel)
library(doParallel)

setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_sCCIgen")

# ------------ simulate data based on parameter file --------
# run_interactive_sCCIgen()

# snRNAseq original data based simulation; 
# gene-gene interaction: Epithelial,Adipocyte,0.06,NULL,NULL,0.1,TRUE,0,0.02
# effect size 0,0.1,0.2,0.3,0.4,0.5


input_pool=paste0("Github/sCCIgen_data/sample_parameter_file/", c("snRNAseq_n500_CCI3_0_param.yml","snRNAseq_n500_CCI3_0.1_param.yml",
                                                                  "snRNAseq_n500_CCI3_0.2_param.yml", "snRNAseq_n500_CCI3_0.3_param.yml",
                                                                  "snRNAseq_n500_CCI3_0.4_param.yml","snRNAseq_n500_CCI3_0.5_param.yml"))
mclapply(1:length(input_pool), function(f) ParaSimulation(input=input_pool[f]))


# ------------ Fig5_C1: demonstrate expression-distance interaction (large effect) ------------
 
expr=fread("R1_outputs/snRNAseq_n500_CCI3_0.5_count_1.tsv") %>% as.data.frame %>% column_to_rownames("GeneName")
pattern=fread("R1_outputs/snRNAseq_n500_CCI3_0.5_expr_pattern_1.tsv")
meta=data.frame(fread("R1_outputs/snRNAseq_n500_CCI3_0.5_meta_1.tsv"))

pattern[1:3,]

# find gene for plotting
expr2=expr %>% filter(rownames(expr) %in% pattern$GeneID[1:(nrow(pattern)/2)])
count2=expr2%>%  apply(., 1, sum)
which(count2>2000)
expr3=expr %>% filter(rownames(expr) %in% pattern$AdjGene[1:(nrow(pattern)/2)])
count3=expr3%>%  apply(., 1, sum)
which(count2 >1000 & count3>1000)

# RHOBTB3    SRPK2    STRN3 SERPINF1 
# 133      176      325      376 

# find interaction pairs
pattern[which(pattern$GeneID=="STRN3"),]
# Type Region        CellType GeneID AdjCellType AdjGene beta
# 1: ExprAssoGenes      1 Epithelial cell  STRN3   Adipocyte  EEFSEC  0.1
dat1=meta %>%
  mutate(`Count in Epithelial`=as.numeric(expr[which(rownames(expr)=="STRN3"),]))%>%
  filter(annotation=="Epithelial")
dat2=meta %>%
  mutate(`Count in Adipocyte`=as.numeric(expr[which(rownames(expr)=="EEFSEC"),]))%>% 
  filter(annotation=="Adipocyte")

para=ParaDigest(input_pool[6])
cutoff=para$spatial_int_expr_1_dist_cutoff
distance=crossdist(dat1$x.loc, dat1$y.loc, dat2$x.loc, dat2$y.loc)
pair=which(distance<cutoff, arr.ind=T)

pair2=data.frame(rbind(cbind(x.loc=dat1[pair[,1],"x.loc"],
                             y.loc=dat1[pair[,1],"y.loc"], idx=1:nrow(pair)),
                       cbind(x.loc=dat2[pair[,2],"x.loc"],
                             y.loc=dat2[pair[,2],"y.loc"], idx=1:nrow(pair))))

library(ggnewscale)
window1=data.frame(rbind(c(0,0), c(0,1), c(1,1), c(1,0), c(0,0)))
map <- ggplot() +
  geom_path(aes(x = X1, y = X2), data = window1)+
  coord_equal()+theme_classic() +
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
# map
m1=max(dat1$`Count in Epithelial`)
m2=max(dat2$`Count in Adipocyte`)
p1=map+
  geom_line(data=pair2, colour="grey", aes(x=x.loc, y=y.loc,
                            group=idx))+
  new_scale_colour() +
  geom_point(data=dat1, aes(x=x.loc, y=y.loc,
                            colour=`Count in Epithelial`, shape=annotation))+
  scale_colour_gradientn(colors=c("gray50", "gray90","red"), limits=c(0, m1),
                         values=c(0, 0.1, 1))+  
  new_scale_colour()+
  geom_point(data=dat2, aes(x=x.loc, y=y.loc, colour=`Count in Adipocyte`,  shape=annotation))+
  scale_colour_gradientn(colors=c("gray50", "gray90","blue"), limits=c(0, m2),
                         values=c(0, 0.1, 1))+
  scale_shape_manual("Cell Type", values = c(15, 19))+ 
  theme(plot.title=element_text(hjust = 0.5),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.6, "lines"))+
  labs(title="Small Effect Size (0.5)")
p1

# ------------ Fig5_C1: demonstrate expression-distance interaction (small effect) ------------

expr=fread("R1_outputs/snRNAseq_n500_CCI3_0.1_count_1.tsv") %>% as.data.frame %>% column_to_rownames("GeneName")
pattern=fread("R1_outputs/snRNAseq_n500_CCI3_0.1_expr_pattern_1.tsv")
meta=data.frame(fread("R1_outputs/snRNAseq_n500_CCI3_0.1_meta_1.tsv"))

# find interaction pairs
pattern[which(pattern$GeneID=="SERPINF1"),]
# Type Region        CellType GeneID AdjCellType AdjGene beta
# 1: ExprAssoGenes      1 Epithelial cell  STRN3   Adipocyte  EEFSEC  0.1
dat1=meta %>%
  mutate(Count=as.numeric(expr[which(rownames(expr)=="STRN3"),]))%>%
  filter(annotation=="Epithelial")
dat2=meta %>%
  mutate(Count=as.numeric(expr[which(rownames(expr)=="EEFSEC"),]))%>% 
  filter(annotation=="Adipocyte")

para=ParaDigest(input_pool[2])
cutoff=para$spatial_int_expr_1_dist_cutoff
distance=crossdist(dat1$x.loc, dat1$y.loc, dat2$x.loc, dat2$y.loc)
pair=which(distance<cutoff, arr.ind=T)

pair2=data.frame(rbind(cbind(x.loc=dat1[pair[,1],"x.loc"],
                             y.loc=dat1[pair[,1],"y.loc"], idx=1:nrow(pair)),
                       cbind(x.loc=dat2[pair[,2],"x.loc"],
                             y.loc=dat2[pair[,2],"y.loc"], idx=1:nrow(pair))))

library(ggnewscale)
window1=data.frame(rbind(c(0,0), c(0,1), c(1,1), c(1,0), c(0,0)))
map <- ggplot() +
  geom_path(aes(x = X1, y = X2), data = window1)+
  coord_equal()+theme_classic() +
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
# map
p2=map+
  geom_line(data=pair2, colour="grey", aes(x=x.loc, y=y.loc,
                                           group=idx))+
  new_scale_colour() +
  geom_point(data=dat1, aes(x=x.loc, y=y.loc,
                            colour=Count, shape=annotation))+
  scale_colour_gradientn(colors=c("gray50", "gray90","red"), limits=c(0, m1),
                         values=c(0, 0.1, 1))+  
  new_scale_colour()+
  geom_point(data=dat2, aes(x=x.loc, y=y.loc, colour=Count,  shape=annotation))+
  scale_colour_gradientn(colors=c("gray50", "gray90","blue"), limits=c(0, m2),
                         values=c(0, 0.1, 1))+
  scale_shape_manual("Cell Type", values = c(15, 19))+ 
  theme(plot.title=element_text(hjust = 0.5))+
  labs(title="Large Effect Size (0.1)")
p2

library(patchwork)
p1+p2
# ----------------- ICG analysis pipeline ------------
# correlation within vs outside by effect size 
eff=c(0,0.1,0.2,0.3,0.4,0.5)
parameter_ICG_expr=function(i, cutoff=0.06) {
  expr=fread(paste0("R1_outputs/snRNAseq_n500_CCI3_", eff[i],"_count_1.tsv")) %>% as.data.frame %>% column_to_rownames("GeneName")
  pattern=fread(paste0("R1_outputs/snRNAseq_n500_CCI3_", eff[i],"_expr_pattern_1.tsv"))
  meta=data.frame(fread(paste0("R1_outputs/snRNAseq_n500_CCI3_", eff[i],"_meta_1.tsv")))
    # cell-cell distance
    CellType1=pattern$CellType[1]
    CellType2=pattern$AdjCellType[1]
    
    dist=crossdist(meta[which(meta$annotation==CellType1), 3], 
                   meta[which(meta$annotation==CellType1),4], 
                   meta[which(meta$annotation==CellType2),3], 
                   meta[which(meta$annotation==CellType2),4])
    rownames(dist)=meta[which(meta$annotation==CellType1),"Cell"]
    colnames(dist)=meta[which(meta$annotation==CellType2),"Cell"]
    
    # data
    dlong=dist %>% as.data.frame() %>% rownames_to_column("CellType1")  %>%
      pivot_longer(-CellType1, names_to="CellType2", values_to="Dist") %>%
      mutate(Group=ifelse(Dist<cutoff, 1, 0)) 
    # group==1, gene 1 = correlation
    g1=pattern[1:(nrow(pattern)/2), "GeneID"] %>%
      as.matrix()%>% as.character()
    g2=pattern[1:(nrow(pattern)/2), "AdjGene"] %>%
      as.matrix()%>% as.character()   
    
    idx1=dlong$CellType1[which(dlong$Group==1)]
    idx2=dlong$CellType2[which(dlong$Group==1)]
    expr1=expr[g1, idx1]
    expr2=expr[g2, idx2]
    
    r=sapply(1:nrow(expr1), function(f)
      cor(as.numeric(expr1[f,]), 
          as.numeric(expr2[f,])))
  return(r)
}

res=sapply(1:6, function(f) parameter_ICG_expr(f))
colnames(res)=eff

# - plot 

library(ggplot2)

a=res %>% as.data.frame() %>% pivot_longer(everything()) 

p3=ggplot(a, aes(x=name, y=value, fill=name)) + 
  geom_boxplot()+ 
  theme_classic() +
  labs(x="Effect Size", y = "Correlation")+ 
  theme(legend.position = "none")+
  ggtitle(" ")
p3 

library(cowplot)
legend1 <- get_legend(p1)
plot_grid(p1+ theme(legend.position="none"), 
          p2+ theme(legend.position="none"), legend1,
          p3, rel_widths = c(1, 1, .45, 1), ncol=4)

ggsave("R1/Figures/fig5_c.pdf", width=9, height=3)
