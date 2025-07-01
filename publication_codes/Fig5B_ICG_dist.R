# Identify ICG on simulated data
#rm(list=ls())
library(tidyverse)
library(data.table)
library(raster)
library(spatstat)
library(rlist)
library(parallel)
library(doParallel)
library(ggnewscale)


setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_sCCIgen")

# ------------ simulate data based on parameter file --------
# run_interactive_sCCIgen()

# snRNAseq original data based simulation; 
# gene-distance interaction: Epithelial,Adipocyte,0.06,NULL,0.1,0,0.02
# effect size 0,0.1,0.2,0.3,0.4,0.5


input_pool=paste0("Github/sCCIgen_data/sample_parameter_file/", c("snRNAseq_n500_CCI2_0_param.yml","snRNAseq_n500_CCI2_0.1_param.yml",
                                                                  "snRNAseq_n500_CCI2_0.2_param.yml", "snRNAseq_n500_CCI2_0.3_param.yml",
                                                                  "snRNAseq_n500_CCI2_0.4_param.yml","snRNAseq_n500_CCI2_0.5_param.yml"))
mclapply(1:length(input_pool), function(f) ParaSimulation(input=input_pool[f]))




# ------------ Fig5b_2: demonstrate expression-distance interaction (small effect) ------------

expr=fread("R1_outputs/snRNAseq_n500_CCI2_0.1_count_1.tsv") %>% as.data.frame %>% column_to_rownames("GeneName")
pattern=fread("R1_outputs/snRNAseq_n500_CCI2_0.1_expr_pattern_1.tsv")
meta=data.frame(fread("R1_outputs/snRNAseq_n500_CCI2_0.1_meta_1.tsv"))

# find gene for plotting
expr2=expr %>% filter(rownames(expr) %in% pattern$GeneID)
count=expr2%>%  apply(., 1, sum)
count[which(count>2000)]

g1=names(which(count>2000))[1]
# find interaction pairs
dat1=meta %>%
  mutate("Count in Epithelial"=as.numeric(expr[which(rownames(expr)==g1),]))%>%
  filter(annotation=="Epithelial")
dat2=meta %>% filter(annotation=="Adipocyte")

para=ParaDigest(input_pool[2])
cutoff=para$spatial_int_dist_1_dist_cutoff
distance=crossdist(dat1$x.loc, dat1$y.loc, dat2$x.loc, dat2$y.loc)
pair=which(distance<cutoff, arr.ind=T)

pair2=data.frame(rbind(cbind(x.loc=dat1[pair[,1],"x.loc"],
                             y.loc=dat1[pair[,1],"y.loc"], idx=1:nrow(pair)),
                       cbind(x.loc=dat2[pair[,2],"x.loc"],
                             y.loc=dat2[pair[,2],"y.loc"], idx=1:nrow(pair))))

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

p1=map+
  geom_line(data=pair2, colour="grey", aes(x=x.loc, y=y.loc,
                                           group=idx))+
  new_scale_colour() +
  geom_point(data=dat1, aes(x=x.loc, y=y.loc,
                            colour=`Count in Epithelial`, shape=annotation))+
  scale_colour_gradientn(colors=c("gray50", "gray90","red"), limits=c(0, 30),
                         values=c(0, 0.2, 1))+  
  geom_point(data=dat2, 
             colour="gray", aes(x=x.loc, y=y.loc, shape=annotation))+
  scale_shape_manual("Cell Type", values = c(15, 19))+ 
  theme(plot.title=element_text(hjust = 0.5))+
  labs(title="Small Effect Size (0.1)")
p1


# ------------ Fig5b_1: demonstrate expression-distance interaction (large effect) ------------
expr=fread("R1_outputs/snRNAseq_n500_CCI2_0.5_count_1.tsv") %>% as.data.frame %>% column_to_rownames("GeneName")
pattern=fread("R1_outputs/snRNAseq_n500_CCI2_0.5_expr_pattern_1.tsv")
meta=data.frame(fread("R1_outputs/snRNAseq_n500_CCI2_0.5_meta_1.tsv"))

# find interaction pairs
dat1=meta %>%
  mutate(Count=as.numeric(expr[which(rownames(expr)==g1),]))%>%
  filter(annotation=="Epithelial")
dat2=meta %>% filter(annotation=="Adipocyte")

para=ParaDigest(input_pool[6])
cutoff=para$spatial_int_dist_1_dist_cutoff
distance=crossdist(dat1$x.loc, dat1$y.loc, dat2$x.loc, dat2$y.loc)
pair=which(distance<cutoff, arr.ind=T)

pair2=data.frame(rbind(cbind(x.loc=dat1[pair[,1],"x.loc"],
                             y.loc=dat1[pair[,1],"y.loc"], idx=1:nrow(pair)),
                       cbind(x.loc=dat2[pair[,2],"x.loc"],
                             y.loc=dat2[pair[,2],"y.loc"], idx=1:nrow(pair))))
dat1[which(dat1$Count>30),]$Count=30
p2=map+
  geom_line(data=pair2, colour="grey", aes(x=x.loc, y=y.loc,
                                           group=idx))+
  new_scale_colour() +
  geom_point(data=dat1, aes(x=x.loc, y=y.loc,
                            colour=Count, shape=annotation))+
  scale_colour_gradientn(colors=c("gray50", "gray90","red"), limits=c(0, 30),
                         values=c(0, 0.2, 1))+
  geom_point(data=dat2, 
             colour="gray", aes(x=x.loc, y=y.loc, shape=annotation))+
  scale_shape_manual("Cell Type", values = c(15, 19))+ 
  theme(plot.title=element_text(hjust = 0.5))+
  labs(title="Large Effect Size (0.5)")


p2



# ----------------- Fig5b_3: ICG analysis pipeline and power/type I figure ------------
eff=c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
parameter_ICG=function(i) {
  
  para=ParaDigest(input_pool[i])
  # cutoff
  cutoff=para$spatial_int_dist_1_dist_cutoff
  res=NULL
    for (j in 1:10) {
    
    # load data
    expr=fread(paste0("R1_outputs/snRNAseq_n500_CCI2_", eff[i], "_count_", j, ".tsv")) %>% as.data.frame %>%
      column_to_rownames("GeneName")
    pattern=fread(paste0("R1_outputs/snRNAseq_n500_CCI2_", eff[i], "_expr_pattern_", j, ".tsv"))
    meta=data.frame(fread(paste0("R1_outputs/snRNAseq_n500_CCI2_", eff[i], "_meta_", j, ".tsv")))
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
      mutate(Group=ifelse(Dist<=cutoff, 1, 0)) %>%
      dplyr::select(CellType1, Group) %>%
      group_by(CellType1) %>%
      summarise(max=max(Group))
    
    idx1=dlong$CellType1[which(dlong$max==1)]
    idx0=dlong$CellType1[which(dlong$max==0)]
    
    # t.test
    
    expr1=expr[, idx1]
    expr0=expr[, idx0]
    excludeID=union(which(apply(expr1, 1, function(f) all(is.na(f)))),
    which(apply(expr0, 1, function(f) all(is.na(f)))))
    
    expr1_clean=expr1[-excludeID, ]
    expr0_clean=expr0[-excludeID, ]
    tt=unlist(mclapply(1:nrow(expr1_clean), function(f)
     t.test(as.numeric(expr1_clean[f, ]),
                 as.numeric(expr0_clean[f, ]))$p.value))
     
    
    sig1=which(rownames(expr1_clean) %in% pattern$GeneID)
    sig0=which(!(rownames(expr1_clean) %in% pattern$GeneID))
    power=mean(tt[sig1]<0.05, na.rm=T)
    type1=mean(tt[sig0]<0.05, na.rm=T)
    res=rbind(res, data.frame(power=power, type1=type1))
    
    }
  return(res)
}
r=mclapply(1:6, function(f) parameter_ICG(f))
r2=sapply(r, function(f) apply(f, 2, mean))
r2
# - plot 
a=rbind(eff,matrix(unlist(r2), nrow=2))  %>%t() %>% as.data.frame()


b=a %>% pivot_longer(!"eff")
p3=ggplot(b, aes(x=eff, y=value, col=name)) + 
  geom_point()+ geom_line()+ 
  theme_classic() +
  ylim(0,1)+
  labs(x="Effect Size", y = " ")+ 
  scale_color_discrete(name = " ", labels = c("Power", "Type I Error"))+
  theme(legend.position =   c(0.7, 0.3)) +
  ggtitle(" ")
p3

# merge sub-figures
library(cowplot)
legend1 <- get_legend(p1)
plot_grid(p2+ theme(legend.position="none"), 
          p1+ theme(legend.position="none"), legend1,
          p3, rel_widths = c(1, 1, .47, 1), ncol=4)

ggsave("R1/Figures/fig5_b.pdf", width=9, height=3)

