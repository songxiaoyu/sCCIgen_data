rm(list=ls())
# Identify spatial regions on simulated data
library(spatstat)
library(sCCIgen)

setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_sCCIgen")

# ---- data generation -------------
# run_interactive_sCCIgen()
# simulated from snRNAseq breast
# Adipocyte,Epithelial,3 (cell-cell Attraction)
# Adipocyte,Epithelial,0
# Adipocyte,Epithelial,-3 (cell-cell inhibition)

input_pool=c("Github/sCCIgen_data/sample_parameter_file/snRNAseq_n500_CCI1_3_param.yml", 
             "Github/sCCIgen_data/sample_parameter_file/snRNAseq_n500_CCI1_0_param.yml", 
             "Github/sCCIgen_data/sample_parameter_file/snRNAseq_n500_CCI1_neg3_param.yml")
parameter_coloc=function(input) {
  para=ParaDigest(input);  
  attach(para); 
  p1=ParaCellsNoST(para=para, seed_list=all_seeds[[1]])[[1]][[1]]; 
  detach(para)
  return(p1)
}

d1=parameter_coloc(input_pool[1])
d2=parameter_coloc(input_pool[2])
d3=parameter_coloc(input_pool[3])


# ---- demonstration plots  -------------

# attraction 
a1=split(d1)[[1]]; dat1=data.frame(x=a1$x, y=a1$y, type="Adipocyte")
a3=split(d1)[[3]]; dat3=data.frame(x=a3$x, y=a3$y, type="Epithelial Cell")
dat=rbind(dat1, dat3)
dist=crossdist(a1, a3)
nn1=apply(dist, 1, min)
p1= ggplot(dat, aes(x=x, y=y, col=type))+geom_point(shape=19, alpha=0.7) + 
  coord_fixed() +
  theme( panel.background =element_blank(),
         panel.border = element_rect(colour="black", fill=NA),
         panel.grid.major =element_blank(),  
         panel.grid.minor=element_blank(),
         axis.title=element_blank(),
         axis.text=element_blank(),
         axis.ticks=element_blank(),
         plot.title = element_text(hjust = 0.5),
         legend.title=element_blank()) + 
  ggtitle("Cell-Cell Attraction")


# inhibition
a1=split(d3)[[1]]; dat1=data.frame(x=a1$x, y=a1$y, type="Adipocyte")
a3=split(d3)[[3]]; dat3=data.frame(x=a3$x, y=a3$y, type="Epithelial Cell")
dat=rbind(dat1, dat3)
dist=crossdist(a1, a3)
nn3=apply(dist, 1, min)
p2= ggplot(dat, aes(x=x, y=y, col=type))+geom_point(shape=19, alpha=0.7) + 
  coord_fixed() +
  theme(panel.background =element_blank(),
         panel.border = element_rect(colour="black", fill=NA),
         panel.grid.major =element_blank(),  
         panel.grid.minor=element_blank(),
         axis.title=element_blank(),
         axis.text=element_blank(),
         axis.ticks=element_blank(),
         plot.title = element_text(hjust = 0.5),
         legend.title=element_blank()) + 
  ggtitle("Cell-Cell Inhibition")

# plot 3
a1=split(d2)[[1]]; a3=split(d2)[[3]];
dist=crossdist(a1, a3)
nn2=apply(dist, 1, min)
temp= rbind(data.frame(dist=nn1, type="Attraction"),
            data.frame(dist=nn2, type="No Effect"),
            data.frame(dist=nn3, type="Inhibition"))
temp$type <- factor(temp$type, levels=c("Attraction",  "No Effect", "Inhibition"))

p3=ggplot(temp, aes(x=type, y=dist, fill=type))+
  geom_boxplot()+
  labs(x=" ", y="Nearest Distance")+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", linewidth = rel(1)),
        panel.background = element_blank(),
        legend.title =element_blank(),
        legend.key = element_blank(),
        legend.position = "none") + 
  ggtitle(" ")

library(cowplot)
legend1 <- get_legend(p1)
plot_grid(p1+ theme(legend.position="none"), 
          p2+ theme(legend.position="none"), legend1,
          p3, rel_widths = c(1, 1, .45, 1), ncol=4)

ggsave("R1/Figures/fig5_a.pdf", width=9, height=3)

