rm(list=ls())
library(tidyverse)
library(edgeR)
library(cowplot)

library(FNN)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(sCCIgen)
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_sCCIgen")

# # run sCCIgen 511 cells to generate SeqFISH+ data. 
# run_interactive_sCCIgen()
model_param_path="Github/sCCIgen_data/real_data_est/SeqFishPlus/SeqFishPlusCortex_fit_wo_cor.RDS"

# parameter files are accessible at https://github.com/songxiaoyu/sCCIgen_data/tree/main/sample_parameter_file.
input="Github/sCCIgen_data/sample_parameter_file/SRT/SeqFishPlus_n511_param.yml" 
ParaSimulation(input=input, ModelFitFile=model_param_path)
input="Github/sCCIgen_data/sample_parameter_file/SRT/SeqFishPlus_n511_CCIs_param.yml"
ParaSimulation(input=input, ModelFitFile=model_param_path)

# load and clean data  ----------------
# Real data
load("Github/sCCIgen_data/input_data/SeqFishPlusCortex_2025_expr.Rdata")
load("Github/sCCIgen_data/input_data/SeqFishPlusCortex_2025_spatial.Rdata")
dim(expr)
dim(spatial)
colnames(spatial)= c("annotation", "x.loc", "y.loc")
# sCCIgen
c_count=read_tsv("R1_outputs/SeqFishPlus_n511_count_1.tsv") %>% column_to_rownames("GeneName") %>% as.data.frame()
c_meta=read_tsv("R1_outputs/SeqFishPlus_n511_meta_1.tsv")%>% column_to_rownames("Cell") %>% as.data.frame()
HVG=rownames(c_count)

#SRTsim
SRTsim511 <- readRDS("R1/Shared/Weiping2/simulated_data/SRTsim_loc1.rds")
a_count1 <- SRTsim511@simCounts  
a_count = a_count1[which(rownames(a_count1) %in% HVG),]
a_meta = SRTsim511@simcolData
a_meta=a_meta[,c(3,1,2)]; colnames(a_meta)=c("annotation", "x.loc", "y.loc")
a_count=as.matrix(a_count); colnames(a_count)=paste0("Cell", 1:ncol(a_count)); a_meta=as.data.frame(a_meta)

# scDesign3
scDesign511 <- readRDS("R1/Shared/Weiping2/simulated_data/scDesign3_output_seqfish_n+1.rds")
b1_count=scDesign511$new_count
b_count = b1_count[which(rownames(b1_count) %in% HVG),]
b_meta=scDesign511$new_covariate
# [1] 512   4
b_meta=b_meta[,c(3,1,2)]; colnames(b_meta)=c("annotation", "x.loc", "y.loc")

# sCCIgen CCI
d_count=read_tsv("R1_outputs/SeqFishPlus_n511_CCIs_count_1.tsv") %>% column_to_rownames("GeneName") %>% as.data.frame()
d_meta=read_tsv("R1_outputs/SeqFishPlus_n511_CCIs_meta_1.tsv")%>% column_to_rownames("Cell") %>% as.data.frame()

data_meta=list(spatial, a_meta, b_meta, c_meta, d_meta)
data_expr=list(expr, a_count, b_count, c_count, d_count)

##%######################################################%##
#                                                          #
####            Evaluate  Spatial                       ####
#                                                          #
##%######################################################%##


# ------------------ Figure 2a: Simulated spatial map -------------------------
cell_col=c("#F8766D","#9ecae1","#00BA38","#660099","#FFCC33","#FF61CC")
ylim=range(b_meta$y.loc)
p0=ggplot(spatial, aes(x=x.loc, y=y.loc, colour = annotation) ) + geom_point(size = 0.7) +
  coord_fixed() +theme_minimal()+labs(colour = " ") +
  theme(legend.key.size = unit(0.7, "lines"))+
  guides(colour = guide_legend(override.aes = list(size = 2)))+
  labs(colour = " ", x = NULL, y = NULL)+
  scale_color_manual(values = cell_col)  +
  xlim(xlim)+
  ylim(ylim)
legend=get_legend(p0)
p1=p0+theme(legend.position = "none") +
  ggtitle("Ref Data")
p1
p2=ggplot(a_meta, aes(x=x.loc, y=y.loc, colour = annotation)) + geom_point(size = 0.7) +
  coord_fixed() +theme_minimal()+
  theme(legend.position = "none") +
  ggtitle("SRTsim")+
  scale_color_manual(values = cell_col)  +
  xlim(xlim)+
  ylim(ylim)+  labs(colour = " ", x = NULL, y = NULL)
p2
p3=ggplot(b_meta, aes(x=x.loc, y=y.loc, colour = annotation)) + geom_point(size = 0.7) +
  coord_fixed() +theme_minimal()+
  theme(legend.position = "none") +
  ggtitle("scDesign3")+
  scale_color_manual(values = cell_col)  +
  xlim(xlim)+
  ylim(ylim)+  labs(colour = " ", x = NULL, y = NULL)
p3
p4=ggplot(c_meta, aes(x=x.loc, y=y.loc, colour = annotation) ) + geom_point(size = 0.7) +
  coord_fixed() +theme_minimal()+
  theme(legend.position = "none") +
  ggtitle("sCCIgen w/o CCI")+
  xlim(xlim)+
  scale_color_manual(values = cell_col)  +
  ylim(ylim)+  labs(colour = " ", x = NULL, y = NULL)
p4
p5=ggplot(d_meta, aes(x=x.loc, y=y.loc, colour = annotation)) + geom_point(size = 0.7) +
  coord_fixed() +theme_minimal()+
  theme(legend.position = "none") +
  ggtitle("sCCIgen w. est. CCI")+
  scale_color_manual(values = cell_col)  +
  xlim(xlim)+
  ylim(ylim)+  labs(colour = " ", x = NULL, y = NULL)

p5


library(cowplot)# First row with custom widths
row1 <- plot_grid(p1, p2,p3,  nrow = 1)
row2 <- plot_grid(p4,  p5,legend, nrow = 1)
legend_row <- plot_grid(legend, nrow = 1)

# Combine all rows into one figure
final_plot <- plot_grid(
  row1,row2,
  ncol = 1,
  rel_heights = c(1,1)  # Adjust height ratio for legend
)

final_plot
ggsave("R1/Figures/fig2_a.pdf", width=10, height=4)



# -------- Figure 2b: Region overlaps-------------------
method_col=c("#1f78b4", "#7570b3", "#66a61e", "#e6ab02", "#e7298a")
# Estimate region of the ref data
library(spatstat)
w=lapply(data_meta, function(f) simu.window(PointLoc=f[,2:3],method="network") )
area=sapply(w, function(f) area.owin(f))

area_ref=area[2:5]/area[1]

overlaps=sapply(w, function(f) area.owin(intersect.owin(w[[1]], f))/area[1])[-1]
uniqe_ref=sapply(w, function(f) area.owin(setminus.owin(w[[1]], f))/area[1])[-1]
uniqe_sim= sapply(w, function(f) area.owin(setminus.owin(f,w[[1]]))/area[1])[-1]
dat=rbind(area_ref, overlaps, uniqe_ref, uniqe_sim)
colnames(dat)=c("SRTsim", "scDesign3", "sCCIgen w/o CCI", "sCCIgen w. est. CCI")  
dat=dat%>% as.data.frame() %>% rownames_to_column("Metric")

dat_long <- pivot_longer(dat, cols=-Metric, names_to = "Method", values_to = "Area (%)")%>%
  mutate(Method = factor(Method, levels = c("SRTsim", "scDesign3", "sCCIgen w/o CCI", "sCCIgen w. est. CCI")),
         Metric = factor(Metric,
                         levels = c("area_ref", "overlaps", "uniqe_ref", "uniqe_sim"),
                         labels = c("Coverage Area (Ideal 100)", "Overlap Area (Ideal 100) ", " Ref Unique Area (Ideal 0)", "Sim. Unique Area (Ideal 0)")))

pp3= # Plot with facet_wrap by Metric
  ggplot(dat_long, aes(x = Method, y = `Area (%)` * 100, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(`Area (%)` * 100, 1)),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  facet_wrap(~ Metric, nrow = 1,scales = "free_y") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = method_col[-1])  +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  labs(
    x = "Spatial Window Comparison",
    y = "Relative to Ref (%)",
    fill = NULL,
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

pp3
ggsave("R1/Figures/fig2_b.pdf", width=10, height=3.5)


# -------- Figure 2c: Cell type composition-------------------
prop=sapply(data_meta, function(f) table(f$annotation)/length(f$annotation)) %>% as.data.frame() %>%
  rownames_to_column("annotation")

RB_prop=data.frame(anno=prop$annotation, (prop[3:6]-prop[,2])/prop[,2])
colnames(RB_prop)[-1]=c("SRTsim", "scDesign3", "sCCIgen w/o CCI", "sCCIgen w. est. CCI")

AB_prop=data.frame(anno=prop$annotation, (prop[3:6]-prop[,2]))
colnames(AB_prop)[-1]=c("SRTsim", "scDesign3", "sCCIgen w/o CCI", "sCCIgen w. est. CCI")


RB_long <- pivot_longer(RB_prop, cols=-anno, names_to = "Method", values_to = "Proportion")
AB_long <- pivot_longer(AB_prop,cols=-anno,names_to = "Method",values_to = "Proportion")
BB_long=rbind( data.frame(Bias="Relative Bias", RB_long),
               data.frame(Bias="Absolute Bias", AB_long))
q1= ggplot(BB_long, aes(x = Method, y = Proportion*100, fill = anno)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  facet_wrap(~ Bias, scales = "free_x")+
  # geom_text(aes(label = round(value*100, 1)), 
  #           position = position_dodge(width = 0.8), 
  #           hjust = 1, size = 3.5, fontface = "bold") + 
  theme_minimal(base_size = 14) +
  coord_flip()+ 
  scale_fill_manual(values = cell_col)  +
  labs(
    x = " ",
    y = "Bias in Cell Type Proportion (%)",
    fill = "Interaction"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.title = element_blank()
  )+
  guides(fill = guide_legend(ncol = 1))

q1

ggsave("R1/Figures/fig2_c.pdf", width=10, height=3)


# ------------ Figure 2d: Estimate CCI1 in these datasets -------------------
library(Giotto)
real_CCI1=read_csv("Github/sCCIgen_data/real_data_est/SeqFishPlus/est_CCI_dist_dist_disNet_unfiltered.csv")
sig_CCI1=read_csv("Github/sCCIgen_data/real_data_est/SeqFishPlus/est_CCI_dist_dist_disNet.csv", col_names =F)
dis.cut=200
e_r=data.frame(unified_int=real_CCI1$unified_int, enrichm_r=real_CCI1$enrichm)
p_r=data.frame(unified_int=real_CCI1$unified_int, p_r=pmin(real_CCI1$p.adj_higher, real_CCI1$p.adj_lower))

spatial_network_name="distance_based_network"

# SRTsim
db=preprocessGiotto(expr_data=a_count, spatial_data=a_meta, run_Dist_network=T, dis.cut=dis.cut) 
a_CCI1=cellProximityEnrichment(gobject = db, cluster_column = "anno",
                               spatial_network_name = spatial_network_name,adjust_method = "fdr",
                               number_of_simulations = 2000) 
e_a=data.frame(unified_int=a_CCI1$enrichm_res$unified_int, enrichm_a=a_CCI1$enrichm_res$enrichm)
p_a=data.frame(unified_int=a_CCI1$enrichm_res$unified_int, p_a=pmin(a_CCI1$enrichm_res$p.adj_higher, a_CCI1$enrichm_res$p.adj_lower))
# scDesign3
db=preprocessGiotto(expr_data=b_count, spatial_data=b_meta, run_Dist_network=T, dis.cut=dis.cut) 
b_CCI1=cellProximityEnrichment(gobject = db,cluster_column = "anno",spatial_network_name = spatial_network_name,adjust_method = "fdr",number_of_simulations = 2000)
e_b=data.frame(unified_int=b_CCI1$enrichm_res$unified_int, enrichm_b=b_CCI1$enrichm_res$enrichm)
p_b=data.frame(unified_int=b_CCI1$enrichm_res$unified_int, p_b=pmin(b_CCI1$enrichm_res$p.adj_higher, b_CCI1$enrichm_res$p.adj_lower))
# sCCIgen1
db=preprocessGiotto(expr_data=c_count, spatial_data=c_meta, run_Dist_network=T, dis.cut=dis.cut) 
c_CCI1=cellProximityEnrichment(gobject = db,cluster_column = "anno",spatial_network_name = spatial_network_name,adjust_method = "fdr",number_of_simulations = 2000)
e_c=data.frame(unified_int=c_CCI1$enrichm_res$unified_int, enrichm_c=c_CCI1$enrichm_res$enrichm)
p_c=data.frame(unified_int=c_CCI1$enrichm_res$unified_int, p_c=pmin(c_CCI1$enrichm_res$p.adj_higher, c_CCI1$enrichm_res$p.adj_lower))

# sCCIgen2
db=preprocessGiotto(expr_data=d_count, spatial_data=d_meta, run_Dist_network=T, dis.cut=dis.cut) 
d_CCI1=cellProximityEnrichment(gobject = db,cluster_column = "anno",
                               spatial_network_name =spatial_network_name,
                               adjust_method = "fdr",number_of_simulations = 2000)
e_d=data.frame(unified_int=d_CCI1$enrichm_res$unified_int, enrichm_d=d_CCI1$enrichm_res$enrichm)
p_d=data.frame(unified_int=d_CCI1$enrichm_res$unified_int, p_d=pmin(d_CCI1$enrichm_res$p.adj_higher, d_CCI1$enrichm_res$p.adj_lower))


# Relative and absolute bias
enrich=e_r%>% full_join(e_a)%>% full_join(e_b)%>% full_join(e_c)%>% full_join(e_d)
pvalue=p_a %>% full_join(p_b)%>% full_join(p_c)%>% full_join(p_d)     
log10p=data.frame(unified_int=pvalue$unified_int, log10p=-log10(pvalue[,-1]+0.001))


dat=enrich %>% filter(unified_int %in% paste0(sig_CCI1$X1, "--", sig_CCI1$X2))
log10p=log10p %>% filter(unified_int %in% paste0(sig_CCI1$X1, "--", sig_CCI1$X2))
dat
AB=data.frame(unified_int=dat$unified_int, (dat[,3:6]-dat[,2]))
RB=data.frame(unified_int=dat$unified_int, (dat[,3:6]-dat[,2])/dat[,2])

RB_long=RB %>%pivot_longer(cols = !"unified_int")%>%
  mutate(name = recode(name,
                       "enrichm_a" = "SRTsim",
                       "enrichm_b" = "scDesign3",
                       "enrichm_c" = "sCCIgen w/o CCI",
                       "enrichm_d" = "sCCIgen w. est. CCI"),
         name = factor(name, levels = c("sCCIgen w. est. CCI",
                                        "sCCIgen w/o CCI", 
                                        "scDesign3", "SRTsim")))
AB_long=AB %>%pivot_longer(cols = !"unified_int")%>%
  mutate(name = recode(name,
                       "enrichm_a" = "SRTsim",
                       "enrichm_b" = "scDesign3",
                       "enrichm_c" = "sCCIgen w/o CCI",
                       "enrichm_d" = "sCCIgen w. est. CCI"),
         name = factor(name, levels = c("sCCIgen w. est. CCI",
                                        "sCCIgen w/o CCI", 
                                        "scDesign3", "SRTsim")))

BB_long= rbind(data.frame(Bias="Relative Bias", RB_long),
               data.frame(Bias="Absolute Bias", AB_long))

pp0= ggplot(BB_long[which(BB_long$Bias=="Absolute Bias" & BB_long$unified_int=="Oligodendrocyte--Oligodendrocyte"),], 
            aes(x = name, y = value*100, fill = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_minimal(base_size = 12) +
  coord_flip()+ 
  labs(x = "", title="Oligodendrocyte--Oligodendrocyte Colocalization",
       y = " Absolute Bias (%)"
  ) +
  
  scale_fill_manual(values = method_col[c(5,4,3,2)])  +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1), plot.title.position = "plot",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(ncol = 1))
legend=get_legend(pp0)
pp1=pp0+theme(legend.position = "none")
pp1

pp2= ggplot(BB_long[which(BB_long$Bias=="Relative Bias" & BB_long$unified_int=="Oligodendrocyte--Oligodendrocyte"),], 
            aes(x = name, y = value*100, fill = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_minimal(base_size = 12) +
  coord_flip()+ 
  labs(x = "", y = " Relative Bias (%)") +
  scale_fill_manual(values = method_col[c(5,4,3,2)])  +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.title = element_blank(),
    plot.title.position = "plot",
    legend.position = "none")+
  guides(fill = guide_legend(ncol = 1))

pp2

pp3= ggplot(BB_long[which(BB_long$Bias=="Absolute Bias" & BB_long$unified_int=="ExcitatoryNeuron--Oligodendrocyte"),], 
            aes(x = name, y = value*100, fill = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_minimal(base_size = 12) +
  coord_flip()+ 
  scale_fill_manual(values = method_col[c(5,4,3,2)])  +
  labs(x = "", title="ExcitatoryNeuron--Oligodendrocyte Colocalization", y = " Absolute Bias (%)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title.position = "plot",
        legend.title = element_blank(),legend.position = "none",
        plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(ncol = 1))

pp3

pp4= ggplot(BB_long[which(BB_long$Bias=="Relative Bias" & BB_long$unified_int=="ExcitatoryNeuron--Oligodendrocyte"),], 
            aes(x = name, y = value*100, fill = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = method_col[c(5,4,3,2)])  +
  coord_flip()+ 
  labs(x = "", y = " Relative Bias (%)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        legend.title = element_blank(),
        plot.title.position = "plot",
        legend.position = "none")+
  guides(fill = guide_legend(ncol = 1))

pp4


library(cowplot)# First row with custom widths
col1 <- plot_grid(pp1, pp2,  ncol = 1, rel_heights = c(1.2, 1))
col2 <- plot_grid(pp3, pp4, ncol = 1, rel_heights = c(1.2, 1))

# Combine all rows into one figure
pp_final <- plot_grid(
  col1,col2, 
  nrow = 1,
  rel_widths = c(1,1)  # Adjust height ratio for legend
)

pp_final
ggsave("R1/Figures/fig2_d.pdf", width=12, height=3)


