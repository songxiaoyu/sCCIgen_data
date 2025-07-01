rm(list=ls())
library(tidyverse)
library(edgeR)
library(cowplot)
library(Giotto)
library(FNN)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(rlist)
library(sCCIgen)
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_sCCIgen")

# # run sCCIgen 511 cells to generate SeqFISH+ data. 

# run_interactive_sCCIgen()
# parameter files are accessible at https://github.com/songxiaoyu/sCCIgen_data/tree/main/sample_parameter_file.
input="Github/sCCIgen_data/sample_parameter_file/SeqFishPlus_cor_param.yml"
ParaSimulation(input=input)
input="Github/sCCIgen_data/sample_parameter_file/SeqFishPlus_cor_CCI2_param.yml"
ParaSimulation(input=input)
input="Github/sCCIgen_data/sample_parameter_file/SeqFishPlus_cor_CCI3_param.yml"
ParaSimulation(input=input)
input="Github/sCCIgen_data/sample_parameter_file/SeqFishPlus_cor_CCI23_param.yml"
ParaSimulation(input=input)

# load and clean data  ----------------
# Real data
load("Github/sCCIgen_data/input_data/SeqFishPlusCortex_2025_expr.Rdata")
load("Github/sCCIgen_data/input_data/SeqFishPlusCortex_2025_spatial.Rdata")
dim(expr)
dim(spatial)
colnames(spatial)= c("annotation", "x.loc", "y.loc")

# sCCIgen1
c_count=read_tsv("R1_outputs/SeqFishPlus_cor_count_1.tsv") %>% column_to_rownames("GeneName") %>% as.data.frame()
c_meta=read_tsv("R1_outputs/SeqFishPlus_cor_meta_1.tsv")%>% column_to_rownames("Cell") %>% as.data.frame()
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


##%######################################################%##
#                                                          #
####            Evaluate  Expression                    ####
#                                                          #
##%######################################################%##

# ------- Figure 3a: UMAP ---------------
cell_col=c("#F8766D","#9ecae1","#00BA38","#660099","#FFCC33","#FF61CC")
# sCCIgen2 CCIs
d_count=read_tsv("R1_outputs/SeqFishPlus_cor_CCI23_count_1.tsv") %>% column_to_rownames("GeneName") %>% as.data.frame()
d_meta=read_tsv("R1_outputs/SeqFishPlus_cor_CCI23_meta_1.tsv")%>% column_to_rownames("Cell") %>% as.data.frame()

data_meta=list(spatial, a_meta, b_meta, c_meta, d_meta)
data_expr=list(expr, a_count, b_count, c_count, d_count)

expr_matrix = do.call(cbind, data_expr)
meta_matrix = do.call(rbind, data_meta)

colnames(expr_matrix) = paste0("cell_", 1:ncol(expr_matrix))

method=rep(c("Ref", "SRTsim", "scDesign3", "sCCIgen w/o CCI", "sCCIgen w est. CCI"), times = sapply(data_meta,nrow))

moreMeta=cbind(anno=meta_matrix[,1], method=method)


# create Giotto object
dat <- createGiottoObject(expression = expr_matrix, spatial_locs = meta_matrix[,2:3])
dat <- addCellMetadata(dat,new_metadata = moreMeta)

# filter
dat <- filterGiotto(gobject = dat,
                    expression_threshold = 1,
                    feat_det_in_min_cells = 10,
                    min_det_feats_per_cell = 10,
                    expression_values = "raw",
                    verbose = F)

# normalize
dat <- normalizeGiotto(gobject = dat,
                       scalefactor = 6000,
                       verbose = F)

# add gene & cell statistics
dat <- addStatistics(gobject = dat)

# adjust expression matrix for technical or known variables
dat <- adjustGiottoMatrix(gobject = dat, 
                                 expression_values = "normalized",
                                 covariate_columns = c("nr_feats", "total_expr"),
                                 return_gobject = TRUE,
                                 name = "custom"
)

## highly variable features (HVF)
dat <- calculateHVF(gobject = dat, zscore_threshold=1)

## select genes based on highly variable features and gene statistics, both found in feature (gene) metadata
gene_metadata <- fDataDT(dat)
featgenes <- gene_metadata[hvf=="yes" & perc_cells > 4 & mean_expr_det > 0.5]$feat_ID

## run PCA on expression values (default)
dat <- runPCA(gobject = dat, 
              # feats_to_use = featgenes, 
              scale_unit = F, 
              center = F)
plotPCA(dat)



# Run UMAP (on PCA)
set.seed(123)
dat <- runUMAP(gobject = dat,
                        dimensions_to_use = 1:15,  # adjust based on your data
                        name = "umap")
# Step 3: Plot UMAP
plotUMAP(gobject = dat, cell_color = "anno")
plotUMAP(gobject = dat, cell_color = "method")
umap <- getDimReduction(gobject = dat, reduction_method = "umap")@coordinates
umap_dat=data.frame(moreMeta, umap) %>%
  mutate(method=factor(method, levels=c("Ref", "SRTsim", "scDesign3", "sCCIgen w/o CCI", "sCCIgen w est. CCI")))


pp3= ggplot(umap_dat, aes(x = Dim.1, y = Dim.2, color = anno)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ method, nrow = 3)+
  theme_bw() +coord_fixed()+
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
  ) +
  
  scale_color_manual(values = cell_col)  +
  theme(legend.position = c(0.8, 0.15),
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  guides(fill = guide_legend(ncol = 1))

pp3

ggsave("R1/Figures/fig3_a.pdf", width=4, height=6)

# 


# ------- Figure 3b: mean, variance, correlation metrics --------------
method_col=c("#1f78b4", "#7570b3", "#66a61e", "#e6ab02", "#e7298a")

# Metrics 1: violin plot of mean log expression

mean_matrix=sapply(data_expr, function(f) apply(f, 1, mean))
colnames(mean_matrix) <- c("Ref", "SRTSim", "scDesign3", "sCCIgen w/o CCI", "sCCIgen w. est. CCI")
mean_long <- as.data.frame(mean_matrix) %>%
  pivot_longer(cols = everything(),
               names_to = "Input",
               values_to = "Expression") %>%
  mutate(Input= factor(Input, levels=c(c("Ref", "SRTSim", "scDesign3", "sCCIgen w/o CCI", "sCCIgen w. est. CCI"))))

p2=ggplot(mean_long, aes(x = Input, y = log(Expression+1), fill= Input)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill= "white") +  # optional: add boxplot inside
  theme_minimal() +
  scale_fill_manual(values = method_col)  +
  labs(x = " ", y = "Mean log expression")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
p2

# Metrics 2: violin plot of var log expression

var_matrix=sapply(data_expr, function(f) apply(log(f+1), 1, var)) 
colnames(var_matrix) <- c("Ref", "SRTSim", "scDesign3", "sCCIgen w/o CCI", "sCCIgen w. est. CCI")

var_long = var_matrix%>% as.data.frame() %>% pivot_longer(cols = everything(),
                                                          names_to = "Input",
                                                          values_to = "Expression")%>%
  mutate(Input= factor(Input, levels=c(c("Ref", "SRTSim", "scDesign3", "sCCIgen w/o CCI", "sCCIgen w. est. CCI"))))

p3=ggplot(var_long, aes(x = Input, y = log(Expression+1) , fill= Input)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill= "white") +  # optional: add boxplot inside
  theme_minimal() +
  scale_fill_manual(values = method_col)  +
  labs(x = " ", y = "Var log expression")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

p3

# Metric 3: violin plot for gene-gene correlation 
hvg=mean_matrix[order(mean_matrix[,1]),1] %>% tail(., 200) %>% names(.)
hvg_data=sapply(data_expr, function(f) f[which(rownames(f) %in% hvg),])

cor_matrix=sapply(hvg_data, function(f) cor(t(f))) 
colnames(cor_matrix) <- c("Ref", "SRTSim", "scDesign3", "sCCIgen w/o CCI", "sCCIgen w. est. CCI")
cor_matrix=cor_matrix[which(cor_matrix[,1]!=1 | cor_matrix[,2]!=1),]
cor_long = cor_matrix%>% as.data.frame() %>% pivot_longer(cols = everything(),
                                                          names_to = "Input",
                                                          values_to = "Value")%>%
  mutate(Input= factor(Input, levels=c(c("Ref", "SRTSim", "scDesign3", "sCCIgen w/o CCI", "sCCIgen w. est. CCI"))))

p4=ggplot(cor_long, aes(x = Input, y = Value , fill= Input)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill= "white") +  # optional: add boxplot inside
  theme_minimal() +
  scale_fill_manual(values = method_col)  +
  labs( x = " ", y = "Correlation of HVG")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
p4


library(cowplot)
plot_grid( p2, p3, p4, ncol = 1)
ggsave("R1/Figures/fig3_b.pdf", width=4, height=7)

# ------- Figure 3c: Evaluate CCI2 in these datasets -------------------
d_count=read_tsv("R1_outputs/SeqFishPlus_cor_CCI2_count_1.tsv") %>% column_to_rownames("GeneName") %>% as.data.frame()
d_meta=read_tsv("R1_outputs/SeqFishPlus_cor_CCI2_meta_1.tsv")%>% column_to_rownames("Cell") %>% as.data.frame()


# Estimate
data_expr1=list(a_count, b_count, c_count, d_count)
data_meta1=list(a_meta, b_meta, c_meta, d_meta)

method=c("SRTsim", "scDesign3", "sCCIgen1", "sCCIgen2")
dis.cut=200
for (i in 1:4) {
  print(i)
  db=preprocessGiotto(expr_data=data_expr1[[i]], spatial_data=data_meta1[[i]], 
                      run_hvg=T, run_Dist_network=T,dis.cut=dis.cut) 
  ExprDistanceTable(gobject=db, in_hvg=T, abs_log2fc_ICG=0, p_adj = 1.1, spatial_network_name = "distance_based_network", 
                    seed=123, output_file=file.path(paste0("R1/Results/est_CCI_dist_expr_", method[i], ".csv")))    
}


# Compare CCI2
real_CCI2=read_csv("Github/sCCIgen_data/real_data_est/SeqFishPlus/est_CCI_dist_expr_disNet.csv", col_names =F)
a_CCI2=read_csv("R1/Results/est_CCI_dist_expr_SRTsim.csv", col_names =F)
b_CCI2=read_csv("R1/Results/est_CCI_dist_expr_scDesign3.csv", col_names =F)
c_CCI2=read_csv("R1/Results/est_CCI_dist_expr_sCCIgen1.csv", col_names =F)
d_CCI2=read_csv("R1/Results/est_CCI_dist_expr_sCCIgen2.csv", col_names =F)
colnames(real_CCI2)[7]=c("Real");colnames(a_CCI2)[7]=c("SRTsim");colnames(b_CCI2)[7]=c("scDesign3");colnames(c_CCI2)[7]=c("sCCIgen1");colnames(d_CCI2)[7]=c("sCCIgen2")

CCI2=real_CCI2[,c(2,3,5,7)]%>% left_join(a_CCI2[,c(2,3,5,7)])%>% left_join(b_CCI2[,c(2,3,5,7)])%>% 
  left_join(c_CCI2[,c(2,3,5,7)])%>% left_join(d_CCI2[,c(2,3,5,7)])

CCI2[is.na(CCI2)]=0
CCI2_long <- CCI2 %>% pivot_longer(cols = c(SRTsim, scDesign3, sCCIgen1, sCCIgen2), names_to = "Method", values_to = "Estimate")%>%
  mutate(Method = recode(Method,"sCCIgen1" = "sCCIgen w/o CCI","sCCIgen2" = "sCCIgen w. est. CCI"),
         Method = factor(Method, levels = c("SRTsim", "scDesign3", "sCCIgen w/o CCI", "sCCIgen w. est. CCI")))
cor_mse_labels <- CCI2_long %>%group_by(Method) %>%
  summarise(cor = cor(Real, Estimate, use = "complete.obs"), mse = mean((Real - Estimate)^2, na.rm = TRUE)) %>%
  mutate(label = paste0("r = ", round(cor, 2), "\nMSE = ", round(mse, 2)))
# Step 3: Plot with annotations
pp4=ggplot(CCI2_long, aes(x = Real, y = Estimate)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  facet_wrap(~ Method, nrow = 1) +
  geom_text(data = cor_mse_labels, 
            aes(x = -Inf, y = Inf, label = label), color="red",
            hjust = -0.1, vjust = 1.2, inherit.aes = FALSE, size = 4) +
  theme_bw(base_size = 14) +
  labs(title = "",
       x = "Ref. Effect Size",
       y = "Simu. Effect Size")+
  xlim(-3.3, 3)+
  ylim(-3.3, 3)

pp4

ggsave("R1/Figures/fig3_c.pdf", width=10, height=3)

# ------- Figure 3d: Evaluate CCI3 in these datasets -------------------
# Estimate
d_count=read_tsv("R1_outputs/SeqFishPlus_cor_CCI3_count_1.tsv") %>% column_to_rownames("GeneName") %>% as.data.frame()
d_meta=read_tsv("R1_outputs/SeqFishPlus_cor_CCI3_meta_1.tsv")%>% column_to_rownames("Cell") %>% as.data.frame()

data_meta1=list(a_meta, b_meta, c_meta, d_meta)
data_expr1=list(a_count, b_count, c_count, d_count)

method=c("SRTsim", "scDesign3", "sCCIgen1", "sCCIgen2")
spatial_network_name="distance_based_network"
dis.cut=200

for (i in 1:4) {
  print(i)
  db=preprocessGiotto(expr_data=data_expr1[[i]], spatial_data=data_meta1[[i]], run_hvg=F, run_Dist_network=T, dis.cut=dis.cut) 
  ExprExprTable(gobject=db, database="mouse", p_adj=1.01, abs_log2fc_LR=0, seed=123, spatial_network_name = spatial_network_name, 
                direction="both",
                output_file=file.path(paste0("R1/Results/est_CCI_expr_expr_", method[i], ".csv")))  
}

# Load results 
real_CCI3=read_csv("Github/sCCIgen_data/real_data_est/SeqFishPlus/est_CCI_expr_expr_disNet.csv", col_names =F)
a_CCI3=read_csv("R1/Results/est_CCI_expr_expr_SRTsim.csv", col_names =F)
b_CCI3=read_csv("R1/Results/est_CCI_expr_expr_scDesign3.csv", col_names =F)
c_CCI3=read_csv("R1/Results/est_CCI_expr_expr_sCCIgen1.csv", col_names =F)
d_CCI3=read_csv("R1/Results/est_CCI_expr_expr_sCCIgen2.csv", col_names =F)
colnames(real_CCI3)[9]=c("Real");colnames(a_CCI3)[9]=c("SRTsim");
colnames(b_CCI3)[9]=c("scDesign3");colnames(c_CCI3)[9]=c("sCCIgen1");colnames(d_CCI3)[9]=c("sCCIgen2")


# Compare CCI3 
CCI3=real_CCI3[,c(2,3,5,6,9)]%>% left_join(a_CCI3[,c(2,3,5,6,9)])%>% left_join(b_CCI3[,c(2,3,5,6,9)])%>% left_join(c_CCI3[,c(2,3,5,6,9)])%>% left_join(d_CCI3[,c(2,3,5,6,9)])
View(CCI3)
CCI3[is.na(CCI3)]=0
#CCI3=CCI3[which(CCI3$Real>0),]
CCI3_long <- CCI3 %>% pivot_longer(cols = c(SRTsim, scDesign3, sCCIgen1, sCCIgen2), names_to = "Method", values_to = "Estimate")%>%
  mutate(Method = recode(Method,"sCCIgen1" = "sCCIgen w/o CCI","sCCIgen2" = "sCCIgen w. est. CCI"),
         Method = factor(Method, levels = c("SRTsim", "scDesign3", "sCCIgen w/o CCI", "sCCIgen w. est. CCI")))
cor_mse_labels <- CCI3_long %>%group_by(Method) %>%
  summarise(cor = cor(Real, Estimate, use = "complete.obs"), mse = mean((Real - Estimate)^2, na.rm = TRUE)) %>%
  mutate(label = paste0("r = ", round(cor, 2), "\nMSE = ", round(mse, 2)))

# Step 3: Plot with annotations
pp5=ggplot(CCI3_long, aes(x = Real, y = Estimate)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  facet_wrap(~ Method, nrow = 1) +
  geom_text(data = cor_mse_labels, 
            aes(x = -Inf, y = Inf, label = label), color="red",
            hjust = -0.1, vjust = 1.2, inherit.aes = FALSE, size = 4) +
  theme_bw(base_size = 14) +
  labs(title = "",
       x = "Ref. Effect Size",
       y = "Simu. Effect Size")+
  ylim(-3.5, 1.5)+
  xlim(-3.5, 1.5)
pp5

ggsave("R1/Figures/fig3_d.pdf", width=10, height=3)

