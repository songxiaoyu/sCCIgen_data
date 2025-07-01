rm(list=ls())
library(tidyverse)
library(edgeR)

library(FNN)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(sCCIgen)
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_sCCIgen")

# # run sCCIgen to generate SeqFISH+ data. 
# run_interactive_sCCIgen()
input="Github/sCCIgen_data/sample_parameter_file/snRNAseq_Region1_n5990_cor_param.yml"
ParaSimulation(input=input)

##%######################################################%##
#                                                          #
####            Evaluation Metrics                      ####
#      1. Cell-type proportion                             #
#      2. Averaged gene count                              #
#      3. Variance of gene count                           #
#      4. Gene-gene correlation                            #
#                                                          #
##%######################################################%##

# load data
load("Github/sCCIgen_data/input_data/snRNAseq_breast_2025_expr.Rdata")
anno=colnames(expr)

# sCCIgen
s_count=read_tsv("R1_outputs/snRNAseq_Region1_n5990_cor_count_1.tsv") %>% column_to_rownames("GeneName") %>% as.data.frame()
s_meta=read_tsv("R1_outputs/snRNAseq_Region1_n5990_cor_meta_1.tsv")%>% column_to_rownames("Cell") %>% as.data.frame()
s_anno=s_meta$annotation



# Combine data
s_count[] <- lapply(s_count, function(x) as.numeric(x))

expr_data=list(expr, s_count)
anno_data=list(anno, s_meta$annotation)

# Metrics 1: barplot of cell type composition 

prop_matrix=sapply(anno_data, function(f) table(f)/length(f)) 
colnames(prop_matrix) <- c("Real", "Simulated")
prop_long <- as.data.frame(prop_matrix) %>%
  rownames_to_column("CellType") %>%
  pivot_longer(-CellType, names_to = "Input", values_to = "Proportion")

p1=ggplot(prop_long, aes(x = Input, y = Proportion*100, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Cell type Composition (%)", x = " ") +
  theme_minimal()+
  theme(legend.title = element_blank(),
        legend.position = "bottom")
p1

# Metrics 2: violin plot of mean log expression

mean_matrix=sapply(expr_data, function(f) apply(f, 1, mean))
colnames(mean_matrix) <- c("Real", "Simulated")
mean_long <- as.data.frame(mean_matrix) %>%
  pivot_longer(cols = everything(),
               names_to = "Input",
               values_to = "Expression"
  )
p2=ggplot(mean_long, aes(x = Input, y = log(Expression+1), fill= Input)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill= "white") +  # optional: add boxplot inside
  theme_minimal() +
  labs(x = " ", y = "Mean log expression")+
  theme(legend.position = "none")
p2

# Metrics 3: violin plot of var log expression

var_matrix=sapply(expr_data, function(f) apply(log(f+1), 1, var)) 
colnames(var_matrix)=c("Real", "Simulated" )
var_long = var_matrix%>% as.data.frame() %>% pivot_longer(cols = everything())

p3=ggplot(var_long, aes(x = name, y = log(value+1) , fill= name)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill= "white") +  # optional: add boxplot inside
  theme_minimal() +
  labs(x = " ", y = "Var log expression")+
  theme(legend.position = "none")

p3

# Metric 4: violin plot for gene-gene correlation 
hvg=mean_matrix[order(mean_matrix[,1]),1] %>% tail(., 200) %>% names(.)
hvg_data=sapply(expr_data, function(f) f[which(rownames(f) %in% hvg),])

cor_matrix=sapply(hvg_data, function(f) cor(t(f))) 
colnames(cor_matrix)=c("Real", "Simulated")
cor_matrix=cor_matrix[which(cor_matrix[,1]!=1 | cor_matrix[,2]!=1),]
cor_long = cor_matrix%>% as.data.frame() %>% pivot_longer(cols = everything())

p4=ggplot(cor_long, aes(x = name, y = value , fill= name)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill= "white") +  # optional: add boxplot inside
  theme_minimal() +
  labs( x = " ", y = "Correlation of HVG")+
  theme(legend.position = "none")
p4


library(cowplot)
plot_grid(p1, p2, p3, p4, nrow = 1, rel_widths = c(1.2, 1, 1, 1))
ggsave("R1/Figures/fig4_ab.pdf", width=10, height=3)