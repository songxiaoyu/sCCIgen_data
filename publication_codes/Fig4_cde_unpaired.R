rm(list=ls())
library(sCCIgen)
library(tidyverse)
library(data.table)
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_sCCIgen")


#####################################

# run_interactive_sCCIgen()

input="Github/sCCIgen_data/sample_parameter_file/Unpaired_cor_param.yml"
ParaSimulation(input=input)

# Time difference of 39.89215 secs

# load real SeqFISH+ data -----------------
load("Github/sCCIgen_data/input_data/MERFISH_OV_2025_expr.Rdata")
load("Github/sCCIgen_data/input_data/MERFISH_OV_2025_spatial.Rdata")
dim(expr)
# [1] 550   248065
dim(spatial)
# [1] 209173   3
anno=colnames(expr)

# Also load simulated data ------------------
# 
expr2=fread(paste0("R1_outputs/Unpaired_cor_count_1.tsv")) %>% as.data.frame %>%
  column_to_rownames("GeneName")
spatial2=as.data.frame(fread(paste0("R1_outputs/Unpaired_cor_meta_1.tsv")))
dim(expr2)
# [1] 550   209173
dim(spatial2)
# [1] 209173   4
anno2=colnames(expr2)=spatial2$annotation
anno_data=list(colnames(expr), spatial$Labels, spatial2$annotation)
expr_data=list(expr, expr2)
########################
# Describe input  data  
########################


# Metrics 1: barplot of cell type composition 
prop_matrix=sapply(anno_data, function(f) table(f)/length(f)) 
colnames(prop_matrix) <- c("Expr. Ref", "Spatial Ref", "Simulated")
prop_long <- as.data.frame(prop_matrix) %>%
  rownames_to_column("CellType") %>%
  pivot_longer(-CellType, names_to = "Input", values_to = "Proportion")
prop_long$Input <- factor(prop_long$Input, levels = c("Expr. Ref", "Spatial Ref", "Simulated"))

p1=ggplot(prop_long, aes(x = Input, y = Proportion*100, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Cell type Composition (%)", x = " ") +
  theme_minimal()+
  theme(legend.title = element_blank(),
        legend.position = "bottom",      # Smaller text
        legend.key.size = unit(0.3, "cm"),            # Smaller boxes
        legend.spacing.y = unit(0.5, "cm"),           # Less vertical space
        plot.margin = margin(5, 5, 5, 5),   
        legend.margin = margin(t = -15, unit = "pt"),         # Tighter around plot
        axis.text.x = element_text(angle = 20, hjust = 1))+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
p1



########################
# spatial map   
########################
# spatial map ------------------
p22=ggplot(data=spatial, aes(x = center_x, y = center_y, colour=Labels))+
  geom_point(size=0.1)+
  ggtitle("Spatial Ref")+
  theme(legend.position = "bottom")+
  theme_void()+
  coord_fixed()+
  labs(colour = NULL)+
  guides(colour = guide_legend(override.aes = list(size = 2)))
p22

# Same as before 
p33=ggplot(data=spatial2, aes(x = x.loc, y = y.loc, colour=annotation))+
  geom_point(size=0.1)+
  ggtitle("Simulated")+
  theme(legend.position = "bottom")+
  theme_void()+
  coord_fixed()+
  labs(colour = NULL)+
  guides(colour = guide_legend(override.aes = list(size = 2)))
p33


# Metrics 2: violin plot of mean log expression

mean_matrix=sapply(expr_data, function(f) apply(f, 1, mean))
colnames(mean_matrix) <- c("Expr. Ref", "Simulated")
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


# 
# mean1_matrix=sapply(1:nrow(expr), function(f) tapply(expr[f,], anno, mean))
# mean2_matrix=sapply(1:nrow(expr2), function(f) tapply(expr2[f,], anno2, mean))
# 
# 
# 
# 
# mean_matrix=sapply(1:nrow(expr2), function(f) tapply(expr[f,], anno, mean))
# 
# colnames(mean_matrix) <- c("Real", "Simulated")
# mean_long <- as.data.frame(mean_matrix) %>%
#   pivot_longer(cols = everything(),
#                names_to = "Input",
#                values_to = "Expression"
#   )
# p2=ggplot(mean_long, aes(x = Input, y = log(Expression+1), fill= Input)) +
#   geom_violin(trim = FALSE) +
#   geom_boxplot(width = 0.1, outlier.shape = NA, fill= "white") +  # optional: add boxplot inside
#   theme_minimal() +
#   labs(x = " ", y = "Mean log expression")+
#   theme(legend.position = "none")
# p2

# Metrics 3: violin plot of var log expression

var_matrix=sapply(expr_data, function(f) apply(log(f+1), 1, var)) 
colnames(var_matrix)=c("Expr. Ref", "Simulated" )
var_long = var_matrix%>% as.data.frame() %>% pivot_longer(cols = everything())

p3=ggplot(var_long, aes(x = name, y = log(value+1) , fill= name)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill= "white") +  # optional: add boxplot inside
  theme_minimal() +
  labs(x = " ", y = "Var log expression")+
  theme(legend.position = "none")

p3

# Metric 4: violin plot for gene-gene correlation 

hvg=mean_matrix[order(mean_matrix[,1]),1] %>% tail(., 50) %>% names(.)
hvg_data=sapply(expr_data, function(f) f[which(rownames(f) %in% hvg), ])


cor_matrix=sapply(hvg_data, function(f) cor(t(f))) 
colnames(cor_matrix)=c("Expr. Ref", "Simulated")
cor_matrix=cor_matrix[which(cor_matrix[,1]!=1 | cor_matrix[,2]!=1),]
cor_long = cor_matrix%>% as.data.frame() %>% pivot_longer(cols = everything())

p4=ggplot(cor_long, aes(x = name, y = value , fill= name)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill= "white") +  # optional: add boxplot inside
  theme_minimal() +
  labs( x = " ", y = "Gene-gene correlation")+
  theme(legend.position = "none")
p4


library(cowplot)# First row with custom widths
row1 <- plot_grid(p22, p33, nrow = 1)

# Second row with default or custom widths
row2 <- plot_grid(p1, p2, p3, p4,nrow = 1, rel_widths = c(1.5, 1, 1, 1)  )

# Combine both rows
final_plot <- plot_grid(row2,row1, ncol = 1)
final_plot
ggsave("R1/Figures/fig4_def.pdf", width=10, height=6)

dev.off()

