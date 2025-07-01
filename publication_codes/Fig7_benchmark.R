
library(tidyverse)
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_sCCIgen")


run_interactive_sCCIgen()

input="Github/sCCIgen_data/sample_parameter_file/SeqFishPlus_n511_m10_cor_CCIs_param.yml"
ParaSimulation(input=input)



d_count=read_tsv("R1_outputs/SeqFishPlus_n511_cor_CCIs_count_1.tsv") %>% column_to_rownames("GeneName") %>% as.data.frame()
d_meta=read_tsv("R1_outputs/SeqFishPlus_n511_cor_CCIs_meta_1.tsv")%>% column_to_rownames("Cell") %>% as.data.frame()





