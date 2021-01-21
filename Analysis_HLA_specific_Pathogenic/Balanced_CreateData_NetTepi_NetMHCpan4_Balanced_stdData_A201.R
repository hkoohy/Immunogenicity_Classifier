# Script altered to account for multi alleles.
library(pROC)
library(ggpubr)
library(Biostrings)
library(data.table)
library(dplyr)
library(PepTools)
library(caret)
library(purrr)
library(readxl)
library(PepTools)
library(ggseqlogo)
library(DiffLogo)
library(tidyverse)
# Requires a working version of netMHCpan 4.0
 # setwd("../Cleaned_Analysis_Full/Analysis_HLA_specific_Pathogenic/")
setwd("../Analysis_HLA_specific_Pathogenic/")
options(scipen = 999)

# Function to provide a closest match. Used to match HLA Alleles across mixed output styles.
ClosestMatch2 = function(string, stringVector){
  
  stringVector[amatch(string, stringVector, maxDist=Inf)]
  
}

FullDataset= readRDS("A201_standardData_forAnalysis_BALANCED.rds")
FullDataset=FullDataset %>% as.data.table %>% ungroup

folds=readRDS(file = "A201_standardData_Folds_BALANCED.rds")

# Run NetMHCpan
TEST_DATA_LOCATION="NetTepi_netMHCpan4_Balanced_stdData_Analysis_A201/"
for(allele_i in 1:length(unique(FullDataset$HLA_Allele))){
  HLA_ALLELE_FOR_TESTING = gsub(x=unique(FullDataset$HLA_Allele)[allele_i],pattern=":",replacement = "")
  
  testdata=paste0(TEST_DATA_LOCATION,"Allele_",HLA_ALLELE_FOR_TESTING,"_NetMHC_test_data.txt")
  write.table(FullDataset %>% filter(HLA_Allele %in% unique(FullDataset$HLA_Allele)[allele_i]) %>% select(Peptide) %>% pull,file=testdata,sep="\n",col.names = F,row.names = F,quote=F)
  # Run model
  RESULTS_OUTPUT = paste0(TEST_DATA_LOCATION,"Allele_",HLA_ALLELE_FOR_TESTING,"_NetMHC","_Results.csv")
  
  system(paste0("/Applications/netMHCpan-4.0/netMHCpan -p ",testdata," -a ",unique(FullDataset$HLA_Allele)[allele_i]," -xls -xlsfile ", RESULTS_OUTPUT))
}

data_path <- "NetTepi_netMHCpan4_Balanced_stdData_Analysis_A201/"
files <- dir(data_path, pattern = "NetMHC_Results.csv") 

data3 <- data_frame(file = files) %>%
  mutate(file_contents = map(file,      
                             ~ fread(file.path(data_path, .))) 
  )

Netmhcpanres <- unnest(data3)

Netmhcpanres=Netmhcpanres %>% mutate(HLA_Allele = gsub(x=Netmhcpanres$file,pattern="Allele_|_NetMHC_Results.csv",replacement = ""))

Netmhcpanres$HLA_Allele = ClosestMatch2(Netmhcpanres$HLA_Allele,unique(FullDataset$HLA_Allele))
Netmhcpanres=Netmhcpanres %>% select(!file)%>% dplyr::rename(Score=`1-log50k`)


data_path <- "NETTEPI/"   
files <- dir(data_path, pattern = "_Results.csv") 

data3 <- data_frame(Fold = files) %>%
  mutate(file_contents = map(Fold,      
                             ~ fread(file.path(data_path, .))) 
  )

NetTepi_CV_RESULTS <- unnest(data3) %>% as.data.table  

NetTepi_CV_RESULTS=NetTepi_CV_RESULTS %>% dplyr::rename(HLA_Allele=Allele) 


NetTepi_CV_RESULTS$Fold= parse_number(str_extract(NetTepi_CV_RESULTS$Fold,"Fold_[:digit:]*"))
NetTepi_CV_RESULTS=NetTepi_CV_RESULTS %>% inner_join(FullDataset)

NetTepi_netMHCResults_FULL <- NetTepi_CV_RESULTS
NetTepi_netMHCResults_FULL = as.data.table(NetTepi_netMHCResults_FULL)

NetTepi_netMHCResults_FULL=NetTepi_netMHCResults_FULL %>% inner_join(Netmhcpanres)

NetTepi.MHCPAN.process_t_calculate_combi_score = function(t) {
  return(NetTepi_netMHCResults_FULL %>% mutate(CombinedScore = (t * Tcell + (1-t) * Score  )) %>% mutate(tWeight=t))
}


tvector=seq(0,1,by=0.001)

NetTepi_netMHCResults_FULL_T_EXPLORATION_DATA=lapply(tvector, NetTepi.MHCPAN.process_t_calculate_combi_score) %>% bind_rows

NetTepi_netMHCResults_FULL_T_EXPLORATION_DATA=NetTepi_netMHCResults_FULL_T_EXPLORATION_DATA %>% inner_join(FullDataset %>% select(Peptide,Immunogenicity))

weightsAnalysis=NetTepi_netMHCResults_FULL_T_EXPLORATION_DATA %>% group_by(tWeight) %>% dplyr::summarise(AUC=as.numeric(roc(Immunogenicity ~ CombinedScore)$auc)) %>% arrange(desc(AUC))

write.csv(weightsAnalysis,file="NetTepi_netMHCpan4_Balanced_stdData_Analysis_A201/NetTepi_netMHCResults_FULL_T_EXPLORATION_WEIGHTS.csv",quote=F,row.names = F)
NetTepi_netMHCResults_FULL_T_EXPLORATION_DATA=NetTepi_netMHCResults_FULL_T_EXPLORATION_DATA%>% filter(tWeight == weightsAnalysis$tWeight[1])
write.csv(NetTepi_netMHCResults_FULL_T_EXPLORATION_DATA,file="NetTepi_netMHCpan4_Balanced_stdData_Analysis_A201/NetTepi_netMHCResults_FULL_T_EXPLORATION_DATA.csv",quote=F,row.names = F)

