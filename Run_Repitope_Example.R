library(Biostrings)
library(tidyverse)
library(data.table)
library(dplyr)
library(tidyr)
options(java.parameters="-Xmx50G")  ## allow JAVA to use large memory space
#dyn.load('/Library/Java/JavaVirtualMachines/adoptopenjdk-8.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
require(Repitope)

warning("This script is unlikely to run as is, and serves as an example of how Repitope is trained and tested. It will require some tweaking")

# Prepare data for feature computation
HLA_SPECIFIC_DATA = readRDS("Analysis_HLA_specific_Pathogenic/A201_standardData_forAnalysis_BALANCED.rds")
HLA_SPECIFIC_DATA=HLA_SPECIFIC_DATA  %>% mutate(Dataset="UNBALANCED_WIDE_IEDB_standardData") %>% dplyr::rename(MHC=HLA_Allele) %>% select(Peptide,Immunogenicity,MHC,Dataset)
write.csv(HLA_SPECIFIC_DATA,row.names = FALSE,quote=FALSE,file = "GS_REPITOPEOUT.csv")
HLA_SPECIFIC_DATA= Epitope_Import(OtherFileNames = list("GS_REPITOPEOUT.csv"))
Data_ToCompute=HLA_SPECIFIC_DATA

# Compute fragmennt library, requires .fst file from Repitope repository.
fragLibDT <- CPP_FragmentLibrary(TCRSet_Public, fragLenSet=3:11, maxFragDepth=100000, seedSet=1:1)
fst::write_fst(fragLibDT, "FragmentLibrary.fst", compress=0)

# Compute Features [MHC-I]
featureDFList_MHCI <- Features(
  peptideSet=unique(c(Data_ToCompute$Peptide[])),
  fragLib="FragmentLibrary.fst",
  aaIndexIDSet="all",
  fragLenSet=3:8,
  fragDepth=10000,
  fragLibType="Weighted",
  seedSet=1:1,                                   ## must be the same random seeds used for preparing the fragment library
  coreN=parallel::detectCores()-4,        ## parallelization
  tmpDir="temp/"   ## where intermediate files are stored
)


saveFeatureDFList(featureDFList_MHCI, "HLA_Specific_")
featureDF_Combi = fst::read_fst("HLA_Specific_FeatureDF_Weighted.10000.fst", as.data.table=T)

# # Feature selection [MHC-I]
# minFeatureSet_MHCI_Human_HLA_SPECIFIC <- Features_MinimumFeatures(
#   featureDFList=list(featureDF_Combi[Peptide%in%Data_ToCompute$Peptide,]),
#   metadataDF=Data_ToCompute[,.(Peptide,Immunogenicity)][,Cluster:=.I],
#   seedSet=1:1,
#   corThreshold=0.75,
#   featureNSet=100,
#   criteria="intersect",
#   returnImpDFList=T
# )

# Train and test models based on folds.

FullDataset = HLA_SPECIFIC_DATA %>% as.data.table

folds = readRDS("Analysis_HLA_specific_Pathogenic/A201_standardData_Folds_BALANCED.rds")

featureDF_MHCI <- featureDF_Combi

for(i in 1:length(folds)) {
  print(c("Folds : ",i))
  # TRAIN AND TEST THE MODEL
  trainingData = FullDataset[-folds[[i]]]
  testData = FullDataset[folds[[i]]]
  # feature set if working, eitherwise or 'all' works fine.
  IEDBGS_TRAINED_MODEL=Immunogenicity_TrainModels(
    featureDF=featureDF_MHCI[Peptide%in%trainingData$Peptide,],
    metadataDF=trainingData[,.(Peptide, Immunogenicity)],
    featureSet=minFeatureSet_MHCI_Human_HLA_SPECIFIC,
    seedSet = 1:5,
    coreN = parallel::detectCores()-4
  )
  
  FOLD_CV_PREDICTIONS = Immunogenicity_Predict(list(featureDF_MHCI[Peptide%in%testData$Peptide,]),IEDBGS_TRAINED_MODEL)
  write.csv(FOLD_CV_PREDICTIONS$ScoreDT_1, paste0("FOLD_CV_PREDICTIONS_",i,"_data.csv"),quote=F,row.names = F)
}




