#Example code showing the pipeline for running a basic meta-analysis of Log2FC and sampling variance values using our previously generated objects MetaAnalysis_FoldChanges & MetaAnalysis_SV
#Megan Hagenauer
#Original version: July 25 2024
#In response to reviewers' comments, this function has been updated to include heterogeneity statistics, publication bias statistics, and robustness statistics
#Updated version: March 10, 2026
#Updated again: July 28, 2026 (to make the code more generalizable for the full 2026 cohort)

######################

#Grabbing input data and setting the working directory:

######################

#Installing and loading relevant code packages:

if (!require("metafor", quietly = TRUE)){
  install.packages("metafor")
}

library(metafor)

if (!require("plyr", quietly = TRUE)){
  install.packages("plyr")
}

library(plyr)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("multtest")                                 

library(multtest)

######################

#Download the two necessary functions from our Github repository:

#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_RunBasicMetaAnalysis.R

#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_FalseDiscoveryCorrection.R

######################

#Read in the functions:

source("Function_RunBasicMetaAnalysis.R")

source("Function_FalseDiscoveryCorrection.R")

######################

#Example Usage:

#Figure out which columns contain differential expression output:
colnames(MetaAnalysis_FoldChanges)

#Then make a variable that is a vector of with the column numbers containing the DE output:
#Example code for making a variable representing the columns 
Columns_DE<-c(29:34)

#Then make a variable that is the shorthand annotation that we will use to label the genes in our output, e.g. MouseRat_GeneSymbol:
Column_GeneName<-35

NumberOfComparisons=6
CutOffForNAs=2
#I have 6 comparisons
#2 NA is too many

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV, Columns_DE, Column_GeneName)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

#Take a peek at the output:

str(metaOutput)
head(metaOutput)
tail(metaOutput)

write.csv(metaOutput, "metaOutput_wHeterogeneityPubBiasRobustMeasures.csv")
write.csv(MetaAnalysis_Annotation, "MetaAnalysis_Annotation_for_metaOutput_wHeterogeneityPubBiasRobustMeasures.csv")

colnames(metaOutput)

write.csv(influence_dfbs, "influence_dfbs.csv")
write.csv(influence_cookd, "influence_cookd.csv" )
write.csv(influence_TF, "influence_TF.csv")

###############

#FDR correction:
  
FalseDiscoveryCorrection(metaOutput, MetaAnalysis_Annotation)

