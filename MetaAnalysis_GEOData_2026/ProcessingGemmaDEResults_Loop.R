#The fast and dirty version of processing the Gemma DE files
#Megan Hagenauer 
#July 23 2026

#Note: by looping this instead of running it individually for each dataset the output will...
#have columns in the Log2FC,Tstat, SE, and SV output with stupid names (either very long or uninterpretable)
#But we can fix those later...
#The loop may also crash on some datasets, in which case just jump to the next dataset (iteration)

#####################

#Install and load the necessary code packages

#should already be installed from previous steps:

library(gemma.R)

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

library(tidyverse)

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

library(dplyr)

########################

#Next, download the functions from our repository and put them in your working directory:

#Function FilteringDEResults_GoodAnnotation:
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDEResults_2024/Function_FilteringDEResults_GoodAnnotation.R

#Function ExtractingDEResultsForContrasts (this one was updated/debugged for 2026):
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_ExtractingDEResultsForContrasts.R

#Function_CollapsingDEResults_OneResultPerGene:
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_CollapsingDEResults_OneResultPerGene.R

##########################

#Then source the functions from their files in your working directory:

source("Function_FilteringDEResults_GoodAnnotation.R")

source("Function_ExtractingDEResultsForContrasts.R")

source("Function_CollapsingDEResults_OneResultPerGene.R")

########################

#And then apply the functions to your differentials object:

YourWorkingDirectory<-getwd()
  
for(i in c(1:length(differentials)) ){

  print(i)
  
  DE_Results<-differentials[[i]]
  CurrentResultSet<-names(differentials)[i]
  
  FilteringDEResults_GoodAnnotation(DE_Results)

  ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts, CurrentResultSet)
  
  CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns, CurrentResultSet)
  
  setwd(YourWorkingDirectory)
}

########################

#Save your workspace!  (under R session)
#Save your code! (under file)
