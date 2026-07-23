
#The fast and dirty version of processing the Gemma DE files
#Megan Hagenauer 
#July 23 2026

#Note: by looping this instead of running it individually for each dataset the output will...
#have columns in the Log2FC,Tstat, SE, and SV output with stupid names (either very long or uninterpretable)
#But we can fix those later...
#The loop may also crash on some datasets, in which case just jump to the next dataset (iteration)


#First step, download the functions from the 2024 BDA github repository and put them in your working directory
https://github.com/hagenaue/BrainDataAlchemy/tree/main/MetaAnalysis_GemmaDEResults_2024

#source the functions from their files:

source("Function_FilteringDEResults_GoodAnnotation.R")

source("Function_ExtractingDEResultsForContrasts.R")

source("Function_CollapsingDEResults_OneResultPerGene.R")


#And then apply the functions to your differentials object:

YourWorkingDirectory<-getwd()
  
for(i in c(1:length(differentials)) ){

  print(i)
  
  DE_Results<-differentials[[i]]
  FilteringDEResults_GoodAnnotation(DE_Results)

  ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)
  
  CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)
  
  setwd(YourWorkingDirectory)
}


#Save your workspace!  (under R session)
#Save your code! (under file)
