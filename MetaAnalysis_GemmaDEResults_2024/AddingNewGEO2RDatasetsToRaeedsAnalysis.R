#The goal of this code is to add the results from GEO2R for two antipsychotic datasets to Ra'eed's meta-analysis

#Challenges: 

#1) The GEO2R output only has Affymetrix probe id and Gene Symbol annotation
#### so we'll need to join by Gene Symbol and then add the results to the rest of our output (which is currently joined by Entrez ID)

#2) The GEO2R output is ordered by p-value, therefore even the output from the same dataset is in a different order for different group comparisons
#####so we'll need a join or merge function


######

#Before getting started:
#1) Open the workspace that contains your meta-analysis R objects/data
#2) Place all of the GEO2R differential expression results in the same folder
#### This means giving each output its own unambiguous name e.g., GSEblah_haloperidolvscontrol_lowdose.tsv
##### This should only be drug vs. control comparisons, no drug vs drug
#3) Set your working directory to where your files are located (with the GUI or with code) setwd()

#############

#We will eventually need to install and load these code libraries:
library(tidyverse)
library(dplyr)

########

#The first step is getting each differential expression output into the correct format 
#to add to the joined differential expression results (log2fc, sv) that are fed into the meta-analysis


#I'm adapting this code from recent work that we just published:
#https://github.com/hagenaue/NIDA_bLRvsbHR_F2Cross_HC_RNASeq/blob/main/NON_HC/GEO2R_Code_GSE88874.R
#This code is adapted from Brain Data Alchemy Project (v.2023 and 2024)

#This line would need to be adapted to reflect each of your datasets and group comparisons (file name, object name):
GSE88874_DE<-read.delim("GSE88874.top.table.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)

str(GSE88874_DE)
#for this analysis we didn't have gene symbol
#so I had to add it
#I have removed that code
#but to make downstream code work we'll need to rename the object


#########################

#Preparing the results to be fed into a meta-analysis:

#Making the column names match up with previous coding:

#This part will need to be applied to your new object
DE_Results<-GSE88874_DE
colnames(DE_Results)

#There won't be an NCBIid column
#Depending on which column # contains Gene Symbol, you will need to rename it:
#So 2 would need to be changed to the column # for Gene Symbol:
colnames(DE_Results)[2]<-"GeneSymbol"

#Adapting the Brain Data Alchemy Function (v.2024):

FilteringDEResults_GoodAnnotation<-function(DE_Results){
  
  print("# of rows in results")
  print(nrow(DE_Results))
  
  # print("# of rows with missing NCBI annotation:")
  # print(sum(DE_Results$NCBIid==""|DE_Results$NCBIid=="null"))
  
  # print("# of rows with NA NCBI annotation:")
  # print(sum(is.na(DE_Results$NCBIid)))
  
  print("# of rows with missing Gene Symbol annotation:")
  print(sum(DE_Results$GeneSymbol==""|DE_Results$GeneSymbol=="null"))
  
  # print("# of rows mapped to multiple NCBI_IDs:")
  # print(length(grep('\\|', DE_Results$NCBIid)))
  
  print("# of rows mapped to multiple Gene Symbols:")
  print(length(grep('\\|', DE_Results$GeneSymbol)))
  
  
  #I originally only wanted the subset of data which contains rows that do not contain an NCBI EntrezID of ""
  #I altered this code to be GeneSymbol instead of NCBIid
  DE_Results_NoNA<-DE_Results[(DE_Results$GeneSymbol==""|DE_Results$GeneSymbol=="null")==FALSE & is.na(DE_Results$GeneSymbol)==FALSE,]
  
  #I also originally only wanted the subset of data that is annotated with a single gene (not ambiguously mapped to more than one gene)
  #I altered this code to be GeneSymbol instead of NCBIid
  #This code block is probably irrelevant because we're working with GEO2R instead of Gemma output
  if(length(grep('\\|', DE_Results_NoNA$GeneSymbol))==0){
    DE_Results_GoodAnnotation<<-DE_Results_NoNA
  }else{
    #I only want rows annotated with a single Gene Symbol (no pipe):
    DE_Results_GoodAnnotation<<-DE_Results_NoNA[-(grep('\\|', DE_Results_NoNA$GeneSymbol)),]
  }
  #I used a double arrow in that conditional to place DE_Results_GoodAnnotation back out in the environment outside the function 
  
  print("# of rows with good annotation")
  print(nrow(DE_Results_GoodAnnotation))
  
  #For record keeping (sometimes useful for troubleshooting later)
  write.csv(DE_Results_GoodAnnotation, "DE_Results_GoodAnnotation.csv")
  
  rm(DE_Results_NoNA, DE_Results)
  
  print("Outputted object: DE_Results_GoodAnnotation")
}

FilteringDEResults_GoodAnnotation(DE_Results)


colnames(DE_Results_GoodAnnotation)

#We may need to change this code to indicate the column names for the Log2FC and t statistic in your object
# NamesOfFoldChangeColumns<-c("logFC")
# NamesOfTstatColumns<-c("t")

#We need to rename this to reflect the actual dataset and group comparison in this object
ComparisonsOfInterest<-c("GSEblah_HaloperidolvsControl_lowdose")

#We already double-checked direction of effect in the differential expression output and don't need to do that now.


#The Brain Data Alchemy (v.2024) function CollapsingDEResults_OneResultPerGene seems to not be working anymore due to issues with the select() function
#I reverted back to the simpler code from v.2023

#Find out the column names:
colnames(DE_Results_GoodAnnotation)

#This code will be dataset specific
#We need to extract the Log2FC ("Coef") and T-statistic ("t.") columns for the statistical contrasts relevant to our meta-analysis and place them into their own matrix
#We may need to change this code to indicate the column names for the Log2FC and t statistic in your object:
FoldChanges<-cbind(DE_Results_GoodAnnotation$logFC)
Tstats<-cbind(DE_Results_GoodAnnotation$t)

#Originally this code made the row names for the Log2FC and Tstat matrices the Entrez ID gene annotation
#Now we'll use GeneSymbol:
row.names(FoldChanges)<-DE_Results_GoodAnnotation$GeneSymbol
row.names(Tstats)<-DE_Results_GoodAnnotation$GeneSymbol

#Let's rename our columns to something nicer describing the effect of interest:
#Note - we later discovered that this name needs to include the dataset identifier (GSEID#) for later joining and plotting purposes
ComparisonsOfInterest<-c("GSE88874_HaloperidolVsControl_LowDose")
colnames(FoldChanges)<-ComparisonsOfInterest
colnames(Tstats)<-ComparisonsOfInterest

#############


#for this function to work with the current GEO2R output we'll need more than just GSE_ID in the GSE_ID object, e.g.
#GSE88874_HaloperidolVsControl_LowDose

CollapsingDEResults_OneResultPerGene<-function(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest){
  
  print("Double check that the vectors containing the fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:")
  
  # print("# of rows with unique NCBI IDs:")
  # print(length(unique(DE_Results_GoodAnnotation$NCBIid)))
  
  print("# of rows with unique Gene Symbols:")
  print(length(unique(DE_Results_GoodAnnotation$GeneSymbol)))
  
  #This makes a folder named after the dataset ID to store your results
  dir.create(paste("./", GSE_ID, sep=""))
  
  #And then sets the working directory to be that folder:
  setwd(paste("./", GSE_ID, sep=""))
  
  #Calculating the average Log2FC and T-stat for each gene:
  
  #Creating an empty list to store our results from each of our statistical contrasts:
  DE_Results_GoodAnnotation_FoldChange_Average<-list()
  DE_Results_GoodAnnotation_Tstat_Average<-list()
  DE_Results_GoodAnnotation_SE_Average<-list()
  
  #For each of the columns containing fold change and t-stat information for our statistical contrasts of interest:
  
  #I changed this code for this application because we only have one result per object (one log2fc and one tstat) and the select function isn't working anymore
  #That will eventually need to be more fully debugged, but not right now:
  
  #for(i in c(1:length(NamesOfFoldChangeColumns))){
    
    #Grab our Log2FC column of interest:
    FoldChangeColumn<-DE_Results_GoodAnnotation$logFC
    
    #Grab our Tstat column of interest:
    TstatColumn<-DE_Results_GoodAnnotation$t
    
    #Calculate the SE
    #Says FoldChangeColumn[[1]] and TstatColumn[[1]] - why?
    DE_Results_GoodAnnotation_SE<-FoldChangeColumn[[1]]/TstatColumn[[1]]
    
    #Calculate the average Log2FC per gene:
    DE_Results_GoodAnnotation_FoldChange_Average[[1]]<-tapply(FoldChangeColumn[[1]], DE_Results_GoodAnnotation$GeneSymbol, mean)
    
    #Calculate the average Tstat per gene:
    DE_Results_GoodAnnotation_Tstat_Average[[1]]<-tapply(TstatColumn[[1]], DE_Results_GoodAnnotation$GeneSymbol, mean)
    
    #Calculate the average SE per gene:
    DE_Results_GoodAnnotation_SE_Average[[1]]<-tapply(DE_Results_GoodAnnotation_SE, DE_Results_GoodAnnotation$GeneSymbol, mean)
    
  #}
  
  #Currently we have all of this averaged Log2FC information stored in a list, with one entry for each statistical contrast
  #Let's convert this to a data frame 
  DE_Results_GoodAnnotation_FoldChange_AveragedByGene<-do.call(cbind, DE_Results_GoodAnnotation_FoldChange_Average)
  
  print("Dimensions of Fold Change matrix, averaged by gene symbol:")
  print(dim(DE_Results_GoodAnnotation_FoldChange_AveragedByGene))
  
  #Name the columns in the dataframe in a manner that describes the dataset and factors for each statistic contrast
  colnames(DE_Results_GoodAnnotation_FoldChange_AveragedByGene)<-ComparisonsOfInterest
  
  #and then write it out to save it:
  write.csv(DE_Results_GoodAnnotation_FoldChange_AveragedByGene, "DE_Results_GoodAnnotation_FoldChange_AveragedByGene.csv")
  
  #And do the same for the T-stats:
  DE_Results_GoodAnnotation_Tstat_AveragedByGene<-do.call(cbind, DE_Results_GoodAnnotation_Tstat_Average)
  
  colnames(DE_Results_GoodAnnotation_Tstat_AveragedByGene)<-ComparisonsOfInterest
  
  write.csv(DE_Results_GoodAnnotation_Tstat_AveragedByGene, "DE_Results_GoodAnnotation_Tstat_AveragedByGene.csv")
  
  #And the same for the SE:
  DE_Results_GoodAnnotation_SE_AveragedByGene<-do.call(cbind, DE_Results_GoodAnnotation_SE_Average)
  
  colnames(DE_Results_GoodAnnotation_SE_AveragedByGene)<-ComparisonsOfInterest
  
  write.csv(DE_Results_GoodAnnotation_SE_AveragedByGene, "DE_Results_GoodAnnotation_SE_AveragedByGene.csv")
  
  #For running our meta-analysis, we are actually going to need the sampling variance instead of the standard error
  #The sampling variance is just the standard error squared.
  
  DE_Results_GoodAnnotation_SV<-(DE_Results_GoodAnnotation_SE_AveragedByGene)^2
  
  #Writing that information out to store it:
  write.csv(DE_Results_GoodAnnotation_SV, "DE_Results_GoodAnnotation_SV.csv")
  
  #Compiling all of our results into a single, big object 
  TempMasterResults<-list(Log2FC=DE_Results_GoodAnnotation_FoldChange_AveragedByGene, Tstat=DE_Results_GoodAnnotation_Tstat_AveragedByGene, SE=DE_Results_GoodAnnotation_SE_AveragedByGene, SV=DE_Results_GoodAnnotation_SV)
  
  #And then sending that object out into our environment outside the function
  #With a name that we can read and understand (DEResults for a particular GSEID)
  assign(paste("DEResults", GSE_ID, sep="_"), TempMasterResults, envir = as.environment(1))
  
  #Letting us know what that name is:
  print(paste("Output: Named DEResults", GSE_ID, sep="_"))
  
  #And then cleaning up our environment of temporary results:
  #Took out Names of fold change and tstat columns
  rm(TempMasterResults, DE_Results_GoodAnnotation, DE_Results_GoodAnnotation_SV, DE_Results_GoodAnnotation_SE, DE_Results_GoodAnnotation_FoldChange_AveragedByGene, DE_Results_GoodAnnotation_FoldChange_Average, DE_Results_GoodAnnotation_Tstat_AveragedByGene, DE_Results_GoodAnnotation_Tstat_Average, DE_Results_GoodAnnotation_SE_Average, FoldChangeColumn, TstatColumn, GSE_ID, ComparisonsOfInterest)
  
  #This sets the working directory back to the parent directory
  setwd("../")
  
}




###########################

#This function will pull out the columns containing the information that we are interested in

#for this function to work with the current GEO2R output we'll need more than just GSE_ID in the GSE_ID object, e.g.
#GSE88874_HaloperidolVsControl_LowDose

ExtractingDEResults<-function(GSE_ID, FoldChanges, Tstats){
  
  #We calculate the standard error by dividing the log2FC by the tstat
  StandardErrors<-FoldChanges/Tstats
  str(StandardErrors)
  
  #For running our meta-analysis, we are actually going to need the sampling variance instead of the standard error
  #The sampling variance is just the standard error squared.
  
  SamplingVars<-(StandardErrors)^2
  str(SamplingVars)
  
  TempMasterResults<-list(Log2FC=FoldChanges, Tstat=Tstats, SE=StandardErrors, SV=SamplingVars)
  
  assign(paste("DEResults", GSE_ID, sep="_"), TempMasterResults, envir = as.environment(1))
  
  print(paste("Output: Named DEResults", GSE_ID, sep="_"))
  
  rm(TempMasterResults, SamplingVars, StandardErrors, FoldChanges, Tstats)
  
}

ExtractingDEResults("GSE88874_HaloperidolVsControl_LowDose", FoldChanges, Tstats)

#... then go back to line 37 and run it for the next dataset and result outputs...
#Maybe just make a new script at that point
#Stay in the same workspace so that you are gradually accumulating all of the DEResults objects in your environment

