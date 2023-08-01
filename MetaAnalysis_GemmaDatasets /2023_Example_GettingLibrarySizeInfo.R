#Adding RNA-Seq Library Size to the summarized experiment for a dataset
#Megan Hagenauer
#July 20, 2023

#Set working directory:
setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Teaching/BrainDataAlchemy")

#load useful code packages
library(SummarizedExperiment)
library(gemma.R)
library(plyr)
#library(tidyr)

#Input the summarized experiment object for your dataset:
SummarizedExperiment_Filtered<-gemma.R::get_dataset_object("GSE81672", type = 'se', filter=TRUE, consolidate="average")
SummarizedExperiment_Filtered
# $`13458`
# class: SummarizedExperiment 
# dim: 22782 99 
# metadata(8): title abstract ... GemmaSuitabilityScore taxon
# assays(1): counts
# rownames(22782): 100009600 100017 ... Averaged from 100039542 100993 Averaged from 11641 677884
# rowData names(4): Probe GeneSymbol GeneName NCBIid
# colnames(99): Sample 18: AMY_resilient_saline Sample 20: AMY_susceptible_ketamine_responder ...
# Sample 108: HIP_susceptible_saline Sample 112: AMY_susceptible_imipramine_non_responder
# colData names(5): factorValues organism part block phenotype treatment

#Library Size:
#Within the Gemma Diagnostics tab, download the MultiQC report
#Within the MultiQC report, choose "copy table"
#Paste the table within an excel file and save as a .csv file in your working directory
GSE81672_LibrarySize<-read.csv("GSE81672_LibrarySize.csv", header=TRUE, stringsAsFactors = FALSE)

#The library size information is annotated with sample information using GSM number. 
head(GSE81672_LibrarySize)
# ExternalID X..Aligned M.Aligned X..Dups X..GC M.Seqs
# 1 GSM2166101     85.00%      20.1                   NA
# 2 GSM2166102     86.10%      33.4                   NA

#As far as I can tell, this identifier is not in the Summarized experiment object
#e.g., 
colData(SummarizedExperiment_Filtered[[1]])
#Sample names are given in the format
#"Sample 18: AMY_resilient_saline"
row.names(colData(SummarizedExperiment_Filtered[[1]]))
# [1] "Sample 18: AMY_resilient_saline"                     
# [2] "Sample 20: AMY_susceptible_ketamine_responder"       
# [3] "Sample 8: NAC_resilient_saline"   

#So we need to read in the experimental design info from the Gemma website, which contains GSM#
#To do this, go to the Experimental Design tab -> "show details"->"download design file"
#It is a tab-delimited text file
GSE81672_expdesign<-read.delim("13458_GSE81672_expdesign.data.txt", sep="\t", comment.char="#", header=TRUE, stringsAsFactors = FALSE)
str(GSE81672_expdesign)
#'data.frame':	99 obs. of  6 variables:

#This column has the GSM #
GSE81672_expdesign$ExternalID
# [1] "GSM2166189" "GSM2166195" "GSM2166180"

#This column includes the Name... and a bunch of other stuff:
GSE81672_expdesign$Bioassay
# [1] "GSE81672_Biomat_64___BioAssayId=450411Name=Sample102.NAC_susceptible_ketamine_non_responder"  
# [2] "GSE81672_Biomat_58___BioAssayId=450412Name=Sample108.HIP_susceptible_saline"
#This splits these long names up at the "Name=" and then combines things back into a matrix again
temp<-do.call(rbind.data.frame, strsplit(GSE81672_expdesign$Bioassay, "Name="))
#This adds the name information to the expdesign matrix:
GSE81672_expdesign$SampleName<-temp[,2]
GSE81672_expdesign$SampleName

#For this dataset, The formatting is still a little different from the Summarized Experiment object, e.g.
GSE81672_expdesign$SampleName[1]
#"Sample102.NAC_susceptible_ketamine_non_responder"
row.names(colData(SummarizedExperiment_Filtered[[1]]))[1]
#"Sample 18: AMY_resilient_saline"
#One version uses a ., the other uses a ": "

#the ": " version is straight off of GEO, so I'm going to assume that for some reason Gemma doesn't like : in the design data.frame
#This is likely to be a dataset specific issue - I'm not sure how to generalize this code.
#For many datasets, things may be fine at this point.
#So we essentially need a find and replace - in R, these are often done with sub() or gsub()
#https://www.digitalocean.com/community/tutorials/sub-and-gsub-function-r
#https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/grep
#this is a little funky, because we're not replacing simple alphanumeric characters, but punctuation that also has a meaning for coding
temp<-gsub("[.]", ": ", GSE81672_expdesign$SampleName)

#But it turns out that the names also differ because of an extra space. 
#To take care fo this we just removed all spaces from both vectors because it is easier
GSE81672_expdesign$SampleName_toJoin<-gsub(" ", "", temp)
SampleName_toJoin<-gsub(" ", "", row.names(colData(SummarizedExperiment_Filtered[[1]])))

#We can't just add library size as a column to experimental design because the samples aren't in the same order
#We have to use a "join" or "merge" function to align them instead
#To join the design matrix with the library size, we need to have two data.frames that have columns that have the same name:
str(GSE81672_LibrarySize)
#In this data.frame, the GSM# is called Sample.Name, whereas in the exp. design data.frame it is called ExternalID
colnames(GSE81672_LibrarySize)[1]<-"ExternalID"

#Now we can join by these columns (the function "merge" also works for this):
GSE81672_expdesign_wLibrarySize<-join(GSE81672_expdesign, GSE81672_LibrarySize, by="ExternalID", type="left")
str(GSE81672_expdesign_wLibrarySize)


#Now we have to add this information into the Summarized Experiment object
#Which also has a different sample order
#I'm going to grab that sample order information and name it after the column with Sample Names in the experimental design object

SamplesOrderedLikeSummarizedExperiment<-data.frame(SampleName_toJoin=SampleName_toJoin)
                                
GSE81672_expdesign_wLibrarySize_Ordered<-join(SamplesOrderedLikeSummarizedExperiment, GSE81672_expdesign_wLibrarySize, by="SampleName_toJoin", type="left")                   
str(GSE81672_expdesign_wLibrarySize_Ordered)
#Double checking that things are actually in the same order
cbind(row.names(colData(SummarizedExperiment_Filtered[[1]])), GSE81672_expdesign_wLibrarySize_Ordered$SampleName_toJoin)

#Adding library size to the Summarized Experiment object
colData(SummarizedExperiment_Filtered[[1]])$LibrarySize<-GSE81672_expdesign_wLibrarySize_Ordered$M.Aligned
