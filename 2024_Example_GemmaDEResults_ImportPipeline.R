#Code for processing the differential expression results that we can automatically pipe in from Gemma in preparation for running the meta-analysis

#Goals:
#1. For each of the datasets that have survived our inclusion/exclusion criteria, we will need to identify the relevant result sets and statistical contrasts that are relevant to our research question
#2. We will download the differential expression results for these statistical contrasts
#Consolidate them to one result per gene
#Then pulling out the useful info, calculating sampling variance
#And putting everything in a format that we can easily use in an effect size meta-analysis

#Megan Hagenauer
#July 13, 2024

###############

#Some of the code is adapted from Ogan's Gemma.R vignette
#https://bioconductor.org/packages/release/bioc/vignettes/gemma.R/inst/doc/metanalysis.html

#################

#I tested out this code
#Using the hippocampal LPS datasets

#################

library(gemma.R)
library(tidyr)
library(dplyr)

#################

#Set your working directory

#Sample code:
setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Teaching/BrainDataAlchemy/Summer_2024/Summer2024_Pipeline/Example_Analysis/QueryResults")

list.files()

#Read in your final chosen dataset spreadsheet:

MyDatasets_Screened<-read.csv("MyDatasets_Screened.csv", stringsAsFactors = FALSE, header=TRUE)

str(MyDatasets_Screened)

# 'data.frame':	3 obs. of  24 variables:
#   $ X                         : int  1 2 3
# $ experiment.shortName      : chr  "GSE126678" "GSE181285" "GSE205325"
# $ experiment.name           : chr  "Enduring and sex-specific changes in hippocampal gene expression after subchronic immune challenge" "Expression data from the hippocampus of LPS induced depression  mice model treated with Luteolin" "Modulation of behavioral and hippocampal transcriptomic responses in rat prolonged chronic unpredictable stress"| __truncated__
# $ experiment.ID             : int  15127 21341 24923
# $ experiment.description    : chr  " The goals of this study include examining long-lasting changes in hippocampal gene expression in males and in "| __truncated__ " Gene expression profiling reveals a potential role of Luteolin in LPS induced depression model LPS depression "| __truncated__ " Here, we examine behavioral and brain transcriptomic (RNA-seq) responses in rat prolonged chronic unpredictabl"| __truncated__
# $ experiment.troubled       : logi  FALSE FALSE FALSE
# $ experiment.accession      : chr  "GSE126678" "GSE181285" "GSE205325"
# $ experiment.database       : chr  "GEO" "GEO" "GEO"
# $ experiment.URI            : chr  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126678" "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181285" "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205325"
# $ experiment.sampleCount    : int  24 16 21
# $ experiment.lastUpdated    : chr  "2023-12-16 09:57:45.253" "2023-12-17 13:14:26.142" "2023-12-18 00:16:39.838"
# $ experiment.batchEffectText: chr  "NO_BATCH_INFO" "NO_BATCH_INFO" "NO_BATCH_INFO"
# $ experiment.batchCorrected : logi  FALSE FALSE FALSE
# $ experiment.batchConfound  : int  0 0 0
# $ experiment.batchEffect    : int  0 0 0
# $ experiment.rawData        : int  1 1 1
# $ geeq.qScore               : num  0.282 0.283 0.28
# $ geeq.sScore               : num  1 0.812 0.75
# $ taxon.name                : chr  "mouse" "mouse" "rat"
# $ taxon.scientific          : chr  "Mus musculus" "Mus musculus" "Rattus norvegicus"
# $ taxon.ID                  : int  2 2 3
# $ taxon.NCBI                : int  10090 10090 10116
# $ taxon.database.name       : chr  "mm10" "mm10" "rn6"
# $ taxon.database.ID         : int  81 81 86


#And pull out the GSE #s for your datasets:

ExperimentIDs<-MyDatasets_Screened$experiment.shortName
ExperimentIDs
#[1] "GSE126678" "GSE181285" "GSE205325"

#Alternatively, you can just enter them by hand as a vector:
ExperimentIDs<-c("GSE126678", "GSE181285", "GSE205325")


#################

#For the meta-analysis, we will be extracting the differential expression results from Gemma. 
#The differential expression results for each dataset may include multiple result sets (e.g., one result set for the subset of the data from the hippocampus, one result set for frontal cortex). 
#Each of these result sets may have multiple statistical contrasts (e.g., drug1 vs. vehicle, drug2 vs. vehicle). 
#Therefore, each of the statistical contrasts is labeled with a result id and contrast id within the Gemma database. 
#We will need to know which of these ids are relevant to our project goals to easily extract their results.

#We will also need to double-check that these statistical contrasts are set up in a manner that makes sense for our experiments:

#First, for experiments that include more than one brain region, we will need to double-check that the results have been subsetted by brain region (instead of including brain region ("OrganismPart") as a factor in the model). If they haven't been subsetted by region, we will probably need to re-run the differential expression analysis.

#Depending on the goals of the meta-analysis, we may also need to re-run the differential expression analysis to remove other unwanted subjects (e.g., removing subjects with genotypes that might interfere with our results)

#Second, we will need to double-check that the comparisons include an appropriate reference group - sometimes they are reversed in Gemma (e.g., having the drug treatment set as the baseline, with vehicle as the manipulation). If this is the case, we will need to invert the effects when we input them into our meta-analysis (multiply the effects by -1).


#This code makes a data.frame that includes all of the contrast ids for each dataset with their basic metadata in a format that is easily readable in a spreadsheet program.

source("Function_GettingResultSetInfoForDatasets.R")

#Example usage:
GettingResultSetInfoForDatasets(ExperimentIDs)

#######################

#By hand, sort through the ResultSets_toScreen.csv to identify the statistical contrasts useful to your meta-analysis

#Delete all of the rows that you determine to be not applicable to your meta-analysis

#Save the file as a comma separate variable file named "ResultSets_Screened.csv"

#Afterwards, it would be good to double check your selections with me before proceeding forward.


#######################

#Read in your screened result sets:

ResultSet_contrasts<-read.csv("ResultSets_Screened.csv", header=TRUE, stringsAsFactors = FALSE )
str(ResultSet_contrasts)

source("Function_DownloadingDEResults.R")

#example usage:
DownloadingDEResults(ResultSet_contrasts)

str(differentials)
#Also seems to work
#List of 3 data frames

#Each of those is a result set containing all of the differential expression results for that particular variable. 
#If the variable has more than one level (e.g., 3 LPS dosages), you will have several statistical contrasts included in the result set differential expression output.

#Here is how you can access and review the differential expression results for a particular result set id:
str(differentials[1])

###############################

#For record-keeping purposes, let's save the differential expression results for each result set:

#The reason for this is because we will have a hard time replicating our results in the future with just our code
#Because Gemma updates the annotation on the data every time a new reference genome is released

setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Teaching/BrainDataAlchemy/Summer_2024/Summer2024_Pipeline/Example_Analysis/Gemma_DE_Results")

source("Function_SavingGemmaDEResults_forEachResultSet.R")

#Example usage:

SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)

#################################

#Next we will start working with cleaning up the results for a single result set

#We need to filter down our differential expression results to just the rows with good gene annotation

DE_Results<-differentials[[i]]

str(DE_Results)
#Classes ‘data.table’ and 'data.frame':	21693 obs. of  19 variables:

#Reading in the function

source("Function_FilteringDEResults_GoodAnnotation.R")

#Example of using the function for a dataset:

FilteringDEResults_GoodAnnotation(differentials[[1]])

# [1] "# of rows in results"
# [1] 21693
# [1] "# of rows with missing NCBI annotation:"
# [1] NA
# [1] "# of rows with NA NCBI annotation:"
# [1] 79
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 79
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 0
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 0
# [1] "# of rows with good annotation"
# [1] 21614
# [1] "Outputted object: DE_Results_GoodAnnotation"

str(DE_Results_GoodAnnotation)

#################

#Next we are going to pull out the differential expression for the specific statistical contrasts that we are interested in

source("Function_ExtractingDEResultsForContrasts.R")

#Example usage:
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)
  
# [1] "GSE126678_lipopolysaccharide has modifier Acute immune challenge , vehicle; lipopolysaccharide has modifier Acute immune challenge , vehicle"                       
# [2] "GSE126678_lipopolysaccharide has modifier Long-term subchronic immune challenge + acute immune challenge"                                                                  
# [3] "GSE126678_lipopolysaccharide has modifier long-term subchronic immune challenge , vehicle; lipopolysaccharide has modifier long-term subchronic immune challenge , vehicle"


#Those names are super unwieldy. I'm going to rename them:
ComparisonsOfInterest<-c("GSE126678_LPS_Acute", "GSE126678_LPS_SubchronicPlusAcute", "GSE126678_Subchronic")


#####################

#Sometimes gene expression is measured using multiple probes (microarray)
#Next we need to collapse our differential expression results down to one result per gene
#At the same time, we will calculate the standard error for our effect size (Log2FC) using the t-statistic
#And then square the standard error to get the sampling variance

source("Function_CollapsingDEResults_OneResultPerGene.R")

#Notes about parameters for function CollapsingDEResults_OneResultPerGene()
#GSE_ID is a string indicating the name of the Gemma dataset
#DE_Results_GoodAnnotation is the data frame outputted by our previous function
#ComparisonsOfInterest is a character vector containing the names of the group comparisons of interest within this dataset. Important: These group comparisons should be listed in exactly the same order as the order that you provide the column names for their associated Fold Change and Tstat output.
#NamesOfFoldChangeColumns is a vector containing the names of the columns of DE_Results_GoodAnnotation containing the FoldChange results for your comparisons of interes, in the same order as the ComparisonsOfInterest vector
#NamesOfTstatColumns is a vector containing the names of the columns of DE_Results_GoodAnnotation containing the Tstat results for your comparisons of interes, in the same order as the ComparisonsOfInterest vector

#Example usage:
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

# [1] "Double check that the vectors containing the fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:"
# [1] "# of rows with unique NCBI IDs:"
# [1] 21614
# [1] "# of rows with unique Gene Symbols:"
# [1] 21614
# [1] "Dimensions of Fold Change matrix, averaged by gene symbol:"
# [1] 21614     3
# [1] "Output: Named DEResults_GSE126678"

#################

FilteringDEResults_GoodAnnotation(differentials[[2]])

# [1] "# of rows in results"
# [1] 45100
# [1] "# of rows with missing NCBI annotation:"
# [1] 13489
# [1] "# of rows with NA NCBI annotation:"
# [1] 0
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 13489
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 717
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 717
# [1] "# of rows with good annotation"
# [1] 30894
# [1] "Outputted object: DE_Results_GoodAnnotation"


ComparisonsOfInterest
#[1] "GSE181285_lipopolysaccharide delivered at dose 850 ug/kg"

#That names is super unwieldy. I'm going to rename it:
ComparisonsOfInterest<-c("GSE181285_LPS_Acute")


CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

# [1] "Double check that the vectors containing the fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:"
# [1] "# of rows with unique NCBI IDs:"
# [1] 18563
# [1] "# of rows with unique Gene Symbols:"
# [1] 18563
# [1] "Dimensions of Fold Change matrix, averaged by gene symbol:"
# [1] 18563     1
# [1] "Output: Named DEResults_GSE181285"

#################

FilteringDEResults_GoodAnnotation(differentials[[3]])

# [1] "# of rows in results"
# [1] 18964
# [1] "# of rows with missing NCBI annotation:"
# [1] NA
# [1] "# of rows with NA NCBI annotation:"
# [1] 1768
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 1768
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 0
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 0
# [1] "# of rows with good annotation"
# [1] 17196
# [1] "Outputted object: DE_Results_GoodAnnotation"


#These are the names of the Log2FC Columns for our statistical contrasts of interest within the differential expression results for this particular result set:

NamesOfFoldChangeColumns<-colnames(DE_Results_GoodAnnotation)[colnames(DE_Results_GoodAnnotation)%in%Contrasts_Log2FC]

#These are the names of the T-statistic Columns for our statistical contrasts of interest within the differential expression results for this particular result set:

NamesOfTstatColumns<-colnames(DE_Results_GoodAnnotation)[colnames(DE_Results_GoodAnnotation)%in%Contrasts_Tstat]

#Next we're going to pull out the contrast IDs associated with each result:

GetContrastIDsforResultSet(NamesOfFoldChangeColumns)

ContrastIDs_inCurrentDF
#[1] "204289"

#Next we're going to grab some metadata to go with those statistical contrasts

#I would like the dataset id for our contrasts:
Datasets_inCurrentDF<-ResultSet_contrasts$ExperimentID[ResultSet_contrasts$ContrastIDs%in%ContrastIDs_inCurrentDF]

#And I would like the experimental factor information for our statistical contrast:
Factors_inCurrentDF<-ResultSet_contrasts$ExperimentalFactors[ResultSet_contrasts$ContrastIDs%in%ContrastIDs_inCurrentDF]

#We can combine those to make an interpretable unique identifier for each statistical comparison:
ComparisonsOfInterest<-paste(Datasets_inCurrentDF, Factors_inCurrentDF, sep="_" )

ComparisonsOfInterest
#[1] "GSE205325_lipopolysaccharide delivered at dose 28 x 0.1 mg/kg"

#That names is super unwieldy. I'm going to rename it:
ComparisonsOfInterest<-c("GSE205325_LPS_Chronic")

GSE_ID<-Datasets_inCurrentDF[1]
GSE_ID
#[1] "GSE205325"

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######################

#To do:
#double-check that I successfully functionalized some of the middle code where we pull out column names etc.
#make sure environment is actually cleaned up after each dataset
#add contrast id to items named with just the dataset name - directory, DEResults object
#make notes that print out in the console also print out into a text document
#split up function script and sample pipeline script

#######################

library(plyr)

AligningRatDatasets<-function(ListOfRatDEResults){
  
  Rat_MetaAnalysis_FoldChange_Dfs<-list()
  
  for(i in c(1:length(ListOfRatDEResults))){
    Rat_MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(Rat_EntrezGene.ID=row.names(ListOfRatDEResults[[i]][[1]]),ListOfRatDEResults[[i]][[1]], stringsAsFactors=FALSE)
  }
  
  print("Rat_MetaAnalysis_FoldChange_Dfs:")
  print(str(Rat_MetaAnalysis_FoldChange_Dfs))
  
  Rat_MetaAnalysis_FoldChanges<<-join_all(Rat_MetaAnalysis_FoldChange_Dfs, by="Rat_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("Rat_MetaAnalysis_FoldChanges:")
  print(str(Rat_MetaAnalysis_FoldChanges))
  
  Rat_MetaAnalysis_SV_Dfs<-list()
  
  for(i in c(1:length(ListOfRatDEResults))){
    Rat_MetaAnalysis_SV_Dfs[[i]]<-data.frame(Rat_EntrezGene.ID=row.names(ListOfRatDEResults[[i]][[4]]),ListOfRatDEResults[[i]][[4]], stringsAsFactors=FALSE)
  }
  
  print("Rat_MetaAnalysis_SV_Dfs:")
  print(str(Rat_MetaAnalysis_SV_Dfs))
  
  Rat_MetaAnalysis_SV<<-join_all(Rat_MetaAnalysis_SV_Dfs, by="Rat_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("Rat_MetaAnalysis_SV:")
  print(str(Rat_MetaAnalysis_SV))
  
  rm(Rat_MetaAnalysis_SV_Dfs, Rat_MetaAnalysis_FoldChange_Dfs)
}

#Example Usage;

ListOfRatDEResults<-list(DEResults_GSE205325)

AligningRatDatasets(ListOfRatDEResults)


###########

AligningMouseDatasets<-function(ListOfMouseDEResults){
  
  Mouse_MetaAnalysis_FoldChange_Dfs<-list()
  
  for(i in c(1:length(ListOfMouseDEResults))){
    Mouse_MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(Mouse_EntrezGene.ID=row.names(ListOfMouseDEResults[[i]][[1]]),ListOfMouseDEResults[[i]][[1]], stringsAsFactors=FALSE)
  }
  
  print("Mouse_MetaAnalysis_FoldChange_Dfs:")
  print(str(Mouse_MetaAnalysis_FoldChange_Dfs))
  
  Mouse_MetaAnalysis_FoldChanges<<-join_all(Mouse_MetaAnalysis_FoldChange_Dfs, by="Mouse_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("Mouse_MetaAnalysis_FoldChanges:")
  print(str(Mouse_MetaAnalysis_FoldChanges))
  
  Mouse_MetaAnalysis_SV_Dfs<-list()
  
  for(i in c(1:length(ListOfMouseDEResults))){
    Mouse_MetaAnalysis_SV_Dfs[[i]]<-data.frame(Mouse_EntrezGene.ID=row.names(ListOfMouseDEResults[[i]][[4]]),ListOfMouseDEResults[[i]][[4]], stringsAsFactors=FALSE)
  }
  
  print("Mouse_MetaAnalysis_SV_Dfs:")
  print(str(Mouse_MetaAnalysis_SV_Dfs))
  
  Mouse_MetaAnalysis_SV<<-join_all(Mouse_MetaAnalysis_SV_Dfs, by="Mouse_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("Mouse_MetaAnalysis_SV:")
  print(str(Mouse_MetaAnalysis_SV))
  
  rm(Mouse_MetaAnalysis_SV_Dfs, Mouse_MetaAnalysis_FoldChange_Dfs)
}

#Example Usage;

ListOfMouseDEResults<-list(DEResults_GSE126678, DEResults_GSE181285)

AligningMouseDatasets(ListOfMouseDEResults)


############

#Code for aligning the rat and mice results:

#We have the ortholog database that we downloaded from Jackson Labs on April 25, 2024
#This database was trimmed and formatted using the code "FormattingRatMouseOrthologDatabase_20240425.R"
MouseVsRat_NCBI_Entrez<-read.csv("MouseVsRat_NCBI_Entrez_JacksonLab_20240425.csv", header=TRUE, stringsAsFactors = FALSE, row.names=1, colClasses=c("character", "character", "character"))

Mouse_MetaAnalysis_FoldChanges_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_FoldChanges, by="Mouse_EntrezGene.ID", type="full")
str(Mouse_MetaAnalysis_FoldChanges_wOrthologs)
#'data.frame':	25288 obs. of  6 variables:

#If there are rat datasets:
MetaAnalysis_FoldChanges<-join(Mouse_MetaAnalysis_FoldChanges_wOrthologs, Rat_MetaAnalysis_FoldChanges, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_FoldChanges)
#'data.frame':	28101 obs. of  7 variables:

#If there aren't any rat datasets:
MetaAnalysis_FoldChanges<-Mouse_MetaAnalysis_FoldChanges_wOrthologs
str(MetaAnalysis_FoldChanges)


Mouse_MetaAnalysis_SV_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_SV, by="Mouse_EntrezGene.ID", type="full")
str(Mouse_MetaAnalysis_SV_wOrthologs)
#'data.frame':	25288 obs. of  6 variables:

#If there are rat datasets:
MetaAnalysis_SV<-join(Mouse_MetaAnalysis_SV_wOrthologs, Rat_MetaAnalysis_SV, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_SV)
#'data.frame':	28101 obs. of  7 variables:

#If there aren't any rat datasets:
MetaAnalysis_SV<-Mouse_MetaAnalysis_SV_wOrthologs

#For simplicity's sake, I'm going to replace that Mouse-Rat Entrez annotation
#Because it is missing entries for any genes in the datasets that *don't* have orthologs
MetaAnalysis_FoldChanges$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_FoldChanges$Mouse_EntrezGene.ID, MetaAnalysis_FoldChanges$Rat_EntrezGene.ID, sep="_")
MetaAnalysis_SV$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_SV$Mouse_EntrezGene.ID, MetaAnalysis_SV$Rat_EntrezGene.ID, sep="_")


#Comparing Log2FC across datasets

#Simple scatterplot... not so promising:
colnames(MetaAnalysis_FoldChanges)

# [1] "Rat_EntrezGene.ID"                  "Mouse_EntrezGene.ID"                "MouseVsRat_EntrezGene.ID"          
# [4] "LPS_SubchronicAndAcute_vs_Vehicle"  "LPS_Acute_vs_Vehicle"               "LPS_Subchronic_vs_Vehicle"         
# [7] "LPS_Acute850ugPerKg_vs_Vehicle"     "LPS_Chronic_vs_Vehicle_AllStressed"


plot(MetaAnalysis_FoldChanges$LPS_Subchronic_vs_Vehicle~MetaAnalysis_FoldChanges$LPS_Chronic_vs_Vehicle_AllStressed)

#Note - many people prefer to plot these relationships using RRHOs (Rank rank hypergeometric overlap plots)
#I like using both.
#The code for the RRHOs is a little complicated, but I'm happy to share if folks are interested.

cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman")
#There isn't much similarity across conditions here (outside of comparisons within the same experiment)

#An illustration of the correlation matrix using a hierarchically clustered heatmap, although somewhat pathetic:
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman"))



