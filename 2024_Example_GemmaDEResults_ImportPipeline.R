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

#Set your working directory - 
#A working directory was already defined when you created your R project
#If not:

#Sample code:
setwd("/Users/hagenaue/Documents/BrainDataAlchemy_2024_Example")

#The functions and files that you will need to should be added to your working directory:

list.files()
# [1] "2024_Example_Code_forGemmaSearch_Part1.R"        
# [2] "2024_Example_GemmaDEResults_ImportPipeline.R"    
# [3] "2024_ExampleCode_SettingUpDirectory.R"           
# [4] "BrainDataAlchemy_2024_Example.Rproj"             
# [5] "Function_CollapsingDEResults_OneResultPerGene.R" 
# [6] "Function_DownloadingDEResults.R"                 
# [7] "Function_ExtractingDEResultsForContrasts.R"      
# [8] "Function_FilteringDEResults_GoodAnnotation.R"    
# [9] "Function_GettingResultSetInfoForDatasets.R"      
# [10] "Function_SavingGemmaDEResults_forEachResultSet.R"
# [11] "MyDatasets_Screened.csv"                         
# [12] "ResultSets_Screened.csv"                         
# [13] "ResultSets_toScreen.csv"  

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

# [1] "These are the result sets that you identified as being of interest"
# [1] 553805 570552 556647
# List of 3
# $ :Classes ‘data.table’ and 'data.frame':	21693 obs. of  19 variables:
#   ..$ Probe                      : int [1:21693] 20344 20558 110454 20339 17067 68774 21946 21366 76905 16190 ...
# ..$ NCBIid                     : int [1:21693] 20344 20558 110454 20339 17067 68774 21946 21366 76905 16190 ...
# ..$ GeneSymbol                 : chr [1:21693] "Selp" "Slfn4" "Ly6a" "Sele" ...
# ..$ GeneName                   : chr [1:21693] "selectin, platelet" "schlafen 4" "lymphocyte antigen 6 family member A" "selectin, endothelial cell" ...
# ..$ pvalue                     : num [1:21693] 1.70e-10 9.11e-09 5.38e-08 7.73e-08 1.00e-07 ...
# ..$ corrected_pvalue           : num [1:21693] 3.69e-06 9.88e-05 4.00e-04 4.00e-04 4.00e-04 ...
# ..$ rank                       : num [1:21693] 4.61e-05 9.22e-05 1.00e-04 2.00e-04 2.00e-04 ...
# ..$ contrast_151617_coefficient: num [1:21693] 3.102 3.27 0.763 2.868 0.513 ...
# ..$ contrast_151617_log2fc     : num [1:21693] 3.102 3.27 0.763 2.868 0.513 ...
# ..$ contrast_151617_tstat      : num [1:21693] 2.97 3.62 3.86 2.31 4.17 ...
# ..$ contrast_151617_pvalue     : num [1:21693] 0.0077 0.0017 0.001 0.0318 0.0005 ...
# ..$ contrast_151618_coefficient: num [1:21693] 7.01 5.21 1.35 5.86 0.82 ...
# ..$ contrast_151618_log2fc     : num [1:21693] 7.01 5.21 1.35 5.86 0.82 ...
# ..$ contrast_151618_tstat      : num [1:21693] 7.81 6.14 7.02 5.28 6.78 ...
# ..$ contrast_151618_pvalue     : num [1:21693] 1.86e-07 5.64e-06 9.01e-07 3.78e-05 1.47e-06 ...
# ..$ contrast_151619_coefficient: num [1:21693] 1.223 -0.303 -0.243 0.118 -0.214 ...
# ..$ contrast_151619_log2fc     : num [1:21693] 1.223 -0.303 -0.243 0.118 -0.214 ...
# ..$ contrast_151619_tstat      : num [1:21693] 1.0192 -0.2593 -1.1221 0.0785 -1.6535 ...
# ..$ contrast_151619_pvalue     : num [1:21693] 0.32 0.798 0.275 0.938 0.114 ...
# ..- attr(*, ".internal.selfref")=<externalptr> 
#   ..- attr(*, "call")= chr "https://gemma.msl.ubc.ca/rest/v2/resultSets/553805"
# ..- attr(*, "env")=<environment: 0x7f897aec2838> 
#   $ :Classes ‘data.table’ and 'data.frame':	45100 obs. of  11 variables:
#   ..$ Probe                      : chr [1:45100] "1426114_PM_at" "1452267_PM_at" "1441388_PM_at" "1436221_PM_at" ...
# ..$ NCBIid                     : chr [1:45100] "" "224613" "" "100039795" ...
# ..$ GeneSymbol                 : chr [1:45100] "" "Flywch1" "" "Ildr2" ...
# ..$ GeneName                   : chr [1:45100] "" "FLYWCH-type zinc finger 1" "" "immunoglobulin-like domain containing receptor 2" ...
# ..$ pvalue                     : num [1:45100] 3.74e-09 1.57e-08 1.32e-08 8.87e-09 2.37e-08 ...
# ..$ corrected_pvalue           : num [1:45100] 2e-04 2e-04 2e-04 2e-04 2e-04 3e-04 3e-04 3e-04 4e-04 5e-04 ...
# ..$ rank                       : num [1:45100] 2.22e-05 8.87e-05 6.65e-05 4.44e-05 1.00e-04 ...
# ..$ contrast_186753_coefficient: num [1:45100] 0.1587 0.0967 -0.5049 0.7694 0.1039 ...
# ..$ contrast_186753_log2fc     : num [1:45100] 0.1587 0.0967 -0.5049 0.7694 0.1039 ...
# ..$ contrast_186753_tstat      : num [1:45100] 1.65 2.28 -5.4 10.4 2.11 ...
# ..$ contrast_186753_pvalue     : num [1:45100] 1.20e-01 3.83e-02 8.18e-05 3.93e-08 5.26e-02 ...
# ..- attr(*, ".internal.selfref")=<externalptr> 
#   ..- attr(*, "call")= chr "https://gemma.msl.ubc.ca/rest/v2/resultSets/570552"
# ..- attr(*, "env")=<environment: 0x7f897ab59f78> 
#   $ :Classes ‘data.table’ and 'data.frame':	18964 obs. of  11 variables:
#   ..$ Probe                      : int [1:18964] 497918 303604 501110 498386 299052 25636 65166 102548697 100364561 24183 ...
# ..$ NCBIid                     : int [1:18964] 497918 303604 501110 498386 299052 25636 65166 102548697 100364561 24183 ...
# ..$ GeneSymbol                 : chr [1:18964] "Smcr8" "Map3k3" "Gsta6" "Cc2d2a" ...
# ..$ GeneName                   : chr [1:18964] "SMCR8-C9orf72 complex subunit" "mitogen activated protein kinase kinase kinase 3" "glutathione S-transferase alpha 6" "coiled-coil and C2 domain containing 2A" ...
# ..$ pvalue                     : num [1:18964] 0.869 0.885 0.23 0.597 0.518 ...
# ..$ corrected_pvalue           : num [1:18964] 1 1 1 1 1 ...
# ..$ rank                       : num [1:18964] 0.855 0.875 0.191 0.569 0.486 ...
# ..$ contrast_204289_coefficient: num [1:18964] 0.0982 -0.0239 -0.2088 -0.1256 -0.3315 ...
# ..$ contrast_204289_log2fc     : num [1:18964] 0.0982 -0.0239 -0.2088 -0.1256 -0.3315 ...
# ..$ contrast_204289_tstat      : num [1:18964] 0.167 -0.146 -1.242 -0.538 -0.659 ...
# ..$ contrast_204289_pvalue     : num [1:18964] 0.869 0.885 0.23 0.597 0.518 ...
# ..- attr(*, ".internal.selfref")=<externalptr> 
#   ..- attr(*, "call")= chr "https://gemma.msl.ubc.ca/rest/v2/resultSets/556647"
# ..- attr(*, "env")=<environment: 0x7f88feede760> 
#   [1] "These are the result sets that had differential expression results:"
# [1] 553805 570552 556647
# [1] "Your differential expression results for each of your result sets are stored in the object named differentials. This object is structured as a list of data frames. Each element in the list represetns a result set, with the data frame containing the differential expression results"
# [1] "These are the columns for the effect sizes for our statistical contrasts of interest (Log(2) Fold Changes"
# [1] "contrast_151618_log2fc" "contrast_151617_log2fc" "contrast_151619_log2fc"
# [4] "contrast_186753_log2fc" "contrast_204289_log2fc"
# [1] "these are the columns for the T-statistics for our statistical contrasts of interest - we will use that information to derive the sampling variances"
# [1] "contrast_151618_tstat" "contrast_151617_tstat" "contrast_151619_tstat"
# [4] "contrast_186753_tstat" "contrast_204289_tstat"


#####

str(differentials)
#List of 3

#Each of those is a result set containing all of the differential expression results for that particular variable. 
#If the variable has more than one level (e.g., 3 LPS dosages), you will have several statistical contrasts included in the result set differential expression output.

#Here is how you can access and review the differential expression results for a particular result set id:
str(differentials[1])
# List of 1
# $ :Classes ‘data.table’ and 'data.frame':	21693 obs. of  19 variables:
#   ..$ Probe                      : int [1:21693] 20344 20558 110454 20339 17067 68774 21946 21366 76905 16190 ...
# ..$ NCBIid                     : int [1:21693] 20344 20558 110454 20339 17067 68774 21946 21366 76905 16190 ...
# ..$ GeneSymbol                 : chr [1:21693] "Selp" "Slfn4" "Ly6a" "Sele" ...
# ..$ GeneName                   : chr [1:21693] "selectin, platelet" "schlafen 4" "lymphocyte antigen 6 family member A" "selectin, endothelial cell" ...
# ..$ pvalue                     : num [1:21693] 1.70e-10 9.11e-09 5.38e-08 7.73e-08 1.00e-07 ...
# ..$ corrected_pvalue           : num [1:21693] 3.69e-06 9.88e-05 4.00e-04 4.00e-04 4.00e-04 ...
# ..$ rank                       : num [1:21693] 4.61e-05 9.22e-05 1.00e-04 2.00e-04 2.00e-04 ...
# ..$ contrast_151617_coefficient: num [1:21693] 3.102 3.27 0.763 2.868 0.513 ...
# ..$ contrast_151617_log2fc     : num [1:21693] 3.102 3.27 0.763 2.868 0.513 ...
# ..$ contrast_151617_tstat      : num [1:21693] 2.97 3.62 3.86 2.31 4.17 ...
# ..$ contrast_151617_pvalue     : num [1:21693] 0.0077 0.0017 0.001 0.0318 0.0005 ...
# ..$ contrast_151618_coefficient: num [1:21693] 7.01 5.21 1.35 5.86 0.82 ...
# ..$ contrast_151618_log2fc     : num [1:21693] 7.01 5.21 1.35 5.86 0.82 ...
# ..$ contrast_151618_tstat      : num [1:21693] 7.81 6.14 7.02 5.28 6.78 ...
# ..$ contrast_151618_pvalue     : num [1:21693] 1.86e-07 5.64e-06 9.01e-07 3.78e-05 1.47e-06 ...
# ..$ contrast_151619_coefficient: num [1:21693] 1.223 -0.303 -0.243 0.118 -0.214 ...
# ..$ contrast_151619_log2fc     : num [1:21693] 1.223 -0.303 -0.243 0.118 -0.214 ...
# ..$ contrast_151619_tstat      : num [1:21693] 1.0192 -0.2593 -1.1221 0.0785 -1.6535 ...
# ..$ contrast_151619_pvalue     : num [1:21693] 0.32 0.798 0.275 0.938 0.114 ...
# ..- attr(*, ".internal.selfref")=<externalptr> 
#   ..- attr(*, "call")= chr "https://gemma.msl.ubc.ca/rest/v2/resultSets/553805"
# ..- attr(*, "env")=<environment: 0x7f89a3033660> 

###############################

#For record-keeping purposes, let's save the differential expression results for each result set:

#The reason for this is because we will have a hard time replicating our results in the future with just our code
#Because Gemma updates the annotation on the data every time a new reference genome is released

source("Function_SavingGemmaDEResults_forEachResultSet.R")

#Example usage:

SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)

#################################

#Next we will start working with cleaning up the results for a single result set

#We need to filter down our differential expression results to just the rows with good gene annotation

DE_Results<-differentials[[1]]

str(DE_Results)
#Classes ‘data.table’ and 'data.frame':	21693 obs. of  19 variables:

#Reading in the function

source("Function_FilteringDEResults_GoodAnnotation.R")

#Example of using the function for a dataset:

FilteringDEResults_GoodAnnotation(DE_Results)

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
# Classes ‘data.table’ and 'data.frame':	21614 obs. of  19 variables:
#   $ Probe                      : int  20344 20558 110454 20339 17067 68774 21946 21366 76905 16190 ...
# $ NCBIid                     : int  20344 20558 110454 20339 17067 68774 21946 21366 76905 16190 ...
# $ GeneSymbol                 : chr  "Selp" "Slfn4" "Ly6a" "Sele" ...
# $ GeneName                   : chr  "selectin, platelet" "schlafen 4" "lymphocyte antigen 6 family member A" "selectin, endothelial cell" ...
# $ pvalue                     : num  1.70e-10 9.11e-09 5.38e-08 7.73e-08 1.00e-07 ...
# $ corrected_pvalue           : num  3.69e-06 9.88e-05 4.00e-04 4.00e-04 4.00e-04 ...
# $ rank                       : num  4.61e-05 9.22e-05 1.00e-04 2.00e-04 2.00e-04 ...
# $ contrast_151617_coefficient: num  3.102 3.27 0.763 2.868 0.513 ...
# $ contrast_151617_log2fc     : num  3.102 3.27 0.763 2.868 0.513 ...
# $ contrast_151617_tstat      : num  2.97 3.62 3.86 2.31 4.17 ...
# $ contrast_151617_pvalue     : num  0.0077 0.0017 0.001 0.0318 0.0005 ...
# $ contrast_151618_coefficient: num  7.01 5.21 1.35 5.86 0.82 ...
# $ contrast_151618_log2fc     : num  7.01 5.21 1.35 5.86 0.82 ...
# $ contrast_151618_tstat      : num  7.81 6.14 7.02 5.28 6.78 ...
# $ contrast_151618_pvalue     : num  1.86e-07 5.64e-06 9.01e-07 3.78e-05 1.47e-06 ...
# $ contrast_151619_coefficient: num  1.223 -0.303 -0.243 0.118 -0.214 ...
# $ contrast_151619_log2fc     : num  1.223 -0.303 -0.243 0.118 -0.214 ...
# $ contrast_151619_tstat      : num  1.0192 -0.2593 -1.1221 0.0785 -1.6535 ...
# $ contrast_151619_pvalue     : num  0.32 0.798 0.275 0.938 0.114 ...
# - attr(*, ".internal.selfref")=<externalptr> 
#   - attr(*, "call")= chr "https://gemma.msl.ubc.ca/rest/v2/resultSets/553805"
# - attr(*, "env")=<environment: 0x7f89998ffea8> 

#################

#Next we are going to pull out the differential expression for the specific statistical contrasts that we are interested in

source("Function_ExtractingDEResultsForContrasts.R")

#Example usage:
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)
  
# [1] "These are all of the columns in the differential expression results for our current result set:"
# [1] "Probe"                       "NCBIid"                     
# [3] "GeneSymbol"                  "GeneName"                   
# [5] "pvalue"                      "corrected_pvalue"           
# [7] "rank"                        "contrast_151617_coefficient"
# [9] "contrast_151617_log2fc"      "contrast_151617_tstat"      
# [11] "contrast_151617_pvalue"      "contrast_151618_coefficient"
# [13] "contrast_151618_log2fc"      "contrast_151618_tstat"      
# [15] "contrast_151618_pvalue"      "contrast_151619_coefficient"
# [17] "contrast_151619_log2fc"      "contrast_151619_tstat"      
# [19] "contrast_151619_pvalue"     
# [1] "These are the names of the Log(2) Fold Change Columns for our statistical contrasts of interest within the differential expression results for this particular result set:"
# [1] "contrast_151617_log2fc" "contrast_151618_log2fc" "contrast_151619_log2fc"
# [1] "These are the names of the T-statistic Columns for our statistical contrasts of interest within the differential expression results for this particular result set:"
# [1] "contrast_151617_tstat" "contrast_151618_tstat" "contrast_151619_tstat"
# [1] "These are the contrast ids for the statistical contrasts of interest within your current result set:"
# [1] "151617" "151618" "151619"
# [1] "This is the dataset id for the result set and statistical contrasts:"
# [1] "GSE126678"
# [1] "These are the current names for your statistical contrasts of interest - if they are unwieldy, you may want to change them"
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


#... and then we would just move on to our next two contrasts:
#differentials[[2]] and differentials [[3]]

DE_Results<-differentials[[2]]

FilteringDEResults_GoodAnnotation(DE_Results)

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

ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

# [1] "These are all of the columns in the differential expression results for our current result set:"
# [1] "Probe"                       "NCBIid"                     
# [3] "GeneSymbol"                  "GeneName"                   
# [5] "pvalue"                      "corrected_pvalue"           
# [7] "rank"                        "contrast_186753_coefficient"
# [9] "contrast_186753_log2fc"      "contrast_186753_tstat"      
# [11] "contrast_186753_pvalue"     
# [1] "These are the names of the Log(2) Fold Change Columns for our statistical contrasts of interest within the differential expression results for this particular result set:"
# [1] "contrast_186753_log2fc"
# [1] "These are the names of the T-statistic Columns for our statistical contrasts of interest within the differential expression results for this particular result set:"
# [1] "contrast_186753_tstat"
# [1] "These are the contrast ids for the statistical contrasts of interest within your current result set:"
# [1] "186753"
# [1] "This is the dataset id for the result set and statistical contrasts:"
# [1] "GSE181285"
# [1] "These are the current names for your statistical contrasts of interest - if they are unwieldy, you may want to change them"
# [1] "GSE181285_lipopolysaccharide delivered at dose 850 ug/kg"

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

#########################

#On to ResultSet #3:

DE_Results<-differentials[[3]]

FilteringDEResults_GoodAnnotation(DE_Results)
# [1] "# of rows in results"
# [1] 18964
# [1] "# of rows with missing NCBI annotation:"
# [1] 0
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

ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

# [1] "These are all of the columns in the differential expression results for our current result set:"
# [1] "Probe"                       "NCBIid"                     
# [3] "GeneSymbol"                  "GeneName"                   
# [5] "pvalue"                      "corrected_pvalue"           
# [7] "rank"                        "contrast_204289_coefficient"
# [9] "contrast_204289_log2fc"      "contrast_204289_tstat"      
# [11] "contrast_204289_pvalue"     
# [1] "These are the names of the Log(2) Fold Change Columns for our statistical contrasts of interest within the differential expression results for this particular result set:"
# [1] "contrast_204289_log2fc"
# [1] "These are the names of the T-statistic Columns for our statistical contrasts of interest within the differential expression results for this particular result set:"
# [1] "contrast_204289_tstat"
# [1] "These are the contrast ids for the statistical contrasts of interest within your current result set:"
# [1] "204289"
# [1] "This is the dataset id for the result set and statistical contrasts:"
# [1] "GSE205325"
# [1] "These are the current names for your statistical contrasts of interest - if they are unwieldy, you may want to change them"
# [1] "GSE205325_lipopolysaccharide delivered at dose 28 x 0.1 mg/kg"

ComparisonsOfInterest<-c("GSE205325_LPS_Chronic")

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

# [1] "Double check that the vectors containing the fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:"
# [1] "# of rows with unique NCBI IDs:"
# [1] 17196
# [1] "# of rows with unique Gene Symbols:"
# [1] 17195
# [1] "Dimensions of Fold Change matrix, averaged by gene symbol:"
# [1] 17196     1
# [1] "Output: Named DEResults_GSE205325"


###################

#My differential expression results are now saved in three objects:
#DEResults_GSE126678
#DEResults_GSE181285
#DEResults_GSE205325

#Since I was debugging functions, I ended up with a lot of other extra, useless objects leftover in my environment
#This code removes them:
rm(differentials, FoldChangeColumn, MatrixOfColumnNames_BrokenUp, TempMasterResults, TstatColumn, ComparisonsOfInterest, ContrastIDs_inCurrentDF, Contrasts_Log2FC, Contrasts_Tstat, Datasets_inCurrentDF, DE_Results, DE_Results_GoodAnnotation, Factors_inCurrentDF, GSE_ID, i, NamesOfFoldChangeColumns, NamesOfTstatColumns, DE_Results_GoodAnnotation_FoldChange_Average, DE_Results_GoodAnnotation_FoldChange_AveragedByGene, DE_Results_GoodAnnotation_SE, DE_Results_GoodAnnotation_SE_Average, DE_Results_GoodAnnotation_SE_AveragedByGene, DE_Results_GoodAnnotation_SV, DE_Results_GoodAnnotation_Tstat_Average, DE_Results_GoodAnnotation_Tstat_AveragedByGene, ColumnNames_BrokenUp)

###################

#Save the code and workspace! 
#Based on your project settings, this may happen automatically
#Many times you will hear advice not to save the workspace because it should be able to be completely recreated using your code.
#In this situation, we definitely want to save the workspace because it can take time to import all of the results from Gemma's API and process them.
#Also, that makes it easier to go back and recreate figures, results etc when we need to for publications, etc
#Because if we import the results from Gemma again, they may be slightly different
#Due to updates to the genome annotation.

#The code file is saved using:
#"File"->"Save"
#The code file will have a file extension ".R"

#Whereas the workspace is saved using:
#"Session"->"Save workspace as"
#The workspace file will have a file extension ".Rdata"


