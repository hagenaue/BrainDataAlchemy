#Messing around with meta-analyses of Gemma output
#This file includes the code for comparing the differential expression output from different datasets and runing a meta-analysis (steps #6-10 mentioned below).  
#This is the fourth version of this file that includes one of the HDRF datasets as well (because the previous meta-analysis was being dominated by the suspicious CMS dataset with it's dubiously small sampling variances) as well as DG datasets (because I wanted to play with the code for adding dissection as a predictor to the model)
#At this point, I also started to functionalize things.
#Megan Hagenauer
#07-26-2022
#Associated R workspace: 
#save.image("~/Documents/Teaching/Grinnell Interns 2022/MetaAnalysis_HC_Stress_wHDRF_andDG_Fixed/Workspace_MetaAnalysis_HC_Stress_HDRFandDG_BugFixed.RData")


####################

#Overview of general coding steps:

#0) Reading in & visualizing the Log2 Expression data and Metadata for the individual samples - an illustration of where the differential expression results come from. 

#1) Read in the results from Gemma for each dataset

#2) Identify the results for the variables of interest

#3) Remove rows of data that have missing or unambiguous gene annotation

#### Question: Which gene annotation should we use? I was originally planning to use NCBI ID (because it is a little more stable than gene symbol), but if we use gene symbol we can run a half-way decent meta-analysis of data from both rats and mice without having to add in a step where we decode which genes in rats and mice are orthologous using an orthology database, as many mice and rat genes have the same gene symbol (76%, last time I checked).

#4) Collapse results (average) if there is more than row of data representing a single gene 

#5) Extract out the information needed for running a meta-analysis: Use the Log2FC and Tstat to calculate the standard error for the Log2FC, and then use the standard error to calculate the sampling variance.

#6) Combine together the relevant results from different studies into a single data-frame for the effect sizes (Log2FC) and a single data.frame for the sampling variances. 

#7) Make a correlation matrix to compare the overall results from the different studies. Further visualize the comparison using a hierarchically-clustered heatmap of the correlation matrix. Which studies seem to have similar results? 

#8) Run a meta-analysis using all of the effect sizes for each gene that has data from at least 2 studies. 

#9) Correct the meta-analysis output to take into account the fact that we are running the statistical calculations many times and therefore have a heightened risk of false discovery (false discovery rate correction) 

#10) Determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. 

#####################################

#Code packages used (may require installation & installation of dependencies):

library(plyr)
library(metafor)
library(reshape)
library(multtest)

######################################

#Set working directory:
setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Teaching/BrainDataAlchemy/Grinnell Interns 2022/2022Cohort_GoogleDrive/GoogleDriveDownload_20230508/ErinHernandez_Antidepressants_Hippocampus/Gemma_ReDownloaded_20230509")

########################################

#FUNCTIONS for reading in the datasets:

#######################################

#1) Read in the differential expression results from Gemma for each dataset:

#To start out with, I downloaded the Gemma differential expression output for the studies of interest.
#I put each study is in its own folder.

#I set working directory to where I downloaded the Gemma differential expression output for this particular study:
setwd("~/Documents/SideProjects/BrainGMT/Gemma/Hippocampus/11627_GSE59070_diffExpAnalysis_94802")
list.files()

#This procedure for reading in the differential expression results was pretty easily generalizable, so I functionalized it:


ReadingInGemmaDE<-function(ResultSetFileNames){
  
  #Reading in the analysis results file:
  TempAnalysisResults<-read.delim("analysis.results.txt", sep="\t", stringsAsFactors = FALSE, comment.char = "#")
  
  #I placed these results in a list format:
  TempResultsToJoin<-list(TempAnalysisResults)
  
  for(i in c(1:length(ResultSetFileNames))){
    #The result set files are the detailed differential expression results for a particular variable:
    TempResultsToJoin[[i]]<-read.delim(ResultSetFileNames[i], sep="\t", stringsAsFactors = FALSE, comment.char = "#")
  }
  
  TempResultsJoined<<-join_all(TempResultsToJoin, by="Element_Name")
  #Note: I've heard from other students that the join() and join_all() functions in the plyr package can cause problems in newer versions of R - you may need to replace this with merge and merge_all
  
  #Saving the joined results:
  write.csv(TempResultsJoined, "TempResultsJoined.csv")
  
  rm(TempAnalysisResults, TempResultsToJoin)
  
  print("Outputted object: TempResultsJoined")
}


#Notes about parameters for function: ReadingInGemmaDE()
#ResultSetFileNames should be a character vector containing the names of the result files

#Example Usage of the function:

#ReadingInGemmaDE(ResultSetFileNames=c("resultset_ID478782.data.txt"))
#[1] "Outputted object: TempResultsJoined"


#######################

#3) Remove rows of data that have missing or unambiguous gene annotation

#Notably, using modern probe annotation, many of these probes now map to more than one gene symbol, and when this happens Gemma separates provides a list of all of the gene symbols that the probe maps to using a pipe. 
#There are also probes that no longer map to any gene symbol (annotated with ""). 

#This data processing step was pretty generalizable, so I functionalized it:

FilteringDEResults_GoodAnnotation<-function(TempResultsJoined){
  
  print("# of rows with missing NCBI annotation:")
  print(sum(TempResultsJoined$NCBI_ID==""|TempResultsJoined$NCBI_ID=="null"))
  
  print("# of rows with missing Gene Symbol annotation:")
  print(sum(TempResultsJoined$Gene_Symbol==""|TempResultsJoined$Gene_Symbol=="null"))
  
  print("# of rows mapped to multiple NCBI_IDs:")
  print(length(grep('\\|', TempResultsJoined$NCBI_ID)))
  
  print("# of rows mapped to multiple Gene Symbols:")
  print(length(grep('\\|', TempResultsJoined$Gene_Symbol)))
  
  #I only want the subset of data which contains rows that do not contain a Gene Symbol of ""
  TempResultsJoined_NoNA<-TempResultsJoined[(TempResultsJoined$Gene_Symbol==""|TempResultsJoined$Gene_Symbol=="null")==FALSE,]
  
  if(length(grep('\\|', TempResultsJoined_NoNA$Gene_Symbol))==0){
    TempResultsJoined_NoNA_NoMultimapped<<-TempResultsJoined_NoNA
  }else{
    #I only want rows annotated with a single Gene Symbol (no pipe):
    TempResultsJoined_NoNA_NoMultimapped<<-TempResultsJoined_NoNA[-(grep('\\|', TempResultsJoined_NoNA$Gene_Symbol)),]
  }
  
  print("# of rows with good annotation")
  print(nrow(TempResultsJoined_NoNA_NoMultimapped))
  
  #For record keeping (sometimes useful for troubleshooting later)
  write.csv(TempResultsJoined_NoNA_NoMultimapped, "TempResultsJoined_NoNA_NoMultimapped.csv")
  
  rm(TempResultsJoined_NoNA, TempResultsJoined_NoNA_NoMultimapped)
  
  print("Outputted object: TempResultsJoined_NoNA_NoMultimapped")
}

#Example usage of function:

#FilteringDEResults_GoodAnnotation(TempResultsJoined)

#############################

#2) Identifying the results for the variables of interest for the meta-analysis:

#At this point, we should double-check:

#A) Whether Gemma has used a reference (baseline) group that makes sense for our analysis 

#B) Which group comparisons are actually of interest to us.


CollapsingDEResults_OneResultPerGene<-function(GSE_ID, TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns){
  
  print("Double check that the vectors containing the two fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:")
  
  print("# of rows with unique NCBI IDs:")
  print(length(unique(TempResultsJoined_NoNA_NoMultimapped$NCBI_ID)))
  
  print("# of rows with unique Gene Symbols:")
  print(length(unique(TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol)))
  
  #We will need both the Log2FC and T-stats averaged:
  
  TempResultsJoined_NoNA_NoMultimapped_FoldChange_Average<-list()
  
  for(i in c(1:length(NamesOfFoldChangeColumns))){
    
    TempResultsJoined_NoNA_NoMultimapped_FoldChange_Average[[i]]<-tapply(NamesOfFoldChangeColumns[i][[1]], TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, mean)
    
  }
  
  TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol<-do.call(cbind, TempResultsJoined_NoNA_NoMultimapped_FoldChange_Average)
  
  print("Dimensions of Fold Change matrix, averaged by gene symbol:")
  print(dim(TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol))
  
  colnames(TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol)<-ComparisonsOfInterest
  
  write.csv(TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol, "TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol.csv")
  
  
  TempResultsJoined_NoNA_NoMultimapped_Tstat_Average<-list()
  
  for(i in c(1:length(NamesOfFoldChangeColumns))){
    
    TempResultsJoined_NoNA_NoMultimapped_Tstat_Average[[i]]<-tapply(NamesOfTstatColumns[i][[1]], TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, mean)
    
  }
  
  TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol<-do.call(cbind, TempResultsJoined_NoNA_NoMultimapped_Tstat_Average)
  
  colnames(TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol)<-ComparisonsOfInterest
  
  write.csv(TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol, "TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol.csv")
  
  
  
  TempResultsJoined_NoNA_NoMultimapped_SE<-list()
  
  for(i in c(1:length(NamesOfFoldChangeColumns))){
    TempResultsJoined_NoNA_NoMultimapped_SE[[i]]<-NamesOfFoldChangeColumns[i][[1]]/NamesOfTstatColumns[i][[1]]
  }
  
  
  TempResultsJoined_NoNA_NoMultimapped_SE_Average<-list()
  
  for(i in c(1:length(NamesOfFoldChangeColumns))){
    
    TempResultsJoined_NoNA_NoMultimapped_SE_Average[[i]]<-tapply(TempResultsJoined_NoNA_NoMultimapped_SE[[i]], TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, mean)
    
  }
  
  TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol<-do.call(cbind, TempResultsJoined_NoNA_NoMultimapped_SE_Average)
  
  colnames(TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol)<-ComparisonsOfInterest
  
  write.csv(TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol, "TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol.csv")
  
  #For running our meta-analysis, we are actually going to need the sampling variance instead of the standard error
  #The sampling variance is just the standard error squared.
  
  TempResultsJoined_NoNA_NoMultimapped_SV<-(TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol)^2
  
  write.csv(TempResultsJoined_NoNA_NoMultimapped_SV, "TempResultsJoined_NoNA_NoMultimapped_SV.csv")
  
  TempMasterResults<-list(Log2FC=TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol, Tstat=TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol, SE=TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol, SV=TempResultsJoined_NoNA_NoMultimapped_SV)
  
  assign(paste("DEResults", GSE_ID, sep="_"), TempMasterResults, envir = as.environment(1))
  
  print(paste("Output: Named DEResults", GSE_ID, sep="_"))
  
  rm(TempMasterResults, TempResultsJoined_NoNA_NoMultimapped_SV, TempResultsJoined_NoNA_NoMultimapped_SE, TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol, TempResultsJoined_NoNA_NoMultimapped_FoldChange_Average, TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol, TempResultsJoined_NoNA_NoMultimapped_Tstat_Average)
  
}


#Notes about parameters for function CollapsingDEResults_OneResultPerGene()
#GSE_ID is a string indicating the name of the Gemma dataset
#TempResultsJoined_NoNA_NoMultimapped is the data frame outputted by our previous function
#ComparisonsOfInterest is a character vector containing a list of group comparisons of interest within this dataset. Important: These group comparisons should be listed in exactly the same order as the order that you provide the column names for their associated Fold Change and Tstat output.
#NamesOfFoldChangeColumns is a list containing the columns of TempResultsJoined_NoNA_NoMultimapped containing the FoldChange results for your comparisons of interes, in the same order as the ComparisonsOfInterest vector
#NamesOfTstatColumns is a list containing the columns of TempResultsJoined_NoNA_NoMultimapped containing the Tstat results for your comparisons of interes, in the same order as the ComparisonsOfInterest vector

#Example function usage:

#CollapsingDEResults_OneResultPerGene(GSE_ID="GSE59070", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest=c("Stress8days_vs_Acute", "Stress13days_vs_Acute"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_8.d, TempResultsJoined_NoNA_NoMultimapped$FoldChange_13.d), NamesOfTstatColumns=list(TempResultsJoined_NoNA_NoMultimapped$Tstat_8.d, TempResultsJoined_NoNA_NoMultimapped$Tstat_13.d))

#######################

#CODE for reading in individual datasets

#######################

#GSE56028

setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Teaching/BrainDataAlchemy/Grinnell Interns 2022/2022Cohort_GoogleDrive/GoogleDriveDownload_20230508/ErinHernandez_Antidepressants_Hippocampus/Gemma_ReDownloaded_20230509/11625_GSE56028_diffExpAnalysis_264081")

list.files()
#[1] "analysis.results.txt"        "resultset_ID528627.data.txt" "resultset_ID528628.data.txt"

ReadingInGemmaDE(ResultSetFileNames = c("resultset_ID528628.data.txt"))
#[1] "Outputted object: TempResultsJoined"

str(TempResultsJoined)
# 'data.frame':	29156 obs. of  16 variables:
# $ Element_Name                                                                                      : int  10870762 10733944 10731195 10776207 10932813 10724593 10799397 10830050 10769441 10804710 ...
# $ Gene_Symbol                                                                                       : chr  "Echdc2" "Trim11" "Sfxn4" "" ...
# $ Gene_Name                                                                                         : chr  "enoyl CoA hydratase domain containing 2" "tripartite motif-containing 11" "sideroflexin 4" "" ...
# $ NCBI_ID                                                                                           : chr  "298381" "360534" "361778" "" ...
# $ FoldChange_7.3.chloro.6.methyl.5.5.dioxo.11H.benzo.c.2.1.benzothiazepin.11.yl.amino.heptanoic.acid: num  0.02819 -0.00414 0.1125 0.00704 0.1256 ...
# $ Tstat_7.3.chloro.6.methyl.5.5.dioxo.11H.benzo.c.2.1.benzothiazepin.11.yl.amino.heptanoic.acid     : num  0.289 -0.0401 1.515 0.0945 1.492 ...
# $ PValue_7.3.chloro.6.methyl.5.5.dioxo.11H.benzo.c.2.1.benzothiazepin.11.yl.amino.heptanoic.acid    : num  0.776 0.969 0.146 0.926 0.152 ...
# $ FoldChange_agomelatine                                                                            : num  -0.1182 -0.2393 -0.0266 -0.0471 -0.0744 ...
# $ Tstat_agomelatine                                                                                 : num  -1.212 -2.317 -0.358 -0.632 -0.883 ...
# $ PValue_agomelatine                                                                                : num  0.2404 0.0318 0.7244 0.5346 0.3881 ...
# $ FoldChange_fluoxetine                                                                             : num  0.0719 -0.1761 0.0745 -0.0827 0.0767 ...
# $ Tstat_fluoxetine                                                                                  : num  0.738 -1.705 1.003 -1.11 0.911 ...
# $ PValue_fluoxetine                                                                                 : num  0.47 0.104 0.328 0.281 0.374 ...
# $ FoldChange_imipramine                                                                             : num  -0.0861 -0.0851 0.1303 -0.0572 0.0903 ...
# $ Tstat_imipramine                                                                                  : num  -0.883 -0.824 1.755 -0.768 1.073 ...
# $ PValue_imipramine                                                                                 : num  0.3883 0.42 0.0953 0.452 0.2967 ...

FilteringDEResults_GoodAnnotation(TempResultsJoined)
# [1] "# of rows with missing NCBI annotation:"
# [1] 11803
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 11803
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 827
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 827
# [1] "# of rows with good annotation"
# [1] 16526
# [1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"

CollapsingDEResults_OneResultPerGene(GSE_ID="GSE56028", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest=c("Tianeptine_10mgPerKg_vs_Ctrl","Agomelatine_40mgPerKg_vs_Ctrl", "Fluoxetine_10mgPerKg_vs_Ctrl", "Imipramine_10mgPerKg_vs_Ctrl"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_7.3.chloro.6.methyl.5.5.dioxo.11H.benzo.c.2.1.benzothiazepin.11.yl.amino.heptanoic.acid, TempResultsJoined_NoNA_NoMultimapped$FoldChange_agomelatine, TempResultsJoined_NoNA_NoMultimapped$FoldChange_fluoxetine, TempResultsJoined_NoNA_NoMultimapped$FoldChange_imipramine), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_7.3.chloro.6.methyl.5.5.dioxo.11H.benzo.c.2.1.benzothiazepin.11.yl.amino.heptanoic.acid, TempResultsJoined_NoNA_NoMultimapped$Tstat_agomelatine, TempResultsJoined_NoNA_NoMultimapped$Tstat_fluoxetine, TempResultsJoined_NoNA_NoMultimapped$Tstat_imipramine)) 
# [1] "Double check that the vectors containing the two fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:"
# [1] "# of rows with unique NCBI IDs:"
# [1] 15724
# [1] "# of rows with unique Gene Symbols:"
# [1] 15724
# [1] "Dimensions of Fold Change matrix, averaged by gene symbol:"
# [1] 15724     4
# [1] "Output: Named DEResults_GSE56028"

#############

#GSE63469

setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Teaching/BrainDataAlchemy/Grinnell Interns 2022/2022Cohort_GoogleDrive/GoogleDriveDownload_20230508/ErinHernandez_Antidepressants_Hippocampus/Gemma_ReDownloaded_20230509/14882_GSE63469_diffExpAnalysis_263154")

list.files()
#[1] "analysis.results.txt"        "resultset_ID526947.data.txt"

ReadingInGemmaDE(ResultSetFileNames = c("resultset_ID526947.data.txt"))
#[1] "Outputted object: TempResultsJoined"

str(TempResultsJoined)
# 'data.frame':	45017 obs. of  10 variables:
# $ Element_Name                      : chr  "1419455_at" "1426606_at" "1450888_at" "1447575_at" ...
# $ Gene_Symbol                       : chr  "Il10rb" "Crtac1" "Napb" "" ...
# $ Gene_Name                         : chr  "interleukin 10 receptor, beta" "cartilage acidic protein 1" "N-ethylmaleimide sensitive fusion protein attachment protein beta" "" ...
# $ NCBI_ID                           : chr  "16155" "72832" "17957" "" ...
# $ FoldChange_30.kg.d.kg_venlafaxine : num  0.1156 -0.0263 -0.1049 0.011 -0.1363 ...
# $ Tstat_30.kg.d.kg_venlafaxine      : num  2.997 -0.426 -1.362 0.162 -1.712 ...
# $ PValue_30.kg.d.kg_venlafaxine     : num  0.028 0.687 0.228 0.877 0.144 ...
# $ FoldChange_100.kg.d.kg_venlafaxine: num  0.0399 0.0109 -0.0718 -0.2167 -0.2859 ...
# $ Tstat_100.kg.d.kg_venlafaxine     : num  1.035 0.177 -0.932 -3.204 -3.592 ...
# $ PValue_100.kg.d.kg_venlafaxine    : num  0.3454 0.8663 0.3916 0.022 0.0142 ...

FilteringDEResults_GoodAnnotation(TempResultsJoined)
# [1] "# of rows with missing NCBI annotation:"
# [1] 13487
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 13487
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 709
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 709
# [1] "# of rows with good annotation"
# [1] 30821
# [1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"

CollapsingDEResults_OneResultPerGene(GSE_ID="GSE63469", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest=c("Venlafaxine_30kgDkg_vs_Ctrl", "Venlafaxine_100kgDkg_vs_Ctrl"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_30.kg.d.kg_venlafaxine, TempResultsJoined_NoNA_NoMultimapped$FoldChange_100.kg.d.kg_venlafaxine), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_30.kg.d.kg_venlafaxine, TempResultsJoined_NoNA_NoMultimapped$Tstat_100.kg.d.kg_venlafaxine)) 
# [1] "Double check that the vectors containing the two fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:"
# [1] "# of rows with unique NCBI IDs:"
# [1] 18532
# [1] "# of rows with unique Gene Symbols:"
# [1] 18532
# [1] "Dimensions of Fold Change matrix, averaged by gene symbol:"
# [1] 18532     2
# [1] "Output: Named DEResults_GSE63469"

########################

#GSE27532

setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Teaching/BrainDataAlchemy/Grinnell Interns 2022/2022Cohort_GoogleDrive/GoogleDriveDownload_20230508/ErinHernandez_Antidepressants_Hippocampus/Gemma_ReDownloaded_20230509/2749_GSE27532_diffExpAnalysis_275286")

list.files()
#[1] "analysis.results.txt"        "resultset_ID546193.data.txt" "resultset_ID546194.data.txt"

ReadingInGemmaDE(ResultSetFileNames = c("resultset_ID546193.data.txt"))
#[1] "Outputted object: TempResultsJoined"

str(TempResultsJoined)
# 'data.frame':	25659 obs. of  7 variables:
# $ Element_Name          : chr  "ILMN_2616509" "ILMN_3120335" "ILMN_1225966" "ILMN_3132050" ...
# $ Gene_Symbol           : chr  "Arnt2" "Gabrb3" "Gipc1" "Rasgrf1" ...
# $ Gene_Name             : chr  "aryl hydrocarbon receptor nuclear translocator 2" "gamma-aminobutyric acid (GABA) A receptor, subunit beta 3" "GIPC PDZ domain containing family, member 1" "RAS protein-specific guanine nucleotide-releasing factor 1" ...
# $ NCBI_ID               : chr  "11864" "14402" "67903" "19417" ...
# $ FoldChange_desipramine: num  -0.1765 -0.0512 0.2352 0.059 0.0154 ...
# $ Tstat_desipramine     : num  -1.605 -0.523 4.344 0.444 0.289 ...
# $ PValue_desipramine    : num  0.1264 0.6078 0.000416 0.6623 0.7756 ...

FilteringDEResults_GoodAnnotation(TempResultsJoined)
# [1] "# of rows with missing NCBI annotation:"
# [1] 1007
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 1007
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 1025
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 1025
# [1] "# of rows with good annotation"
# [1] 23627
# [1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"

CollapsingDEResults_OneResultPerGene(GSE_ID="GSE27532", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest=c("Desipramine_7mgPerKg_vs_Ctrl"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_desipramine), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_desipramine)) 
# [1] "Double check that the vectors containing the two fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:"
# [1] "# of rows with unique NCBI IDs:"
# [1] 16378
# [1] "# of rows with unique Gene Symbols:"
# [1] 16378
# [1] "Dimensions of Fold Change matrix, averaged by gene symbol:"
# [1] 16378     1
# [1] "Output: Named DEResults_GSE27532"

########################

#FUNCTIONS and CODE for meta-analysis

setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Teaching/BrainDataAlchemy/Grinnell Interns 2022/2022Cohort_GoogleDrive/GoogleDriveDownload_20230508/ErinHernandez_Antidepressants_Hippocampus/RESULTS_MH_moreExclusions")

########################
#6) Combine together the relevant results from different studies into a single data-frame for the effect sizes (Log2FC) and a single data.frame for the sampling variances. 

#The Log2FC values are the first element in the differential expression result object for each study
#The gene symbols are the row.names - 75% of gene symbols for rats and mice are the same, so for simplicity sake, instead of using a gene orthology database to inform equivalency, we're just going to align the datasets by gene symbol.


AligningDatasets<-function(ListOfDEResults){
  
  MetaAnalysis_FoldChange_Dfs<-list()
  
  for(i in c(1:length(ListOfDEResults))){
    MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(x=row.names(ListOfDEResults[[i]][[1]]),ListOfDEResults[[i]][[1]], stringsAsFactors=FALSE)
  }
  
  print("MetaAnalysis_FoldChange_Dfs:")
  print(str(MetaAnalysis_FoldChange_Dfs))
  
  MetaAnalysis_FoldChanges<<-join_all(MetaAnalysis_FoldChange_Dfs, by="x", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("MetaAnalysis_FoldChanges:")
  print(str(MetaAnalysis_FoldChanges))
  
  MetaAnalysis_SV_Dfs<-list()
  
  for(i in c(1:length(ListOfDEResults))){
    MetaAnalysis_SV_Dfs[[i]]<-data.frame(x=row.names(ListOfDEResults[[i]][[4]]),ListOfDEResults[[i]][[4]], stringsAsFactors=FALSE)
  }
  
  print("MetaAnalysis_SV_Dfs:")
  print(str(MetaAnalysis_SV_Dfs))
  
  MetaAnalysis_SV<<-join_all(MetaAnalysis_SV_Dfs, by="x", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("MetaAnalysis_SV:")
  print(str(MetaAnalysis_SV))
  
  rm(MetaAnalysis_SV_Dfs, MetaAnalysis_FoldChange_Dfs)
}

#Example Usage;

ListOfDEResults<-list(DEResults_GSE27532, DEResults_GSE56028, DEResults_GSE63469)

#I just discovered that if I have datasets with *comparisons with the same name*, when I join the datasets some of those comparisons disappear. :(

AligningDatasets(ListOfDEResults)
# [1] "MetaAnalysis_FoldChange_Dfs:"
# List of 3
# $ :'data.frame':	16378 obs. of  2 variables:
#   ..$ x                           : chr [1:16378] "0610009B22Rik" "0610010K14Rik" "0610012G03Rik" "0610040J01Rik" ...
# ..$ Desipramine_7mgPerKg_vs_Ctrl: num [1:16378] 0.03548 -0.00465 0.1887 -0.09332 -0.00764 ...
# $ :'data.frame':	15724 obs. of  5 variables:
#   ..$ x                            : chr [1:15724] "A1bg" "A1cf" "A2m" "A3galt2" ...
# ..$ Tianeptine_10mgPerKg_vs_Ctrl : num [1:15724] 0.1322 0.0188 -0.1332 0.0397 0.0291 ...
# ..$ Agomelatine_40mgPerKg_vs_Ctrl: num [1:15724] 0.121 -0.0832 -0.088 -0.052 -0.0643 ...
# ..$ Fluoxetine_10mgPerKg_vs_Ctrl : num [1:15724] -0.01397 0.1506 -0.04246 0.00956 -0.06863 ...
# ..$ Imipramine_10mgPerKg_vs_Ctrl : num [1:15724] 0.0515 0.0602 -0.2672 0.0789 -0.0039 ...
# $ :'data.frame':	18532 obs. of  3 variables:
#   ..$ x                           : chr [1:18532] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010K14Rik" ...
# ..$ Venlafaxine_30kgDkg_vs_Ctrl : num [1:18532] 0.2384 0.0695 0.0605 0.0921 0.0248 ...
# ..$ Venlafaxine_100kgDkg_vs_Ctrl: num [1:18532] 0.0502 0.0408 0.1478 0.1317 0.0407 ...
# NULL
# [1] "MetaAnalysis_FoldChanges:"
# 'data.frame':	22772 obs. of  8 variables:
#   $ x                            : chr  "0610009B22Rik" "0610010K14Rik" "0610012G03Rik" "0610040J01Rik" ...
# $ Desipramine_7mgPerKg_vs_Ctrl : num  0.03548 -0.00465 0.1887 -0.09332 -0.00764 ...
# $ Tianeptine_10mgPerKg_vs_Ctrl : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Agomelatine_40mgPerKg_vs_Ctrl: num  NA NA NA NA NA NA NA NA NA NA ...
# $ Fluoxetine_10mgPerKg_vs_Ctrl : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Imipramine_10mgPerKg_vs_Ctrl : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Venlafaxine_30kgDkg_vs_Ctrl  : num  0.0695 0.0921 0.0248 -0.0467 0.1309 ...
# $ Venlafaxine_100kgDkg_vs_Ctrl : num  0.0408 0.1317 0.0407 0.1255 0.1296 ...
# NULL
# [1] "MetaAnalysis_SV_Dfs:"
# List of 3
# $ :'data.frame':	16378 obs. of  2 variables:
#   ..$ x                           : chr [1:16378] "0610009B22Rik" "0610010K14Rik" "0610012G03Rik" "0610040J01Rik" ...
# ..$ Desipramine_7mgPerKg_vs_Ctrl: num [1:16378] 0.00518 0.01285 0.01194 0.00369 0.00201 ...
# $ :'data.frame':	15724 obs. of  5 variables:
#   ..$ x                            : chr [1:15724] "A1bg" "A1cf" "A2m" "A3galt2" ...
# ..$ Tianeptine_10mgPerKg_vs_Ctrl : num [1:15724] 0.01222 0.01078 0.02587 0.00668 0.00968 ...
# ..$ Agomelatine_40mgPerKg_vs_Ctrl: num [1:15724] 0.01221 0.01077 0.02588 0.00668 0.00968 ...
# ..$ Fluoxetine_10mgPerKg_vs_Ctrl : num [1:15724] 0.01222 0.01077 0.02587 0.00668 0.00968 ...
# ..$ Imipramine_10mgPerKg_vs_Ctrl : num [1:15724] 0.01222 0.01077 0.02588 0.00668 0.00968 ...
# $ :'data.frame':	18532 obs. of  3 variables:
#   ..$ x                           : chr [1:18532] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010K14Rik" ...
# ..$ Venlafaxine_30kgDkg_vs_Ctrl : num [1:18532] 0.00133 0.00117 0.00114 0.00315 0.00244 ...
# ..$ Venlafaxine_100kgDkg_vs_Ctrl: num [1:18532] 0.00133 0.00117 0.00114 0.00315 0.00244 ...
# NULL
# [1] "MetaAnalysis_SV:"
# 'data.frame':	22772 obs. of  8 variables:
#   $ x                            : chr  "0610009B22Rik" "0610010K14Rik" "0610012G03Rik" "0610040J01Rik" ...
# $ Desipramine_7mgPerKg_vs_Ctrl : num  0.00518 0.01285 0.01194 0.00369 0.00201 ...
# $ Tianeptine_10mgPerKg_vs_Ctrl : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Agomelatine_40mgPerKg_vs_Ctrl: num  NA NA NA NA NA NA NA NA NA NA ...
# $ Fluoxetine_10mgPerKg_vs_Ctrl : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Imipramine_10mgPerKg_vs_Ctrl : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Venlafaxine_30kgDkg_vs_Ctrl  : num  0.00117 0.00315 0.00244 0.00192 0.00174 ...
# $ Venlafaxine_100kgDkg_vs_Ctrl : num  0.00117 0.00315 0.00244 0.00192 0.00174 ...
# NULL

#I'm going to double-check that the row.names are still in the same order:
cbind(MetaAnalysis_FoldChanges$x, MetaAnalysis_SV$x)[c(1:100),]
#Looks good
sum(MetaAnalysis_FoldChanges$x==MetaAnalysis_SV$x)
#[1] 22772
sum(MetaAnalysis_FoldChanges$x!=MetaAnalysis_SV$x)
#[1] 0
#Looks good.


####################################

#7) Make a correlation matrix to compare the overall results from the different studies. Further visualize the comparison using a hierarchically-clustered heatmap of the correlation matrix. Which studies seem to have similar results? 

#We can generally compare the differential expression associated with different datasets or variable comparisons using a scatterplot and correlation analysis.

#The code for making scatterplots is very similar to the boxplot code that we used earlier. It uses a y~x formula, and can include many of the same parameters (e.g., x and y labels, color)

#plot(MetaAnalysis_FoldChanges$GSE109315_StressResilient_Vs_Ctrl~MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Ctrl, ylab="Stress Resilient Log2FC", xlab="Stress Susceptible Log2FC" )

#Within this plot, each data point represents the differential expression results (log2FC) for a particular gene for two different comparisons: stress susceptible vs. no stress (x-axis) and stress resilient vs. no stress (y axis)

#From looking at this plot, we can see that, in general, if a gene shows a positive log2FC for the stress susceptible vs. no stress comparison (i.e., the stress susceptible group has greater log2 expression for that gene than the no stress group), then that gene is also likely to have a positive log2FC for the stress resilient vs. no stress comparison.

#Similarly, if a gene shows a negative log2FC for the stress susceptible vs. no stress comparison (i.e., the stress susceptible group has lower log2 expression for that gene than the no stress group), then that gene is also likely to have a negative log2FC for the stress resilient vs. no stress comparison.

#This means that the differential expression results associated with the stress susceptible and stress resilient comparisons are positively correlated - they show a similar direction of effect.

#We can illustrate this by adding a trendline to the plot:

#We use a linear regression model to calculate the intercept and slope for the linear relationship between the two variables using the y~x formula above:
#Trendline<-lm(MetaAnalysis_FoldChanges$GSE109315_StressResilient_Vs_Ctrl~MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Ctrl)
#And then add that line to our scatterplot:
#abline(Trendline, col=2, lwd=3)

#If we want to know whether that linear relationship is stronger than we might expect due to random chance:
#summary.lm(Trendline)
# Call:
#   lm(formula = MetaAnalysis_FoldChanges$GSE109315_StressResilient_Vs_Ctrl ~ 
#        MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Ctrl)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.2960 -0.2212 -0.0068  0.2179  4.7673 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                  0.008552   0.002710   3.156   0.0016 ** 
#   MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Ctrl 0.739991   0.006558 112.842   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.507 on 35178 degrees of freedom
# (3878 observations deleted due to missingness)
# Multiple R-squared:  0.2658,	Adjusted R-squared:  0.2657 
# F-statistic: 1.273e+04 on 1 and 35178 DF,  p-value: < 2.2e-16

#The estimate for intercept is where the trend line crosses the y-axis
#The estimate for "Stress Susceptible vs. Control" is the slope - how much of an increase in Log2FC you should expect in the "Stress Resilient vs. Control" if you see a one unit increase in Log2FC for "Stress Susceptible vs. Control" 
#The Pr(>|t|) is the p-value for that relationship, in this case it is smaller than R is willing to display (<2e-16)

#If we want to gene a normalized correlation coefficient for this relationship (ranging between -1 to 1, with -1 being a perfect negative correlation and +1 being a perfect positive correlation), we can run a correlation analysis:
#While running this correlation analysis, we have to tell R to ignore any rows of differential expression output that don't have Log2FC for one of our variables (use "pairwise complete observations")
#cor(MetaAnalysis_FoldChanges$GSE109315_StressResilient_Vs_Ctrl,MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Ctrl, use="pairwise.complete.obs")
#[1] 0.5155259


#So here's where things get particularly cool: 
#We can actually run a correlation analysis comparing the Log2FC for every set of differential expression results in our meta-analysis using a single line of code.
#This is called a correlation matrix.
#cor(as.matrix(MetaAnalysis_FoldChanges[,-1]), use="pairwise.complete.obs")

#I find that these correlation matrices can be a little easier to look at in Excel, so I often output them into a file:
write.csv(cor(as.matrix(MetaAnalysis_FoldChanges[,-1]), use="pairwise.complete.obs"), "CorMatrix_Log2FC.csv")

#In the output, each cell includes the correlation coefficient reflecting the similarity of the effects for the variable in the row and column.  
#So at the intersection of the row for "StressSusceptible_Vs_Control" and column "StressResilient_Vs_Control" we can see the correlation coefficient that we calculated earlier (0.5155)
#The intersection of a variable with itself creates a coefficient of 1 (i.e., identical)

#Looking at the full matrix, most of the correlation coefficients are very close to 0.
#The only correlation coefficients that are larger, positive numbers are comparisons that come from the same dataset originally. e.g., "StressSusceptible_Vs_Control" and  "StressResilient_Vs_Control"

#Disappointing, but not surprising.  The comparisons that come from the same dataset (e.g.) often reflect comparisons with the same reference group. Therefore, any random variation in the reference group is going to be artificially shared in the differential expression results for both comparisons.


#You can also visualize that correlation matrix using a hierarchically-clustered heatmap:
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[,-1]), use="pairwise.complete.obs"))
#In this heatmap, white/yellow indicates a more positive correlation
#The groups are placed in order by similarity, as determined by hierarchical clustering.
#The lines ("tree branches") on the left and top illustrate that similarity (clustering) using a "dendrogram"


#################################


#8) Run a meta-analysis using all of the effect sizes for each gene that has data from at least 2 studies. 


#We can only run a meta-analysis if there are differential expression results from more than one comparison.
#Since I know that the differential expression results from the same study (dataset) are artificially correlated, I would actually prefer that there are results from more than one dataset.

#How many genes satisfy this criteria?

#This code caculates the number of NAs in each row:
MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges, 1, function(y) sum(is.na(y)))

#I'm going to make a histogram of the results because I'm curious to see how they are distributed
hist(MetaAnalysis_FoldChanges_NAsPerRow)

#Or, since there are a limited number of integer answers (0-11), I could make a table of the results:
table(MetaAnalysis_FoldChanges_NAsPerRow)
#MetaAnalysis_FoldChanges_NAsPerRow
MetaAnalysis_FoldChanges_NAsPerRow
# 0     1     2     3     4     5     6 
# 11562  1356   496  2310  2886  2728  1434  

#Let's try running a meta-analysis using genes that were found in at least 7 sets of differential expression results
#Since there are 7 sets of differential expression results, that means that the genes that we are including need to have 0 or fewer NAs in their results
#I set this conservatively, because there are so few studies in this meta-analysis.
#1 NAs is too many

RunBasicMetaAnalysis<-function(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV){
  
  MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges, 1, function(y) sum(is.na(y)))
  
  print("Table of # of NAs per Row (Gene):")
  print(table(MetaAnalysis_FoldChanges_NAsPerRow))
  
  MetaAnalysis_FoldChanges_ForMeta<<-MetaAnalysis_FoldChanges[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  MetaAnalysis_SV_ForMeta<<-MetaAnalysis_SV[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  
  print("MetaAnalysis_FoldChanges_ForMeta:")
  print(str(MetaAnalysis_FoldChanges_ForMeta))
  
  #I'm going to make an empty matrix to store the results of my meta-analysis:
  metaOutput<-matrix(NA, length(MetaAnalysis_FoldChanges_ForMeta$x), 6)
  
  #And then run a loop that run's a meta-analysis on the differential expression results (columns 2-10) for each gene (row):
  for(i in c(1:length(MetaAnalysis_FoldChanges_ForMeta$x))){
    
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-1])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-1])
    
    #I added a function tryCatch that double-checks that the meta-analysis function (rma) doesn't produce errors (which breaks the loop):
    skip_to_next <- FALSE
    tryCatch(TempMeta<-rma(effect, var), error = function(e) {skip_to_next <<- TRUE})
    
    if(skip_to_next){}else{
      TempMeta<-rma(effect, var)
      metaOutput[i, 1]<-TempMeta$b #gives estimate Log2FC
      metaOutput[i, 2]<-TempMeta$se #gives standard error
      metaOutput[i, 3]<-TempMeta$pval #gives pval
      metaOutput[i, 4]<-TempMeta$ci.lb #gives confidence interval lower bound
      metaOutput[i, 5]<-TempMeta$ci.ub #gives confidence interval upper bound
      metaOutput[i, 6]<-NumberOfComparisons-sum(is.na(effect))#Number of comparisons with data
      rm(TempMeta)
    }
    rm(effect, var)
  }
  
  colnames(metaOutput)<-c("Log2FC_estimate", "SE", "pval", "CI_lb", "CI_ub", "Number_Of_Comparisons")
  row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta$x
  
  metaOutput<<-metaOutput
  return(metaOutput)
  
  print("metaOutput:")
  print(str(metaOutput))
  
  print("Top of metaOutput:")
  print(head(metaOutput))
  
  print("Bottom of metaOutput")
  print(tail(metaOutput))
  
}

#Example Usage:
NumberOfComparisons=7
CutOffForNAs=1
#I want at least 7 comparisons
#1 NA is too many

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

# [1] "Table of # of NAs per Row (Gene):"
# MetaAnalysis_FoldChanges_NAsPerRow
# 0     1     2     3     4     5     6 
# 11562  1356   496  2310  2886  2728  1434 
# [1] "MetaAnalysis_FoldChanges_ForMeta:"
# 'data.frame':	11562 obs. of  8 variables:
#   $ x                            : chr  "A4galt" "Aaas" "Aacs" "Aadac" ...
# $ Desipramine_7mgPerKg_vs_Ctrl : num  0.00882 0.03177 -0.02216 -0.06205 0.08798 ...
# $ Tianeptine_10mgPerKg_vs_Ctrl : num  0.0291 -0.0445 -0.0257 -0.0961 0.2043 ...
# $ Agomelatine_40mgPerKg_vs_Ctrl: num  -0.0643 -0.2562 0.0177 -0.1124 -0.0331 ...
# $ Fluoxetine_10mgPerKg_vs_Ctrl : num  -0.06863 -0.00253 -0.01577 -0.1021 0.1074 ...
# $ Imipramine_10mgPerKg_vs_Ctrl : num  -0.0039 -0.0874 -0.1441 -0.1393 0.00518 ...
# $ Venlafaxine_30kgDkg_vs_Ctrl  : num  0.2085 0.1111 -0.0447 -0.1497 -0.0331 ...
# $ Venlafaxine_100kgDkg_vs_Ctrl : num  0.2465 0.1304 -0.1482 0.0866 -0.1226 ...

str(metaOutput)
# num [1:11562, 1:6] 0.06022 -0.01543 -0.06567 -0.07498 -0.00186 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:11562] "A4galt" "Aaas" "Aacs" "Aadac" ...
# ..$ : chr [1:6] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

head(metaOutput)
# Log2FC_estimate         SE       pval       CI_lb       CI_ub Number_Of_Comparisons
# A4galt     0.060215268 0.04928678 0.22180877 -0.03638505  0.15681559                     7
# Aaas      -0.015428086 0.04761888 0.74594492 -0.10875938  0.07790321                     7
# Aacs      -0.065673071 0.02722289 0.01584698 -0.11902896 -0.01231718                     7
# Aadac     -0.074981686 0.03104035 0.01570846 -0.13581966 -0.01414371                     7
# Aadat     -0.001855932 0.04559224 0.96752932 -0.09121508  0.08750321                     7
# Aagab      0.014365668 0.02539447 0.57159737 -0.03540659  0.06413792                     7

tail(metaOutput)
# Log2FC_estimate         SE        pval       CI_lb        CI_ub Number_Of_Comparisons
# Zswim4      0.03188732 0.02941755 0.278384515 -0.02577003  0.089544665                     7
# Zswim5     -0.02784781 0.02089780 0.182672524 -0.06880675  0.013111128                     7
# Zswim8     -0.09537509 0.03564146 0.007451569 -0.16523106 -0.025519117                     7
# Zup1       -0.05477842 0.03025440 0.070203843 -0.11407596  0.004519123                     7
# Zw10        0.01668070 0.02166933 0.441428173 -0.02579041  0.059151812                     7
# Zyx         0.02567283 0.04867763 0.597912326 -0.06973357  0.121079233                     7

########################################

## Multiple Comparison corrections
#The following code applies two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli) 
#Meta-analysis output with adjusted p-values is then outputted along with effect size information.

#9) Correct the meta-analysis output to take into account the fact that we are running the statistical calculations many times and therefore have a heightened risk of false discovery (false discovery rate correction) 


#Let's functionalize it!
FalseDiscoveryCorrection<-function(metaOutput){
  
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,3], proc=c("BH"))
  
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  colnames(metaOutputFDR)[7]<-"FDR"
  
  metaOutputFDR<<-metaOutputFDR
  
  print("metaOutputFDR:")
  print(str(metaOutputFDR))
  
  write.csv(metaOutputFDR, "metaOutputFDR.csv")
  
  #a version of the output in order by p-value:
  metaOutputFDR_OrderbyPval<<-metaOutputFDR[order(metaOutputFDR[,3]),]
  
  #Let's write out a version of the output in order by p-value:
  write.csv(metaOutputFDR_OrderbyPval, "metaOutputFDR_orderedByPval_wHDRFData.csv")
  
  print("Do we have any genes that are statistically significant following false discovery rate correction?")
  print(sum(metaOutputFDR[,7]<0.10, na.rm=TRUE))
  
  print("What are the top results?")
  print(head(metaOutputFDR[order(metaOutputFDR[,3]),]))
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
}

#Example usage:

FalseDiscoveryCorrection(metaOutput)
# [1] "metaOutputFDR:"
# num [1:11562, 1:7] 0.06022 -0.01543 -0.06567 -0.07498 -0.00186 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:11562] "A4galt" "Aaas" "Aacs" "Aadac" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# NULL
# [1] "Do we have any genes that are statistically significant following false discovery rate correction?"
# [1] 797
# [1] "What are the top results?"
# Log2FC_estimate         SE         pval      CI_lb       CI_ub Number_Of_Comparisons          FDR
# Parm1       -0.1798498 0.02389052 5.149049e-14 -0.2266743 -0.13302524                     7 5.953331e-10
# Lpl          0.2138084 0.03015666 1.341870e-12  0.1547024  0.27291434                     7 7.757348e-09
# Gabrr2       0.1965106 0.02859386 6.309657e-12  0.1404677  0.25255358                     7 2.431742e-08
# Pla2g5       0.1593702 0.02375177 1.948623e-11  0.1128176  0.20592280                     7 5.632494e-08
# Akap8l      -0.1174562 0.01848935 2.116479e-10 -0.1536947 -0.08121775                     7 4.894145e-07
# Dusp1       -0.4415588 0.07101670 5.045755e-10 -0.5807489 -0.30236858                     7 9.723171e-07

############################

#10) Determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. 

row.names(metaOutputFDR_OrderbyPval)[c(1:100)]
# [1] "Parm1"    "Lpl"      "Gabrr2"   "Pla2g5"   "Akap8l"   "Dusp1"    "Grm8"     "Pla2g12a" "Smarca1"  "Flcn"     "Tmem150b" "Tpbg"     "Themis2"  "Krt73"    "Unc13c"   "Wnt4"    
# [17] "Kdm1b"    "Ccdc117"  "Creb5"    "Lancl3"   "Vgf"      "Sall2"    "Emilin1"  "Fmnl1"    "Omg"      "Ypel2"    "Clcn2"    "Atf4"     "Chst10"   "Zfp691"   "Acot1"    "Ucp1"    
# [33] "Isy1"     "Kansl3"   "Rap1gap2" "Ppp4r1"   "Isg20"    "Mrps28"   "Prkd2"    "Cbfa2t3"  "Hcar1"    "Rgs14"    "Bbs4"     "Rfc5"     "Ppp2r1b"  "Ptrh2"    "Zcchc7"   "Surf2"   
# [49] "Insrr"    "Ubxn2a"   "Adar"     "Snhg11"   "Hoxb3"    "Twist1"   "Il2rg"    "Ier3ip1"  "Mtm1"     "Slc2a4"   "Abcd4"    "Aard"     "Soat2"    "Pkn1"     "Myl12a"   "Lpar6"   
# [65] "Krt82"    "Trpm4"    "Cyb561"   "Atp2b2"   "Mpped1"   "Pptc7"    "Zfp830"   "Tgm4"     "Pkd1"     "Mcpt2"    "Il10ra"   "Il1r1"    "Celsr3"   "Fbxo9"    "Sptlc1"   "Nccrp1"  
# [81] "Zfp697"   "Ccdc59"   "Paox"     "Grk3"     "Madcam1"  "Zfp35"    "Gpr135"   "Mrpl10"   "Aqp9"     "Mertk"    "Tnni1"    "Polg"     "Mta1"     "Aldh1l2"  "Csrnp2"   "Lrp1"    
# [97] "Evl"      "Ufsp2"    "Gfra4"    "Ccl20"        

#Let's plot some of those top results!

#Quickly looking at the range of Log2FC values to figure out the limits for the x-axis for the forest plots:
hist(metaOutputFDR[,1], breaks=40)
#Range is mostly -0.2 to 0.2, but there are a few with Log2FC as big as -0.5-0.5

MakeForestPlots<-function(GeneSymbol){
  
  pdf(paste("ForestPlot_", GeneSymbol, ".pdf", sep=""), height=5, width=8)
  
  effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
  var<-as.numeric(MetaAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
  
  forest.rma(rma(effect, var),slab=colnames(MetaAnalysis_FoldChanges_ForMeta)[-1],  xlim=c(-3, 3))
  
  mtext(paste(GeneSymbol), line=-1.5, cex=2)
  dev.off()
}


MakeForestPlots("Parm1") #believable
MakeForestPlots("Lpl") #The venlafaxine study with the pooled samples has disproportionate effect on the results
MakeForestPlots("Gabrr2") #very believable
MakeForestPlots("Pla2g5")#The venlafaxine study with the pooled samples has disproportionate effect on the results
MakeForestPlots("Akap8l") #Very believable
MakeForestPlots("Dusp1")#Very believable
MakeForestPlots("Vgf")#believable
MakeForestPlots("Omg")#believable

#The fact that I added more datasets/comparisons seems to have made it so that the top results are more believable 
#Maybe the debugging fo the SE calculations helped too.

#In general, it seems to me like the results that have both a significant FDR & large effect size tend to be the most convincing 

#Here's a summary of the distribution for the Log2FCs
summary(metaOutputFDR[,1])
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
# -0.441559 -0.032601  0.000023 -0.000460  0.032454  0.311504        17 

#Let's see how many results are both statistically significant (using a more conservative FDR) and have a larger estimated Log2FC (>0.10 or <-0.10)
sum(metaOutputFDR[,7]<0.05 & abs(metaOutputFDR[,1])>0.1, na.rm=TRUE)
#[1] 224

#How about bigger?
sum(metaOutputFDR[,7]<0.05 & abs(metaOutputFDR[,1])>0.15, na.rm=TRUE)
#[1] 36

#What are their gene symbols?
row.names(metaOutputFDR_OrderbyPval)[metaOutputFDR_OrderbyPval[,7]<0.05 & abs(metaOutputFDR_OrderbyPval[,1])>0.15 & is.na(metaOutputFDR_OrderbyPval[,1])==FALSE]
# [1] "Parm1"    "Lpl"      "Gabrr2"   "Pla2g5"   "Dusp1"    "Grm8"     "Smarca1"  "Tpbg"     "Krt73"    "Unc13c"   "Ccdc117"  "Lancl3"   "Emilin1"  "Prkd2"    "Rgs14"    "Twist1"  
# [17] "Soat2"    "Cyb561"   "Nccrp1"   "Ccdc59"   "Tnni1"    "Defb1"    "Tac1"     "Dnajb1"   "Irf9"     "Krt31"    "Tnfrsf25" "P4ha1"    "Dusp6"    "Il20"     "Tamalin"  "Kcnq1"   
# [33] "Dpysl5"   "Slc17a6"  "Banp"     "Nxph3"   



#What if we look at the correlation between datasets again, but only using the top genes found in all datasets?
sum(metaOutputFDR[,7]<0.05, na.rm=TRUE)
#[1] 451

TopGenesInAllDatasets<-row.names(metaOutputFDR)[metaOutputFDR[,7]<0.05 & is.na(metaOutputFDR[,6])==FALSE]

cor(as.matrix(MetaAnalysis_FoldChanges[which(MetaAnalysis_FoldChanges$x%in%TopGenesInAllDatasets),-1]), use="pairwise.complete.obs")
# Desipramine_7mgPerKg_vs_Ctrl Tianeptine_10mgPerKg_vs_Ctrl Agomelatine_40mgPerKg_vs_Ctrl Fluoxetine_10mgPerKg_vs_Ctrl Imipramine_10mgPerKg_vs_Ctrl
# Desipramine_7mgPerKg_vs_Ctrl                     1.0000000                    0.4745928                     0.4037031                    0.4770679                    0.4009778
# Tianeptine_10mgPerKg_vs_Ctrl                     0.4745928                    1.0000000                     0.5385285                    0.7782184                    0.7396772
# Agomelatine_40mgPerKg_vs_Ctrl                    0.4037031                    0.5385285                     1.0000000                    0.5742822                    0.5655176
# Fluoxetine_10mgPerKg_vs_Ctrl                     0.4770679                    0.7782184                     0.5742822                    1.0000000                    0.6739490
# Imipramine_10mgPerKg_vs_Ctrl                     0.4009778                    0.7396772                     0.5655176                    0.6739490                    1.0000000
# Venlafaxine_30kgDkg_vs_Ctrl                      0.5229038                    0.6590241                     0.6284662                    0.7054354                    0.5891862
# Venlafaxine_100kgDkg_vs_Ctrl                     0.5756095                    0.6952051                     0.6774385                    0.7296030                    0.6407904
# Venlafaxine_30kgDkg_vs_Ctrl Venlafaxine_100kgDkg_vs_Ctrl
# Desipramine_7mgPerKg_vs_Ctrl                    0.5229038                    0.5756095
# Tianeptine_10mgPerKg_vs_Ctrl                    0.6590241                    0.6952051
# Agomelatine_40mgPerKg_vs_Ctrl                   0.6284662                    0.6774385
# Fluoxetine_10mgPerKg_vs_Ctrl                    0.7054354                    0.7296030
# Imipramine_10mgPerKg_vs_Ctrl                    0.5891862                    0.6407904
# Venlafaxine_30kgDkg_vs_Ctrl                     1.0000000                    0.8447345
# Venlafaxine_100kgDkg_vs_Ctrl                    0.8447345                    1.0000000

################

#Workspace saved as:
save.image("~/Dropbox (University of Michigan)/LaptopBackup_20221123/Teaching/BrainDataAlchemy/Grinnell Interns 2022/2022Cohort_GoogleDrive/GoogleDriveDownload_20230508/ErinHernandez_Antidepressants_Hippocampus/METHODS SECTION/MetaAnalysis_ReDoneMoreExclusions_MH_20230509.RData")
