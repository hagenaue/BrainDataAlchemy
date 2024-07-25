#This is code that provides functions for aligning our differential expression results across datasets
#Megan Hagenauer July 25 2024

#Goals:
#Each dataset has differential expression results from a slightly different list of genes
#Depending on the exact tissue dissected, the sensitivity of the transcriptional profiling platform, the representation on the transcriptional profiling platform (for microarray), and the experimental conditions
#The differential expression results from different datasets will also be in a slightly different order
#We want to align these results so that the differential expression results from each dataset are columns, with each row representing a different gene

#Reading in the plyr code package:
library(plyr)

################

#A function for aligning all of our rat differential expression results from different datasets into a single data frame for Log2FCs and sampling variances (SVs):

AligningRatDatasets<-function(ListOfRatDEResults){
  
  #Making an empty list to hold our results:
  Rat_MetaAnalysis_FoldChange_Dfs<-list()
  
  #Looping over all of the rat differential expression results:
  for(i in c(1:length(ListOfRatDEResults))){
    
    #Placing each of the log2FC results for each dataset into a single list
    #Each element in the list is formatted so that the rownames are Rat Entrez Gene ID and then there are columns containing the Log2FC for the differential expression results:
    Rat_MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(Rat_EntrezGene.ID=row.names(ListOfRatDEResults[[i]][[1]]),ListOfRatDEResults[[i]][[1]], stringsAsFactors=FALSE)
  }
  
  #Letting the user know the structure of the set of Log2FC dataframes that we are starting out with:
  print("Rat_MetaAnalysis_FoldChange_Dfs:")
  print(str(Rat_MetaAnalysis_FoldChange_Dfs))
  
  #Running "join all" on the list to align all of the results by Entrez Gene ID and make them a single data frame:
  Rat_MetaAnalysis_FoldChanges<<-join_all(Rat_MetaAnalysis_FoldChange_Dfs, by="Rat_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  #Letting the user know the structure of the dataframe that we just created:
  print("Rat_MetaAnalysis_FoldChanges:")
  print(str(Rat_MetaAnalysis_FoldChanges))
  
  #Doing the same steps for the sampling variances:
  
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
  
  #Cleaning up our environment to remove unneeded objects:
  rm(Rat_MetaAnalysis_SV_Dfs, Rat_MetaAnalysis_FoldChange_Dfs)
}

#Example Usage;

#ListOfRatDEResults<-list(DEResults_GSE205325)
# 
#AligningRatDatasets(ListOfRatDEResults)
# [1] "Rat_MetaAnalysis_FoldChange_Dfs:"
# List of 1
# $ :'data.frame':	17196 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID    : chr [1:17196] "24153" "24157" "24158" "24159" ...
# ..$ GSE205325_LPS_Chronic: num [1:17196] 0.3485 0.0288 -0.1887 -0.2126 0.1235 ...
# NULL
# [1] "Rat_MetaAnalysis_FoldChanges:"
# 'data.frame':	17196 obs. of  2 variables:
#   $ Rat_EntrezGene.ID    : chr  "24153" "24157" "24158" "24159" ...
# $ GSE205325_LPS_Chronic: num  0.3485 0.0288 -0.1887 -0.2126 0.1235 ...
# NULL
# [1] "Rat_MetaAnalysis_SV_Dfs:"
# List of 1
# $ :'data.frame':	17196 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID    : chr [1:17196] "24153" "24157" "24158" "24159" ...
# ..$ GSE205325_LPS_Chronic: num [1:17196] 0.0745 0.0229 0.0383 0.0257 0.0135 ...
# NULL
# [1] "Rat_MetaAnalysis_SV:"
# 'data.frame':	17196 obs. of  2 variables:
#   $ Rat_EntrezGene.ID    : chr  "24153" "24157" "24158" "24159" ...
# $ GSE205325_LPS_Chronic: num  0.0745 0.0229 0.0383 0.0257 0.0135 ...
# NULL


###########

#A function for aligning all of our mouse differential expression results from different datasets into a single data frame for Log2FCs and sampling variances (SVs):

#This function works the same way as the rat alignment function:

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

#ListOfMouseDEResults<-list(DEResults_GSE126678, DEResults_GSE181285)

#AligningMouseDatasets(ListOfMouseDEResults)

# [1] "Mouse_MetaAnalysis_FoldChange_Dfs:"
# List of 2
# $ :'data.frame':	21614 obs. of  4 variables:
#   ..$ Mouse_EntrezGene.ID              : chr [1:21614] "11287" "11298" "11302" "11303" ...
# ..$ GSE126678_LPS_Acute              : num [1:21614] 1.9397 0.0805 0.0595 0.0306 0.276 ...
# ..$ GSE126678_LPS_SubchronicPlusAcute: num [1:21614] 1.2967 -0.0472 -0.1459 0.1367 1.5651 ...
# ..$ GSE126678_LPS_Subchronic         : num [1:21614] 0.0582 0.203 -0.1144 0.1361 -0.0051 ...
# $ :'data.frame':	18563 obs. of  2 variables:
#   ..$ Mouse_EntrezGene.ID: chr [1:18563] "100008567" "100009600" "100012" "100017" ...
# ..$ GSE181285_LPS_Acute: num [1:18563] 0.0198 0.0225 0.0641 -0.0049 -0.0588 ...
# NULL
# [1] "Mouse_MetaAnalysis_FoldChanges:"
# 'data.frame':	24287 obs. of  5 variables:
#   $ Mouse_EntrezGene.ID              : chr  "11287" "11298" "11302" "11303" ...
# $ GSE126678_LPS_Acute              : num  1.9397 0.0805 0.0595 0.0306 0.276 ...
# $ GSE126678_LPS_SubchronicPlusAcute: num  1.2967 -0.0472 -0.1459 0.1367 1.5651 ...
# $ GSE126678_LPS_Subchronic         : num  0.0582 0.203 -0.1144 0.1361 -0.0051 ...
# $ GSE181285_LPS_Acute              : num  -0.042 -0.0368 -0.0534 0.1067 -0.5258 ...
# NULL
# [1] "Mouse_MetaAnalysis_SV_Dfs:"
# List of 2
# $ :'data.frame':	21614 obs. of  4 variables:
#   ..$ Mouse_EntrezGene.ID              : chr [1:21614] "11287" "11298" "11302" "11303" ...
# ..$ GSE126678_LPS_Acute              : num [1:21614] 0.62127 0.14737 0.00437 0.01624 1.34369 ...
# ..$ GSE126678_LPS_SubchronicPlusAcute: num [1:21614] 0.66434 0.14559 0.00438 0.01567 0.98359 ...
# ..$ GSE126678_LPS_Subchronic         : num [1:21614] 0.84004 0.14211 0.00439 0.01592 1.40671 ...
# $ :'data.frame':	18563 obs. of  2 variables:
#   ..$ Mouse_EntrezGene.ID: chr [1:18563] "100008567" "100009600" "100012" "100017" ...
# ..$ GSE181285_LPS_Acute: num [1:18563] 0.11456 0.01612 0.00329 0.00487 0.00719 ...
# NULL
# [1] "Mouse_MetaAnalysis_SV:"
# 'data.frame':	24287 obs. of  5 variables:
#   $ Mouse_EntrezGene.ID              : chr  "11287" "11298" "11302" "11303" ...
# $ GSE126678_LPS_Acute              : num  0.62127 0.14737 0.00437 0.01624 1.34369 ...
# $ GSE126678_LPS_SubchronicPlusAcute: num  0.66434 0.14559 0.00438 0.01567 0.98359 ...
# $ GSE126678_LPS_Subchronic         : num  0.84004 0.14211 0.00439 0.01592 1.40671 ...
# $ GSE181285_LPS_Acute              : num  0.00391 0.02738 0.00601 0.0101 0.03332 ...
# NULL

#####################
