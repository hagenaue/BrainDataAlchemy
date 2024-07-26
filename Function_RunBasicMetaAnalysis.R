#This code document includes the code for a function that is designed to run a basic meta-analysis of Log2FC and sampling variance values using our previously generated objects MetaAnalysis_FoldChanges & MetaAnalysis_SV
#Megan Hagenauer
#July 25 2024

######################

#Installing and loading relevant code packages:
install.packages("metafor")
library(metafor)
library(plyr) 

######################

#Function:

RunBasicMetaAnalysis<-function(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV){
  
  #The function first provides information about how many of the statistical contrasts have NA values as their differential expression results for each gene:
  MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges[,-c(1:3)], 1, function(y) sum(is.na(y)))
  
  print("Table of # of NAs per Row (Gene):")
  print(table(MetaAnalysis_FoldChanges_NAsPerRow))
  
  #Then any row (gene) that has too many NAs is removed from the analysis:
  MetaAnalysis_FoldChanges_ForMeta<<-MetaAnalysis_FoldChanges[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  MetaAnalysis_SV_ForMeta<<-MetaAnalysis_SV[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  
  print("MetaAnalysis_FoldChanges_ForMeta:")
  print(str(MetaAnalysis_FoldChanges_ForMeta))
  
  #I'm going to make an empty matrix to store the results of my meta-analysis:
  metaOutput<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), 6)
  
  #And then run a loop that run's a meta-analysis on the differential expression results (i.e., the columns that aren't annotation) for each gene (row):
  for(i in c(1:nrow(MetaAnalysis_FoldChanges_ForMeta))){
    
    #When pulling out the log2FC values and sampling variances (SV) for each gene, we use the function as.numeric to make sure they are in numeric matrix format because this is the required input format for the meta-analysis function that we will use:
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-c(1:3)])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-c(1:3)])
    
    #I added a function tryCatch that double-checks that the meta-analysis function (rma) doesn't produce errors (which breaks the loop):
    skip_to_next <- FALSE
    tryCatch(TempMeta<-rma(effect, var), error = function(e) {skip_to_next <<- TRUE})
    
    #If everything looks good, we move on to running the meta-analysis using a model that treats the variation in Log2FC across studies as random effects:
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
  
  #Naming the columns in our output:
  colnames(metaOutput)<-c("Log2FC_estimate", "SE", "pval", "CI_lb", "CI_ub", "Number_Of_Comparisons")
  
  #The row names for our output are the combined mouse-rat entrez ids: 
  row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta[,3]
  
  #We return this output back into our global environment
  metaOutput<<-metaOutput
  MetaAnalysis_Annotation<<-MetaAnalysis_FoldChanges_ForMeta[,c(1:3)]
  return(metaOutput)
  return(MetaAnalysis_Annotation)
  
  #... and provide the user with an update about the newly created object:
  
  print("metaOutput:")
  print(str(metaOutput))
  
  print("Top of metaOutput:")
  print(head(metaOutput))
  
  print("Bottom of metaOutput")
  print(tail(metaOutput))
  
}

######################

#Example Usage:

#NumberOfComparisons=5
#CutOffForNAs=2
#I have 5 comparisons
#2 NA is too many

#metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

######################

#Example Output:

# [1] "Table of # of NAs per Row (Gene):"
# MetaAnalysis_FoldChanges_NAsPerRow
# 0     1     2     3     4     5 
# 13355  3200  5059   277  5293   917 
# [1] "MetaAnalysis_FoldChanges_ForMeta:"
# 'data.frame':	16555 obs. of  8 variables:
#   $ Rat_EntrezGene.ID                 : chr  "114087" "191569" "246307" "65041" ...
# $ Mouse_EntrezGene.ID               : chr  "23825" "18585" "66514" "20480" ...
# $ MouseVsRat_EntrezGene.ID          : chr  "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# $ LPS_SubchronicAndAcute_vs_Vehicle : num  -0.0239 -0.0858 -0.0686 0.0891 0.0376 ...
# $ LPS_Acute_vs_Vehicle              : num  -0.0938 0.2524 -0.109 0.0788 -0.0475 ...
# $ LPS_Subchronic_vs_Vehicle         : num  0.0355 0.2735 0.0415 -0.0116 0.0299 ...
# $ LPS_Acute850ugPerKg_vs_Vehicle    : num  0.1754 0.1651 0.1878 -0.0146 0.103 ...
# $ LPS_Chronic_vs_Vehicle_AllStressed: num  0.223 -0.3824 -0.197 -0.0051 0.0393 ...
# NULL

#str(metaOutput)
# num [1:16555, 1:6] 0.0234 0.1315 -0.023 0.0437 0.0769 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:16555] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:6] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

#head(metaOutput)
# Log2FC_estimate         SE       pval        CI_lb      CI_ub Number_Of_Comparisons
# 23825_114087      0.02335239 0.04876404 0.63202016 -0.072223377 0.11892815                     5
# 18585_191569      0.13152931 0.06603932 0.04640599  0.002094618 0.26096399                     5
# 66514_246307     -0.02304066 0.05315084 0.66465470 -0.127214402 0.08113308                     5
# 20480_65041       0.04374446 0.03359984 0.19294219 -0.022110021 0.10959893                     5
# 13726_25437       0.07691445 0.04407249 0.08095348 -0.009466050 0.16329495                     5
# 16952_25380       0.04465958 0.10253423 0.66315765 -0.156303828 0.24562298                     5                   

#tail(metaOutput)
# Log2FC_estimate         SE       pval       CI_lb      CI_ub Number_Of_Comparisons
# 108168734_NA      0.11566278 0.07742137 0.13519162 -0.03608031 0.26740587                     4
# 108168883_NA      0.16641238 0.09078313 0.06679127 -0.01151929 0.34434405                     4
# 108168923_NA      0.07757892 0.16987291 0.64789533 -0.25536587 0.41052370                     4
# 108168987_NA     -0.13375235 0.07407534 0.07097679 -0.27893734 0.01143264                     4
# 108169023_NA     -0.06106329 0.04794113 0.20276476 -0.15502618 0.03289960                     4
# 113002583_NA      0.02982367 0.07843726 0.70377969 -0.12391053 0.18355788                     4
