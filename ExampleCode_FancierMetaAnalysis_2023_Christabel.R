#Example: Running a fancier meta-analysis - controlling for co-variates/additional predictor variables
#Megan Hagenauer
#2022-08-03

#This code assumes that you have already run the beginning of the code for a basic meta-analysis from the code document:
#Example_MetaAnalysis_GemmaOutput_ComparingDatasets_DebuggedMaybe.R
#..up through aligning datasets, and therefore have two objects available to you named:
#MetaAnalysis_FoldChanges, MetaAnalysis_SV
#This code would then replace this section:
#8) Run a meta-analysis using all of the effect sizes for each gene that has data from at least 2 studies. 

#Current workspace:

load("/Users/christabelmclain/Documents/Mount Sinai/PhD/Results/BrainData/Analysis/DEWorkspace.RData")

####################

#Creating the predictor variables of interest
#This whole section of code would need to be adapted to fit the research question of interest


#First, I need to create a vector that identifies the dissection type for each of my datasets/comparisons:
colnames(MetaAnalysis_FoldChanges)
# [1] "Rat_EntrezGene.ID"             "Mouse_EntrezGene.ID"           "MouseVsRat_EntrezGene.ID"      "GSE85136.Male_StressvsCtrl"   
# [5] "GSE85136.Female_StressvsCtrl"  "GSE90962.Male_StressvsCtrl"    "GSE90962.Female_StressvsCtrl"  "GSE102556.Male_StressvsCtrl"  
# [9] "GSE102556.Female_StressvsCtrl"


Sex<-as.factor(c("M", "F","M", "F","M", "F"))
str(Sex)
#Factor w/ 2 levels "DG","HC": 2 2 2 1 2 1 2 2 2 2 ...
Sex<-relevel(Sex, ref="M")
str(Sex)
#Factor w/ 2 levels "HC","DG": 1 1 1 2 1 2 1 1 1 1 ...

#Double checking that I assigned the correct sex to each dataset:

cbind(colnames(MetaAnalysis_FoldChanges)[-c(1:3)], Sex)
 

#####################

#Adapting the Meta-Analysis code to handle additional predictors:
library(metafor)

#Example Usage:
NumberOfComparisons=6
CutOffForNAs=2

#This code is the same as what used to be within our previous function:

  MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges, 1, function(y) sum(is.na(y)))
  
  print("Table of # of NAs per Row (Gene):")
  print(table(MetaAnalysis_FoldChanges_NAsPerRow))
  
  MetaAnalysis_FoldChanges_ForMeta<<-MetaAnalysis_FoldChanges[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  MetaAnalysis_SV_ForMeta<<-MetaAnalysis_SV[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  
  print("MetaAnalysis_FoldChanges_ForMeta:")
  print(str(MetaAnalysis_FoldChanges_ForMeta))
  
  
#This part of the code needs to be customized to fit the new number of parameters in the model:
  
  #Test run:
  #To figure out how to adapt it, I'm going to try out the meta-analysis with one row of data first (one gene):
  
  i<-2
  effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-c(1:3)])
  var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-c(1:3)])
  TempMeta<-rma(yi=effect~Sex, vi=var)

  TempMeta
  
  #Actual code for loop:
  
  #This could be changed to any two predictors:
  Predictor1<-Sex
  
  #I'm going to make an empty matrix to store the results of my meta-analysis:
  metaOutput<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), 11)
  
  #And then run a loop that run's a meta-analysis on the differential expression results (columns 2-10) for each gene (row):
  for(i in c(1:nrow(MetaAnalysis_FoldChanges_ForMeta))){
    
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-c(1:3)])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-c(1:3)])
    
    #I added a function tryCatch that double-checks that the meta-analysis function (rma) doesn't produce errors (which breaks the loop):
    skip_to_next <- FALSE
    tryCatch(TempMeta<-rma(yi=effect~Predictor1, vi=var), error = function(e) {skip_to_next <<- TRUE})
    
    if(skip_to_next){}else{
      TempMeta<-rma(yi=effect~Predictor1, vi=var)
      metaOutput[i, c(1:2)]<-TempMeta$b #gives estimate Log2FC
      metaOutput[i, c(3:4)]<-TempMeta$se #gives standard error
      metaOutput[i, c(5:6)]<-TempMeta$pval #gives pval
      metaOutput[i, c(7:8)]<-TempMeta$ci.lb #gives confidence interval lower bound
      metaOutput[i, c(9:10)]<-TempMeta$ci.ub #gives confidence interval upper bound
      metaOutput[i, 11]<-NumberOfComparisons-sum(is.na(effect))#Number of comparisons with data
      rm(TempMeta)
    }
    rm(effect, var)
  }
  
  #The column names need updating:
  colnames(metaOutput)<-c("Main_Log2FC_estimate", "Predictor1_Log2FC_estimate", "Main_SE", "Predictor1_SE", "Main_pval", "Predictor1_pval","Main_CI_lb","Predictor1_CI_lb", "Main_CI_ub", "Predictor1_CI_ub","Number_Of_Comparisons")
  row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID
  
  print(str(metaOutput))
  
  print(head(metaOutput))
  
  print(tail(metaOutput))
  
  
  #################################
  
  
colnames(metaOutput)
  # [1] "Main_Log2FC_estimate"       "Predictor1_Log2FC_estimate" "Predictor2_Log2FC_estimate"
  # [4] "Main_SE"                    "Predictor1_SE"              "Predictor2_SE"             
  # [7] "Main_pval"                  "Predictor1_pval"            "Predictor2_pval"           
  # [10] "Main_CI_lb"                 "Predictor1_CI_lb"           "Predictor2_CI_lb"          
  # [13] "Main_CI_ub"                 "Predictor1_CI_ub"           "Predictor2_CI_ub"          
  # [16] "Number_Of_Comparisons"    
  
  #FDR correction:
  
  #For main effect (column 1??)
  
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,5], proc=c("BH"))
  
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  colnames(metaOutputFDR)[12]<-"Main_FDR"
  
 #For predictor 1 (column 2??)
  
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,6], proc=c("BH"))
  
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  metaOutputFDR<-cbind(metaOutputFDR, metaPvalAdj[,2])
  
  colnames(metaOutputFDR)[13]<-"Predictor1_FDR"
  
  
  print("metaOutputFDR:")
  print(str(metaOutputFDR))
  
  write.csv(metaOutputFDR, "metaOutputFDR.csv")
  
  #a version of the output in order by p-value for the main effect:
  metaOutputFDR_OrderbyPval<-metaOutputFDR[order(metaOutputFDR[,5]),]
  
  #Let's write out a version of the output in order by p-value:
  write.csv(metaOutputFDR_OrderbyPval, "metaOutputFDR_orderedByPval_wHDRFData.csv")
  

  ## sex log2FC correlate with interaction term ?
  ## forrest plots for some genes
  
  MakeForestPlots<-function(EntrezIDAsCharacter){
    
    pdf(paste("ForestPlot_", EntrezIDAsCharacter, ".pdf", sep=""), height=5, width=8)
    
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID==EntrezIDAsCharacter,-c(1:3)])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID==EntrezIDAsCharacter,-c(1:3)])
    
    forest.rma(rma(effect, var),slab=colnames(MetaAnalysis_FoldChanges_ForMeta)[-c(1:3)],  xlim=c(-3, 3))
    
    mtext(paste(EntrezIDAsCharacter), line=-1.5, cex=2)
    dev.off()
  }
  
  
  

