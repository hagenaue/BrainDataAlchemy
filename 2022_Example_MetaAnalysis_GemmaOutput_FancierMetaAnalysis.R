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
save.image("~/Documents/Teaching/Grinnell Interns 2022/MetaAnalysis_HC_Stress_wHDRF_andDG_BugFixed/Workspace_MetaAnalysis_ChronicStress_HC_DG_wCovariates.RData")

####################

#Creating the predictor variables of interest
#This whole section of code would need to be adapted to fit the research question of interest

#In my meta-analysis, I have included data from whole hippocampal dissections as well as dentate gyrus dissections
#Dissection tends to be important, so I would like to add it as a predictor to my meta-analysis model

#First, I need to create a vector that identifies the dissection type for each of my datasets/comparisons:
colnames(MetaAnalysis_FoldChanges)
# [1] "x"                                   "GSE109315_StressResilient_Vs_Ctrl"  
# [3] "GSE109315_StressSusceptible_Vs_Ctrl" "GSE116009_CUMS_vs_Control"          
# [5] "GSE132819_CSDS_vs_Ctrl"              "GSE151807_CMS_vs_Ctrl"              
# [7] "GSE56028_CMS_vs_Ctrl"                "GSE59070_Stress8days_vs_Acute"      
# [9] "GSE59070_Stress13days_vs_Acute"      "GSE81672_StressResistent_vs_Ctrl"   
# [11] "GSE81672_StressSusceptible_vs_Ctrl"  "GSE84183_CMS_vs_Ctrl"   

#These are the DG datasets:
#GSE56028 - unpredictable chronic mild stress (CUMS), Dentate Gyrus (DG)
#GSE132819 - CSDS, DG
#GSE84183 - CUMS, DG

Dissection<-as.factor(c("HC", "HC", "HC", "DG", "HC", "DG", "HC", "HC", "HC", "HC", "DG"))
str(Dissection)
#Factor w/ 2 levels "DG","HC": 2 2 2 1 2 1 2 2 2 2 ...
Dissection<-relevel(Dissection, ref="HC")
str(Dissection)
#Factor w/ 2 levels "HC","DG": 1 1 1 2 1 2 1 1 1 1 ...

#Double checking that I assigned the correct dissection to each dataset:

cbind(colnames(MetaAnalysis_FoldChanges)[-1], Dissection)
#Dissection
# [1,] "GSE109315_StressResilient_Vs_Ctrl"   "1"       
# [2,] "GSE109315_StressSusceptible_Vs_Ctrl" "1"       
# [3,] "GSE116009_CUMS_vs_Control"           "1"       
# [4,] "GSE132819_CSDS_vs_Ctrl"              "2"       
# [5,] "GSE151807_CMS_vs_Ctrl"               "1"       
# [6,] "GSE56028_CMS_vs_Ctrl"                "2"       
# [7,] "GSE59070_Stress8days_vs_Acute"       "1"       
# [8,] "GSE59070_Stress13days_vs_Acute"      "1"       
# [9,] "GSE81672_StressResistent_vs_Ctrl"    "1"       
# [10,] "GSE81672_StressSusceptible_vs_Ctrl"  "1"       
# [11,] "GSE84183_CMS_vs_Ctrl"                "2"   

#Another variable that might be an interesting source of heterogeneity in the data is type of stress.
TypeOfStress<-as.factor(c("CSDS", "CSDS", "CMS", "CSDS", "CMS", "CMS", "CSDS", "CSDS", "CSDS", "CSDS", "CMS"))
str(TypeOfStress)
#Factor w/ 2 levels "CMS","CSDS": 2 2 1 2 1 1 2 2 2 2 ...

cbind(colnames(MetaAnalysis_FoldChanges)[-1], Dissection, TypeOfStress)

                                      #Dissection TypeOfStress
# [1,] "GSE109315_StressResilient_Vs_Ctrl"   "1"        "2"         
# [2,] "GSE109315_StressSusceptible_Vs_Ctrl" "1"        "2"         
# [3,] "GSE116009_CUMS_vs_Control"           "1"        "1"         
# [4,] "GSE132819_CSDS_vs_Ctrl"              "2"        "2"         
# [5,] "GSE151807_CMS_vs_Ctrl"               "1"        "1"         
# [6,] "GSE56028_CMS_vs_Ctrl"                "2"        "1"         
# [7,] "GSE59070_Stress8days_vs_Acute"       "1"        "2"         
# [8,] "GSE59070_Stress13days_vs_Acute"      "1"        "2"         
# [9,] "GSE81672_StressResistent_vs_Ctrl"    "1"        "2"         
# [10,] "GSE81672_StressSusceptible_vs_Ctrl"  "1"        "2"         
# [11,] "GSE84183_CMS_vs_Ctrl"                "2"        "1"  

#Example of numeric predictors:

DurationOfSleepDep<-c(12,6,5,3,6,6,10,6,18)

#The intercept will be artificially set at 0 for this variable, which can make interpreting results confusing.
#So sometimes people subtract the "center" of the variable to make so that your values span the y-axis

#center of your data:
mean(DurationOfSleepDep)
#[1] 8

DurationOfSleepDep_Centered<-DurationOfSleepDep-8

DurationOfSleepDep_Centered
#[1]  4 -2 -3 -5 -2 -2  2 -2 10
#So all values that were less than 8 are negative, all values that were more than 8 are positive

DurationOfRecovery<-c(0,1,4,0,0,6,12,0,2)

mean(DurationOfRecovery)
#[1] 2.777778

#I might round just to make this more interpretable

DurationOfRecovery_Centered<-DurationOfRecovery-3

DurationOfRecovery_Centered
#[1] -3 -2  1 -3 -3  3  9 -3 -1

#####################

#Adapting the Meta-Analysis code to handle additional predictors:
library(metafor)

#Example Usage:
NumberOfComparisons=11
CutOffForNAs=5

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
  effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-1])
  var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-1])
  TempMeta<-rma(yi=effect~Dissection+TypeOfStress, vi=var)
  # Warning messages:
  #   1: In rma(yi = effect ~ Dissection + TypeOfStress, vi = var) :
  #   Studies with NAs omitted from model fitting.
  # 2: In rma(yi = effect ~ Dissection + TypeOfStress, vi = var) :
  #   Redundant predictors dropped from the model.
  TempMeta
  # Mixed-Effects Model (k = 7; tau^2 estimator: REML)
  # 
  # tau^2 (estimated amount of residual heterogeneity):     0 (SE = 0.0416)
  # tau (square root of estimated tau^2 value):             0
  # I^2 (residual heterogeneity / unaccounted variability): 0.00%
  # H^2 (unaccounted variability / sampling variability):   1.00
  # R^2 (amount of heterogeneity accounted for):            100.00%
  # 
  # Test for Residual Heterogeneity: 
  #   QE(df = 5) = 2.8510, p-val = 0.7229
  # 
  # Test of Moderators (coefficient(s) 2): 
  #   QM(df = 1) = 5.2288, p-val = 0.0222
  # 
  # Model Results:
  #   
  #                   estimate      se     zval    pval    ci.lb    ci.ub    
  #   intrcpt            -0.2051  0.0645  -3.1810  0.0015  -0.3315  -0.0787  **
  #   TypeOfStressCSDS    0.3048  0.1333   2.2867  0.0222   0.0435   0.5661   *
  #   
  #   ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
  
  #Alternative coding: TempMeta<-rma(yi=effect,vi=var, mods = cbind(Dissection, TypeOfStress))
  
  
  #hmmm.... that's an interesting little problem:  
  #This only works if the NA data.points don't contain all of the data for one of the comparison groups (e.g., all of the DG experiments)
  

  #Soooo....As long as there is guaranteed enough data to represent each predictor, we can loop things like before, but with more things in our storage matrix:
  
  #Actual code for loop:
  
  #This could be changed to any two predictors:
  Predictor1<-Dissection
  Predictor2<-TypeOfStress
  
  #I'm going to make an empty matrix to store the results of my meta-analysis:
  metaOutput<-matrix(NA, length(MetaAnalysis_FoldChanges_ForMeta$x), 16)
  
  #And then run a loop that run's a meta-analysis on the differential expression results (columns 2-10) for each gene (row):
  for(i in c(1:length(MetaAnalysis_FoldChanges_ForMeta$x))){
    
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-1])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-1])
    
    #I added a function tryCatch that double-checks that the meta-analysis function (rma) doesn't produce errors (which breaks the loop):
    skip_to_next <- FALSE
    tryCatch(TempMeta<-rma(yi=effect~Predictor1+Predictor2, vi=var), error = function(e) {skip_to_next <<- TRUE})
    
    if(skip_to_next){}else{
      TempMeta<-rma(yi=effect~Predictor1+Predictor2, vi=var)
      metaOutput[i, c(1:3)]<-TempMeta$b #gives estimate Log2FC
      metaOutput[i, c(4:6)]<-TempMeta$se #gives standard error
      metaOutput[i, c(7:9)]<-TempMeta$pval #gives pval
      metaOutput[i, c(10:12)]<-TempMeta$ci.lb #gives confidence interval lower bound
      metaOutput[i, c(13:15)]<-TempMeta$ci.ub #gives confidence interval upper bound
      metaOutput[i, 16]<-NumberOfComparisons-sum(is.na(effect))#Number of comparisons with data
      rm(TempMeta)
    }
    rm(effect, var)
  }
  
  #The column names need updating:
  colnames(metaOutput)<-c("Main_Log2FC_estimate", "Predictor1_Log2FC_estimate", "Predictor2_Log2FC_estimate", "Main_SE", "Predictor1_SE", "Predictor2_SE",  "Main_pval", "Predictor1_pval","Predictor2_pval", "Main_CI_lb","Predictor1_CI_lb","Predictor2_CI_lb", "Main_CI_ub", "Predictor1_CI_ub","Predictor2_CI_ub","Number_Of_Comparisons")
  row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta$x
  
  metaOutput<<-metaOutput
  return(metaOutput)
  
  print("metaOutput:")
  print(str(metaOutput))
  
  print("Top of metaOutput:")
  print(head(metaOutput))
  
  print("Bottom of metaOutput")
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
  
  #For main effect (column 3)
  
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,7], proc=c("BH"))
  
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  colnames(metaOutputFDR)[17]<-"Main_FDR"
  
 #For predictor 1 (column4)
  
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,8], proc=c("BH"))
  
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  metaOutputFDR<-cbind(metaOutputFDR, metaPvalAdj[,2])
  
  colnames(metaOutputFDR)[18]<-"Predictor1_FDR"
  
  #For predictor 2 (column5)
  
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,9], proc=c("BH"))
  
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  metaOutputFDR<-cbind(metaOutputFDR, metaPvalAdj[,2])
  
  colnames(metaOutputFDR)[19]<-"Predictor1_FDR"
  
  
  print("metaOutputFDR:")
  print(str(metaOutputFDR))
  
  write.csv(metaOutputFDR, "metaOutputFDR.csv")
  
  #a version of the output in order by p-value for the main effect:
  metaOutputFDR_OrderbyPval<-metaOutputFDR[order(metaOutputFDR[,7]),]
  
  #Let's write out a version of the output in order by p-value:
  write.csv(metaOutputFDR_OrderbyPval, "metaOutputFDR_orderedByPval_wHDRFData.csv")
  

  

