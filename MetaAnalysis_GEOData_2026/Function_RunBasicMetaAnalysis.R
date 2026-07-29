#This code document includes the code for a function that is designed to run a basic meta-analysis of Log2FC and sampling variance values using our previously generated objects MetaAnalysis_FoldChanges & MetaAnalysis_SV
#Megan Hagenauer
#Original version: July 25 2024
#In response to reviewers' comments, this function has been updated to include heterogeneity statistics, publication bias statistics, and robustness statistics
#Updated version: March 10, 2026
#Updated again: July 28, 2026 (to make the code more generalizable for the full 2026 cohort)

#Function:

RunBasicMetaAnalysis<-function(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV, Columns_DE, Column_GeneName){
  
  #The function first provides information about how many of the statistical contrasts have NA values as their differential expression results for each gene:
  MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges[,Columns_DE], 1, function(y) sum(is.na(y)))
  
  print("Table of # of NAs per Row (Gene):")
  print(table(MetaAnalysis_FoldChanges_NAsPerRow))
  
  #Then any row (gene) that has too many NAs is removed from the analysis:
  MetaAnalysis_FoldChanges_ForMeta<<-MetaAnalysis_FoldChanges[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  MetaAnalysis_SV_ForMeta<<-MetaAnalysis_SV[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  
  print("MetaAnalysis_FoldChanges_ForMeta:")
  print(str(MetaAnalysis_FoldChanges_ForMeta))
  
  #I'm going to make an empty matrix to store the results of my meta-analysis:
  #This matrix was originally 6 columns
  #It was made larger to incorporate heterogeneity statistics (8 stats)
  #Then it was made larger to incorporate publication bias (3 stats) 
  #...And even larger to incorporate robustness stats (4 stats)
  metaOutput<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), 21)
  
  influence_dfbs<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), ncol(MetaAnalysis_FoldChanges_ForMeta[Columns_DE]))
  influence_cookd<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), ncol(MetaAnalysis_FoldChanges_ForMeta[Columns_DE]))
  influence_TF<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), ncol(MetaAnalysis_FoldChanges_ForMeta[Columns_DE]))
  
  colnames(influence_dfbs)<-colnames(MetaAnalysis_FoldChanges_ForMeta)[Columns_DE]
  colnames(influence_cookd)<-colnames(MetaAnalysis_FoldChanges_ForMeta)[Columns_DE]
  colnames(influence_TF)<-colnames(MetaAnalysis_FoldChanges_ForMeta)[Columns_DE]
  
  row.names(influence_dfbs)<-MetaAnalysis_FoldChanges_ForMeta[,Column_GeneName]
  row.names(influence_cookd)<-MetaAnalysis_FoldChanges_ForMeta[,Column_GeneName]
  row.names(influence_TF)<-MetaAnalysis_FoldChanges_ForMeta[,Column_GeneName]
  
  #And then run a loop that run's a meta-analysis on the differential expression results (i.e., the columns that aren't annotation) for each gene (row):
  for(i in c(1:nrow(MetaAnalysis_FoldChanges_ForMeta))){
    
    print(i)
    
    #When pulling out the log2FC values and sampling variances (SV) for each gene, we use the function as.numeric to make sure they are in numeric matrix format because this is the required input format for the meta-analysis function that we will use:
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,Columns_DE])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,Columns_DE])
    
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
      metaOutput[i, 7]<-TempMeta$k #Metafor output: number of studies (contrasts) included in the analysis - sanity check, should be the same as column 6
      metaOutput[i, 8]<-TempMeta$p #Metafor output: number of coefficients in model
      metaOutput[i, 9]<-TempMeta$tau2 #estimated amount of (residual) heterogeneity
      metaOutput[i, 10]<-TempMeta$se.tau2 #SE of the estimated amount of (residual) heterogeneity
      metaOutput[i, 11]<-TempMeta$QE #test statistic of the test for (residual) heterogeneity (Cochran’s Q-test)
      metaOutput[i, 12]<-TempMeta$QEp #p-value of the test for (residual) heterogeneity (Cochran’s Q-test)
      metaOutput[i, 13]<-TempMeta$I2 #the I 2 statistic, which estimates (in percent) how much of the total variability in the observed effect sizes or outcomes can be attributed to heterogeneity among the true effects
      metaOutput[i, 14]<-TempMeta$H2 #the H2 statistic, which estimates the ratio of the total amount of variability in the observed effect sizes or outcomes to the amount of sampling variability
      
      #Testing for evidence of publication bias (funnel plot asymmetry)
      skip_to_next <- FALSE
      tryCatch(PubBias<-regtest(TempMeta), error = function(e) {skip_to_next <<- TRUE})
      if(skip_to_next){}else{
        
      PubBias<-regtest(TempMeta) #The regression test by Egger et al. (1997) 
      
      metaOutput[i, 15]<-PubBias$zval #the value of the Egger test statistic
      metaOutput[i, 16]<-PubBias$pval #the corresponding Egger test p-value 
      metaOutput[i, 17]<-PubBias$dfs #the corresponding Egger test degrees of freedom
      rm(PubBias)
      }
      
      #Testing for robustness (sensitivity of the meta-analysis to the results of individual contrasts)
      skip_to_next <- FALSE
      tryCatch(Robustness<-leave1out(TempMeta), error = function(e) {skip_to_next <<- TRUE})
      if(skip_to_next){}else{
        
      Robustness<-leave1out(TempMeta) #In a leave-one-out analysis, the same model is repeatedly fitted, leaving out one study at a time. By doing so, we can assess how much the results are influenced by each individual study
      
      metaOutput[i, 18]<-min(Robustness$estimate, na.rm=TRUE) #minimum estimated effect in leave-one-out analysis
      metaOutput[i, 19]<-max(Robustness$estimate, na.rm=TRUE) #maximum estimated effect in leave-one-out analysis
      metaOutput[i, 20]<-min(Robustness$pval, na.rm=TRUE) #minimum pval in leave-one-out analysis
      metaOutput[i, 21]<-max(Robustness$pval, na.rm=TRUE) #maximum pval in leave-one-out analysis 
      rm(Robustness)
      }
      
      #Testing for influential datapoints:
      skip_to_next <- FALSE
      tryCatch(InfluenceStats<-influence(TempMeta), error = function(e) {skip_to_next <<- TRUE})
      if(skip_to_next){}else{
        
      InfluenceStats<-influence(TempMeta) #Pulling out some stats for influential datapoints (contrasts)
      
      influence_dfbs[i,InfluenceStats$ids]<-InfluenceStats$dfbs$intrcpt
  
      influence_cookd[i,InfluenceStats$ids]<-InfluenceStats$inf$cook.d
      
      influence_TF[i,InfluenceStats$ids]<-InfluenceStats$is.infl
      
      rm(InfluenceStats)
      }
      
      rm(TempMeta)
    }
    rm(effect, var)
  }
  
  #Naming the columns in our output:
  colnames(metaOutput)<-c("Log2FC_estimate", "SE", "pval", "CI_lb", "CI_ub", "Number_Of_Comparisons", "Number_of_ Contrasts", "Number_of_Coefficients", "tau2_ResidualHeterogeneity", "SE_tau2_ResidualHeterogeneity", "QE_CochransQ_Teststat", "QEp_CochransQ_pval", "I2_PercentVar_TrueHeterogeneity", "H2_Ratio_EffectHetero_overSamplVar", "PubBias_Egger_Zstat", "PubBias_Egger_pval", "PubBias_Egger_DF", "Leave1Out_Min_Log2FC","Leave1Out_Max_Log2FC", "Leave1Out_Min_Pval","Leave1Out_Max_Pval")
  
  #The row names for our output are the combined mouse-rat entrez ids: 
  row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta[,Column_GeneName]
  
  #We return this output back into our global environment
  metaOutput<<-metaOutput
  MetaAnalysis_Annotation<<-MetaAnalysis_FoldChanges_ForMeta[,-Columns_DE]
  influence_dfbs<<-influence_dfbs
  influence_cookd<<-influence_cookd
  influence_TF<<-influence_TF
  
  return(metaOutput)
  return(MetaAnalysis_Annotation)
  return(influence_dfbs)
  return(influence_cookd)
  return(influence_TF)
  
  #... and provide the user with an update about the newly created object:
  
  print("metaOutput:")
  print(str(metaOutput))
  
  print("Top of metaOutput:")
  print(head(metaOutput))
  
  print("Bottom of metaOutput")
  print(tail(metaOutput))
  
}
