# This code document includes a function that corrects the meta-analysis output to take into account the fact that we are running the statistical calculations thousands of times and therefore have a heightened risk of false discovery (false discovery rate correction) 
#Megan Hagenauer
#July 25, 2024

####################

#Installing and loading relevant code packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("multtest")                                 
library(multtest)

####################

#Function:

FalseDiscoveryCorrection<-function(metaOutput, HOM_MouseVsRat, MetaAnalysis_Annotation){
  
  #This calculates the false discovery rate, or q-value, for each of our p-values using the Benjamini-Hochberg procedure:
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,3], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[7]<-"FDR"
  
  #These results are returned to our global environment:
  metaOutputFDR<<-metaOutputFDR
  
  #We let the user know the basic structure of the meta-analysis output with FDR added to it (just to make sure everything still looks good)
  print("metaOutputFDR:")
  print(str(metaOutputFDR))
  
  #Then we make a dataframe that adds the annotation to that output:
  TempDF<-cbind.data.frame(metaOutputFDR, MetaAnalysis_Annotation)
  #And then adds even more detailed gene annotation:
  
  #First the detailed annotation for the mouse genes:
  TempDF2<-join(TempDF, HOM_MouseVsRat[,c(4:5,9:11)], by="Mouse_EntrezGene.ID", type="left", match="first")
  
  #Next the annnotation for the rat genes:
  TempDF3<-join(TempDF2, HOM_MouseVsRat[,c(15:16,20:22)], by="Rat_EntrezGene.ID", type="left", match="first")
  
  #This is renamed and returned to our global environment:
  metaOutputFDR_annotated<-TempDF3
  metaOutputFDR_annotated<<-metaOutputFDR_annotated
  
  #And written out into our working directory:
  write.csv(metaOutputFDR_annotated, "metaOutputFDR_annotated.csv")
  
  #Then we make a version of the output in order by p-value:
  metaOutputFDR_OrderbyPval<<-metaOutputFDR_annotated[order(metaOutputFDR_annotated[,5]),]
  
  #Let's write out a version of the output in order by p-value:
  write.csv(metaOutputFDR_OrderbyPval, "metaOutputFDR_orderedByPval.csv")
  
  #And give the user some information about their results:
  
  print("Do we have any genes that are statistically significant following loose false discovery rate correction (FDR<0.10)?")
  print(sum(metaOutputFDR_annotated[,9]<0.10, na.rm=TRUE))
  
  print("Do we have any genes that are statistically significant following traditional false discovery rate correction (FDR<0.05)?")
  print(sum(metaOutputFDR_annotated[,9]<0.05, na.rm=TRUE))
  
  print("What are the top results?")
  print(head(metaOutputFDR_annotated[order(metaOutputFDR_annotated[,5]),]))
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
}

#####################

#Example usage:

#FalseDiscoveryCorrection(metaOutput, HOM_MouseVsRat, #MetaAnalysis_Annotation)

#####################

#Example output:

# [1] "metaOutputFDR:"
# num [1:16555, 1:7] 0.0234 0.1315 -0.023 0.0437 0.0769 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:16555] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# NULL
# [1] "Do we have any genes that are statistically significant following loose false discovery rate correction (FDR<0.10)?"
# [1] 253
# [1] "Do we have any genes that are statistically significant following traditional false discovery rate correction (FDR<0.05)?"
# [1] 136
# [1] "What are the top results?"
# Rat_EntrezGene.ID Mouse_EntrezGene.ID Log2FC_estimate         SE         pval      CI_lb       CI_ub
# 13170_24309               24309               13170      -0.4293013 0.07014204 9.330764e-10 -0.5667771 -0.29182542
# 74772_362645             362645               74772      -0.1313814 0.02330647 1.729185e-08 -0.1770612 -0.08570151
# 192188_282580            282580              192188       0.3562365 0.06392332 2.505735e-08  0.2309491  0.48152389
# 54006_83632               83632               54006      -0.1635032 0.02996583 4.861038e-08 -0.2222351 -0.10477121
# 216565_305556            305556              216565       0.1764604 0.03363489 1.551430e-07  0.1105372  0.24238352
# 268445_360575            360575              268445      -0.1584491 0.03081248 2.712959e-07 -0.2188405 -0.09805779
# Number_Of_Comparisons          FDR MouseVsRat_EntrezGene.ID Mouse_Symbol Mouse_Genetic.Location
# 13170_24309                       5 1.544708e-05              13170_24309          Dbp               Chr7  cM
# 74772_362645                      5 1.382748e-04             74772_362645      Atp13a2               Chr4  cM
# 192188_282580                     5 1.382748e-04            192188_282580        Stab2              Chr10  cM
# 54006_83632                       5 2.011862e-04              54006_83632        Deaf1               Chr7  cM
# 216565_305556                     5 5.136784e-04            216565_305556        Ehbp1              Chr11  cM
# 268445_360575                     5 7.485505e-04            268445_360575     Ankrd13b              Chr11  cM
# Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38.                              Mouse_Name Rat_Symbol
# 13170_24309                                Chr7:45354658-45359579(+) D site albumin promoter binding protein        Dbp
# 74772_362645                             Chr4:140714184-140734641(+)                        ATPase type 13A2    Atp13a2
# 192188_282580                             Chr10:86677062-86843889(-)                              stabilin 2      Stab2
# 54006_83632                              Chr7:140877093-140907603(-)             DEAF1, transcription factor      Deaf1
# 216565_305556                             Chr11:21955825-22237086(-)             EH domain binding protein 1      Ehbp1
# 268445_360575                             Chr11:77361311-77380504(-)               ankyrin repeat domain 13b   Ankrd13b
# Rat_Genetic.Location Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38.                                    Rat_Name
# 13170_24309               Chr1 q22                                                   NA D-box binding PAR bZIP transcription factor
# 74772_362645              Chr5 q36                                                   NA             ATPase cation transporting 13A2
# 192188_282580                 Chr7                                                   NA                                  stabilin 2
# 54006_83632               Chr1 q41                                                   NA                  DEAF1 transcription factor
# 216565_305556            Chr14 q22                                                   NA                 EH domain binding protein 1
# 268445_360575            Chr10 q24                                                   NA                   ankyrin repeat domain 13B



