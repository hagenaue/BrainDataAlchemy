################

#Adapting the old meta-analysis code for our new objects:



#8) Run a meta-analysis using all of the effect sizes for each gene that has data from at least 2 studies. 


#We can only run a meta-analysis if there are differential expression results from more than one comparison.
#Since I know that the differential expression results from the same study (dataset) are artificially correlated, I would actually prefer that there are results from more than one dataset.

#How many genes satisfy this criteria?

#This code caculates the number of NAs in each row:
MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges[,-c(1:3)], 1, function(y) sum(is.na(y)))

#I'm going to make a histogram of the results because I'm curious to see how they are distributed
hist(MetaAnalysis_FoldChanges_NAsPerRow)

#Or, since there are a limited number of integer answers (0-3), I could make a table of the results:
table(MetaAnalysis_FoldChanges_NAsPerRow)
# 0     1     2     3 
# 11043  4823 13346   875

#Let's try running a meta-analysis using genes that were found in at least 2 sets of differential expression results
#Since there are 3 sets of differential expression results, that means that the genes that we are including need to have 1 or fewer NAs in their results
#I set this conservatively, because there are so few studies in this meta-analysis.
#2 NA is too many

library(metafor)

RunBasicMetaAnalysis<-function(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV){
  
  MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges[,-c(1:3)], 1, function(y) sum(is.na(y)))
  
  print("Table of # of NAs per Row (Gene):")
  print(table(MetaAnalysis_FoldChanges_NAsPerRow))
  
  MetaAnalysis_FoldChanges_ForMeta<<-MetaAnalysis_FoldChanges[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  MetaAnalysis_SV_ForMeta<<-MetaAnalysis_SV[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  
  print("MetaAnalysis_FoldChanges_ForMeta:")
  print(str(MetaAnalysis_FoldChanges_ForMeta))
  
  #I'm going to make an empty matrix to store the results of my meta-analysis:
  metaOutput<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), 6)
  
  #And then run a loop that run's a meta-analysis on the differential expression results (columns 2-10) for each gene (row):
  for(i in c(1:nrow(MetaAnalysis_FoldChanges_ForMeta))){
    
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-c(1:3)])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-c(1:3)])
    
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
  row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta[,3]
  
  metaOutput<<-metaOutput
  MetaAnalysis_Annotation<<-MetaAnalysis_FoldChanges_ForMeta[,c(1:3)]
  return(metaOutput)
  return(MetaAnalysis_Annotation)
  
  print("metaOutput:")
  print(str(metaOutput))
  
  print("Top of metaOutput:")
  print(head(metaOutput))
  
  print("Bottom of metaOutput")
  print(tail(metaOutput))
  
}

#Example Usage:
NumberOfComparisons=3
CutOffForNAs=2
#I have 3 comparisons
#2 NA is too many

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

# [1] "Table of # of NAs per Row (Gene):"
# MetaAnalysis_FoldChanges_NAsPerRow
# 0     1     2     3 
# 11043  4823 13346   875 
# [1] "MetaAnalysis_FoldChanges_ForMeta:"
# 'data.frame':	15866 obs. of  6 variables:
# $ Rat_EntrezGene.ID            : chr  "498097" "114521" "360576" "24628" ...
# $ Mouse_EntrezGene.ID          : chr  "68980" "21665" "237858" "18591" ...
# $ MouseVsRat_EntrezGene.ID     : chr  "68980_498097" "21665_114521" "237858_360576" "18591_24628" ...
# $ StressEffects_InteractionWSex: num  0.0831 -0.0504 -0.7557 -0.0782 0.161 ...
# $ CuffOperation_vs_Ctrl        : num  0.01656 0.00964 -0.0427 0.01155 NA ...
# $ ChronicConstriction_vs_Ctrl  : num  -0.0608 0.0631 -0.4516 0.4982 1.4642 ...
# NULL

str(metaOutput)
# num [1:15866, 1:6] 0.012 0.0143 -0.1497 0.157 1.3608 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:15866] "68980_498097" "21665_114521" "237858_360576" "18591_24628" ...
# ..$ : chr [1:6] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

head(metaOutput)
# Log2FC_estimate         SE       pval       CI_lb      CI_ub
# 68980_498097       0.01203600 0.05083559 0.81284048 -0.08759993 0.11167192
# 21665_114521       0.01426912 0.03385763 0.67343034 -0.05209062 0.08062887
# 237858_360576     -0.14968156 0.17385300 0.38925664 -0.49042717 0.19106405
# 18591_24628        0.15695379 0.17760171 0.37683643 -0.19113918 0.50504675
# 229323_688737      1.36081008 0.45349519 0.00269346  0.47197584 2.24964432
# 240549_309227      0.04938634 0.14315287 0.73010174 -0.23118814 0.32996081
# Number_Of_Comparisons
# 68980_498097                      3
# 21665_114521                      3
# 237858_360576                     3
# 18591_24628                       3
# 229323_688737                     2
# 240549_309227                     2                    

tail(metaOutput)
# Log2FC_estimate         SE       pval       CI_lb      CI_ub
# 97423_NA      0.08335277 0.22999011 0.71703901 -0.36741957 0.53412510
# 97775_NA      0.19686310 0.17990823 0.27384896 -0.15575056 0.54947676
# 98303_NA      0.06781550 0.04093082 0.09755343 -0.01240744 0.14803843
# 98736_NA      0.03495936 0.04096462 0.39343552 -0.04532983 0.11524854
# 98870_NA      0.01442392 0.03195141 0.65167762 -0.04819970 0.07704753
# 99100_NA      0.11964317 0.11449473 0.29603839 -0.10476237 0.34404872
# Number_Of_Comparisons
# 97423_NA                     2
# 97775_NA                     2
# 98303_NA                     2
# 98736_NA                     2
# 98870_NA                     2
# 99100_NA                     2

########################################

## Multiple Comparison corrections
#The following code applies two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli) 
#Meta-analysis output with adjusted p-values is then outputted along with effect size information.

#We can add some additional gene annotation at this point too:
HOM_MouseVsRat <- read.csv("HOM_MouseVsRat.csv", header = TRUE, row.names = 1)
HOM_MouseVsRat$Mouse_EntrezGene.ID <- as.character(HOM_MouseVsRat$Mouse_EntrezGene.ID)
HOM_MouseVsRat$Rat_EntrezGene.ID <- as.character(HOM_MouseVsRat$Rat_EntrezGene.ID)

                                            
#9) Correct the meta-analysis output to take into account the fact that we are running the statistical calculations many times and therefore have a heightened risk of false discovery (false discovery rate correction) 

library(multtest)

#Let's functionalize it!
FalseDiscoveryCorrection<-function(metaOutput, HOM_MouseVsRat, MetaAnalysis_Annotation){
  
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,3], proc=c("BH"))
  
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  colnames(metaOutputFDR)[7]<-"FDR"
  
  metaOutputFDR<<-metaOutputFDR
  
  print("metaOutputFDR:")
  print(str(metaOutputFDR))
  
  TempDF<-cbind.data.frame(metaOutputFDR, MetaAnalysis_Annotation)
  
  TempDF2<-join(TempDF, HOM_MouseVsRat[,c(4:5,9:11)], by="Mouse_EntrezGene.ID", type="left", match="first")
  
  TempDF3<-join(TempDF2, HOM_MouseVsRat[,c(15:16,20:22)], by="Rat_EntrezGene.ID", type="left", match="first")
  
  metaOutputFDR_annotated<-TempDF3
  metaOutputFDR_annotated<<-metaOutputFDR_annotated
  
  write.csv(metaOutputFDR_annotated, "metaOutputFDR_annotated.csv")
  
  #a version of the output in order by p-value:
  metaOutputFDR_OrderbyPval<<-metaOutputFDR_annotated[order(metaOutputFDR_annotated[,5]),]
  
  #Let's write out a version of the output in order by p-value:
  write.csv(metaOutputFDR_OrderbyPval, "metaOutputFDR_orderedByPval.csv")
  
  print("Do we have any genes that are statistically significant following false discovery rate correction?")
  print(sum(metaOutputFDR_annotated[,9]<0.10, na.rm=TRUE))
  
  print("What are the top results?")
  print(head(metaOutputFDR_annotated[order(metaOutputFDR_annotated[,5]),]))
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
}

#Example usage:

FalseDiscoveryCorrection(metaOutput, HOM_MouseVsRat)
# [1] "metaOutputFDR:"
# num [1:15866, 1:7] 0.012 0.0143 -0.1497 0.157 1.3608 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:15866] "68980_498097" "21665_114521" "237858_360576" "18591_24628" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# NULL
# [1] "Do we have any genes that are statistically significant following false discovery rate correction?"
# [1] 17
# [1] "What are the top results?"
# Rat_EntrezGene.ID Mouse_EntrezGene.ID Log2FC_estimate         SE
# X215418_363165               363165              215418       0.6358098 0.08742879
# X55985_498335                498335               55985       2.9953890 0.56341937
# X271944_100363111         100363111              271944      -0.8210441 0.16711368
# X14281_314322                314322               14281       0.9824057 0.20007763
# X269593_79428                 79428              269593       0.1943280 0.04299957
# X233332_293004               293004              233332      -0.7318502 0.16338717
# pval      CI_lb      CI_ub Number_Of_Comparisons          FDR
# X215418_363165    3.533738e-13  0.4644526  0.8071671                     3 5.606628e-09
# X55985_498335     1.058133e-07  1.8911073  4.0996707                     2 8.394166e-04
# X271944_100363111 8.965341e-07 -1.1485808 -0.4935073                     2 3.610292e-03
# X14281_314322     9.101958e-07  0.5902607  1.3745506                     3 3.610292e-03
# X269593_79428     6.204449e-06  0.1100503  0.2786056                     2 1.877459e-02
# X233332_293004    7.490953e-06 -1.0520832 -0.4116172                     2 1.877459e-02
# MouseVsRat_EntrezGene.ID Mouse_Symbol Mouse_Genetic.Location
# X215418_363165               215418_363165       Csrnp1               Chr9  cM
# X55985_498335                 55985_498335       Cxcl13               Chr5  cM
# X271944_100363111         271944_100363111       C2cd4d               Chr3  cM
# X14281_314322                 14281_314322          Fos              Chr12  cM
# X269593_79428                 269593_79428        Luzp1               Chr4  cM
# X233332_293004               233332_293004     Adamts17               Chr7  cM
# Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38.p7.
# X215418_363165                                  Chr9:119800229-119813724(-)
# X55985_498335                                     Chr5:96104810-96108927(+)
# X271944_100363111                                 Chr3:94269751-94271874(+)
# X14281_314322                                    Chr12:85520664-85524047(+)
# X269593_79428                                   Chr4:136197072-136282091(+)
# X233332_293004                                    Chr7:66489483-66802919(+)
# Mouse_Name Rat_Symbol
# X215418_363165                       cysteine-serine-rich nuclear protein 1     Csrnp1
# X55985_498335                               C-X-C motif chemokine ligand 13     Cxcl13
# X271944_100363111                 C2 calcium-dependent domain containing 4D     C2cd4d
# X14281_314322                                     FBJ osteosarcoma oncogene        Fos
# X269593_79428                                      leucine zipper protein 1      Luzp1
# X233332_293004    ADAM metallopeptidase with thrombospondin type 1 motif 17   Adamts17
# Rat_Genetic.Location
# X215418_363165                Chr8 q32
# X55985_498335                Chr14 p22
# X271944_100363111                 Chr2
# X14281_314322                 Chr6 q31
# X269593_79428                 Chr5 q36
# X233332_293004                Chr1 q22
# Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38.p7.
# X215418_363165                                                           
# X55985_498335                                                            
# X271944_100363111                                                        
# X14281_314322                                                            
# X269593_79428                                                            
# X233332_293004                                                           
# Rat_Name
# X215418_363165                    cysteine and serine rich nuclear protein 1
# X55985_498335                                C-X-C motif chemokine ligand 13
# X271944_100363111                  C2 calcium-dependent domain containing 4D
# X14281_314322          Fos proto-oncogene, AP-1 transcription factor subunit
# X269593_79428                                       leucine zipper protein 1
# X233332_293004    ADAM metallopeptidase with thrombospondin type 1 motif, 17

#Note- I need to fix the annotation issues here
#It probably would have made sense early on to create a combo mouse-rat entrez id
#Maybe a combined mouse-rat gene symbol too

############################

#10) Determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. 

metaOutputFDR_OrderbyPval$Mouse_Symbol[c(1:17)]
# [1] "Csrnp1"   "Cxcl13"   "C2cd4d"   "Fos"      "Luzp1"    "Adamts17" "Morf4l2" 
# [8] "Wdr36"    "Cdc42ep3" "Tbc1d4"   "Rbmx2"    "Usp34"    "Phyh"     "Sema3c"  
# [15] "Cit"      "Gigyf1"   "Klf2"     

metaOutputFDR_OrderbyPval$Mouse_EntrezGene.ID[c(1:17)]
# [1] "215418" "55985"  "271944" "14281"  "269593" "233332" "56397"  "225348" "260409"
# [10] "210789" "209003" "17847"  "16922"  "20348"  "12704"  "57330"  "16598" 

#Let's plot some of those top results!

#Quickly looking at the range of Log2FC values to figure out the limits for the x-axis for the forest plots:
hist(metaOutputFDR[,1], breaks=40)
#Range is mostly -1 to 1, but there are a few with Log2FC as big as -3-3

MakeForestPlots<-function(EntrezIDAsCharacter){
  
  pdf(paste("ForestPlot_", EntrezIDAsCharacter, ".pdf", sep=""), height=5, width=8)
  
  effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID==EntrezIDAsCharacter,-c(1:3)])
  var<-as.numeric(MetaAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID==EntrezIDAsCharacter,-c(1:3)])
  
  forest.rma(rma(effect, var),slab=colnames(MetaAnalysis_FoldChanges_ForMeta)[-c(1:3)],  xlim=c(-3, 3))
  
  mtext(paste(EntrezIDAsCharacter), line=-1.5, cex=2)
  dev.off()
}


#Example Usage:
MakeForestPlots("215418") 

