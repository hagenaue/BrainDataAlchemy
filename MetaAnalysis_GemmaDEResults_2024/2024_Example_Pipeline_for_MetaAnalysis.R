#This code illustrates the analysis pipeline for running a meta-analysis using the effect sizes and sampling variances for each gene stored in the objects MetaAnalysis_FoldChanges and MetaAnalysis_SV
#Megan Hagenauer
#July 25 2024

##############

#Specifying the minimum amount of results necessary to run a meta-anlaysis for a gene:

#For any particular gene, it is likely that some datasets may be missing differential expression results. 

#This is especially true for genes that have low levels of expression and may not be detected by less sensitive assays

#It is also likely to be true for genes that were discovered more recently (i.e., not targeted by older microarray platforms) or that lack a clear ortholog in rats/mice

#We can only run a meta-analysis if there are differential expression results from more than one statistical contrast.

#Since the differential expression results from the same study (dataset) are often artificially correlated (especially if they use the same control group as the comparison), I would actually prefer that there are results from *more than one dataset.* (not just more than one statistical contrast)

#And ideally, we would like to have even more than that...

#Before moving forward, we need to decide what we want to be the minimum number of differential expression results allowable for a gene to be included in our meta-analysis.

#This code caculates the number of NAs (i.e., the number of statistical contrasts lacking real differential expression results) in each row (i.e., for each gene):
MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges[,-c(1:3)], 1, function(y) sum(is.na(y)))

#I'm going to make a histogram of the results because I'm curious to see how they are distributed
hist(MetaAnalysis_FoldChanges_NAsPerRow)

#Or, since there are a limited number of integer answers, I could make a table of the results:
table(MetaAnalysis_FoldChanges_NAsPerRow)
# 0     1     2     3     4     5 
# 13355  3200  5059   277  5293   917

#For this dataset, I'm going to try running a meta-analysis using genes that were found in at least 4 sets of differential expression results
#Since there are 5 sets of differential expression results, that means that the genes that we are including need to have 1 or fewer NAs in their results
#I set this conservatively, because there are so few studies in this meta-analysis.
#2 NA is too many

NumberOfComparisons=5
CutOffForNAs=2
#I have 5 statistical contrasts total (comparisons)
#and 2 NA is too many


########################

#Running a basic meta-analysis:

#This function  is designed to run a basic meta-analysis of Log2FC and sampling variance values using our previously generated objects MetaAnalysis_FoldChanges & MetaAnalysis_SV

source("Function_RunBasicMetaAnalysis.R")                         
#Example usage:

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...


#Example output:

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

str(metaOutput)
# num [1:16555, 1:6] 0.0234 0.1315 -0.023 0.0437 0.0769 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:16555] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:6] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

head(metaOutput)
# Log2FC_estimate         SE       pval        CI_lb      CI_ub Number_Of_Comparisons
# 23825_114087      0.02335239 0.04876404 0.63202016 -0.072223377 0.11892815                     5
# 18585_191569      0.13152931 0.06603932 0.04640599  0.002094618 0.26096399                     5
# 66514_246307     -0.02304066 0.05315084 0.66465470 -0.127214402 0.08113308                     5
# 20480_65041       0.04374446 0.03359984 0.19294219 -0.022110021 0.10959893                     5
# 13726_25437       0.07691445 0.04407249 0.08095348 -0.009466050 0.16329495                     5
# 16952_25380       0.04465958 0.10253423 0.66315765 -0.156303828 0.24562298                     5                   

tail(metaOutput)
# Log2FC_estimate         SE       pval       CI_lb      CI_ub Number_Of_Comparisons
# 108168734_NA      0.11566278 0.07742137 0.13519162 -0.03608031 0.26740587                     4
# 108168883_NA      0.16641238 0.09078313 0.06679127 -0.01151929 0.34434405                     4
# 108168923_NA      0.07757892 0.16987291 0.64789533 -0.25536587 0.41052370                     4
# 108168987_NA     -0.13375235 0.07407534 0.07097679 -0.27893734 0.01143264                     4
# 108169023_NA     -0.06106329 0.04794113 0.20276476 -0.15502618 0.03289960                     4
# 113002583_NA      0.02982367 0.07843726 0.70377969 -0.12391053 0.18355788                     4

########################################

#Accessing some additional gene annotation:

#Reading in a database containing more detailed gene annotation:
HOM_MouseVsRat <- read.csv("HOM_MouseVsRat_20240425.csv", header = TRUE, row.names = 1)

colnames(HOM_MouseVsRat)
# [1] "DB.Class.Key"                                          
# [2] "Mouse_Common.Organism.Name"                            
# [3] "Mouse_NCBI.Taxon.ID"                                   
# [4] "Mouse_Symbol"                                          
# [5] "Mouse_EntrezGene.ID"                                   
# [6] "Mouse_Mouse.MGI.ID"                                    
# [7] "Mouse_HGNC.ID"                                         
# [8] "Mouse_OMIM.Gene.ID"                                    
# [9] "Mouse_Genetic.Location"                                
# [10] "Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38."
# [11] "Mouse_Name"                                            
# [12] "Mouse_Synonyms"                                        
# [13] "Rat_Common.Organism.Name"                              
# [14] "Rat_NCBI.Taxon.ID"                                     
# [15] "Rat_Symbol"                                            
# [16] "Rat_EntrezGene.ID"                                     
# [17] "Rat_Mouse.MGI.ID"                                      
# [18] "Rat_HGNC.ID"                                           
# [19] "Rat_OMIM.Gene.ID"                                      
# [20] "Rat_Genetic.Location"                                  
# [21] "Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38."  
# [22] "Rat_Name"                                              
# [23] "Rat_Synonyms" 

#Renaming the columns so that we can easily join the annotation to our meta-analysis results:
HOM_MouseVsRat$Mouse_EntrezGene.ID <- as.character(HOM_MouseVsRat$Mouse_EntrezGene.ID)

HOM_MouseVsRat$Rat_EntrezGene.ID <- as.character(HOM_MouseVsRat$Rat_EntrezGene.ID)

#################

## Multiple comparison corrections

#This code runs a function that corrects the meta-analysis output to take into account the fact that we are running the statistical calculations thousands of times and therefore have a heightened risk of false discovery (false discovery rate correction) 

source("Function_FalseDiscoveryCorrection.R")

#Example usage:

FalseDiscoveryCorrection(metaOutput, HOM_MouseVsRat, MetaAnalysis_Annotation)

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

############################

#Next we will determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. 

#Here are the top results as listed by mouse gene symbol:
metaOutputFDR_OrderbyPval$Mouse_Symbol[c(1:20)]
# [1] "Dbp"      "Atp13a2"  "Stab2"    "Deaf1"    "Ehbp1"    "Ankrd13b" "Pbdc1"    "Sall1"    "Hnrnpc"   "Bhlhb9"   "Snx14"   
# [12] "Lrrc4b"   "Ccdc85b"  "Rbms3"    "Gas2l1"   "Zc3h15"   "Cdc27"    "Garem2"   "Dnajc7"   "Strn4"     

#Or as listed by mouse entrez id:
metaOutputFDR_OrderbyPval$Mouse_EntrezGene.ID[c(1:20)]
# [1] "13170"  "74772"  "192188" "54006"  "216565" "268445" "67683"  "58198"  "15381"  "70237"  "244962" "272381" "240514" "207181"
# [15] "78926"  "69082"  "217232" "242915" "56354"  "97387" 


#Let's plot some of those top results!

#Quickly looking at the range of Log2FC values to figure out the limits for the x-axis for the forest plots:
hist(metaOutputFDR[,1], breaks=40)
#Range is mostly -1 to 1, but there are a few with Log2FC as big as -2-3

source("Function_MakeForestPlots.R")
#Note - this function currently uses Entrez ID (NCBI ID) as it's input
#It needs to be formatted as a character (not an integer) to work
#It can accept either mouse annotation (Entrez ID) or rat annotation (Entrez ID) as its input
                                            
#Example Usage:

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="13170", species="Mouse") #Dbp

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="74772", species="Mouse")#Atp13a2

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="282580", species="Rat") #Stab2

######################

#This would be a good time to double check again that you have saved your code and your workspace

#It is very common to decide that you want to come back to these results and make more forest plots for top genes


