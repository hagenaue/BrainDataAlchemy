#Example pipeline for aligning our results across datasets:
#Megan Hagenauer
#July 27 2026

############

#Goals:
#Each dataset has differential expression results from a slightly different list of genes
#Depending on the exact tissue dissected, the sensitivity of the transcriptional profiling platform, the representation on the transcriptional profiling platform (for microarray), and the experimental conditions
#The differential expression results from different datasets will also be in a slightly different order
#We want to align these results so that the differential expression results from each dataset are columns, with each row representing a different gene

############

if (!requireNamespace("plyr", quietly = TRUE)) {
  install.packages("plyr")
}

############

#Download the necessary functions from our Github repository into your working directory:
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDEResults_2024/Function_AligningDEResults.R

############

#Reading in the functions:

source("Function_AligningDEResults.R")

###########

#Aligning the mouse datasets with each other:

#Example Usage:

ListOfMouseDEResults<-list(DEResults_GSE111212_515554, DEResults_GSE123582_519225)

AligningMouseDatasets(ListOfMouseDEResults)

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

#################

#Code for aligning the rat and mice results:

#This code isn't nicely functionalized yet
#It also assumes that there are mouse datasets
#It will break if there are only rat datasets - this needs to be fixed.

################

#First: What are gene orthologs?

# Homology refers to biological features including genes and their products that are descended from a feature present in a common ancestor.

# Homologous genes become separated in evolution in two different ways: separation of two populations with the ancestral gene into two species or gene duplication of the ancestral gene within a lineage:

### Genes separated by speciation are called orthologs.
### Genes separated by gene duplication events are called paralogs.

#This definition came from NCBI (https://www.nlm.nih.gov/ncbi/workshops/2023-08_BLAST_evol/ortho_para.html)


#We have the ortholog database that we downloaded from Jackson Labs on April 25, 2024
#This database was trimmed and formatted using the code "FormattingRatMouseOrthologDatabase_20240425.R"

MouseVsRat_NCBI_Entrez<-read.csv("HOM_MouseVsRat_EntrezEnsemblAgree_NoMultimapped_20260511.csv", header=TRUE, stringsAsFactors = FALSE, row.names=1, colClasses=c("character", "character", "character"))

#We want to join this ortholog database to our mouse results (Log2FC and SV):

Mouse_MetaAnalysis_FoldChanges_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_FoldChanges, by="Mouse_EntrezGene.ID", type="full")

str(Mouse_MetaAnalysis_FoldChanges_wOrthologs)
#'data.frame':	28920 obs. of  30 variables:

Mouse_MetaAnalysis_SV_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_SV, by="Mouse_EntrezGene.ID", type="full")

str(Mouse_MetaAnalysis_SV_wOrthologs)
#'data.frame':	28920 obs. of  30 variables:


#*If there are rat datasets*, we then want to join our mouse Log2FC and SV results to the rat results using the ortholog information:
MetaAnalysis_FoldChanges<-join(Mouse_MetaAnalysis_FoldChanges_wOrthologs, Rat_MetaAnalysis_FoldChanges, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_FoldChanges)
#'data.frame':	36792 obs. of  32 variables:

MetaAnalysis_SV<-join(Mouse_MetaAnalysis_SV_wOrthologs, Rat_MetaAnalysis_SV, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_SV)
#'data.frame':	36792 obs. of  32 variables:


#*If there aren't any rat datasets*, we just rename the dataframes so that our downstream code works:
MetaAnalysis_FoldChanges<-Mouse_MetaAnalysis_FoldChanges_wOrthologs
str(MetaAnalysis_FoldChanges)

MetaAnalysis_SV<-Mouse_MetaAnalysis_SV_wOrthologs
str(MetaAnalysis_SV)

###############################

#Not all of these have annotation...
#Those would be the genes that have Entrez IDs that aren't 1:1 with Ensembl
#This matters more this year because we will be joining across datasets with different annotation

sum(is.na(MetaAnalysis_FoldChanges$Mouse_ENSEMBLGene.ID) & is.na(MetaAnalysis_FoldChanges$Rat_Ensembl))
#[1] 14097

MetaAnalysis_FoldChanges<-MetaAnalysis_FoldChanges[(is.na(MetaAnalysis_FoldChanges$Mouse_ENSEMBLGene.ID) & is.na(MetaAnalysis_FoldChanges$Rat_Ensembl))==FALSE,]

dim(MetaAnalysis_FoldChanges)
#[1] 22695    32

MetaAnalysis_SV<-MetaAnalysis_SV[(is.na(MetaAnalysis_SV$Mouse_ENSEMBLGene.ID) & is.na(MetaAnalysis_SV$Rat_Ensembl))==FALSE,]

dim(MetaAnalysis_SV)
#[1] 22695    32

#For simplicity's sake for labeling later charts, let's make a combo annotation with both mouse and rat gene symbol:
MetaAnalysis_FoldChanges$MouseRat_GeneSymbol<-paste(MetaAnalysis_FoldChanges$Mouse_Symbol, MetaAnalysis_FoldChanges$Rat_Symbol, sep="_")

MetaAnalysis_SV$MouseVsRat_GeneSymbol<-paste(MetaAnalysis_SV$Mouse_Symbol, MetaAnalysis_SV$Rat_Symbol, sep="_")


#######################

#We should probably spend some time renaming the comparisons included in our MetaAnalysis_FoldChanges and MetaAnalysis_SV objects now...

colnames(MetaAnalysis_FoldChanges)
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
# [13] "Mouse_ENSEMBLGene.ID"                                  
# [14] "Ensembl_Entrez"                                        
# [15] "Rat_Common.Organism.Name"                              
# [16] "Rat_NCBI.Taxon.ID"                                     
# [17] "Rat_Symbol"                                            
# [18] "Rat_EntrezGene.ID"                                     
# [19] "Rat_Mouse.MGI.ID"                                      
# [20] "Rat_HGNC.ID"                                           
# [21] "Rat_OMIM.Gene.ID"                                      
# [22] "Rat_Genetic.Location"                                  
# [23] "Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38."  
# [24] "Rat_Name"                                              
# [25] "Rat_Synonyms"                                          
# [26] "Rat_Ensembl"                                           
# [27] "Rat_DBSymbol"                                          
# [28] "Rat_DBName"                                            
# [29] "GSE111212_exercise"                                    
# [30] "GSE123582_F1..sedentary...F1..sedentary."              
# [31] "GSE270831_social.isolation"                            
# [32] "GSE299436_exercise"                                    
# [33] "MouseRat_GeneSymbol" 

###################

#Comparing Log2FC across datasets

#Simple scatterplot... not so promising:

#Example scatter plot comparing two datasets:
plot(MetaAnalysis_FoldChanges$GSE299436_exercise~MetaAnalysis_FoldChanges$GSE111212_exercise)

#Note - many people prefer to plot these relationships using RRHOs (Rank rank hypergeometric overlap plots)
#I like using both.
#The code for the RRHOs is a little complicated, but I'm happy to share if folks are interested.

#Here's code for looking at the correlation of all of our log2FC results with all of our other log2FC results
#This is called a correlation matrix:

cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:28,33)]), use="pairwise.complete.obs", method="spearman")
#There isn't much similarity across conditions here (outside of comparisons within the same experiment)

#An illustration of the correlation matrix using a hierarchically clustered heatmap, although somewhat pathetic:
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:28,33)]), use="pairwise.complete.obs", method="spearman"))
