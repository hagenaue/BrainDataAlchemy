#Example pipeline for aligning our results across datasets:
#Megan Hagenauer
#July 25 2024

############

#Goals:
#Each dataset has differential expression results from a slightly different list of genes
#Depending on the exact tissue dissected, the sensitivity of the transcriptional profiling platform, the representation on the transcriptional profiling platform (for microarray), and the experimental conditions
#The differential expression results from different datasets will also be in a slightly different order
#We want to align these results so that the differential expression results from each dataset are columns, with each row representing a different gene

############

#Reading in the functions:

source("Function_AligningDEResults.R")

###########

#Aligning the mouse datasets with each other:

#Example Usage:

ListOfMouseDEResults<-list(DEResults_GSE126678, DEResults_GSE181285)

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

#Aligning the rat datasets with each other:

#Example Usage:
#We only have one rat dataset, but we still run it through this code to get it in the right format for our next steps:

ListOfRatDEResults<-list(DEResults_GSE205325)

AligningRatDatasets(ListOfRatDEResults)
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

################

############

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

MouseVsRat_NCBI_Entrez<-read.csv("MouseVsRat_NCBI_Entrez_JacksonLab_20240425.csv", header=TRUE, stringsAsFactors = FALSE, row.names=1, colClasses=c("character", "character", "character"))

#We want to join this ortholog database to our mouse results (Log2FC and SV):

Mouse_MetaAnalysis_FoldChanges_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_FoldChanges, by="Mouse_EntrezGene.ID", type="full")

str(Mouse_MetaAnalysis_FoldChanges_wOrthologs)
#'data.frame':	25288 obs. of  6 variables:

Mouse_MetaAnalysis_SV_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_SV, by="Mouse_EntrezGene.ID", type="full")

str(Mouse_MetaAnalysis_SV_wOrthologs)
#'data.frame':	25288 obs. of  6 variables:


#*If there are rat datasets*, we then want to join our mouse Log2FC and SV results to the rat results using the ortholog information:
MetaAnalysis_FoldChanges<-join(Mouse_MetaAnalysis_FoldChanges_wOrthologs, Rat_MetaAnalysis_FoldChanges, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_FoldChanges)
#'data.frame':	28101 obs. of  7 variables:

MetaAnalysis_SV<-join(Mouse_MetaAnalysis_SV_wOrthologs, Rat_MetaAnalysis_SV, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_SV)
#'data.frame':	28101 obs. of  7 variables:

#*If there aren't any rat datasets*, we just rename the dataframes so that our downstream code works:
MetaAnalysis_FoldChanges<-Mouse_MetaAnalysis_FoldChanges_wOrthologs
str(MetaAnalysis_FoldChanges)

MetaAnalysis_SV<-Mouse_MetaAnalysis_SV_wOrthologs
str(MetaAnalysis_SV)

#For simplicity's sake, I'm going to replace that Mouse-Rat Entrez annotation
#Because it is missing entries for any genes in the datasets that *don't* have orthologs
MetaAnalysis_FoldChanges$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_FoldChanges$Mouse_EntrezGene.ID, MetaAnalysis_FoldChanges$Rat_EntrezGene.ID, sep="_")

MetaAnalysis_SV$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_SV$Mouse_EntrezGene.ID, MetaAnalysis_SV$Rat_EntrezGene.ID, sep="_")


#Comparing Log2FC across datasets

#Simple scatterplot... not so promising:
colnames(MetaAnalysis_FoldChanges)

# [1] "Rat_EntrezGene.ID"                 "Mouse_EntrezGene.ID"              
# [3] "MouseVsRat_EntrezGene.ID"          "GSE126678_LPS_Acute"              
# [5] "GSE126678_LPS_SubchronicPlusAcute" "GSE126678_LPS_Subchronic"         
# [7] "GSE181285_LPS_Acute"               "GSE205325_LPS_Chronic" 

#Example scatter plot comparing two datasets:
plot(MetaAnalysis_FoldChanges$GSE126678_LPS_Acute~MetaAnalysis_FoldChanges$GSE181285_LPS_Acute)

#Note - many people prefer to plot these relationships using RRHOs (Rank rank hypergeometric overlap plots)
#I like using both.
#The code for the RRHOs is a little complicated, but I'm happy to share if folks are interested.

#Here's code for looking at the correlation of all of our log2FC results with all of our other log2FC results
#This is called a correlation matrix:

cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman")
#There isn't much similarity across conditions here (outside of comparisons within the same experiment)

#An illustration of the correlation matrix using a hierarchically clustered heatmap, although somewhat pathetic:
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman"))

