#Example code for making a Heatmap illustrating the Log2FCs from the individual studies for each of the top meta-analysis DEGs
#Megan Hagenauer & Jinglin Xiong
#June 6 2025

#This code is not functionalized and beautiful, but it is well annotated...

################

#Set the working directory to where you are storing the R workspaces for your project
#This can be done using coding or using the drop-down menu in R studio ("Session" ->"Set Working Directory").
#If you use the drop down menu, I recommend grabbing the code that appears in the console and saving it in your code so that you have documentation as to where your input file came from.

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_SophieMensch_ADHD/R_Code_And_Workspaces/Workspaces")

#This will give you the file names in that folder:
list.files()

#############

#Read in the final workspace for the meta-analysis
#This can be done with code or using the drop-down menu ("Session" ->"Load Workspace").
#If you use the drop down menu, I recommend grabbing the code that appears in the console and saving it in your code so that you have documentation as to where your input file came from.

load("Workspace_20240730.RData")

#############

#Loading and installing the relevant code libraries:

install.packages("pheatmap")
library(pheatmap)

install.packages("dichromat")
library(dichromat)

install.packages("plyr")
library(plyr)

#############

#There should be two objects in the workspace that are of interest for this task. We will need to locate each of these objects and determine which column contains the gene identifier information.

#1. The meta-analysis output ordered by p-value, so that the most significant results are at the top of the data-frame:

#In this case, the object is named metaOutputFDR_OrderbyPval

#This shows a glimpse of the object:
head(metaOutputFDR_OrderbyPval)

#This provides the column names:
colnames(metaOutputFDR_OrderbyPval)
# [1] "Rat_EntrezGene.ID"                                     
# [2] "Mouse_EntrezGene.ID"                                   
# [3] "Log2FC_estimate"                                       
# [4] "SE"                                                    
# [5] "pval"                                                  
# [6] "CI_lb"                                                 
# [7] "CI_ub"                                                 
# [8] "Number_Of_Comparisons"                                 
# [9] "FDR"                                                   
# [10] "MouseVsRat_EntrezGene.ID"                              
# [11] "Mouse_Symbol"                                          
# [12] "Mouse_Genetic.Location"                                
# [13] "Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38."
# [14] "Mouse_Name"                                            
# [15] "Rat_Symbol"                                            
# [16] "Rat_Genetic.Location"                                  
# [17] "Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38."  
# [18] "Rat_Name"          

#2. The log2 fold changes from the individual studies that were used to run the meta-analysis

#In this case, the object is named MetaAnalysis_FoldChanges_ForMeta

head(MetaAnalysis_FoldChanges_ForMeta)

colnames(MetaAnalysis_FoldChanges_ForMeta)
# [1] "Rat_EntrezGene.ID"                  "Mouse_EntrezGene.ID"               
# [3] "MouseVsRat_EntrezGene.ID"           "GSE87700_ethanol"                  
# [5] "GSE66468_Dup16p11.2"                "GSE66468_DupPlusDel16p11.2"        
# [7] "GSE97181_Foxp1"                     "GSE171191_Scn1aHeterozygous"       
# [9] "GSE171191_Scn1aHeterozyPlusOverExp" "GSE62385_hypoxia"                  
# [11] "GSE60000_ethanol"                   "GSE29417_MSN"                    
# [13] "GSE117357_Adgrl3"                   "GSE240887_Kdm5b"                 
# [15] "GSE226730_SHR"                      "GSE12457_PCBs"                   
# [17] "GSE140598_bHRvbLR"                  "GSE140598_bHRvbIR"               
# [19] "GSE140596_bHR"                      "GSE140287_bHR"  

##################

#In this case, it looks like the column of gene identifiers that we will want to use is MouseVsRat_EntrezGene.ID - that gene identifier is located in both data frames 

#Let's join the two data frames using that column
metaOutputFDR_wFoldChanges<-join(metaOutputFDR_OrderbyPval, MetaAnalysis_FoldChanges_ForMeta, by="MouseVsRat_EntrezGene.ID", type="inner")

#Let's grab the top 50 genes:

metaOutputFDR_wFoldChanges_Top50<-metaOutputFDR_wFoldChanges[c(1:50),]

head(metaOutputFDR_wFoldChanges_Top50)

#And then re-order the results by Log2FC estimate from the meta-analysis

metaOutputFDR_wFoldChanges_Top50_byLog2FC<-metaOutputFDR_wFoldChanges_Top50[order(metaOutputFDR_wFoldChanges_Top50$Log2FC_estimate),]

#We will eventually need a name for each gene that can be used in the heatmap.
#Gene symbols are most recognizable, but some of these genes don't have gene symbols, just entrez ids...
#So we could make some sort of combination label like this:

metaOutputFDR_wFoldChanges_Top50_byLog2FC$geneSymbols_top_50 <- paste("Mm", metaOutputFDR_wFoldChanges_Top50_byLog2FC$Mouse_EntrezGene.ID, metaOutputFDR_wFoldChanges_Top50_byLog2FC$Mouse_Symbol,"; Rn", metaOutputFDR_wFoldChanges_Top50_byLog2FC$Rat_EntrezGene.ID, metaOutputFDR_wFoldChanges_Top50_byLog2FC$Rat_Symbol, sep=" ")

#Depending on the meta-analysis cohort (2022, 2023, 2024), the gene identifier column may just be the row.names, named "X", or "GeneSymbol", in which case, our code might look like this instead:
#metaOutputFDR_wFoldChanges_Top50_byLog2FC$gene_top_50 <- rownames(metaOutputFDR_wFoldChanges_Top50_byLog2FC)[1:50]
#or:
#metaOutputFDR_wFoldChanges_Top50_byLog2FC$gene_top_50 <- metaOutputFDR_wFoldChanges_Top50_byLog2FC$X[1:50]
#or:
#metaOutputFDR_wFoldChanges_Top50_byLog2FC$gene_top_50 <- metaOutputFDR_wFoldChanges_Top50_byLog2FC$GeneSymbol[1:50]

#############################################

#Grabbing the rows of fold change information for those genes

colnames(metaOutputFDR_wFoldChanges_Top50_byLog2FC)
# [1] "Rat_EntrezGene.ID"                                     
# [2] "Mouse_EntrezGene.ID"                                   
# [3] "Log2FC_estimate"                                       
# [4] "SE"                                                    
# [5] "pval"                                                  
# [6] "CI_lb"                                                 
# [7] "CI_ub"                                                 
# [8] "Number_Of_Comparisons"                                 
# [9] "FDR"                                                   
# [10] "MouseVsRat_EntrezGene.ID"                              
# [11] "Mouse_Symbol"                                          
# [12] "Mouse_Genetic.Location"                                
# [13] "Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38."
# [14] "Mouse_Name"                                            
# [15] "Rat_Symbol"                                            
# [16] "Rat_Genetic.Location"                                  
# [17] "Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38."  
# [18] "Rat_Name"                                              
# [19] "Rat_EntrezGene.ID"                                     
# [20] "Mouse_EntrezGene.ID"                                   
# [21] "GSE87700_ethanol"                                      
# [22] "GSE66468_Dup16p11.2"                                   
# [23] "GSE66468_DupPlusDel16p11.2"                            
# [24] "GSE97181_Foxp1"                                        
# [25] "GSE171191_Scn1aHeterozygous"                           
# [26] "GSE171191_Scn1aHeterozyPlusOverExp"                    
# [27] "GSE62385_hypoxia"                                      
# [28] "GSE60000_ethanol"                                      
# [29] "GSE29417_MSN"                                          
# [30] "GSE117357_Adgrl3"                                      
# [31] "GSE240887_Kdm5b"                                       
# [32] "GSE226730_SHR"                                         
# [33] "GSE12457_PCBs"                                         
# [34] "GSE140598_bHRvbLR"                                     
# [35] "GSE140598_bHRvbIR"                                     
# [36] "GSE140596_bHR"                                         
# [37] "GSE140287_bHR"                                         
# [38] "geneSymbols_top_50"  

#At this point, we need to get rid of the columns containing the gene identifier information to make a numerical matrix of Log2FC information from each study:
#The columns that are removed will depend on which columns in the object contain gene information (string) instead of Log2FC:
Log2FC_Subsetted_Matrix <- as.matrix(metaOutputFDR_wFoldChanges_Top50_byLog2FC[,c(21:37)])

#For making a visualization, most readers will recognize the gene symbols but not the entrez ids, so lets make that the gene identifiers:

row.names(Log2FC_Subsetted_Matrix) <-metaOutputFDR_wFoldChanges_Top50_byLog2FC$geneSymbols_top_50

str(Log2FC_Subsetted_Matrix)
# num [1:50, 1:17] NA NA NA NA -0.123 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:50] "Mm 66737 NA ; Rn NA NA" "Mm NA NA ; Rn 498276 Fcgr2al1" "Mm 431706 Zfp457 ; Rn NA NA" "Mm NA NA ; Rn 102552664 NA" ...
# ..$ : chr [1:17] "GSE87700_ethanol" "GSE66468_Dup16p11.2" "GSE66468_DupPlusDel16p11.2" "GSE97181_Foxp1" ...


########################

#Set the working directory to the folder where you are storing R output:

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_SophieMensch_ADHD/ROutput_And_Results")

########################

#Making a hierarchically clustered heatmap:

#Unlike Jinglin's project, I discovered issues with clustering by row because it is impossible to cluster the genes that are only found in mice or only found in rats because they aren't found in the same studies.

#... so I turned off clustering by row and just ordered the rows by the estimated Log2FC in the meta-analysis

pdf("Heatmap_TopMetaGenes_50_color_scale_narrow.pdf",
    height = 11, width = 8.5)

pheatmap(Log2FC_Subsetted_Matrix,
         color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         scale = "none",
         breaks = seq(-1.5, 1.5, length.out = 101),
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         fontsize_col = 8,  
         width = 8.5,
         height = 11,
         border_color = NA)

dev.off()
