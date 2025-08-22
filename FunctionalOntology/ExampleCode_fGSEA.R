#Example usage of Brain.GMT:
#Example Code written for R 4.3.3. using fgsea v.1.2.1 
#Toni Duan and Megan Hagenauer
#July 8th, 2024
#Updated Sept 27 2024 when we realized that we had missed a dataset in our meta-analysis

#Link to the fast gene set enrichment analysis (fGSEA) documentation:
# https://bioconductor.org/packages/release/bioc/html/fgsea.html

#installing the R package fast Gene Set Enrichment Analysis (fGSEA):

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea", force = TRUE)

library(fgsea)

#This analysis assumes a differential expression (DE) output file structure similar to that produced by the Limma or EdgeR pipelines 
#Rows=all genes included in the DE analysis, columns=gene annotation and DE statistical output
#At least one of the annotation columns must be official gene symbol
#At least one of the columns of differential statistics must include DE effect size (e.g., Log2 Fold Change)

#Read in the full DE results for a condition from the working directory 
#Replace "DEResults.csv" in the code with your file name
DEResults<-read.csv("metaOutputFDR_orderedByPval.csv", header=TRUE, stringsAsFactors = FALSE)

#Remove rows of DE results that are missing gene symbol annotation or effect size information
#Replace $gene_symbol in the code with the column name containing gene symbols in your DE output
#Replace $Log2FC in the code with the column name containing effect sizes in your DE output
DEResults_noNA<-DEResults[is.na(DEResults$Mouse_Symbol)==FALSE & is.na(DEResults$Log2FC_estimate)==FALSE,]

#The analysis only works if there is one effect size (e.g., log2 fold change or Log2FC) per gene symbol.
#One way to deal with multiple effect sizes mapping to the same gene (e.g., multiple transcripts or probes) is to average them:
#Replace $Log2FC in the code with the column name containing effect sizes in your DE output
#Replace $gene_symbol in the code with the column name containing gene symbols in your DE output
DEResults_Log2FC_forGSEA<-tapply(X=DEResults_noNA$Log2FC_estimate, INDEX=DEResults_noNA$Mouse_Symbol, FUN=mean)
names(DEResults_Log2FC_forGSEA)<-names(table(DEResults_noNA$Mouse_Symbol))

#The effect sizes should be ordered from smallest to largest:
#Replace $Log2FC in the code with the column name containing effect sizes in your DE output
DEResults_Log2FC_forGSEA_Ranked<-DEResults_Log2FC_forGSEA[order(DEResults_Log2FC_forGSEA)]

str(DEResults_Log2FC_forGSEA_Ranked)
# num [1:11885(1d)] -1.388 -0.801 -0.781 -0.726 -0.678 ...
# - attr(*, "dimnames")=List of 1
# ..$ : chr [1:11885] "Fam50a" "Steap1" "Cartpt" "Tfap2d" ...

#Read in Brain.GMT for your species of interest (this example uses rat)
#If you get a warning about an incomplete line in the .gmt file, just ignore it
BrainGMT<-gmtPathways("BrainGMTv2_wGO_MouseOrthologs_TrimmedToCortex.gmt.txt")

#Run fast fGSEA on your ranked, averaged effect sizes:
#This code should be compatible with updated fgsea packages - if you have an updated package, this code will run as fgseaSimple()
GSEA_Results<-fgsea(BrainGMT, DEResults_Log2FC_forGSEA_Ranked, nperm=10000, minSize = 10, maxSize = 1000)

#Pull out the names for the genes that are driving the enrichment of differential expression in each gene set:
GSEA_Results$leadingEdge<-vapply(GSEA_Results$leadingEdge, paste, collapse= ",", character(1L))

#Write out the results:
write.csv(GSEA_Results, "GSEA_Results.csv")

#You can easily view these results in Excel
# Sort by p-value
# padj: false discovery rate (FDR) corrected p-value. This value is normally used to set the threshold for significance (FDR<0.05) 
# ES & NES: Enrichment Score and Normalized Enrichment Score for each gene set. 
# Positive ES & NES values mean that the gene set is enriched with upregulation in response to your variable of interest
# Negative ES & NES values mean that the gene set is enriched with downregulation in response to your variable of interest

# Other aspects of the output can be deciphered by referencing the original GSEA publication: Subramanian et al. 2005
# https://www.pnas.org/doi/10.1073/pnas.0506580102
