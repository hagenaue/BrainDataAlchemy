#Example code: Exploring a Dataset (with more detail)
#Megan Hagenauer
#July 13, 2023

##################################

#This code document is set up to provide an example of analyzing a single Gemma dataset
#It combines many elements that we've already talked about.
#To keep things neat and tidy, I'm only analyzing one dataset in this document & workspace

#Example Dataset:
#GSE81672 is an RNA-Seq experiment examining the effects of antidepressants in stressed mice.

##################################

#Setting the working directory
setwd("/Users/hagenaue/Documents/Example_GSE81672")

##################################

#I was previously demonstrating analyses using just the expression data

#Gemma also has the ability to import some larger objects that include sample and gene/probe metadata
#These objects are called Summarized Experiments and Expression Sets
#Gemma recommends using Summarized Experiments, and peeking at these objects they definitely contain more info
#... but I hadn't worked with them before, so I had to read a little about them:
https://bioconductor.org/packages/release/bioc/manuals/SummarizedExperiment/man/SummarizedExperiment.pdf
https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html

#One of the benefits of these two objects is that they also allow us to read in the filtered expression data
#Gemma applies a standardized filter across all platforms (RNA-Seq/microarray) to remove genes with minimal variance
#These genes are likely to suffer from severe floor effects (very low expression/unmeasurable)
#Or not actually be expressed at all and the measurements that we have are just noise.

#I heard back from Paul about the actual cut offs that they use:
#"There are two stages to our variance filter. 
# First is just to remove genes with zero variance. 
# The second is to remove genes that have too many identical values (within a small tolerance), 
# where "too many" is currently defined as <70% distinct values, but we have tweaked this over time (slightly - we're never perfectly happy with anything)."
#I call this a variance filter as a shorthand, of course it's not that exactly.
#But the latter filter was originally intended to address cases where microarray data had been clipped by the submitter.
#But it also has the effect of removing RNA-seq data that has very low counts (like, 1 count in a couple of samples).
#It's quite permissive - it filters out data that causes problems for the analysis, but otherwise lets it go.
#We used to have an expression level filter but we got rid of it long ago as being too contentious."

#I am satisfied by this answer. 
#Excluding low level expressed genes within an individual dataset reduces false positives due to noise within low-variance data
#It also decreases the severity of the multiple comparisons correction (because there are fewer genes in the dataset)
#Within a meta-analysis, we do not need to worry about effects caused by noise in individual datasets quite as much
#So I would prefer to err on the side of including as many genes as possible.
#So Paul's liberal filter to only toss out genes that are likely to cause problems with existing algorithms sounds good.

#Instructions for function 
gemma.R::get_dataset_object
https://rdrr.io/github/PavlidisLab/Gemma-API/man/get_dataset_object.html
#looks like we want to add filter=TRUE to use Gemma's filtered results (dropping genes with low variability)
#looks like we can add consolidate="average" to consolidate the data for probes representing the same gene.
#That won't matter as much for RNA-Seq data
#E.g. it makes a small difference in this RNA-Seq dataset (~600 genes less in the consolidated version)
#If we need to revert back to an RNA-Seq count matrix at some point, we may need to use the version with no consolidation

#Let's compare the distributions of the unfiltered and filtered objects:
#We won't want to do this for every dataset, this is just for demonstration purposes

SummarizedExperiment<-gemma.R::get_dataset_object("GSE81672", type = 'se', consolidate="average")
str(SummarizedExperiment)

library(SummarizedExperiment)

#Getting a matrix of expression data
ExpressionData<-assay(SummarizedExperiment[[1]])
str(ExpressionData)
#num [1:34986, 1:99] -2.87 -5.86 -5.86 -5.86 -5.86 ...

#Making a histogram of the log2 cpm values:
hist(ExpressionData)
#Big floor effect in the unfiltered data

#The minimum log2 cpm value:
min(ExpressionData)
#[1] -5.860108

#Let's look at the mean vs. variance curve:
ExpressionData_MeanPerGene<-apply(ExpressionData, 1, mean)
ExpressionData_SDPerGene<-apply(ExpressionData, 1, sd)

plot(ExpressionData_SDPerGene~ExpressionData_MeanPerGene)
#This curve is definitely not flat - the log2 transformation did not solve heteroskedasticity issues, especially at the low end

#How many genes have zero variance?
sum(ExpressionData_SDPerGene==0)
#[1] 5032

#Reading in a summarized experiment object for the filtered data:

SummarizedExperiment_Filtered<-gemma.R::get_dataset_object("GSE81672", type = 'se', filter=TRUE, consolidate="average")
# SummarizedExperiment_Filtered
# $`13458`
# class: SummarizedExperiment 
# dim: 22782 99 
# metadata(8): title abstract ... GemmaSuitabilityScore taxon
# assays(1): counts
# rownames(22782): 100009600 100017 ... Averaged from 100039542 100993 Averaged from 11641 677884
# rowData names(4): Probe GeneSymbol GeneName NCBIid
# colnames(99): Sample 18: AMY_resilient_saline Sample 20: AMY_susceptible_ketamine_responder ...
# Sample 108: HIP_susceptible_saline Sample 112: AMY_susceptible_imipramine_non_responder
# colData names(5): factorValues organism part block phenotype treatment

ExpressionData_Filtered<-assay(SummarizedExperiment_Filtered[[1]])
str(ExpressionData_Filtered)
#num [1:22782, 1:99] -2.87 1.73 6.71 -5.09 -2.42 ...

hist(ExpressionData_Filtered)
#The floor effect is much smaller now.

min(ExpressionData_Filtered)
#[1] -5.860108

#Let's look at the mean vs. variance curve:
ExpressionData_Filtered_MeanPerGene<-apply(ExpressionData_Filtered, 1, mean)
ExpressionData_Filtered_SDPerGene<-apply(ExpressionData_Filtered, 1, sd)

plot(ExpressionData_Filtered_SDPerGene~ExpressionData_Filtered_MeanPerGene)
#This curve is more flat, but still has heteroskedasticity issues, especially at the low end
#The voom function will help this data still be useable in regression equations

#How many genes have zero variance?
sum(ExpressionData_Filtered_SDPerGene==0)
#[1] 0
#good - those genes were filtered out as expected.

#Unfortunately, if we subset the data (e.g., to focus on a region or specific groups) we may still have genes with zero variance 
#So we'll have to come back and filter some more later

###########################

#Dataset Subsetting:

#Before we do much more with the dataset, let's subset down to the samples that we actually plan to use:

#First, we need to know what we have:

#How to access different parts of the Summarized Experiment object:

colData(SummarizedExperiment_Filtered[[1]])
#Sample data - organism part, treatment, batch, etc

rowData(SummarizedExperiment_Filtered[[1]])
#All of the annotation - Probe, GeneSymbol, GeneName, NCBIid

#The distribution of samples from different brain regions
table(SummarizedExperiment_Filtered[[1]]$`organism part`)
# basolateral amygdaloid nuclear complex                      nucleus accumbens 
# 25                                     25 
# prefrontal cortex                   ventral,Ammon's horn 
#                                     25                                     24

#The brain region that I'm interested is the hippocampus, so I will subset down to the "ventral,Ammon's horn" 

#The distribution of the different phenotypes in this experiment
table(SummarizedExperiment_Filtered[[1]]$phenotype)
# Non-responder,susceptible toward                    resistant to               susceptible toward 
# 28                               16                               12 
# susceptible toward,Responder               wild type genotype 
# 24                               19 

#This experiment has a weird design
#I'm interested in antidepressant effects, and it looks like only the stress susceptible subjects were treated
#So I am going to subset down to just the stress susceptible groups
#I'm going to keep both antidepressant responders and non-responders 
#because I'm looking at antidepressant effects (not antidepressant effectiveness)
#and most of my other datasets will not screen for antidepressant effectiveness

#When combining together criteria "&" means "AND", "|" (also called the pipe) means "OR"
#"==" means equals, whereas "!=" means doesn't equal
SampleFilter<-
  SummarizedExperiment_Filtered[[1]]$`organism part`=="ventral,Ammon's horn" &
  (SummarizedExperiment_Filtered[[1]]$phenotype=="susceptible toward,Responder"|
     SummarizedExperiment_Filtered[[1]]$phenotype=="Non-responder,susceptible toward"|
     SummarizedExperiment_Filtered[[1]]$phenotype=="susceptible toward")
   
#Subsetting the data to have only hippocampus:
SummarizedExperiment_Subset<-SummarizedExperiment_Filtered[[1]][,SampleFilter]

SummarizedExperiment_Subset
# class: SummarizedExperiment 
# dim: 22782 16 
#This dataset is much smaller now

#Sanity Check: Double-checking that the subsetting worked properly:

table(SummarizedExperiment_Subset$`organism part`)
#ventral,Ammon's horn 
#16 

table(SummarizedExperiment_Subset$phenotype, SummarizedExperiment_Subset$treatment)
#                                     imipramine ketamine reference substance role,saline
# Non-responder,susceptible toward          4        3                               0
# susceptible toward                        0        0                               3
# susceptible toward,Responder              3        3                               0

#The sample sizes are now pretty sad. :(
#Good thing we are running a meta-analysis.

#Pulling out a matrix of gene expression data for this subset (to use in functions that require matrices)
ExpressionData_Subset<-assay(SummarizedExperiment_Subset)
str(ExpressionData_Subset)
# num [1:22782, 1:16] -2.39 2.22 6.36 -4.27 -3.54 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:22782] "100009600" "100017" "100019" "100033459" ...
# ..$ : chr [1:16] "Sample 15: HIP_susceptible_imipramine_non_responder" "Sample 37: HIP_susceptible_ketamine_non_responder" "Sample 13: HIP_susceptible_saline" "Sample 84: HIP_susceptible_imipramine_responder" ...

###########################

#Outlier Removal:

#This would be a good time to check for outliers and remove them if they are present.

#Creating an example sample-sample correlation scatterplot (data for all genes for 1 sample versus the data for all genes for the second sample)
plot(ExpressionData_Subset[,1]~ExpressionData_Subset[,2])

#Creating a matrix showing how each sample correlates with every other sample:
CorMatrix<-cor(ExpressionData_Subset)

#Writing that matrix out to your working directory to save it:
write.csv(CorMatrix, "CorMatrix.csv")

#Creating a hierarchically clustered heatmap illustrating the sample-sample correlation matrix:
heatmap(CorMatrix)
#Cool - there is actually some clustering amongst the samples that are from the antidepressant vs. saline groups

#Creating a boxplot illustrating the sample-sample correlations for each sample. Outliers should be obvious at this point.
boxplot(CorMatrix)
#None of these samples look like clear outliers

#If there was an outlier sample, you could remove it using subsetting similar to above by identifying it's column name
#E.g.
OutlierFilter<-colnames(ExpressionData_Subset)!="Sample 15: HIP_susceptible_imipramine_non_responder"
OutlierFilter
#[1] FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE

SummarizedExperiment_Subset_noBad<-SummarizedExperiment_Subset[,OutlierFilter]
SummarizedExperiment_Subset_noBad
# class: SummarizedExperiment 
# dim: 22782 15 

#If we don't have any outliers to remove, we can just rename our object so it works with the downstream code:
SummarizedExperiment_Subset_noBad<-SummarizedExperiment_Subset

#And then we will need to recreate the ExpressionData Matrix as well:
ExpressionData_Subset_noBad<-assay(SummarizedExperiment_Subset_noBad)
str(ExpressionData_Subset_noBad)
#num [1:22782, 1:16] -2.39 2.22 6.36 -4.27 -3.54 ...


###########################


#Filtering Genes...Again:

#Now that we have subsetted our data, we have a new problem:
#Gemma filtered out genes that lacked variability in the full dataset
#..but that doesn't mean all of the remaining genes have variability in this particular subset of samples

ExpressionData_Subset_noBad_SDperGene<-apply(ExpressionData_Subset_noBad, 1, sd)
min(ExpressionData_Subset_noBad_SDperGene)
#[1] 0

sum(ExpressionData_Subset_noBad_SDperGene==0)
#[1] 121
#121 genes have zero variability in the HC dataset.

#These genes are going to cause problems - you can't run stats on a variable with no variability! - let's get rid of them.

#...But there is more:
#We may still have issues with genes that lack any variability *for any particular subgroup of interest* as well.

#This function calculates the sd for each treatment group for a particular gene (row of data):
tapply(ExpressionData_Subset_noBad[1,], SummarizedExperiment_Subset_noBad$treatment, sd)
# imipramine                        ketamine reference substance role,saline 
# 0.6547312                       0.9273463                       0.7106250

#We want the minimum sd for all of our groups to not be zero:
min(tapply(ExpressionData_Subset_noBad[1,], SummarizedExperiment_Subset_noBad$treatment, sd))!=0
#[1] TRUE

#... and we need to know that for all rows.
GenesToFilter<-apply(ExpressionData_Subset_noBad, 1, function(y) (min(tapply(y, SummarizedExperiment_Subset_noBad$treatment, sd))!=0))

#How many genes we'll end up keeping
sum(GenesToFilter)
#[1] 21686
#So we're not losing that many:
22782-21686
#[1] 1096

SummarizedExperiment_Subset_noBad_Filtered<-SummarizedExperiment_Subset_noBad[GenesToFilter,]
SummarizedExperiment_Subset_noBad_Filtered
#class: SummarizedExperiment 
#dim: 21686 16 

#And again, remaking the expression set:
ExpressionData_Subset_noBad_Filtered<-assay(SummarizedExperiment_Subset_noBad_Filtered)
str(ExpressionData_Subset_noBad_Filtered)
#num [1:21686, 1:16] -2.39 2.22 6.36 -4.27 -3.54 ...

###########################

#Checking for batch confound:

#Since we are down to the subset of samples that we plan to use, this would be a good time to check for confounds
#Unfortunately, processing batches are often unbalanced in regards to variables of interest

#Currently, Gemma has all of the processing batch information lumped into one variable

table(SummarizedExperiment_Subset_noBad_Filtered$block)
# Device=11V6WR1:Run=206:Flowcell=C3HBLACXX:Lane=1    Device=11V6WR1:Run=206:Flowcell=C3HBLACXX:Lane=3 
# 1                                                   3 
# Device=HWI-D00147:Run=59:Flowcell=C2GWMACXX:Lane=2  Device=HWI-D00147:Run=59:Flowcell=C2GWMACXX:Lane=4 
# 3                                                   2 
# Device=HWI-D00147:Run=59:Flowcell=C2GWMACXX:Lane=5  Device=HWI-D00147:Run=59:Flowcell=C2GWMACXX:Lane=6 
# 1                                                   1 
# Device=HWI-D00147:Run=59:Flowcell=C2GWMACXX:Lane=7  Device=HWI-D00147:Run=59:Flowcell=C2GWMACXX:Lane=8 
# 1                                                   1 
# Device=HWI-ST1147:Run=158:Flowcell=C34A7ACXX:Lane=3 Device=HWI-ST1147:Run=158:Flowcell=C34A7ACXX:Lane=7 
# 1                                                   1 
# Device=HWI-ST276:Run=345:Flowcell=C3H70ACXX:Lane=1 
# 1 

#this function breaks apart these "blocks" into specific batch-related variables:
strsplit(SummarizedExperiment_Subset_noBad_Filtered$block, ":")

#To make that into an easier-to-use data.frame
BatchVariables<-do.call(rbind.data.frame, strsplit(SummarizedExperiment_Subset_noBad_Filtered$block, ":")) 
str(BatchVariables)
#I'm going to rename these (note - this will be dataset specific)
colnames(BatchVariables)<-c("Device", "Run", "FlowCell", "Lane")

table(SummarizedExperiment_Subset_noBad_Filtered$treatment, BatchVariables$Device)
#                                   Device=11V6WR1 Device=HWI-D00147 Device=HWI-ST1147 Device=HWI-ST276
# imipramine                                   2                 4                 1                0
# ketamine                                     1                 4                 1                0
# reference substance role,saline              1                 1                 0                1

table(SummarizedExperiment_Subset_noBad_Filtered$treatment, BatchVariables$Run)
#                                   Run=158 Run=206 Run=345 Run=59
# imipramine                            1       2       0      4
# ketamine                              1       1       0      4
# reference substance role,saline       0       1       1      1

#Looks like run and device may be redundant
table(BatchVariables$Device, BatchVariables$Run)
#                     Run=158 Run=206 Run=345 Run=59
# Device=11V6WR1          0       4       0      0
# Device=HWI-D00147       0       0       0      9
# Device=HWI-ST1147       2       0       0      0
# Device=HWI-ST276        0       0       1      0

table(SummarizedExperiment_Subset_noBad_Filtered$treatment, BatchVariables$FlowCell)
#Looks also potentially redundant

table(BatchVariables$FlowCell, BatchVariables$Run)
#                       Run=158 Run=206 Run=345 Run=59
# Flowcell=C2GWMACXX       0       0       0      9
# Flowcell=C34A7ACXX       2       0       0      0
# Flowcell=C3H70ACXX       0       0       1      0
# Flowcell=C3HBLACXX       0       4       0      0
#yep

table(SummarizedExperiment_Subset_noBad_Filtered$treatment, BatchVariables$Lane)
#                                   Lane=1 Lane=2 Lane=3 Lane=4 Lane=5 Lane=6 Lane=7 Lane=8
# imipramine                           1      1      1      1      0      1      1      1
# ketamine                             0      1      2      1      1      0      1      0
# reference substance role,saline      1      1      1      0      0      0      0      0

#Too many Lanes to take into consideration

###########################

#Renaming conditions
#Some of the levels for my variables are a mouthful. 
#I'm going to rename them:

table(SummarizedExperiment_Subset_noBad_Filtered$treatment)
#imipramine               ketamine reference substance role,saline 
#7                               6                              11 

#I'm going to rename the control condition because it is a mouthful
SummarizedExperiment_Subset_noBad_Filtered$treatment[SummarizedExperiment_Subset_noBad_Filtered$treatment=="reference substance role,saline"]<-"saline"

table(SummarizedExperiment_Subset_noBad_Filtered$treatment)
# imipramine   ketamine     saline 
# 7          6         11

SummarizedExperiment_Subset_noBad_Filtered$phenotype[SummarizedExperiment_Subset_noBad_Filtered$phenotype=="Non-responder,susceptible toward"]<-"non-responder"
SummarizedExperiment_Subset_noBad_Filtered$phenotype[SummarizedExperiment_Subset_noBad_Filtered$phenotype=="susceptible toward,Responder"]<-"responder"
SummarizedExperiment_Subset_noBad_Filtered$phenotype[SummarizedExperiment_Subset_noBad_Filtered$phenotype=="susceptible toward"]<-"untreated"

table(SummarizedExperiment_Subset_noBad_Filtered$phenotype)
# non-responder     responder     untreated 
# 7             6             3

###########################

#Principal components analysis:

library(stats)
pca_output<-prcomp(t(ExpressionData_Subset_noBad_Filtered), scale=TRUE)

PCeigenvectors<-pca_output$rotation[ ,c(1:4)]
PCeigenvectors2<-cbind(PCeigenvectors, rowData(SummarizedExperiment_Subset_noBad_Filtered))
write.csv(PCeigenvectors2, "PCeigenvectors.csv")

PC1<-pca_output$x[,1]
PC2<-pca_output$x[,2]
PC3<-pca_output$x[,3]
PC4<-pca_output$x[,4]

#Output a scree plot for the PCA (no outliers):
#This plot illustrates the proportion of variance explained by each principal component (PC):
png("10 PCA Scree Plot1.png")
plot(summary(pca_output)$importance[2,]~(c(1:length(summary(pca_output)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC1 vs PC2.png")
plot(PC1~PC2, main="Principal Components Analysis")
dev.off()

#You can color these plots in using different variables in the dataset 
#This can help explain the main sources of variation in the data
#Technical variables (brain region, batch) tend to be the main culprits
plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(BatchVariables$Run))
#PC1 and 2, don't seem to be related to our main batch variable
plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Subset_noBad_Filtered$treatment))
#But they do seem related to antidepressant treatment, esp. PC1
plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Subset_noBad_Filtered$phenotype))
#Less related to phenotype

#If we want to zoom in on the relationship between PC1 and treatment we can make a boxplot
boxplot(PC1~SummarizedExperiment_Subset_noBad_Filtered$treatment, las=2, xlab="")
#Oh interesting - it is actually the ketamine samples that look different (from imipramine and saline)

#############################

#Plotting data for a particular gene:

#Getting the expression data for a particular gene:
Pvalb<-assay(SummarizedExperiment_Subset_noBad_Filtered)[rowData(SummarizedExperiment_Subset_noBad_Filtered)$GeneSymbol=="Pvalb",]

#Plotting a boxplot of the expression of Pvalb Log2 CPM across antidepressant groups:
boxplot(Pvalb~SummarizedExperiment_Subset_noBad_Filtered$treatment, xlab="", ylab="Log2 CPM", main="Pvalb", las=2)

#We can add jittered data points to our graph so that we can see the values for the individual samples in each group.
#To do this, we use the function stripchart and the exact same y~x formula that we used for the boxplot
#The parameter pch determines the size of the data points.
#The parameter "add" places the points on top of the boxplot that we already created.
stripchart(Pvalb~SummarizedExperiment_Subset_noBad_Filtered$treatment, pch = 19, method = "jitter", jitter = 0.2, vertical = TRUE, add=TRUE)

#It looks like antidepressant treatment might decrease PVALB
#We would need to run inferential statistics to learn whether this effect is significant

#If we were looking at the data from a single gene, and our data was truly numeric (which RNA-Seq is not quite...)
#We could run a simple linear regression

#First we would just set up our variables so that they have intuitive reference groups:
SummarizedExperiment_Subset_noBad_Filtered$treatment_factor<-as.factor(SummarizedExperiment_Subset_noBad_Filtered$treatment)
levels(SummarizedExperiment_Subset_noBad_Filtered$treatment_factor)
#[1] "imipramine" "ketamine"   "saline" 
#We want saline to be our control (reference) group:
SummarizedExperiment_Subset_noBad_Filtered$treatment_factor<-relevel(SummarizedExperiment_Subset_noBad_Filtered$treatment_factor, ref="saline")

#And then we could just run linear regression:
summary.lm(lm(Pvalb~SummarizedExperiment_Subset_noBad_Filtered$treatment_factor))
# Call:
#   lm(formula = Pvalb ~ SummarizedExperiment_Subset_noBad_Filtered$treatment_factor)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.28155 -0.10604 -0.04275  0.03161  0.76254 
# 
# Coefficients:
#                                                                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                            4.67023    0.14493  32.225 8.67e-14 ***
# SummarizedExperiment_Subset_noBad_Filtered$treatment_factorimipramine -0.11417    0.17322  -0.659    0.521    
# SummarizedExperiment_Subset_noBad_Filtered$treatment_factorketamine   -0.03725    0.17750  -0.210    0.837    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.251 on 13 degrees of freedom
# Multiple R-squared:  0.04008,	Adjusted R-squared:  -0.1076 
# F-statistic: 0.2714 on 2 and 13 DF,  p-value: 0.7665

#The Intercept is the saline group (reference)
#The Estimate is the Log2FC (the difference between the group of interest and the saline group)
#The Pr(>|t|) is the p-value showing the probability of this effect arising under conditions of random chance.

######################

#But, of course, running a differential expression analysis for an entire dataset is a little more complicated than that...
#And for RNA-Seq it is even more complicated...



