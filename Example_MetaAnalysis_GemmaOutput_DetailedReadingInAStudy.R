#Messing around with meta-analyses of Gemma output
#This file includes the detailed version of how to read in and extract the preprocessed data and statistical output for a single study (steps #0-5 mentioned below).  I'll put a more streamlined version in a separate file.
#Megan Hagenauer
#07-21-2022
#Associated R workspace: Workspace_Example_MetaAnalysis_GemmaOutput.RData

#At the top of the code document, I recommend giving the code document:
#1) a title that will help you figure out what it contains in the future
#2) your name, in case you share the code
#3) the date (or general dates) that you worked on it.
#4) the R workspace associated with the code document.

#Sometimes folks also include a general overview of what their goals were, as well as the steps used to reach those goals.
#Outlining the coding steps that you plan to use to reach your goals is called "pseudo code". It is a good place to start for any project

####################


#Goal:

#For practice, I decided to focus on hippocampal chronic stress datasets that weren't included in Yusra's HDRF meta-analysis

#I piggybacked on Jinglin Xiong's efforts triaging chronic stress datasets, using her spreadsheet: "ChronicStressData(after exclusion).xlsx"

#These were the hippocampal datasets that seemed to meet criteria (chronic adult stress):
#GSE59070 - chronic social stress (CSDS), "Ammon's horn", mouse, Agilent platform
#GSE86392 - chronic restraint stress, HC, rat, RNA-Seq
#GSE116009 - CUMS, HC, mouse, RNA-Seq
#GSE109315 -CSDS, ventral (V) HC - mouse, RNA-Seq
#GSE151807 - CMS, HC, mice, Affymetrix 430 2.0

#These datasets were hippocampal, but focused on specific subregions:
#GSE11211 - CSDS, CA1
#GSE56028 - unpredictable chronic mild stress (CUMS), Dentate Gyrus (DG)
#GSE84185 - CUMS, DG
#GSE132819 - CSDS, DG
#GSE84183 - CUMS, DG

#Excluded, included in our HDRF meta-analysis:
#GSE81672 - CSDS, HC - Bagot et al.
#GSE72343 - CSDS, "Ammon's Horn" - Bagot et al.

#Excluded, no associated publication:
#GSE109445 - CUMS, "Ammon's horn"

#Excluded - focused on stress during adolescence instead of adulthood:
#GSE172451- chronic social instability stress, vHC

#Excluded:
#GSE102965 - Learned helplessness, HC - this one may not actually work for the simple type of analysis that we're doing - the Gemma results aren't specific to a single brain area, big frowny face on Gemma

####################

#Overview of general coding steps:

#0) Reading in & visualizing the Log2 Expression data and Metadata for the individual samples - an illustration of where the differential expression results come from. 

#1) Read in the results from Gemma for each dataset

#2) Identify the results for the variables of interest

#3) Remove rows of data that have missing or unambiguous gene annotation

#### Question: Which gene annotation should we use? I was originally planning to use NCBI ID (because it is a little more stable than gene symbol), but if we use gene symbol we can run a half-way decent meta-analysis of data from both rats and mice without having to add in a step where we decode which genes in rats and mice are orthologous using an orthology database, as many mice and rat genes have the same gene symbol (76%, last time I checked).

#4) Collapse results (average) if there is more than row of data representing a single gene 

#5) Extract out the information needed for running a meta-analysis: Use the Log2FC and Tstat to calculate the standard error for the Log2FC, and then use the standard error to calculate the sampling variance.

#6) Combine together the relevant results from different studies into a single data-frame for the effect sizes (Log2FC) and a single data.frame for the sampling variances.  - *This will go in a separate code file.*

#7) Make a correlation matrix to compare the overall results from the different studies. Further visualize the comparison using a hierarchically-clustered heatmap of the correlation matrix. Which studies seem to have similar results? - *This will go in a separate code file.*

#8) Run a meta-analysis using all of the effect sizes for each gene that has data from at least 2 studies. - *This will go in a separate code file.*

#9) Correct the meta-analysis output to take into account the fact that we are running the statistical calculations many times and therefore have a heightened risk of false discovery (false discovery rate correction) - *This will go in a separate code file.*

#10) Determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. - *This will go in a separate code file.*

#####################################

#Code packages used (may require installation & installation of dependencies):

library(plyr)
library(metafor)


######################################

#0) Reading in & visualizing the Log2 Expression data and Metadata for the individual samples - an illustration of where the differential expression results come from. 


#After writing the code necessary for the meta-analysis of differential expression results, I realized that it might be helpful to start out with an illustration of what the Log2 gene expression data from the individual samples looks like.


#To play with the Log2 Expression data and Metadata for the individual samples, the easiest thing to do is to go to the Gemma website for the dataset, e.g.:
https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=11627

#And then download one of the "profiles" into an intuitively named folder/location:
###the unfiltered version is the Log2 gene expression data for all genes/probes represented on the transcriptional profiling platform for all samples
###the filtered version is the Log2 gene expression data for all genes/probes that meet some quality control criteria (typically sufficient expression levels to be considered reliably "measurable") for all samples

#The file containing the Log2 gene expression data for all samples (in this case, 11627_GSE59070_expmat.data.txt.gz) is in a compressed format and will probably need to be "unzipped" to work with it. Typically double-clicking on the file in your folder will take care of this.

#Then set the working directory for R to the folder where you stashed the downloaded data:
setwd("~/Documents/SideProjects/BrainGMT/Gemma/Hippocampus/11627_GSE59070_diffExpAnalysis_94802")

#The data is read in as a tab ("\t") delimited text file.
#We set stringsAsFactors as false, so that R doesn't get confused and think that every gene symbol is a different "level" in a factor variable
#This data file contains metadata at the top. If we try to read that into R, R will get confused about the file structure and not know that the rest of the file is a giant data frame. Therefore, we tell R to ignore lines of data that are "commented" using a pound sign.
TempDataSet_Filtered<-read.delim("11627_GSE59070_expmat.data.txt", sep="\t", stringsAsFactors = FALSE, comment.char = "#")

#Here's the structure of the dataset:
#Each of the 37,371 rows represents the Log2-transformed gene expression data measured using a particular probe on this microarray platform. 
#Several of the columns contain probe annotation information - e.g., the Gene Symbol, NCBI ID, etc
#All of the other columns represent the Log2 expression data for a sample, and the columns contain all of the defining characterstics of that sample.
str(TempDataSet_Filtered)
# 'data.frame':	37371 obs. of  30 variables:
#   $ Probe                                                                                              : chr  "A_51_P100021" "A_51_P100034" "A_51_P100052" "A_51_P100063" ...
# $ Sequence                                                                                           : chr  "2290_sequence" "1644_sequence" "7021_sequence" "20775_sequence" ...
# $ GeneSymbol                                                                                         : chr  "Hivep3" "Mif4gd|Mrps7" "Slitrk2" "Fip1l1|Lnx1" ...
# $ GeneName                                                                                           : chr  "human immunodeficiency virus type I enhancer binding protein 3" "MIF4G domain containing|mitchondrial ribosomal protein S7" "SLIT and NTRK-like family, member 2" "FIP1 like 1 (S. cerevisiae)|ligand of numb-protein X 1" ...
# $ GemmaId                                                                                            : chr  "532267" "682279|609073" "838885" "652416|535115" ...
# $ NCBIid                                                                                             : chr  "16656" "69674|50529" "245450" "66899|16924" ...
# $ GSE59070_Biomat_3___BioAssayId.424074Name.Hippocampus.Stress.acute.Biol.Rep.2.dyeswap              : num  0.0367 -0.0554 -0.038 0.2107 -0.0744 ...
# $ GSE59070_Biomat_24___BioAssayId.424075Name.Hippocampus.Stress.acute.Biol.Rep.2                     : num  -0.1154 -0.1459 0.2439 -0.0535 0.1614 ...
# $ GSE59070_Biomat_22___BioAssayId.424076Name.Hippocampus.Stress.acute.Biol.Rep.1.dyeswap             : num  0.3253 0.04 -0.1892 0.3261 0.0868 ...
# $ GSE59070_Biomat_19___BioAssayId.424077Name.Hippocampus.Stress.acute.Biol.Rep.1                     : num  0.0322 -0.0571 -0.395 -0.1841 0.0409 ...
# $ GSE59070_Biomat_7___BioAssayId.424072Name.Hippocampus.Stress.acute.Biol.Rep.3.dyeswap              : num  0.4487 -0.0363 0.2131 0.0795 0.0295 ...
# $ GSE59070_Biomat_5___BioAssayId.424073Name.Hippocampus.Stress.acute.Biol.Rep.3                      : num  -0.3191 -0.1093 -0.0815 -0.0926 0.1974 ...
# $ GSE59070_Biomat_15___BioAssayId.424068Name.Hippocampus.Stress.8days.Biol.Rep.2.dyeswap             : num  -0.086 0.1455 0.0631 0.0593 0.0278 ...
# $ GSE59070_Biomat_13___BioAssayId.424069Name.Hippocampus.Stress.8days.Biol.Rep.2                     : num  -0.4003 0.1001 -0.0373 -0.1439 0.0532 ...
# $ GSE59070_Biomat_11___BioAssayId.424070Name.Hippocampus.Stress.8days.Biol.Rep.1.dyeswap             : num  0.0173 0.315 -0.4984 0.4893 -0.0259 ...
# $ GSE59070_Biomat_9___BioAssayId.424071Name.Hippocampus.Stress.8days.Biol.Rep.1                      : num  -0.4401 -0.1497 0.1392 -0.0563 0.1525 ...
# $ GSE59070_Biomat_20___BioAssayId.424066Name.Hippocampus.Stress.8days.Biol.Rep.3.dyeswap             : num  0.1635 -0.2634 -0.1439 0.3025 -0.0341 ...
# $ GSE59070_Biomat_17___BioAssayId.424067Name.Hippocampus.Stress.8days.Biol.Rep.3                     : num  -0.39151 -0.03517 -0.12069 0.00596 -0.12932 ...
# $ GSE59070_Biomat_6___BioAssayId.424062Name.Hippocampus.Stress.13days.Biol.Rep.2.dyeswap             : num  0.1622 -0.0533 -0.1609 0.021 -0.1497 ...
# $ GSE59070_Biomat_4___BioAssayId.424063Name.Hippocampus.Stress.13days.Biol.Rep.2                     : num  0.00892 0.01282 0.15891 -0.26541 0.02349 ...
# $ GSE59070_Biomat_2___BioAssayId.424064Name.Hippocampus.Stress.13days.Biol.Rep.1.dyeswap             : num  0.4176 -0.0159 0.0915 0.2216 0.213 ...
# $ GSE59070_Biomat_23___BioAssayId.424065Name.Hippocampus.Stress.13days.Biol.Rep.1                    : num  -0.418 0.12 -0.296 -0.331 0.11 ...
# $ GSE59070_Biomat_10___BioAssayId.424060Name.Hippocampus.Stress.13days.Biol.Rep.3.dyeswap            : num  0.2164 -0.0456 -0.1167 0.0968 -0.1678 ...
# $ GSE59070_Biomat_8___BioAssayId.424061Name.Hippocampus.Stress.13days.Biol.Rep.3                     : num  -0.0669 -0.0611 0.6093 -0.1309 0.0111 ...
# $ GSE59070_Biomat_18___BioAssayId.424056Name.Hippocampus.Stress.13days.5daysofrest.Biol.Rep.2.dyeswap: num  0.2435 0.1271 -0.0674 0.1538 -0.1622 ...
# $ GSE59070_Biomat_16___BioAssayId.424057Name.Hippocampus.Stress.13days.5daysofrest.Biol.Rep.2        : num  -0.1247 0.0767 -0.0903 -0.0619 -0.0987 ...
# $ GSE59070_Biomat_14___BioAssayId.424058Name.Hippocampus.Stress.13days.5daysofrest.Biol.Rep.1.dyeswap: num  0.7005 -0.0712 -0.3043 0.1158 -0.0723 ...
# $ GSE59070_Biomat_12___BioAssayId.424059Name.Hippocampus.Stress.13days.5daysofrest.Biol.Rep.1        : num  -0.3158 -0.1181 -0.2702 -0.195 0.0385 ...
# $ GSE59070_Biomat_1___BioAssayId.424054Name.Hippocampus.Stress.13days.5daysofrest.Biol.Rep.3.dyeswap : num  -0.352 -0.2524 0.1796 0.0943 0.0924 ...
# $ GSE59070_Biomat_21___BioAssayId.424055Name.Hippocampus.Stress.13days.5daysofrest.Biol.Rep.3        : num  -0.855 0.167 0.547 -0.119 -0.104 ...

#The first 6 columns are probe annotation - Let's separate them out:
TempDataSet_Filtered_Annotation<-TempDataSet_Filtered[,c(1:6)]
TempDataSet_Log2Expression<-TempDataSet_Filtered[,-c(1:6)]

#The sample meta-data is located in the column names, with individual variables separated by "."
#I would like to separate this data out into a matrix where I can grab the values associated with each of the variables for each sample easily 
#I use the function string-split (strsplit) to do this, because I'm splitting up a long string of characters (the column name)
#Because "." is a special character, I have to put two "\\" in front of it to let R know that I want it to recognize the "." as a character instead of part of the coding.
#Example of splitting up the column names:
Temp_Colnames_SplitIntoList<-strsplit(colnames(TempDataSet_Log2Expression), "\\.")

#Now each sample (previous column name) is treated as an element in a list, and within that element is all of the sample information. 
#It would be much easier to work with this information if it was in a data.frame format.
#The function do.call() applies a function to each element of a list (in this case, each sample or previous column name)
#The function that I'm applying to each element in the list is to make the content within that element (all of the variables that I split apart from the column name) a row in the data.frame.
#Note that this code only works if each element in the list has the same number of components (in this case, values for each variable)
TempDataSet_MetaData<-do.call(rbind.data.frame, Temp_Colnames_SplitIntoList)

#Let's look at the structure for the data.frame of sample metadata that we created:
str(TempDataSet_MetaData)
#The columns don't have variable names yet. If you click on the object in the overview of the Environment (on the right), it will open it up in a manner like a spreadsheet so that you can view it easily 

#From looking at it, the information that defines our main treatment groups is found in these two columns:
TempDataSet_MetaData[,5]
# [1] acute  acute  acute  acute  acute  acute  8days  8days  8days  8days  8days  8days  13days 13days
# [15] 13days 13days 13days 13days 13days 13days 13days 13days 13days 13days
# Levels: 13days 8days acute
TempDataSet_MetaData[,6]
# [1] Biol        Biol        Biol        Biol        Biol        Biol        Biol        Biol       
# [9] Biol        Biol        Biol        Biol        Biol        Biol        Biol        Biol       
# [17] Biol        Biol        5daysofrest 5daysofrest 5daysofrest 5daysofrest 5daysofrest 5daysofrest
# Levels: 5daysofrest Biol

#I'm going to combine the two columns into a vector containing a single grouping variable:
#I'm using the paste function, and separating the two pieces of information for each sample with a "_"
TempDataSet_StressVar<-paste(TempDataSet_MetaData[,5], TempDataSet_MetaData[,6], sep="_")
str(TempDataSet_StressVar)
#chr [1:24] "acute_Biol" "acute_Biol" "acute_Biol" "acute_Biol" "acute_Biol" "acute_Biol" ...

#Alright, let's practice visualizing the log2 expression data for the samples:
#Let's play with the Log2 gene expression data measured using the first probe listed in the dataset.

#This is the annotation that accompanies the probe for the first row of data:
TempDataSet_Filtered_Annotation[1,]
# Probe      Sequence GeneSymbol                                                       GeneName
# 1 A_51_P100021 2290_sequence     Hivep3 human immunodeficiency virus type I enhancer binding protein 3
# GemmaId NCBIid
# 1  532267  16656

#So this probe was targeting the Hivep3 gene

#And this is the first row of data:
TempDataSet_Log2Expression[1,]

#We can visualize the Log2 expression in the different groups using a boxplot
#Boxplots are better than commonly used bar charts because they illustrate the distribution of the data (not just the mean value for each group)
#To make a boxplot, I have to force R to recognize that the data from different columns (samples) from the same row (probe) can be treated as a numeric vector using as.numeric()
#The ~ is used to write out the formula for the plot as y~x. The y (dependent variable) in this case is Log2 expression data for the first probe in the dataset (targeting Hivep3), where as the x (predictor or independent variable) is the stress treatment group
#The parameters ylab and xlab allow you to provide labels for the x and y axes of the chart.
#The parameter col allows you to color in the chart. There are many ways to specify colors in R. In this case, I'm going to use the easy 0-9 numbering system. If you would like a wider variety of colors, google the topic and you can gain access to a wide variety of named colors or just the straight up RGB hexadecimal system.
boxplot(as.numeric(TempDataSet_Log2Expression[1,])~TempDataSet_StressVar, ylab="Hivep3 Log2 Expression", xlab="Stress Condition", col=8)
#When reading a boxplot, the middle black line represents the median of the data for each group, the box represents the interquartile range (25the percentile-75th percentile), and the "whiskers" (dotted lines) represent either the full range of the data for the group or 1.5x the interquartile range (with "outlier" data points beyond that range shown as open circles)

#We can add jittered data points to our graph so that we can see the values for the individual samples in each group.
#To do this, we use the function stripchart and the exact same y~x formula that we used for the boxplot
#The parameter pch determines the size of the data points.
#The parameter "add" places the points on top of the boxplot that we already created.
stripchart(as.numeric(TempDataSet_Log2Expression[1,])~TempDataSet_StressVar, pch = 19, method = "jitter", jitter = 0.2, vertical = TRUE, add=TRUE)
#Note that if you accidentally run this code twice, you will add a second set of data points that have the same height on the vertical y-axis (i.e., the same y values) but are in a slightly different spot on the x-axis within the box for their treatment group - that is because the "jitter" assigning variation to the x-axis location is random. 

#Without adding jitter, all of the data points would be in a straight vertical line (because there isn't actually any true variation within a categorical group!) and often inconveniently overlap with each other, making data visualization difficult. Here's an example of what it looks like without jitter where I've colored the data-points red so that you can see the contrast:
stripchart(as.numeric(TempDataSet_Log2Expression[1,])~TempDataSet_StressVar, pch = 19, method = "jitter", jitter = 0.0, vertical = TRUE, add=TRUE, col=2)

#If we want to save our boxplot (or output it into a file that has better proportions so that we can read the x-axis labels better!), we can write it out into a file:
#The .pdf part of this code opens a connection to a .pdf file of a particular size. For this chart, I made the width bigger than the height.
#The dev.off part of this code closes that connection.
pdf("Boxplot_Hivep3Expression_vs_StressCondition.pdf", height=5, width=7)
boxplot(as.numeric(TempDataSet_Log2Expression[1,])~TempDataSet_StressVar, ylab="Hivep3 Log2 Expression", xlab="Stress Condition", col=8)
stripchart(as.numeric(TempDataSet_Log2Expression[1,])~TempDataSet_StressVar, pch = 19, method = "jitter", jitter = 0.2, vertical = TRUE, add=TRUE)
dev.off()

#Btw - for future reference, instead of using a row number, you can make charts like this for a particular/gene probe using the gene or probe name too. e.g.,
#This is all of the Log2 gene expression data for Hivep3 in the dataset
TempDataSet_Log2Expression[TempDataSet_Filtered_Annotation$GeneSymbol=="Hivep3",]
#It looks like there are three probes targeting this gene, so we can grab the data for just one of those probes by specifying a row of the data matching the Hivep3 gene:
TempDataSet_Log2Expression[TempDataSet_Filtered_Annotation$GeneSymbol=="Hivep3",][1,]

#And then can toss that information into the boxplot code instead: 
pdf("Boxplot_Hivep3Expression_vs_StressCondition.pdf", height=5, width=7)
boxplot(as.numeric(TempDataSet_Log2Expression[TempDataSet_Filtered_Annotation$GeneSymbol=="Hivep3",][1,])~TempDataSet_StressVar, ylab="Hivep3 Log2 Expression", xlab="Stress Condition", col=8)
stripchart(as.numeric(TempDataSet_Log2Expression[TempDataSet_Filtered_Annotation$GeneSymbol=="Hivep3",][1,])~TempDataSet_StressVar, pch = 19, method = "jitter", jitter = 0.2, vertical = TRUE, add=TRUE)
dev.off()

#Alright, so what sort of pattern of group differences do we see for this probe/gene?

#Within this study, the "reference group" - or the group that all other groups are compared to - is acute stress (i.e., there is no true control no-stress condition - bah!). 
#When discussing differential expression, the Log2 expression in the reference group is subtracted from the other treatment groups, e.g., 8days_Biol - acute_Biol.  
#Therefore, if the condition of interest (e.g., 8days_Biol) has lower Log2 gene expression that the reference group, the difference between the groups (Log2FoldChange or Log2FC) will be negative
#If the condition of interest has higher log2 gene expression than the reference group, the difference between the groups (Log2FC) will be positive

#Here's an example of the calculations done for this single probe/gene:

#We can use tapply to apply an average to each level in a grouping variable (e.g., our stress variable):
tapply(as.numeric(TempDataSet_Log2Expression[1,]), TempDataSet_StressVar, mean)
# 13days_5daysofrest  13days_Biol         8days_Biol         acute_Biol 
# -0.11721922         0.05331872        -0.18953243         0.06805238 

#To calculate the Log2FC for 8days_Biol vs. acute_Biol:
-0.18953243- 0.06805238 
#[1] -0.2575848

#To calculate the Log2FC for 13days_Biol:
0.05331872-0.06805238
#[1] -0.01473366
#Note that the Log2FC for 13days_Biol vs. acute_Biol is negative, even though in the boxplot it looks like 13days_Biol has slightly higher Log2 expression - that is because the Log2FC is calculated using the mean of the groups, whereas the boxplots illustrate the medians and interquartile ranges.

#Sometimes when discussing group differences, we use a standardized "effect size" called Cohen's d. This the difference between the groups (Log2FC), but in units representing the average amount of variability within the groups.

#The average variation within the groups is called standard deviation (sd):
tapply(as.numeric(TempDataSet_Log2Expression[1,]), TempDataSet_StressVar, sd)
# 13days_5daysofrest  13days_Biol         8days_Biol         acute_Biol 
# 0.5362083          0.2863431          0.2554044          0.2817382 
#So on average, a data point (sample) in the the 8days_Biol group tends to vary by 0.255 Log2 Expression units from the mean for that group (-0.189)

#When working with two groups we often "pool" the standard deviation of the two groups to come up with a value that is pretty representative for both of them.
#If we have an equal number of samples in each group (in this case, n=6/group), we basically average the standard deviation in the two groups, e.g. for 8days_Biol and acute_Biol:
(0.2554044+0.2817382)/2 
#[1] 0.2685713

#The standardized effect size (Cohen's d) for this group difference (Log2FC) would then be:
#Log2FC/pooled standard deviation = Cohen's d
-0.2575848/0.2685713
#[1] -0.9590928
#So basically the mean Log2 gene expression for the 8day_Biol stress treatment group is almost one full standard deviation away from the mean for the Log2 gene expression for the acute stress group.


#Alright, so we now know how to calculate the group differences (Log2FC) and standardized effect sizes, but how do we know how likely those group differences are likely to arise just due to random chance? (vs. due to a real effect of treatment group)

#One way is to individually compare each of the treatment groups to the reference group using a t-test.
#Within a t-test, it is assumed that the probability of finding a particular magnitude of difference between two groups (e.g., a Log2FC of -0.2575848 between 8day_Biol vs. acute_Biol) just by random chance is a property of:
#1) The amount of variability present in the data. More variability in the data, means more possiblity of seeing a large group difference just by random chance.
#2) The size of the samples. Larger samples will give you a better sense of the true mean Log2 expression for each group, whereas small samples are more likely to show group differences just due to random variability within the groups.

#These two values are combined into a single metric, which is called the standard error (SE). 
#For any particular group, the standard error is the standard deviation for that group divided by the square root of the sample size, e.g., for the 8day_Biol group (n=6):
#SE=sd/sqrt(n)
0.2554044/sqrt(6)
#[1] 0.1042684
#or for the reference (acute stress group):
0.2817382/sqrt(6)
#[1] 0.1150191

#The standard error for the group difference is very similar to adding the standard errors for the two groups being compared (almost)
0.1042684+0.1150191
#[1] 0.2192875
#... but not quite, because it needs to reflect a joint probability. I.e., if you have a mean for the 8day_Biol group that by random chance just happens to be pretty different from the true mean for that group, it is unlikely that the Acute_Biol group will also have a group mean that also just happens to be pretty different from the true mean for that group. 
#So the standard error for the group difference is smaller than the sum of the standard errors for each group.
#The formula for calculating it (if we are assuming that the amount of variability in the two groups isn't equivalent)
#sqrt(((sd1^2)/n1)+((sd2^2)/n2))
sqrt(((0.2554044^2)/6)+((0.2817382^2)/6))
#[1] 0.1552459

#Therefore, if we were to run this experiment over and over again, on average we would find group differences of +/-0.155 Log2 Expression units just by random chance.

#So how much bigger is our group difference (Log2FC) than the amount of variability in group differences that would be expected just due to random chance (SE)?
-0.2575848/0.1552459
#[1] -1.659205
#This is what is called the t-statistic.
#The fact that it is negative represents the fact that the treatment group of interest (8days_Biol) has lower Log2 Expression than the reference group (acute_Biol).

#The probability of finding a group difference that is that much larger (1.65x) than what you would find by random chance depends a little bit on the sample size...sort of. 
#What it actually depends on is the amount of information that we have that is independent of all of our calculations - this is called the degrees of freedom (df)
#Degrees of freedom is a weird concept, so I think examples help.  When running a simple two sample t-test, calculating degrees of freedom is relatively straightforward. For this experiment, we have 12 samples (6 samples per group) and we are running two calculations using that information to run the t-test - we are calculating the group difference (Log2FC) and the standard deviation (or average amount of variability within the groups). 
#12 data points - 2 calculations = 10 df left over.  
#To understand why this matters, if you imagine running a t-test with a sample size of 2 (n=1 in each group), you could calculate the group difference (Log2FC), but you would not be able to calculate the standard deviation within the group (because there is none!): 2-2= 0 df (not enough information!). If you only had a sample size of 1 (n=1 in 1 group and n=0 in another group), you couldn't even calculate the group difference!
#If we are doing a t-test that doesn't assume equal variance (e.g., a Welch's t-test), calculating degrees of freedom is not so straightforward, and you often end up with values that aren't whole numbers (e.g., df=9.9).

#So a t-test calculates the probability of finding a group difference that is that much more extreme (the t-statistic - in this case, 1.65x) than what you would find by random chance, using a probability distribution based on the degrees of freedom for your experiment (in this case, something close to df=10).
#This probability is called a p-value. A low probability (p-value) means that you would be unlikely to encounter a group difference that extreme due to random chance, and is often used as evidence that supports the hypothesis that a group difference exists.  In truth, it is actually a little more complicated than that, because if the probability that a group difference exists is *extremely improbable*, even a low possibility of encountering the group difference that you are observing based on random chance may not be sufficient evidence to convince you. 
#Likewise, if you run a statistical test on many samples, you will find low p-values some of the time just due to random chance (e.g., a p<0.05 will pop up something like 1 out of every 20 tests just by random chance). So when you are running lots and lots of statistical tests (e.g., for every gene in a gene expression dataset with 17,000 rows of data!) you will find lots and lots of low p-values.


#So let's show an example of running a t-test on this data:

#This the data for only the 8 day stress samples for that probe:
Temp_8dayData<-as.numeric(TempDataSet_Log2Expression[1,TempDataSet_StressVar=="8days_Biol"])

#This the data for only the acute stress samples for that probe:
Temp_AcuteData<-as.numeric(TempDataSet_Log2Expression[1,TempDataSet_StressVar=="acute_Biol"])

#In Intro to Statistics you were probably taught about basic two sample T-tests. One of the assumptions of these tests is that the variance in the two groups is very similar (var.equal=TRUE):
t.test(Temp_8dayData, Temp_AcuteData, var.equal=TRUE)

# Two Sample t-test
# 
# data:  Temp_8dayData and Temp_AcuteData
# t = -1.6592, df = 10, p-value = 0.1281
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.60349431  0.08832469
# sample estimates:
#   mean of x   mean of y 
# -0.18953243  0.06805238 

#However, when running a t-test, I tend to *not* assume that the variation in the two groups is similar because I often find that statistical assumption is violated in neuroscience data, so I say var.equal=FALSE and R runs a "Welch" two sample t-test:
t.test(Temp_8dayData, Temp_AcuteData, var.equal=FALSE)

# Welch Two Sample t-test
# 
# data:  Temp_8dayData and Temp_AcuteData
# t = -1.6592, df = 9.9052, p-value = 0.1284
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.60394346  0.08877384
# sample estimates:
#   mean of x   mean of y 
# -0.18953243  0.06805238 

#In both cases, our p-value is 0.1284.  This means that there is a decent probability (>1/10) that we might encounter this group difference just due to random chance, so we would say that our group difference isn't "statistically significant" and our hypothesis that there is a group difference is not supported.
#However, this p-value *should not* be used as evidence that the group difference doesn't exist, because if you think about it, there is a still a decent probability that we *didn't* encounter the group difference due to random chance, right? We just can't be sure.


#Alright, to add a layer of complication, most of our experimental designs are not just comparing two groups (or even just one predictor variable).  In this case, we might include all group comparisons (or all predictor variables) in a single big model predicting the dependent variable (Log2 gene expression).
#If we do this, our predictions tend to be more accurate, because we can use the information about variability within the data from all of the different groups to get a better sense of how much the Log2 expression for the gene might vary just by random chance.

#There are many ways to do this, but most are based on a general method called linear regression. Here's an example of running a linear regression analysis:

#First I need to make the stress variable a factor so that I can set the reference level as acute stress:
TempDataSet_StressVar_AsFactor<-as.factor(TempDataSet_StressVar)
levels(TempDataSet_StressVar_AsFactor)
#[1] "13days_5daysofrest" "13days_Biol"        "8days_Biol"         "acute_Biol"  
#the first "level" of the variable that is listed is the reference level - i.e., the group that all other groups are compared to

#Then I set the reference level to be the acute stress condition:
TempDataSet_StressVar_AsFactor<-relevel(TempDataSet_StressVar_AsFactor, ref="acute_Biol")
levels(TempDataSet_StressVar_AsFactor)
#[1] "acute_Biol"         "13days_5daysofrest" "13days_Biol"        "8days_Biol" 

#The formula used for running the linear regression model is similar to the one used for making a boxplot (y~x):
summary.lm(lm(as.numeric(TempDataSet_Log2Expression[1,])~TempDataSet_StressVar_AsFactor))
# Call:
#   lm(formula = as.numeric(TempDataSet_Log2Expression[1, ]) ~ TempDataSet_StressVar_AsFactor)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.73764 -0.20418 -0.03362  0.21940  0.81776 
# 
# Coefficients:
#                                                   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                       0.06805    0.14636   0.465    0.647
# TempDataSet_StressVar_AsFactor13days_5daysofrest -0.18527    0.20699  -0.895    0.381
# TempDataSet_StressVar_AsFactor13days_Biol        -0.01473    0.20699  -0.071    0.944
# TempDataSet_StressVar_AsFactor8days_Biol         -0.25758    0.20699  -1.244    0.228
# 
# Residual standard error: 0.3585 on 20 degrees of freedom
# Multiple R-squared:  0.1018,	Adjusted R-squared:  -0.03296 
# F-statistic: 0.7554 on 3 and 20 DF,  p-value: 0.5322

#In the output, the "estimate" is the Log2 Fold Change - the difference between the group listed (e.g., 8days_Biol) in that row of the output and the reference group (acute_Biol).
#The Std.Error is the standard error (SE) associated with that Log2FC
#The t value is the t-statistic associated with that Log2FC
#And the Pr(>|t|) is the probability of encountering that t-stat by random chance - i.e., the p-value.
#Note that the SE, t-stat, and p-value are different from what we encountered when running a simple t-test group comparison - that is because this tests uses information about variability from all of the groups to estimate the amount of random variability in the data.


#Gemma uses a larger linear regression model similar to this: 

#When we calculate the effect of each stress condition within a larger linear regression model containing all data, we get a t-statistic as output that matches what Gemma has as output (-1.244) - so it looks like this was their method (see demonstration below after we read in Gemma's output).


#For running our meta-analysis:

#We can use the Log2FC and the standard error (SE) for the group comparisons for our meta-analysis

#... but unfortunately Gemma doesn't provide the SE for the coefficient (Log2FoldChange), only the t-value, so we will have to calculate backwards from the t-statistic:

#As a reminder, the t-statistic in this output is just the estimated group difference (Log2 Fold Change) divided by the Standard Error for this group difference, e.g., for the 8day_Biol vs. acute_Biol group comparison:
#Estimate (Log2FoldChange)=-0.25758
#Standard error (SE)=0.20699 
-0.25758/0.20699 
#[1] -1.244408

#Therefore, if we divide the estiate(Log2FoldChange) by the t-statistic, we get the standard error (SE):
-0.25758/-1.244
#[1] 0.2070579
#very close to the standard error - a little off due to the fact that the output from the statistical model was rounded)

#######################################

#1) Read in the differential expression results from Gemma for each dataset:

#To start out with, I downloaded the Gemma differential expression output for the studies of interest.
#I put each study is in its own folder.

#I set working directory to where the downloaded Gemma differential expression  output is at:
setwd("~/Documents/SideProjects/BrainGMT/Gemma/Hippocampus")

#Then inputted the results from the first study, GSE59070:
#This is an agilent dataset from mice:

setwd("~/Documents/SideProjects/BrainGMT/Gemma/Hippocampus/11627_GSE59070_diffExpAnalysis_94802")
list.files()
#"resultset_ID478782.data.txt" 
#"analysis.results.txt" 

#Reading in the analysis results file. This file is an overview of the differential expression statistical output encoded as tab-delimited text, with additional comments/metadata at the top that is indicated with a #
TempAnalysisResults<-read.delim("analysis.results.txt", sep="\t", stringsAsFactors = FALSE, comment.char = "#")

#Each row represents the results for a particular probe/gene:
str(TempAnalysisResults)
# 'data.frame':	41264 obs. of  6 variables:
# $ Element_Name    : chr  "A_52_P278538" "A_51_P481930" "A_52_P208681" "A_51_P374476" ...
# $ Gene_Symbol     : chr  "Hba-a1|Hba-a2" "Cdh15" "Hba-a1|Hba-a2" "Hbb-b1|Hbb-b2|Hbb-bt|Hbb-bs" ...
# $ Gene_Name       : chr  "hemoglobin alpha, adult chain 1|hemoglobin alpha, adult chain 2" "cadherin 15" "hemoglobin alpha, adult chain 1|hemoglobin alpha, adult chain 2" "hemoglobin, beta adult major chain|hemoglobin, beta adult minor chain|hemoglobin, beta adult t chain|hemoglobin"| __truncated__ ...
# $ NCBI_ID         : chr  "15122|110257" "12555" "15122|110257" "15129|15130|101488143|100503605" ...
# $ QValue_timepoint: num  0.02 0.02 0.02 0.02 0.02 ...
# $ PValue_timepoint: num  1.53e-06 2.07e-06 2.12e-06 2.52e-06 2.74e-06 ...

#I placed these results in a list format, because I'm going to join them with the specific statistical output associated with each predictor variable (found in the result set files)
TempResultsToJoin<-list(TempAnalysisResults)
str(TempResultsToJoin)

#The result set files are the detailed differential expression results for a particular variable (e.g., stress time point):
TempResultsToJoin[[2]]<-read.delim("resultset_ID478782.data.txt", sep="\t", stringsAsFactors = FALSE, comment.char = "#")

str(TempResultsToJoin[[2]])
# 'data.frame':	41264 obs. of  13 variables:
#   $ Element_Name               : chr  "A_52_P266365" "A_52_P762911" "A_52_P200244" "A_52_P1029978" ...
# $ Gene_Symbol                : chr  "Rbm27" "" "Sptlc2" "Lasp1" ...
# $ Gene_Name                  : chr  "RNA binding motif protein 27" "" "serine palmitoyltransferase, long chain base subunit 2" "LIM and SH3 protein 1" ...
# $ NCBI_ID                    : chr  "225432" "" "20773" "16796" ...
# $ FoldChange_13.d.5.d.of.rest: num  0.03501 -0.1365 -0.000978 0.02177 0.09924 ...
# $ Tstat_13.d.5.d.of.rest     : num  0.371 -0.5066 -0.0117 0.1818 0.2917 ...
# $ PValue_13.d.5.d.of.rest    : num  0.715 0.618 0.991 0.858 0.773 ...
# $ FoldChange_8.d             : num  -0.169 -0.09793 0.06481 -0.000861 -0.1216 ...
# $ Tstat_8.d                  : num  -1.791 -0.3635 0.7737 -0.00719 -0.3576 ...
# $ PValue_8.d                 : num  0.0885 0.7201 0.4482 0.9943 0.7244 ...
# $ FoldChange_13.d            : num  -0.07923 -0.1743 -0.01604 -0.00859 0.04578 ...
# $ Tstat_13.d                 : num  -0.8396 -0.647 -0.1915 -0.0717 0.1346 ...
# $ PValue_13.d                : num  0.411 0.525 0.85 0.944 0.894 ...

#Note that even though the analysis.results and result sets are statistical output from the same dataset and thus have the same number of rows of data, you can't just combine the two data.frames together using cbind() because they list the results for each of the probes in a different order

#Also note: 
#This example had only one result set file, but other data sets have multiple (for each predictor variable). 
#You can add other result sets [i] to our list object by reading the result sets into TempResultsToJoin[[i]]

# So... what are these stats?

#The FoldChange in the Gemma output is the same thing as what we were referring to as Log2 Fold Change or Log2FC  - it is the difference in the Log2 gene expression levels for a gene between the two groups (e.g., 8 days of stress vs. acute stress)

#The Tstat is the test statistic (in this case, t-statistic) provided for a particular group comparison (e.g., 8 days of stress vs. acute stress) within the output from a larger linear model that includes the full sample and all of other the relevant variables for the experiment.

#The PValue is the p-value provided for a particular group comparison (e.g., 8 days of stress vs. acute stress) within the output from a large linear model that includes the full sample and all of other the relevant variables for the experiment. A p-value tells you the probability of observing this large of a group difference by random chance if you were to run the experiment over and over again. 
#A small p-value is typically interpreted as support that there may be a real group difference, but note again that in a large dataset we observe a large number of small p-values (false discovery) just because we are running statistical tests thousands of times (for each probe/gene) - i.e., lots of opportunities to observe the full range of values that can be produced by random chance!

#The false discovery rate (FDR) or q-value provides something like a p-value that is corrected for false discovery rate. It gives you a sense of how many more results you have at a particular level of statistical significance (p-value) than you would expect simply due to running a large number of statistical tests. For example, if you have 10,000 genes in your dataset, you would expect that just by random chance you would have at least 10 genes with a p-value<0.001 (10,000*0.001=10). If you actually observe 100 genes with p<0.001, then you have 10x more than you would expect by random chance due to running a statistical test 10,000 times. Those genes would have an FDR or q-value<0.10 (or 1/10).

#For easy of use, I join the results overview and the specific results for my variable of interest together into a single data-frame by aligning them using the Gemma annotation for the individual genes/transcripts/probes (rows):
TempResultsJoined<-join_all(TempResultsToJoin, by="Element_Name")
#Note: I've heard from other students that the join() and join_all() functions in the plyr package can cause problems in newer versions of R - you may need to replace this with merge and merge_all

#For joining: Since these are both results from the same dataset, they should contain the same elements (genes/transcripts/probes) and the joining process should not increase the number of rows in the data frame:
str(TempResultsJoined)
# 'data.frame':	41264 obs. of  18 variables:
#   $ Element_Name               : chr  "A_52_P278538" "A_51_P481930" "A_52_P208681" "A_51_P374476" ...
# $ Gene_Symbol                : chr  "Hba-a1|Hba-a2" "Cdh15" "Hba-a1|Hba-a2" "Hbb-b1|Hbb-b2|Hbb-bt|Hbb-bs" ...
# $ Gene_Name                  : chr  "hemoglobin alpha, adult chain 1|hemoglobin alpha, adult chain 2" "cadherin 15" "hemoglobin alpha, adult chain 1|hemoglobin alpha, adult chain 2" "hemoglobin, beta adult major chain|hemoglobin, beta adult minor chain|hemoglobin, beta adult t chain|hemoglobin"| __truncated__ ...
# $ NCBI_ID                    : chr  "15122|110257" "12555" "15122|110257" "15129|15130|101488143|100503605" ...
# $ QValue_timepoint           : num  0.02 0.02 0.02 0.02 0.02 ...
# $ PValue_timepoint           : num  1.53e-06 2.07e-06 2.12e-06 2.52e-06 2.74e-06 ...
# $ Gene_Symbol                : chr  "Hba-a1|Hba-a2" "Cdh15" "Hba-a1|Hba-a2" "Hbb-b1|Hbb-b2|Hbb-bt|Hbb-bs" ...
# $ Gene_Name                  : chr  "hemoglobin alpha, adult chain 1|hemoglobin alpha, adult chain 2" "cadherin 15" "hemoglobin alpha, adult chain 1|hemoglobin alpha, adult chain 2" "hemoglobin, beta adult major chain|hemoglobin, beta adult minor chain|hemoglobin, beta adult t chain|hemoglobin"| __truncated__ ...
# $ NCBI_ID                    : chr  "15122|110257" "12555" "15122|110257" "15129|15130|101488143|100503605" ...
# $ FoldChange_13.d.5.d.of.rest: num  0.456 -0.137 0.441 0.538 0.387 ...
# $ Tstat_13.d.5.d.of.rest     : num  2.17 -2.09 2.12 2.1 1.89 ...
# $ PValue_13.d.5.d.of.rest    : num  0.0426 0.0499 0.0471 0.0489 0.0733 ...
# $ FoldChange_8.d             : num  1.141 0.355 1.081 1.314 0.974 ...
# $ Tstat_8.d                  : num  5.42 5.42 5.18 5.12 4.76 ...
# $ PValue_8.d                 : num  2.67e-05 2.65e-05 4.50e-05 5.16e-05 1.20e-04 ...
# $ FoldChange_13.d            : num  1.567 -0.0247 1.528 1.859 1.484 ...
# $ Tstat_13.d                 : num  7.44 -0.378 7.329 7.25 7.253 ...
# $ PValue_13.d                : num  3.51e-07 7.10e-01 4.39e-07 5.16e-07 5.12e-07 ...


#Alright, so how do Gemma's results compare to the results that we calculated using the linear regression model above?

#Grabbing t-stat results for the same probe for the Hivep3 gene:
TempResultsJoined$Tstat_8.d[TempResultsJoined$Element_Name=="A_51_P100021"]
#[1] -1.244
#Matches the results above

#Grabbing a few more stats from the Gemma output to show that they match our calculations:
TempResultsJoined$PValue_8.d[TempResultsJoined$Element_Name=="A_51_P100021"]
#[1] 0.2277

TempResultsJoined$FoldChange_8.d[TempResultsJoined$Element_Name=="A_51_P100021"]
#[1] -0.2576


#Also useful for navigating the Gemma output:

#You can view the distribution for any variable using a histogram:
hist(TempResultsJoined$FoldChange_8.d)

#You can see how many results satisfied a traditional statistical cut-off for significance by determining the number of results that have a false discovery rate (FDR) corrected p-value (also called q-value) less than 0.1:
#For this experiment, Gemma calculated these q-values using an ANOVA (similar to a big linear regression) analysis that determined whether there was any effect of stress time point (i.e., using all stress time point groups)

sum(TempResultsJoined$QValue_timepoint<0.10) 
#[1] 46

#And you can take a peek at some of those top differential expression results:
TempResultsJoined$Gene_Symbol[TempResultsJoined$QValue_timepoint<0.10]
# [1] "Hba-a1|Hba-a2"               "Cdh15"                       "Hba-a1|Hba-a2"              
# [4] "Hbb-b1|Hbb-b2|Hbb-bt|Hbb-bs" "Spp1"                        "Hbb-bt|Hbb-b2"              
# [7] "Churc1"                      "Aldh1a1"                     "Hbb-b1|Hbb-b2|Hbb-bt|Hbb-bs"
# [10] "Cxxc4"                       ""                            "Cstad"                      
# [13] "Scx"                         "Banp"                        "Tent5c"                     
# [16] "Akr1b10"                     "Alas2"                       ""                           
# [19] "Vcam1"                       "Klhl3"                       "Hbb-b1|Hbb-bt|Hbb-bs|Hbb-b2"
# [22] "Dbi"                         ""                            "Magel2"                     
# [25] "Thbs4"                       "Htr2c"                       "Sphk2|Dbp"                  
# [28] ""                            "Mt1"                         "S100a8"                     
# [31] "Hc"                          "Cadps2"                      "Fabp7"                      
# [34] "Myh7"                        "Zic1"                        "Hba-ps4"                    
# [37] "S100a9"                      "Sln"                         "Hs3st2"                     
# [40] "Zic4"                        "Mgp"                         "Cdhr4"                      
# [43] ""                            "Vat1l"                       "Zic1"                       
# [46] "Ctsh"  

#So the top results are related to hemoglobin (Hba) - i.e., a major component of red blood cells.
#Is this a real result or an artifact of some sort of (small) dissection variability across the treatment groups?
#This is where meta-analysis comes in handy.

####################

#2) Identifying the results for the variables of interest for the meta-analysis:

#At this point, we should double-check:

#A) Whether Gemma has used a reference (baseline) group that makes sense for our analysis - in this case, the unstressed controls.
######### And it turns out that this experiment *does not include unstressed controls*, only acute stress ("initial timepoint"). This is probably a good reason to exclude this study, but I'm going to pretend that I'm planning on including it so that you can see the rest of the coding involved in processing this dataset.
######### According to Gemma, this reference ("initial timepoint") group included 6 samples, no outliers removed

#B) Which group comparisons are actually of interest to us.
######I am most interested in the chronic stress time points: 8 days (vs. acute stress) and 13 days (vs. acute stress)
######### According to Gemma, each of these groups included 6 samples, no outliers removed



#############################

#3) Remove rows of data that have missing or unambiguous gene annotation

#Notably, using modern probe annotation, many of these probes now map to more than one gene symbol, and when this happens Gemma separates provides a list of all of the gene symbols that the probe maps to using a pipe (e.g., Hba-a1|Hba-a2)
#There are also probes that no longer map to any gene symbol (annotated with ""). To double-check that these probes are genuinely not mapping to an annotated part of the genome (vs. simply a part that doesn't have an official gene symbol) we can also check on the NCBI (Entrez) ID, which should be a little more reliable/stable:

TempResultsJoined$NCBI_ID[TempResultsJoined$QValue_timepoint<0.10]
# [1] "15122|110257"                    "12555"                          
# [3] "15122|110257"                    "15129|15130|101488143|100503605"
# [5] "20750"                           "101488143|15130"                
# [7] "211151"                          "11668"                          
# [9] "15129|15130|101488143|100503605" "319478"                         
# [11] ""                                "78617"                          
# [13] "20289"                           "53325"                          
# [15] "74645"                           "67861"                          
# [17] "11656"                           ""                               
# [19] "22329"                           "100503085"                      
# [21] "15129|101488143|100503605|15130" "13167"                          
# [23] ""                                "27385"                          
# [25] "21828"                           "15560"                          
# [27] "56632|13170"                     ""                               
# [29] "17748"                           "20201"                          
# [31] "15139"                           "320405"                         
# [33] "12140"                           "140781"                         
# [35] "22771"                           "383229"                         
# [37] "20202"                           "66402"                          
# [39] "195646"                          "22774"                          
# [41] "17313"                           "69398"                          
# [43] ""                                "270097"                         
# [45] "22771"                           "13036"   


#Let's see how many of the results do not map unambigously to a single gene:

#To locate all of the rows of data with no NCBI Entrez ID annotation:
str(TempResultsJoined$NCBI_ID=="")
#logi [1:41264] FALSE FALSE FALSE FALSE FALSE FALSE ...

#If you apply the function "sum" to a logical vector, it counts the TRUEs:
sum(TempResultsJoined$NCBI_ID=="")
#[1] 11443
#11,443 that lack NCBI Entrez IDs

#Let's see if we get similar results for gene symbol:
sum(TempResultsJoined$Gene_Symbol=="")
#[1] 11443
#The same. Interesting, and also quite handy.

#To locate all of the rows of data that have multiple NCBI Entrez IDs, I searched for pipes (|).
#Normally, we could just use the grep function to search for the character of interest, but the pipe is a special character (not just alphanumeric), so I had to add two \\ in front of it to tell R to treat the pipe as a character instead of part of the formula.
str(grep('\\|', TempResultsJoined$NCBI_ID))
#int [1:1838] 1 3 4 6 9 21 27 62 100 140 ...
#These integers indicate the rows containing multiple NCBI Entrez IDs separated by pipes

#What if I used gene symbol annotation instead?
str(grep('\\|', TempResultsJoined$Gene_Symbol))
#int [1:1838] 1 3 4 6 9 21 27 62 100 140 ...
#Interesting - identical results again. Also quite handy.

#I would like to throw out the results from any probes that do not unambiguously map to a single gene.

#For the sake of simple coding, I'm going to do this stepwise:

#I only want the subset of data which contains rows that do not contain an NCBI ID of ""
TempResultsJoined_NoNA<-TempResultsJoined[TempResultsJoined$Gene_Symbol!="",]

str(TempResultsJoined_NoNA)
#'data.frame':	29821 obs. of  18 variables:

#... and I also only want the subset of that data which contains the rows that do not contain more than one NCBI ID separated by a pipe:
#Adding the negative symbol before the indices for the rows means "not those rows"
TempResultsJoined_NoNA_NoMultimapped<-TempResultsJoined_NoNA[-(grep('\\|', TempResultsJoined_NoNA$Gene_Symbol)),]

str(TempResultsJoined_NoNA_NoMultimapped)
#'data.frame':	27983 obs. of  18 variables

#Let's do a sanity check and make sure that the rows that are left actually only seem to contain annotation with single Gene Symbol:
#I'm going to peek at the first 100 entries:
TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol[c(1:100)]
# [1] "Cdh15"         "Spp1"          "Churc1"        "Aldh1a1"       "Cxxc4"         "Cstad"        
# [7] "Scx"           "Banp"          "Tent5c"        "Akr1b10"       "Alas2"         "Vcam1"        
# [13] "Klhl3"         "Dbi"           "Magel2"        "Thbs4"         "Htr2c"         "Mt1"          
# [19] "S100a8"        "Hc"            "Cadps2"        "Fabp7"         "Myh7"          "Zic1"         
# [25] "Hba-ps4"       "S100a9"        "Sln"           "Hs3st2"        "Zic4"          "Mgp"          
# [31] "Cdhr4"         "Vat1l"         "Zic1"          "Ctsh"          "Kirrel2"       "Myh7"         
# [37] "Cpne9"         "Nipal2"        "Mt2"           "Asb4"          "Zfp36l1"       "Ebf3"         
# [43] "Avp"           "Dusp11"        "Tac2"          "Zfhx3"         "Col18a1"       "Ccnl2"        
# [49] "Zic3"          "Alpl"          "Gkn3"          "Isg20"         "Nuf2"          "Csrp2"        
# [55] "Tgif2"         "Ly9"           "Cep135"        "Tcf7l2"        "Gjc3"          "Saxo2"        
# [61] "Col2a1"        "Dnah10"        "Gstm1"         "St8sia3"       "Acoxl"         "Sox11"        
# [67] "Tjp2"          "Pcp4l1"        "Pdlim7"        "Nmb"           "Eloa"          "Mitd1"        
# [73] "1810058I24Rik" "Tcf7l2"        "Dleu2"         "Tgfbi"         "Cd93"          "Ciart"        
# [79] "2410022M11Rik" "Tmem51"        "Echdc2"        "Mkks"          "Plscr1"        "Erap1"        
# [85] "Ano1"          "Plscr1"        "Ptpmt1"        "Nid1"          "Igfbp7"        "Mpp7"         
# [91] "Slc7a11"       "Smim11"        "Calb2"         "Nrbp2"         "Tnfrsf25"      "Ramp2"        
# [97] "Enpp4"         "Thbs4"         "Krtcap3"       "Cd59a" 


#For record keeping (sometimes useful for troubleshooting later)
write.csv(TempResultsJoined_NoNA_NoMultimapped, "TempResultsJoined_NoNA_NoMultimapped.csv")

#Another thing worth checking: How many of these rows of data represent probes targeting the *same* gene?

#This code shows all of the unique NCBI Ids in the dataset:
unique(TempResultsJoined_NoNA_NoMultimapped$NCBI_ID)

#If I look at the length of that vector, I can determine how many unique NCBI Ids there are:
length(unique(TempResultsJoined_NoNA_NoMultimapped$NCBI_ID))
#[1] 18503

#So 18503 out of 27983 represent unique NCBI IDs.
#Let's see if those also have the same gene symbols:
length(unique(TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol))
#[1] 18503
#The same as NCBI ID again. Interesting.

#Let's see if those also have the same gemma IDs:
#I'm going to plug in a gene that I suspect has more than one row of data:
TempResultsJoined_NoNA_NoMultimapped[TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol=="Bdnf",]
#       Element_Name Gene_Symbol                         Gene_Name NCBI_ID QValue_timepoint
# 2077  A_52_P384100        Bdnf brain derived neurotrophic factor   12064           0.6537
# 12841 A_51_P261991        Bdnf brain derived neurotrophic factor   12064           0.8828
#So separate element names, but the same gene symbol and NCBI_ID
#I would be curious to know if they target distinct transcripts, or are just Gemma's way of representing different probes.


#To easily move forward, we are likely to need to figure out some way to collapse our results so that there is one result per Gene Symbol. 
#There are many ways to do this. Since we are working with log2FC and t-statistics (both of which are continuous numeric variables that tend to be normally distributed, with both negative and positive values), I think that a simple average is probably sufficient.
#Although note that if we were running the meta-analysis starting from the raw data, we probably would have collapsed the rows representing the same gene prior to differential expression analysis.

#We will need both the Log2FC and T-stats averaged:

TempResultsJoined_NoNA_NoMultimapped_FoldChange_8d_Average<-tapply(TempResultsJoined_NoNA_NoMultimapped$FoldChange_8.d, TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, mean)

TempResultsJoined_NoNA_NoMultimapped_FoldChange_13d_Average<-tapply(TempResultsJoined_NoNA_NoMultimapped$FoldChange_13.d, TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, mean)

TempResultsJoined_NoNA_NoMultimapped_Tstat_8d_Average<-tapply(TempResultsJoined_NoNA_NoMultimapped$Tstat_8.d, TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, mean)

TempResultsJoined_NoNA_NoMultimapped_Tstat_13d_Average<-tapply(TempResultsJoined_NoNA_NoMultimapped$Tstat_13.d, TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, mean)

#Here's what these vectors look like now:
str(TempResultsJoined_NoNA_NoMultimapped_Tstat_8d_Average)
# num [1:18503(1d)] -0.0319 1.6845 -1.311 -0.3 -0.3194 ...
# - attr(*, "dimnames")=List of 1
# ..$ : chr [1:18503] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
#As an aside: note that the order of the rows is now based on "alphabetical order" of the gene symbols (even though some of those symbols are partially numeric). Tapply's output is always in "alphabetical order".

#Comparing the distributions for the t-stats following averaging:
hist(TempResultsJoined_NoNA_NoMultimapped$Tstat_8.d)
hist(TempResultsJoined_NoNA_NoMultimapped_Tstat_8d_Average)
#Averaging makes the distribution for the t-stats a little tighter (i.e., there are fewer extreme values)

#########################

#5) Extract out the information needed for running a meta-analysis: Use the Log2FC and Tstat to calculate the standard error for the Log2FC, and then use the standard error to calculate the sampling variance.

#Putting the averaged Log2FC into the same matrix (rows=Genes, columns=the 8d and 13d stress comparisons)
TempResultsJoined_NoNA_NoMultimapped_FoldChangeAveragedByGeneSymbol<-cbind(TempResultsJoined_NoNA_NoMultimapped_FoldChange_8d_Average, TempResultsJoined_NoNA_NoMultimapped_FoldChange_13d_Average)

str(TempResultsJoined_NoNA_NoMultimapped_FoldChangeAveragedByGeneSymbol)
# num [1:18503, 1:2] -0.00551 0.16795 -0.09224 -0.05196 -0.02993 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:18503] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# ..$ : chr [1:2] "TempResultsJoined_NoNA_NoMultimapped_FoldChange_8d_Average" "TempResultsJoined_NoNA_NoMultimapped_FoldChange_13d_Average"

write.csv(TempResultsJoined_NoNA_NoMultimapped_FoldChangeAveragedByGeneSymbol, "TempResultsJoined_NoNA_NoMultimapped_FoldChangeAveragedByGeneSymbol.csv")

#Now we need the to calculate the standard error (SE) that would normally accompany these fold change value.
#In general, the standard error (in this case, for an estimated group difference) indicates how much you expect that the estimated group difference might vary if you were to run the experiment over and over again.
#In its simplest form, it is calculated using the standard deviation - which is the average amount of variation associated with a variable (e.g., Log2 expression of a particular gene) - divided by the square root of the sample size. 
#On an intuitive level, that means that the larger the sample size for a study, the less you would expect that your estimated group difference might vary if you were to run the experiment over and over again.
#Similarly, the smaller the amount of variation associated with a variable (e.g., gene with very similar Log2 expression in all samples), the less you would expect that your estimated group difference might vary if you were to run the experiment over and over again.

#Unfortunately, Gemma does not provide us with the standard error (SE) that accompanies their estimates of Log2FoldChange, so we need to calculate it using the information that we have.
#The t-statistic is the group difference (Log2FoldChange) divided by the SE
#Therefore, if we divide Log2FoldChange by the t-statistic, we get the SE.

#In R, you can run calculations with entire vectors of values, and it is the same as running the calculation individually for each index of values.
#e.g., 
TempResultsJoined_NoNA_NoMultimapped_FoldChange_8d_Average[1]
# 0610005C13Rik 
# -0.005506 
#this is the FoldChange for this gene for this comparison (group difference)
TempResultsJoined_NoNA_NoMultimapped_Tstat_8d_Average[1]
# 0610005C13Rik 
# -0.03191 
#This is the Tstat for this gene for this comparison
TempResultsJoined_NoNA_NoMultimapped_FoldChange_8d_Average[1]/TempResultsJoined_NoNA_NoMultimapped_Tstat_8d_Average[1]
# 0610005C13Rik 
# 0.1725478 
#This is the SE for this gene for this comparison

#If I run the calculations on the entire vector:
TempResultsJoined_NoNA_NoMultimapped_SE_8d<-TempResultsJoined_NoNA_NoMultimapped_FoldChange_8d_Average/TempResultsJoined_NoNA_NoMultimapped_Tstat_8d_Average

str(TempResultsJoined_NoNA_NoMultimapped_SE_8d)

TempResultsJoined_NoNA_NoMultimapped_SE_8d[1]
# 0610005C13Rik 
# 0.1725478
#I get the same results. :)

TempResultsJoined_NoNA_NoMultimapped_SE_13d<-TempResultsJoined_NoNA_NoMultimapped_FoldChange_13d_Average/TempResultsJoined_NoNA_NoMultimapped_Tstat_13d_Average

str(TempResultsJoined_NoNA_NoMultimapped_SE_13d)

#For running our meta-analysis, we are actually going to need the sampling variance instead of the standard error (also called the sampling standard deviation - i.e. the average deviation of the estimated group difference across different samples)
#The sampling variance is just the standard error squared.

TempResultsJoined_NoNA_NoMultimapped_SE_8d[1]^2
# 0610005C13Rik 
# 0.02977274 

#We can apply this calculation easily to the entire vector too:

TempResultsJoined_NoNA_NoMultimapped_SV_8d<-(TempResultsJoined_NoNA_NoMultimapped_SE_8d)^2

str(TempResultsJoined_NoNA_NoMultimapped_SV_8d)

TempResultsJoined_NoNA_NoMultimapped_SV_8d[1]
# 0610005C13Rik 
# 0.02977274 
# :)

TempResultsJoined_NoNA_NoMultimapped_SV_13d<-(TempResultsJoined_NoNA_NoMultimapped_SE_13d)^2

#I'm going to stash the Log2FC and SV from this study into single data.frames:
#In the process, I gave the output simpler names so that it would be easier to decipher later on, where there is output from many studies.

GSE59070_Log2FC<-cbind.data.frame(GSE59070_Stress_8days_vsAcute=TempResultsJoined_NoNA_NoMultimapped_FoldChange_8d_Average, GSE59070_Stress_13days_vsAcute=TempResultsJoined_NoNA_NoMultimapped_FoldChange_13d_Average)
str(GSE59070_Log2FC)
# 'data.frame':	18503 obs. of  2 variables:
#   $ GSE59070_Stress_8days_vsAcute : num [1:18503(1d)] -0.00551 0.16795 -0.09224 -0.05196 -0.02993 ...
# ..- attr(*, "dimnames")=List of 1
# .. ..$ : chr  "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# $ GSE59070_Stress_13days_vsAcute: num [1:18503(1d)] 0.1983 0.0759 -0.1685 -0.036 -0.0333 ...
# ..- attr(*, "dimnames")=List of 1
# .. ..$ : chr  "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...

GSE59070_SV<-cbind.data.frame(GSE59070_Stress_8days_vsAcute=TempResultsJoined_NoNA_NoMultimapped_SV_8d, GSE59070_Stress_13days_vsAcute=TempResultsJoined_NoNA_NoMultimapped_SV_13d)
str(GSE59070_SV)
# 'data.frame':	18503 obs. of  2 variables:
#   $ GSE59070_Stress_8days_vsAcute : num [1:18503(1d)] 0.02977 0.00994 0.00495 0.03 0.00878 ...
# ..- attr(*, "dimnames")=List of 1
# .. ..$ : chr  "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# $ GSE59070_Stress_13days_vsAcute: num [1:18503(1d)] 0.02979 0.00991 0.00495 0.01624 0.00878 ...
# ..- attr(*, "dimnames")=List of 1
# .. ..$ : chr  "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...

#Before reading in the next dataset, we need to remove all of the temporary variables from our environment so nothing gets confused:

#For documentation purposes:
write.csv(GSE59070_Log2FC, "GSE59070_Log2FC.csv")
write.csv(GSE59070_SV, "GSE59070_SV.csv")

rm(TempResultsToJoin,TempAnalysisResults, TempResultsJoined, TempResultsJoined_NoNA, TempResultsJoined_NoNA_NoMultimapped, TempResultsJoined_NoNA_NoMultimapped_FoldChange_13d_Average, TempResultsJoined_NoNA_NoMultimapped_FoldChange_8d_Average, TempResultsJoined_NoNA_NoMultimapped_FoldChangeAveragedByGeneSymbol, TempResultsJoined_NoNA_NoMultimapped_SE_13d, TempResultsJoined_NoNA_NoMultimapped_SE_8d, TempResultsJoined_NoNA_NoMultimapped_SV_13d, TempResultsJoined_NoNA_NoMultimapped_SV_8d, TempResultsJoined_NoNA_NoMultimapped_Tstat_13d_Average, TempResultsJoined_NoNA_NoMultimapped_Tstat_8d_Average, TempResultsJoined_NoNA_NoMultimapped_TstatsAveragedByGeneSymbol)

###########

#Alright, that code worked, but since I will be reading in a lot of datasets, I would really prefer to stream-line it to decrease the probabilty of error - I'll put that in a new file. :)


