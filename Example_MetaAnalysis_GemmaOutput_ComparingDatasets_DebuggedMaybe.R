#Messing around with meta-analyses of Gemma output
#This file includes the code for comparing the differential expression output from different datasets and runing a meta-analysis (steps #6-10 mentioned below).  
#This is the fourth version of this file that includes one of the HDRF datasets as well (because the previous meta-analysis was being dominated by the suspicious CMS dataset with it's dubiously small sampling variances) as well as DG datasets (because I wanted to play with the code for adding dissection as a predictor to the model)
#At this point, I also started to functionalize things.
#Megan Hagenauer
#07-26-2022
#Associated R workspace: 
#save.image("~/Documents/Teaching/Grinnell Interns 2022/MetaAnalysis_HC_Stress_wHDRF_andDG_Fixed/Workspace_MetaAnalysis_HC_Stress_HDRFandDG_BugFixed.RData")


####################

#Goal:

#For practice, I decided to focus on hippocampal chronic stress datasets that weren't included in Yusra's HDRF meta-analysis

#I piggybacked on Jinglin Xiong's efforts triaging chronic stress datasets, using her spreadsheet: "ChronicStressData(after exclusion).xlsx"

#These were the hippocampal datasets that seemed to meet criteria (chronic adult stress):
#GSE59070 - chronic social stress (CSDS), "Ammon's horn", mouse, Agilent platform
#GSE86392 - chronic restraint stress, HC, rat, RNA-Seq - *has annotation issues*
#GSE116009 - CUMS, HC, mouse, RNA-Seq 
#GSE109315 -CSDS, ventral (V) HC - mouse, RNA-Seq
#GSE151807 - CMS, HC, mice, Affymetrix 430 2.0 - this dataset has a small sample size but suspiciously large t-stats - confounded?

#These are the datasets that I tried to add in to increase my sample size:
#included in our HDRF meta-analysis:
  #GSE81672 - CSDS, HC - Bagot et al.
  #GSE72343 - CSDS, "Ammon's Horn" - Bagot et al. - it turned out that this dataset didn't have Gemma Output that we could use because they had analyzed the data for all brain regions in the same model. We would need to re-run the analysis.

#And these are the DG datasets that I added, to experiment with adding dissection as a predictor variable:
#GSE56028 - unpredictable chronic mild stress (CUMS), Dentate Gyrus (DG)
#GSE132819 - CSDS, DG
#GSE84183 - CUMS, DG

#GSE84185 - CUMS, DG - it turns out this one is actually only anterior cingulate cortex


####################

#Overview of general coding steps:

#0) Reading in & visualizing the Log2 Expression data and Metadata for the individual samples - an illustration of where the differential expression results come from. 

#1) Read in the results from Gemma for each dataset

#2) Identify the results for the variables of interest

#3) Remove rows of data that have missing or unambiguous gene annotation

#### Question: Which gene annotation should we use? I was originally planning to use NCBI ID (because it is a little more stable than gene symbol), but if we use gene symbol we can run a half-way decent meta-analysis of data from both rats and mice without having to add in a step where we decode which genes in rats and mice are orthologous using an orthology database, as many mice and rat genes have the same gene symbol (76%, last time I checked).

#4) Collapse results (average) if there is more than row of data representing a single gene 

#5) Extract out the information needed for running a meta-analysis: Use the Log2FC and Tstat to calculate the standard error for the Log2FC, and then use the standard error to calculate the sampling variance.

#6) Combine together the relevant results from different studies into a single data-frame for the effect sizes (Log2FC) and a single data.frame for the sampling variances. 

#7) Make a correlation matrix to compare the overall results from the different studies. Further visualize the comparison using a hierarchically-clustered heatmap of the correlation matrix. Which studies seem to have similar results? 

#8) Run a meta-analysis using all of the effect sizes for each gene that has data from at least 2 studies. 

#9) Correct the meta-analysis output to take into account the fact that we are running the statistical calculations many times and therefore have a heightened risk of false discovery (false discovery rate correction) 

#10) Determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. 

#####################################

#Code packages used (may require installation & installation of dependencies):

library(plyr)
library(metafor)
library(reshape)
library(multtest)

######################################

#Set working directory:
setwd("~/Documents/Teaching/Grinnell Interns 2022/MetaAnalysis_HC_Stress_wHDRF_andDG_Fixed")

########################################

#6) Combine together the relevant results from different studies into a single data-frame for the effect sizes (Log2FC) and a single data.frame for the sampling variances. 

#The Log2FC values are the first element in the differential expression result object for each study
#The gene symbols are the row.names - 75% of gene symbols for rats and mice are the same, so for simplicity sake, instead of using a gene orthology database to inform equivalency, we're just going to align the datasets by gene symbol.


AligningDatasets<-function(ListOfDEResults){

MetaAnalysis_FoldChange_Dfs<-list()

for(i in c(1:length(ListOfDEResults))){
MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(x=row.names(ListOfDEResults[[i]][[1]]),ListOfDEResults[[i]][[1]], stringsAsFactors=FALSE)
}

print("MetaAnalysis_FoldChange_Dfs:")
print(str(MetaAnalysis_FoldChange_Dfs))
  
MetaAnalysis_FoldChanges<<-join_all(MetaAnalysis_FoldChange_Dfs, by="x", type="full")
#This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)

print("MetaAnalysis_FoldChanges:")
print(str(MetaAnalysis_FoldChanges))

MetaAnalysis_SV_Dfs<-list()

for(i in c(1:length(ListOfDEResults))){
MetaAnalysis_SV_Dfs[[i]]<-data.frame(x=row.names(ListOfDEResults[[i]][[4]]),ListOfDEResults[[i]][[4]], stringsAsFactors=FALSE)
}

print("MetaAnalysis_SV_Dfs:")
print(str(MetaAnalysis_SV_Dfs))

MetaAnalysis_SV<<-join_all(MetaAnalysis_SV_Dfs, by="x", type="full")
#This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)

print("MetaAnalysis_SV:")
print(str(MetaAnalysis_SV))

rm(MetaAnalysis_SV_Dfs, MetaAnalysis_FoldChange_Dfs)
}

#Example Usage;

ListOfDEResults<-list(DEResults_GSE109315, DEResults_GSE116009, DEResults_GSE132819, DEResults_GSE151807, DEResults_GSE56028, DEResults_GSE59070, DEResults_GSE81672, DEResults_GSE84183)

#I just discovered that if I have datasets with *comparisons with the same name*, when I join the datasets some of those comparisons disappear. :(

AligningDatasets(ListOfDEResults)
# [1] "MetaAnalysis_FoldChange_Dfs:"
# List of 8
# $ :'data.frame':	35180 obs. of  3 variables:
#   ..$ x                                  : chr [1:35180] "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# ..$ GSE109315_StressResilient_Vs_Ctrl  : num [1:35180] -0.789 0.221 0.392 -0.725 -0.472 ...
# ..$ GSE109315_StressSusceptible_Vs_Ctrl: num [1:35180] 0.0791 0.3545 -0.279 -0.9989 0.0642 ...
# $ :'data.frame':	18540 obs. of  2 variables:
#   ..$ x                        : chr [1:18540] "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" "0610010F05Rik" ...
# ..$ GSE116009_CUMS_vs_Control: num [1:18540] 0.08169 0.2029 -0.3167 -0.0491 0.00719 ...
# $ :'data.frame':	20240 obs. of  2 variables:
#   ..$ x                     : chr [1:20240] "0610009B22Rik" "0610009L18Rik" "0610010K14Rik" "0610012G03Rik" ...
# ..$ GSE132819_CSDS_vs_Ctrl: num [1:20240] -0.09734 -0.00622 0.0192 -0.00389 -0.01279 ...
# $ :'data.frame':	18337 obs. of  2 variables:
#   ..$ x                    : chr [1:18337] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# ..$ GSE151807_CMS_vs_Ctrl: num [1:18337] -0.2051 -0.1937 -0.2101 -0.0245 -0.3498 ...
# $ :'data.frame':	15617 obs. of  2 variables:
#   ..$ x                   : chr [1:15617] "A1bg" "A1cf" "A2m" "A3galt2" ...
# ..$ GSE56028_CMS_vs_Ctrl: num [1:15617] -0.0819 -0.0412 0.2211 -0.041 0.0327 ...
# $ :'data.frame':	18503 obs. of  3 variables:
#   ..$ x                             : chr [1:18503] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# ..$ GSE59070_Stress8days_vs_Acute : num [1:18503] -0.00551 0.16795 -0.09224 -0.05196 -0.02993 ...
# ..$ GSE59070_Stress13days_vs_Acute: num [1:18503] 0.1983 0.0759 -0.1685 -0.036 -0.0333 ...
# $ :'data.frame':	21766 obs. of  3 variables:
#   ..$ x                                 : chr [1:21766] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ GSE81672_StressResistent_vs_Ctrl  : num [1:21766] 0.155 0.1 0.7492 -0.1075 -0.0212 ...
# ..$ GSE81672_StressSusceptible_vs_Ctrl: num [1:21766] 0.9419 -0.0184 -0.2896 -0.1392 -0.0749 ...
# $ :'data.frame':	17336 obs. of  2 variables:
#   ..$ x                   : chr [1:17336] "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" "0610010K14Rik" ...
# ..$ GSE84183_CMS_vs_Ctrl: num [1:17336] -0.0125 0.01127 -0.02818 -0.02382 -0.00972 ...
# NULL
# [1] "MetaAnalysis_FoldChanges:"
# 'data.frame':	39058 obs. of  12 variables:
#   $ x                                  : chr  "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# $ GSE109315_StressResilient_Vs_Ctrl  : num  -0.789 0.221 0.392 -0.725 -0.472 ...
# $ GSE109315_StressSusceptible_Vs_Ctrl: num  0.0791 0.3545 -0.279 -0.9989 0.0642 ...
# $ GSE116009_CUMS_vs_Control          : num  NA NA 0.0817 0.2029 -0.3167 ...
# $ GSE132819_CSDS_vs_Ctrl             : num  NA NA -0.09734 NA -0.00622 ...
# $ GSE151807_CMS_vs_Ctrl              : num  -0.205 NA -0.194 NA -0.21 ...
# $ GSE56028_CMS_vs_Ctrl               : num  NA NA NA NA NA NA NA NA NA NA ...
# $ GSE59070_Stress8days_vs_Acute      : num  -0.00551 NA 0.16795 NA -0.09224 ...
# $ GSE59070_Stress13days_vs_Acute     : num  0.1983 NA 0.0759 NA -0.1685 ...
# $ GSE81672_StressResistent_vs_Ctrl   : num  0.155 NA 0.1 0.749 -0.107 ...
# $ GSE81672_StressSusceptible_vs_Ctrl : num  0.9419 NA -0.0184 -0.2896 -0.1392 ...
# $ GSE84183_CMS_vs_Ctrl               : num  NA NA -0.0125 0.0113 -0.0282 ...
# NULL
# [1] "MetaAnalysis_SV_Dfs:"
# List of 8
# $ :'data.frame':	35180 obs. of  3 variables:
#   ..$ x                                  : chr [1:35180] "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# ..$ GSE109315_StressResilient_Vs_Ctrl  : num [1:35180] 0.7637 0.1475 0.0319 0.516 0.0515 ...
# ..$ GSE109315_StressSusceptible_Vs_Ctrl: num [1:35180] 0.622 0.12 0.026 0.421 0.042 ...
# $ :'data.frame':	18540 obs. of  2 variables:
#   ..$ x                        : chr [1:18540] "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" "0610010F05Rik" ...
# ..$ GSE116009_CUMS_vs_Control: num [1:18540] 0.0188 0.0706 0.1622 0.0175 0.0235 ...
# $ :'data.frame':	20240 obs. of  2 variables:
#   ..$ x                     : chr [1:20240] "0610009B22Rik" "0610009L18Rik" "0610010K14Rik" "0610012G03Rik" ...
# ..$ GSE132819_CSDS_vs_Ctrl: num [1:20240] 0.014067 0.002855 0.000527 0.008743 0.007016 ...
# $ :'data.frame':	18337 obs. of  2 variables:
#   ..$ x                    : chr [1:18337] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# ..$ GSE151807_CMS_vs_Ctrl: num [1:18337] 0.004157 0.000271 0.007074 0.001403 0.005392 ...
# $ :'data.frame':	15617 obs. of  2 variables:
#   ..$ x                   : chr [1:15617] "A1bg" "A1cf" "A2m" "A3galt2" ...
# ..$ GSE56028_CMS_vs_Ctrl: num [1:15617] 0.01006 0.00868 0.02312 0.00477 0.00764 ...
# $ :'data.frame':	18503 obs. of  3 variables:
#   ..$ x                             : chr [1:18503] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# ..$ GSE59070_Stress8days_vs_Acute : num [1:18503] 0.02977 0.0119 0.00495 0.00969 0.00878 ...
# ..$ GSE59070_Stress13days_vs_Acute: num [1:18503] 0.02979 0.0119 0.00495 0.00969 0.00878 ...
# $ :'data.frame':	21766 obs. of  3 variables:
#   ..$ x                                 : chr [1:21766] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ GSE81672_StressResistent_vs_Ctrl  : num [1:21766] 0.54687 0.00951 0.2498 0.01194 0.00673 ...
# ..$ GSE81672_StressSusceptible_vs_Ctrl: num [1:21766] 0.63824 0.01109 0.29149 0.01394 0.00786 ...
# $ :'data.frame':	17336 obs. of  2 variables:
#   ..$ x                   : chr [1:17336] "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" "0610010K14Rik" ...
# ..$ GSE84183_CMS_vs_Ctrl: num [1:17336] 3.49e-05 3.71e-04 1.65e-04 5.44e-04 5.77e-05 ...
# NULL
# [1] "MetaAnalysis_SV:"
# 'data.frame':	39058 obs. of  12 variables:
#   $ x                                  : chr  "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# $ GSE109315_StressResilient_Vs_Ctrl  : num  0.7637 0.1475 0.0319 0.516 0.0515 ...
# $ GSE109315_StressSusceptible_Vs_Ctrl: num  0.622 0.12 0.026 0.421 0.042 ...
# $ GSE116009_CUMS_vs_Control          : num  NA NA 0.0188 0.0706 0.1622 ...
# $ GSE132819_CSDS_vs_Ctrl             : num  NA NA 0.01407 NA 0.00286 ...
# $ GSE151807_CMS_vs_Ctrl              : num  0.004157 NA 0.000271 NA 0.007074 ...
# $ GSE56028_CMS_vs_Ctrl               : num  NA NA NA NA NA NA NA NA NA NA ...
# $ GSE59070_Stress8days_vs_Acute      : num  0.02977 NA 0.0119 NA 0.00495 ...
# $ GSE59070_Stress13days_vs_Acute     : num  0.02979 NA 0.0119 NA 0.00495 ...
# $ GSE81672_StressResistent_vs_Ctrl   : num  0.54687 NA 0.00951 0.2498 0.01194 ...
# $ GSE81672_StressSusceptible_vs_Ctrl : num  0.6382 NA 0.0111 0.2915 0.0139 ...
# $ GSE84183_CMS_vs_Ctrl               : num  NA NA 3.49e-05 3.71e-04 1.65e-04 ...
# NULL

#I'm going to double-check that the row.names are still in the same order:
cbind(MetaAnalysis_FoldChanges$x, MetaAnalysis_SV$x)[c(1:100),]
#Looks good
sum(MetaAnalysis_FoldChanges$x==MetaAnalysis_SV$x)
#[1] 39058
sum(MetaAnalysis_FoldChanges$x!=MetaAnalysis_SV$x)
#[1] 0
#Looks good.


####################################

#7) Make a correlation matrix to compare the overall results from the different studies. Further visualize the comparison using a hierarchically-clustered heatmap of the correlation matrix. Which studies seem to have similar results? 

#We can generally compare the differential expression associated with different datasets or variable comparisons using a scatterplot and correlation analysis.

#The code for making scatterplots is very similar to the boxplot code that we used earlier. It uses a y~x formula, and can include many of the same parameters (e.g., x and y labels, color)

plot(MetaAnalysis_FoldChanges$GSE109315_StressResilient_Vs_Ctrl~MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Ctrl, ylab="Stress Resilient Log2FC", xlab="Stress Susceptible Log2FC" )

#Within this plot, each data point represents the differential expression results (log2FC) for a particular gene for two different comparisons: stress susceptible vs. no stress (x-axis) and stress resilient vs. no stress (y axis)

#From looking at this plot, we can see that, in general, if a gene shows a positive log2FC for the stress susceptible vs. no stress comparison (i.e., the stress susceptible group has greater log2 expression for that gene than the no stress group), then that gene is also likely to have a positive log2FC for the stress resilient vs. no stress comparison.

#Similarly, if a gene shows a negative log2FC for the stress susceptible vs. no stress comparison (i.e., the stress susceptible group has lower log2 expression for that gene than the no stress group), then that gene is also likely to have a negative log2FC for the stress resilient vs. no stress comparison.

#This means that the differential expression results associated with the stress susceptible and stress resilient comparisons are positively correlated - they show a similar direction of effect.

#We can illustrate this by adding a trendline to the plot:

#We use a linear regression model to calculate the intercept and slope for the linear relationship between the two variables using the y~x formula above:
Trendline<-lm(MetaAnalysis_FoldChanges$GSE109315_StressResilient_Vs_Ctrl~MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Ctrl)
#And then add that line to our scatterplot:
abline(Trendline, col=2, lwd=3)

#If we want to know whether that linear relationship is stronger than we might expect due to random chance:
summary.lm(Trendline)
# Call:
#   lm(formula = MetaAnalysis_FoldChanges$GSE109315_StressResilient_Vs_Ctrl ~ 
#        MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Ctrl)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.2960 -0.2212 -0.0068  0.2179  4.7673 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                  0.008552   0.002710   3.156   0.0016 ** 
#   MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Ctrl 0.739991   0.006558 112.842   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.507 on 35178 degrees of freedom
# (3878 observations deleted due to missingness)
# Multiple R-squared:  0.2658,	Adjusted R-squared:  0.2657 
# F-statistic: 1.273e+04 on 1 and 35178 DF,  p-value: < 2.2e-16

#The estimate for intercept is where the trend line crosses the y-axis
#The estimate for "Stress Susceptible vs. Control" is the slope - how much of an increase in Log2FC you should expect in the "Stress Resilient vs. Control" if you see a one unit increase in Log2FC for "Stress Susceptible vs. Control" 
#The Pr(>|t|) is the p-value for that relationship, in this case it is smaller than R is willing to display (<2e-16)

#If we want to gene a normalized correlation coefficient for this relationship (ranging between -1 to 1, with -1 being a perfect negative correlation and +1 being a perfect positive correlation), we can run a correlation analysis:
#While running this correlation analysis, we have to tell R to ignore any rows of differential expression output that don't have Log2FC for one of our variables (use "pairwise complete observations")
cor(MetaAnalysis_FoldChanges$GSE109315_StressResilient_Vs_Ctrl,MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Ctrl, use="pairwise.complete.obs")
#[1] 0.5155259


#So here's where things get particularly cool: 
#We can actually run a correlation analysis comparing the Log2FC for every set of differential expression results in our meta-analysis using a single line of code.
#This is called a correlation matrix.
cor(as.matrix(MetaAnalysis_FoldChanges[,-1]), use="pairwise.complete.obs")

#                                         GSE109315_StressResilient_Vs_Ctrl GSE109315_StressSusceptible_Vs_Ctrl
# GSE109315_StressResilient_Vs_Ctrl                         1.000000000                        0.5155259229
# GSE109315_StressSusceptible_Vs_Ctrl                       0.515525923                        1.0000000000
# GSE116009_CUMS_vs_Control                                 0.017331296                        0.0148853563
# GSE132819_CSDS_vs_Ctrl                                    0.068887778                        0.0771771259
# GSE151807_CMS_vs_Ctrl                                    -0.015895116                       -0.0030642261
# GSE56028_CMS_vs_Ctrl                                      0.014029455                        0.0013679424
# GSE59070_Stress8days_vs_Acute                             0.006698949                       -0.0242257862
# GSE59070_Stress13days_vs_Acute                            0.007615876                       -0.0537225267
# GSE81672_StressResistent_vs_Ctrl                          0.021964570                       -0.0003280983
# GSE81672_StressSusceptible_vs_Ctrl                       -0.024724161                       -0.0073031582
# GSE84183_CMS_vs_Ctrl                                     -0.040624118                       -0.0101767354
#                                       GSE116009_CUMS_vs_Control GSE132819_CSDS_vs_Ctrl GSE151807_CMS_vs_Ctrl
# GSE109315_StressResilient_Vs_Ctrl                  0.01733130            0.068887778          -0.015895116
# GSE109315_StressSusceptible_Vs_Ctrl                0.01488536            0.077177126          -0.003064226
# GSE116009_CUMS_vs_Control                          1.00000000            0.049956253          -0.023006851
# GSE132819_CSDS_vs_Ctrl                             0.04995625            1.000000000           0.017558967
# GSE151807_CMS_vs_Ctrl                             -0.02300685            0.017558967           1.000000000
# GSE56028_CMS_vs_Ctrl                               0.00473424            0.009987690          -0.001404139
# GSE59070_Stress8days_vs_Acute                     -0.08690626            0.009107949           0.040419069
# GSE59070_Stress13days_vs_Acute                    -0.05294394            0.048642563           0.051656470
# GSE81672_StressResistent_vs_Ctrl                   0.01412631            0.027020489          -0.030827798
# GSE81672_StressSusceptible_vs_Ctrl                 0.04936568            0.003994891          -0.008024615
# GSE84183_CMS_vs_Ctrl                              -0.06381986            0.025978290          -0.020355223
#                               GSE56028_CMS_vs_Ctrl GSE59070_Stress8days_vs_Acute GSE59070_Stress13days_vs_Acute
# GSE109315_StressResilient_Vs_Ctrl            0.014029455                   0.006698949                    0.007615876
# GSE109315_StressSusceptible_Vs_Ctrl          0.001367942                  -0.024225786                   -0.053722527
# GSE116009_CUMS_vs_Control                    0.004734240                  -0.086906260                   -0.052943935
# GSE132819_CSDS_vs_Ctrl                       0.009987690                   0.009107949                    0.048642563
# GSE151807_CMS_vs_Ctrl                       -0.001404139                   0.040419069                    0.051656470
# GSE56028_CMS_vs_Ctrl                         1.000000000                   0.003038779                    0.019724309
# GSE59070_Stress8days_vs_Acute                0.003038779                   1.000000000                    0.512418765
# GSE59070_Stress13days_vs_Acute               0.019724309                   0.512418765                    1.000000000
# GSE81672_StressResistent_vs_Ctrl             0.021994566                   0.054188625                    0.102127025
# GSE81672_StressSusceptible_vs_Ctrl          -0.045349672                  -0.041519866                   -0.024496361
# GSE84183_CMS_vs_Ctrl                         0.013386654                   0.076812633                    0.090190090
#                                             GSE81672_StressResistent_vs_Ctrl GSE81672_StressSusceptible_vs_Ctrl
# GSE109315_StressResilient_Vs_Ctrl                       0.0219645704                       -0.024724161
# GSE109315_StressSusceptible_Vs_Ctrl                    -0.0003280983                       -0.007303158
# GSE116009_CUMS_vs_Control                               0.0141263068                        0.049365676
# GSE132819_CSDS_vs_Ctrl                                  0.0270204894                        0.003994891
# GSE151807_CMS_vs_Ctrl                                  -0.0308277982                       -0.008024615
# GSE56028_CMS_vs_Ctrl                                    0.0219945661                       -0.045349672
# GSE59070_Stress8days_vs_Acute                           0.0541886248                       -0.041519866
# GSE59070_Stress13days_vs_Acute                          0.1021270248                       -0.024496361
# GSE81672_StressResistent_vs_Ctrl                        1.0000000000                        0.467915666
# GSE81672_StressSusceptible_vs_Ctrl                      0.4679156660                        1.000000000
# GSE84183_CMS_vs_Ctrl                                    0.0146845362                       -0.034426991
#                                          GSE84183_CMS_vs_Ctrl
# GSE109315_StressResilient_Vs_Ctrl            -0.04062412
# GSE109315_StressSusceptible_Vs_Ctrl          -0.01017674
# GSE116009_CUMS_vs_Control                    -0.06381986
# GSE132819_CSDS_vs_Ctrl                        0.02597829
# GSE151807_CMS_vs_Ctrl                        -0.02035522
# GSE56028_CMS_vs_Ctrl                          0.01338665
# GSE59070_Stress8days_vs_Acute                 0.07681263
# GSE59070_Stress13days_vs_Acute                0.09019009
# GSE81672_StressResistent_vs_Ctrl              0.01468454
# GSE81672_StressSusceptible_vs_Ctrl           -0.03442699
# GSE84183_CMS_vs_Ctrl                          1.00000000

#I find that these correlation matrices can be a little easier to look at in Excel, so I often output them into a file:
write.csv(cor(as.matrix(MetaAnalysis_FoldChanges[,-1]), use="pairwise.complete.obs"), "CorMatrix_HC_StressComparisons_Log2FC.csv")

#In the output, each cell includes the correlation coefficient reflecting the similarity of the effects for the variable in the row and column.  
#So at the intersection of the row for "StressSusceptible_Vs_Control" and column "StressResilient_Vs_Control" we can see the correlation coefficient that we calculated earlier (0.5155)
#The intersection of a variable with itself creates a coefficient of 1 (i.e., identical)

#Looking at the full matrix, most of the correlation coefficients are very close to 0.
#The only correlation coefficients that are larger, positive numbers are comparisons that come from the same dataset originally. e.g., "StressSusceptible_Vs_Control" and  "StressResilient_Vs_Control"

#Disappointing, but not surprising.  The comparisons that come from the same dataset (e.g.) often reflect comparisons with the same reference group. Therefore, any random variation in the reference group is going to be artificially shared in the differential expression results for both comparisons.


#You can also visualize that correlation matrix using a hierarchically-clustered heatmap:
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[,-1]), use="pairwise.complete.obs"))
#In this heatmap, white/yellow indicates a more positive correlation
#The groups are placed in order by similarity, as determined by hierarchical clustering.
#The lines ("tree branches") on the left and top illustrate that similarity (clustering) using a "dendrogram"


#################################


#8) Run a meta-analysis using all of the effect sizes for each gene that has data from at least 2 studies. 


#We can only run a meta-analysis if there are differential expression results from more than one comparison.
#Since I know that the differential expression results from the same study (dataset) are artificially correlated, I would actually prefer that there are results from more than one dataset.

#How many genes satisfy this criteria?

#This code caculates the number of NAs in each row:
MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges, 1, function(y) sum(is.na(y)))

#I'm going to make a histogram of the results because I'm curious to see how they are distributed
hist(MetaAnalysis_FoldChanges_NAsPerRow)


#Or, since there are a limited number of integer answers (0-11), I could make a table of the results:
table(MetaAnalysis_FoldChanges_NAsPerRow)
# MetaAnalysis_FoldChanges_NAsPerRow
# 0    1    2    3    4    5    6    7    8    9   10 
# 9070 4434 1490 1582 1313 1135 3040 3442 1221 9822 2509 

#So approximately 1/4 of the genes are found in all 11 sets of differential expression results (0 NAs)
#Approximately 1/3 of the genes are only found in 2 or 3 of the sets of differential expression results (9-10 NAs)

#Let's try running a meta-analysis using genes that were found in at least 7 sets of differential expression results
#Since there are 11 sets of differential expression results, that means that the genes that we are including need to have 4 or fewer NAs in their results

#5 NAs is too many

RunBasicMetaAnalysis<-function(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV){

MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges, 1, function(y) sum(is.na(y)))

print("Table of # of NAs per Row (Gene):")
print(table(MetaAnalysis_FoldChanges_NAsPerRow))
  
MetaAnalysis_FoldChanges_ForMeta<<-MetaAnalysis_FoldChanges[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
MetaAnalysis_SV_ForMeta<<-MetaAnalysis_SV[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]

print("MetaAnalysis_FoldChanges_ForMeta:")
print(str(MetaAnalysis_FoldChanges_ForMeta))

#I'm going to make an empty matrix to store the results of my meta-analysis:
metaOutput<-matrix(NA, length(MetaAnalysis_FoldChanges_ForMeta$x), 6)

#And then run a loop that run's a meta-analysis on the differential expression results (columns 2-10) for each gene (row):
for(i in c(1:length(MetaAnalysis_FoldChanges_ForMeta$x))){
 
  effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-1])
  var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-1])
  
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
row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta$x

metaOutput<<-metaOutput
return(metaOutput)

print("metaOutput:")
print(str(metaOutput))

print("Top of metaOutput:")
print(head(metaOutput))

print("Bottom of metaOutput")
print(tail(metaOutput))

}

#Example Usage:
NumberOfComparisons=11
CutOffForNAs=5
#I want at least 7 datasets (so 5 na's is too many)

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

# [1] "Table of # of NAs per Row (Gene):"
# MetaAnalysis_FoldChanges_NAsPerRow
# 0    1    2    3    4    5    6    7    8    9   10 
# 9070 4434 1490 1582 1313 1135 3040 3442 1221 9822 2509 
# [1] "MetaAnalysis_FoldChanges_ForMeta:"
# 'data.frame':	17889 obs. of  12 variables:
#   $ x                                  : chr  "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# $ GSE109315_StressResilient_Vs_Ctrl  : num  -0.7895 0.3921 -0.4718 -0.0136 -0.188 ...
# $ GSE109315_StressSusceptible_Vs_Ctrl: num  0.0791 -0.279 0.0642 -0.1109 0.0157 ...
# $ GSE116009_CUMS_vs_Control          : num  NA 0.08169 -0.3167 -0.0491 0.00719 ...
# $ GSE132819_CSDS_vs_Ctrl             : num  NA -0.09734 -0.00622 NA 0.0192 ...
# $ GSE151807_CMS_vs_Ctrl              : num  -0.2051 -0.1937 -0.2101 -0.0245 -0.3498 ...
# $ GSE56028_CMS_vs_Ctrl               : num  NA NA NA NA NA NA NA NA NA NA ...
# $ GSE59070_Stress8days_vs_Acute      : num  -0.00551 0.16795 -0.09224 -0.05196 -0.02993 ...
# $ GSE59070_Stress13days_vs_Acute     : num  0.1983 0.0759 -0.1685 -0.036 -0.0333 ...
# $ GSE81672_StressResistent_vs_Ctrl   : num  0.155 0.1 -0.1075 -0.0212 -0.0457 ...
# $ GSE81672_StressSusceptible_vs_Ctrl : num  0.9419 -0.0184 -0.1392 -0.0749 -0.0587 ...
# $ GSE84183_CMS_vs_Ctrl               : num  NA -0.0125 -0.0282 NA -0.0238 ...
# NULL
# There were 50 or more warnings (use warnings() to see the first 50)

str(metaOutput)
# num [1:17889, 1:6] -0.03251 -0.000632 -0.088493 -0.039134 -0.065414 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:17889] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# ..$ : chr [1:6] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

head(metaOutput)
#                 Log2FC_estimate          SE        pval       CI_lb        CI_ub Number_Of_Comparisons
# 0610005C13Rik   -0.0325103784 0.116647983 0.780471211 -0.26113622  0.196115468                     7
# 0610009B22Rik   -0.0006324904 0.050853168 0.990076502 -0.10030287  0.099037888                    10
# 0610009L18Rik   -0.0884931270 0.032116119 0.005861834 -0.15143956 -0.025546691                    10
# 0610010F05Rik   -0.0391341232 0.025806845 0.129412132 -0.08971461  0.011446363                     8
# 0610010K14Rik   -0.0654136223 0.039271822 0.095780550 -0.14238498  0.011557735                    10
# 0610025J13Rik   -0.0099393092 0.007539784 0.187420449 -0.02471701  0.004838395                     9

tail(metaOutput)
#            Log2FC_estimate         SE      pval        CI_lb       CI_ub Number_Of_Comparisons
# Macf1     -0.002941478 0.01738809 0.8656659 -0.037021502 0.031138546                     9
# Malat1     0.078979811 0.12193120 0.5171526 -0.160000946 0.317960569                     8
# Map2       0.015123787 0.03301830 0.6469219 -0.049590893 0.079838467                     9
# Pcdh9      0.012891547 0.01753582 0.4622450 -0.021478026 0.047261121                     9
# Sptan1     0.021281603 0.01526441 0.1632579 -0.008636086 0.051199293                     9
# Syne1     -0.030869644 0.02062401 0.1344501 -0.071291964 0.009552675                     8

########################################

## Multiple Comparison corrections
#The following code applies two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli) 
#Meta-analysis output with adjusted p-values is then outputted along with effect size information.

#9) Correct the meta-analysis output to take into account the fact that we are running the statistical calculations many times and therefore have a heightened risk of false discovery (false discovery rate correction) 

colnames(metaOutput)
# [1] "Log2FC_estimate"       "SE"                    "pval"                  "CI_lb"                
# [5] "CI_ub"                 "Number_Of_Comparisons"
#The p-values are in column 3

#I'm going to grab just the column with the pvalues and run a multiple comparisons correction using the Benjamini-Hochberg method ("FDR" or "q-value")
tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,3], proc=c("BH"))

#unfortunately, the output from that function re-orders the FDR-corrected p-values in order of "significance"
#We would like the FDR-corrected p-values to be in the their original order (i. the order of the rest of our statistical output!). This order is recorded in the index (row numbers) for the p-values:
metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]

str(metaPvalAdj)

#adjusted pvalue object is in same orientation as metaOutput so can simply be binded together
#We can double-check that it is in the same order by comparing the uncorrected p-values in the metaPvalAdj object to the pvalues in the metaOutput matrix

plot(metaPvalAdj[,1]~metaOutput[,3])
#Perfectly straight line - they are exactly the same! *phew*

metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])

str(metaOutputFDR)
# num [1:18340, 1:7] -0.03251 0.00799 -0.11333 -0.03398 -0.07757 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:18340] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

#Let's functionalize it!
FalseDiscoveryCorrection<-function(metaOutput){
  
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,3], proc=c("BH"))
  
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]

  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  colnames(metaOutputFDR)[7]<-"FDR"
  
  metaOutputFDR<<-metaOutputFDR
  
  print("metaOutputFDR:")
  print(str(metaOutputFDR))

  write.csv(metaOutputFDR, "metaOutputFDR.csv")
 
  #a version of the output in order by p-value:
  metaOutputFDR_OrderbyPval<<-metaOutputFDR[order(metaOutputFDR[,3]),]
  
  #Let's write out a version of the output in order by p-value:
  write.csv(metaOutputFDR_OrderbyPval, "metaOutputFDR_orderedByPval_wHDRFData.csv")
  
  print("Do we have any genes that are statistically significant following false discovery rate correction?")
  print(sum(metaOutputFDR[,7]<0.10, na.rm=TRUE))
  
  print("What are the top results?")
  print(head(metaOutputFDR[order(metaOutputFDR[,3]),]))
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
}
  
#Example usage:

FalseDiscoveryCorrection(metaOutput)
# [1] "metaOutputFDR:"
# num [1:17889, 1:7] -0.03251 -0.000632 -0.088493 -0.039134 -0.065414 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:17889] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# NULL
# [1] "Do we have any genes that are statistically significant following false discovery rate correction?"
# [1] 680
# [1] "What are the top results?"
# Log2FC_estimate          SE         pval       CI_lb       CI_ub Number_Of_Comparisons          FDR
# Slc7a12      0.07674825 0.005395505 6.458292e-46  0.06617325  0.08732324                     7 1.155324e-41
# Rpl35       -0.02592991 0.001866834 7.306435e-44 -0.02958884 -0.02227098                     7 6.535241e-40
# Zfp277      -0.07803360 0.005819842 5.413801e-41 -0.08944028 -0.06662692                    10 3.228250e-37
# Klhdc7a      0.14554013 0.011157316 6.842285e-39  0.12367219  0.16740807                    10 3.060041e-35
# Erich6       0.11640096 0.009066277 9.928527e-38  0.09863138  0.13417053                    11 3.552229e-34
# Lrp5         0.08882563 0.008249790 4.927053e-27  0.07265634  0.10499492                    11 1.469001e-23

############################

#10) Determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. 

row.names(metaOutputFDR_OrderbyPval)[c(1:100)]
# [1] "Slc7a12"       "Rpl35"         "Zfp277"        "Klhdc7a"       "Erich6"        "Lrp5"          "Ppil1"        
# [8] "Fbxo39"        "Gnpnat1"       "Ncoa5"         "Utp4"          "Pdcd2l"        "Kcnj13"        "Zrsr1"        
# [15] "Pla2g12a"      "Anapc16"       "Triobp"        "Dimt1"         "Apobec3"       "Clip1"         "Gabarap"      
# [22] "Ppcs"          "Gm11190"       "Afg3l2"        "1700016H13Rik" "Arglu1"        "Zfp956"        "Exosc8"       
# [29] "4933440N22Rik" "Mettl7a3"      "Nlgn2"         "Kctd8"         "Mt3"           "Syap1"         "Pno1"         
# [36] "Psmb5"         "Slc25a14"      "Selenow"       "4930430F08Rik" "Galnt1"        "Zmat5"         "Usp38"        
# [43] "Cwc15"         "Adrm1"         "Endog"         "Cd209f"        "Arl9"          "Grk2"          "S100a13"      
# [50] "Lyz1"          "Smdt1"         "Mrpl53"        "Atp2b3"        "Dusp12"        "Haus7"         "Slc48a1"      
# [57] "Mirg"          "Slc22a15"      "Frzb"          "Gm11681"       "Lmtk3"         "Myo6"          "Mrpl34"       
# [64] "Fth1"          "2510039O18Rik" "Hmces"         "Uqcrq"         "Snw1"          "Snx15"         "Usp12"        
# [71] "Lao1"          "Lig1"          "Kbtbd4"        "4933412E12Rik" "Pex6"          "Carf"          "Dock4"        
# [78] "Pcca"          "Pkhd1"         "Ywhaq"         "Atp6v1a"       "Tpsab1"        "Rps15"         "Sqor"         
# [85] "Git1"          "Pip5kl1"       "Parp16"        "Abcb4"         "Nmi"           "Zfp941"        "4632428C04Rik"
# [92] "Ecd"           "Prrg4"         "Pcgf1"         "Ypel3"         "2310009B15Rik" "Fam102b"       "Erbin"        
# [99] "Rfc5"          "Csdc2"     

#Let's plot some of those top results!

#Quickly looking at the range of Log2FC values to figure out the limits for the x-axis for the forest plots:
hist(metaOutputFDR[,1], breaks=40)
#Range is mostly -0.5 to 0.5, but there are a few with Log2FC as big as -1 or 1.5

MakeForestPlots<-function(GeneSymbol){

pdf(paste("ForestPlot_", GeneSymbol, ".pdf", sep=""), height=5, width=8)

  effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
  var<-as.numeric(MetaAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
  
 forest.rma(rma(effect, var),slab=colnames(MetaAnalysis_FoldChanges_ForMeta)[-1],  xlim=c(-3, 3))

mtext(paste(GeneSymbol), line=-1.5, cex=2)
dev.off()
}


MakeForestPlots("Slc7a12") 
MakeForestPlots("Rpl35") #completely unbelievable - still driven by that one CMS dataset.
MakeForestPlots("Zfp277")
MakeForestPlots("Klhdc7a")
MakeForestPlots("Erich6")
MakeForestPlots("Lrp5")
MakeForestPlots("Ppil1") #there seems to be a pattern that the 4 comparisons/datasets near the top tend to show different results than the ones near the bottom (DG vs. HC?)
MakeForestPlots("Fbxo39")
MakeForestPlots("Gnpnat1")

#The fact that I added more datasets/comparisons seems to have made it so that the top results are more believable 
#Maybe the debugging fo the SE calculations helped too.

#However, despite all of the debugging, that one CMS dataset (GSE84183_CMS_vs_Ctrl) still has ridiculously small confidence intervals for a dataset with n=3/group :(
#... which I guess makes sense, because it also had bizarrely large t-statistics
#... but I was still hoping that correcting the SE calculation would fix some of it.
# and because the confidence intervals are so crazy small, they have a really disproportionate impact on the results.

#In general, it seems to me like the results that have both a significant FDR & large effect size tend to be the most convincing 

#Here's a summary of the distribution for the Log2FCs
summary(metaOutputFDR[,1])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -0.85640 -0.02916 -0.00116  0.00319  0.03222  1.18583      207 

#Let's see how many results are both statistically significant (using a more conservative FDR) and have a larger estimated Log2FC (>0.10 or <-0.10)
sum(metaOutputFDR[,7]<0.05 & abs(metaOutputFDR[,1])>0.1, na.rm=TRUE)
#[1] 99

#What are their gene symbols?
row.names(metaOutputFDR_OrderbyPval)[metaOutputFDR_OrderbyPval[,7]<0.05 & abs(metaOutputFDR_OrderbyPval[,1])>0.1 & is.na(metaOutputFDR_OrderbyPval[,1])==FALSE]
# [1] "Klhdc7a"       "Erich6"        "Fbxo39"        "Kcnj13"        "Apobec3"       "4933440N22Rik" "Mettl7a3"     
# [8] "Cd209f"        "Arl9"          "Lyz1"          "Mirg"          "Frzb"          "Gm11681"       "Myo6"         
# [15] "Lao1"          "Pkhd1"         "Pip5kl1"       "Erbin"         "Klrg1"         "Oas1a"         "Cbln1"        
# [22] "Zfand4"        "Prlr"          "Rgs20"         "Kcnt1"         "Baiap2l1"      "Rxfp3"         "Tspear"       
# [29] "Tomm6"         "Gk2"           "C1rl"          "Fap"           "Agpat2"        "Cyp2t4"        "Cd74"         
# [36] "Rsph9"         "Slc44a5"       "Csf2rb2"       "Kcnq5"         "Clec4a1"       "Gdf10"         "Popdc3"       
# [43] "Chrna2"        "Slc38a5"       "Sox3"          "Bbs12"         "Ctsh"          "Fat4"          "Rpl39l"       
# [50] "Spn"           "Zfp324"        "Ace"           "Clec2d"        "Serpina11"     "Stac"          "Clec10a"      
# [57] "Slc13a1"       "Ybx2"          "Pdlim3"        "S100a8"        "Tes"           "Gpnmb"         "Polr1b"       
# [64] "Cd83"          "Zfp474"        "Fndc9"         "D130040H23Rik" "Ccl12"         "Milr1"         "Mx1"          
# [71] "Tet1"          "Optc"          "Ptgfrn"        "Epha6"         "Dpp10"         "Filip1"        "S100a11"      
# [78] "Prg4"          "Car4"          "Tpbgl"         "Ube3a"         "Dclre1c"       "Tmem267"       "Stk26"        
# [85] "Folr2"         "1810014B01Rik" "Tbata"         "Dhx58"         "Manea"         "Emp3"          "Zpld1"        
# [92] "Hcar2"         "Vwf"           "Banp"          "Tspan18"       "Glipr1"        "Pcolce2"       "C230035I16Rik"
# [99] "Cybb" 

MakeForestPlots("Kcnj13")
#Interesting - for this gene, two of the CMS datasets have really small confidence intervals (GSE151807_CMS_vs_Ctrl & GSE84183_CMS_vs_Ctrl)
MakeForestPlots("Apobec3")
MakeForestPlots("4933440N22Rik") #This one also has small confidence intervals for both CMS datasets
MakeForestPlots("Mettl7a3") 

#Instead of effect sizes, that if I focus on genes represented in at least 10 datasets?

colnames(metaOutputFDR)

sum(metaOutputFDR[,7]<0.05 & metaOutputFDR[,6]>9, na.rm=TRUE)
#[1] 427

row.names(metaOutputFDR_OrderbyPval)[metaOutputFDR_OrderbyPval[,7]<0.05 & metaOutputFDR_OrderbyPval[,6]>9 & is.na(metaOutputFDR_OrderbyPval[,1])==FALSE]

#Adding in a p-value threshold too:
sum(metaOutputFDR[,7]<0.05 & metaOutputFDR[,6]>9 & metaOutputFDR[,3]<0.00001, na.rm=TRUE)
#[1] 148

row.names(metaOutputFDR_OrderbyPval)[metaOutputFDR_OrderbyPval[,7]<0.05 & metaOutputFDR_OrderbyPval[,6]>9 & metaOutputFDR_OrderbyPval[,3]<0.00001 & is.na(metaOutputFDR_OrderbyPval[,1])==FALSE]
# [1] "Zfp277"        "Klhdc7a"       "Erich6"        "Lrp5"          "Ppil1"         "Ncoa5"         "Utp4"         
# [8] "Pdcd2l"        "Zrsr1"         "Pla2g12a"      "Anapc16"       "Triobp"        "Dimt1"         "Apobec3"      
# [15] "Clip1"         "Gabarap"       "Ppcs"          "Afg3l2"        "Arglu1"        "Zfp956"        "Exosc8"       
# [22] "4933440N22Rik" "Nlgn2"         "Kctd8"         "Mt3"           "Syap1"         "Pno1"          "Psmb5"        
# [29] "Slc25a14"      "Selenow"       "4930430F08Rik" "Galnt1"        "Zmat5"         "Usp38"         "Cwc15"        
# [36] "Adrm1"         "Endog"         "Grk2"          "S100a13"       "Smdt1"         "Mrpl53"        "Atp2b3"       
# [43] "Dusp12"        "Slc48a1"       "Slc22a15"      "Frzb"          "Myo6"          "Mrpl34"        "2510039O18Rik"
# [50] "Hmces"         "Uqcrq"         "Snw1"          "Snx15"         "Usp12"         "Lao1"          "Lig1"         
# [57] "Kbtbd4"        "Pex6"          "Carf"          "Dock4"         "Pcca"          "Ywhaq"         "Atp6v1a"      
# [64] "Rps15"         "Sqor"          "Git1"          "Pip5kl1"       "Parp16"        "Abcb4"         "Nmi"          
# [71] "Zfp941"        "Ecd"           "Prrg4"         "Pcgf1"         "Ypel3"         "Fam102b"       "Rfc5"         
# [78] "Csdc2"         "Rsl1d1"        "Rbfox1"        "Zer1"          "Rpf2"          "Klrg1"         "Mrto4"        
# [85] "Dscaml1"       "Trmt11"        "Gli3"          "Commd2"        "Ctdsp1"        "Fez1"          "Cbln1"        
# [92] "Vxn"           "Ube4b"         "Ift52"         "Amn1"          "Zfand4"        "Tsc1"          "Tmem63c"      
# [99] "Ppp1r37"       "Ccdc92"        "Piezo1"        "Rptor"         "Uqcrfs1"       "Emc6"          "Prlr"         
# [106] "Ddt"           "Kcnd3"         "Mab21l1"       "Ttc23"         "Rnf128"        "Vps29"         "Mas1"         
# [113] "Ywhag"         "Mpv17l2"       "Pld3"          "Ralyl"         "Nat10"         "Phtf1"         "Lmo7"         
# [120] "Pih1d1"        "Rgs20"         "Deaf1"         "Plekho1"       "Camsap3"       "Pofut1"        "Lsg1"         
# [127] "Dennd2a"       "Rhpn2"         "Kcnt1"         "Aqp11"         "Mterf4"        "Mex3a"         "Rgs7"         
# [134] "2700097O09Rik" "Fbrsl1"        "B9d1"          "Baiap2l1"      "Eda"           "Rxfp3"         "Stat5b"       
# [141] "Ccdc96"        "Tipin"         "Nifk"          "Zfp830"        "Mfsd13a"       "Ctsb"          "Arpc4"        
# [148] "Psma7"    



#What if we look at the correlation between datasets again, but only using the top genes found in all 11 datasets?
sum(metaOutputFDR[,7]<0.05 & metaOutputFDR[,6]>10, na.rm=TRUE)
#[1] 326

TopGenesInAll11Datasets<-row.names(metaOutputFDR)[metaOutputFDR[,7]<0.05 & metaOutputFDR[,6]>10 & is.na(metaOutputFDR[,6])==FALSE]

cor(as.matrix(MetaAnalysis_FoldChanges[which(MetaAnalysis_FoldChanges$x%in%TopGenesInAll11Datasets),-1]), use="pairwise.complete.obs")

#                                      GSE109315_StressResilient_Vs_Ctrl GSE109315_StressSusceptible_Vs_Ctrl
# GSE109315_StressResilient_Vs_Ctrl                         1.000000000                          0.77821766
# GSE109315_StressSusceptible_Vs_Ctrl                       0.778217657                          1.00000000
# GSE116009_CUMS_vs_Control                                 0.005140818                          0.09840061
# GSE132819_CSDS_vs_Ctrl                                    0.221374288                          0.32121161
# GSE151807_CMS_vs_Ctrl                                     0.173286428                          0.29989888
# GSE56028_CMS_vs_Ctrl                                      0.098103639                          0.11649547
# GSE59070_Stress8days_vs_Acute                             0.193479354                          0.21682172
# GSE59070_Stress13days_vs_Acute                            0.129568117                          0.12236028
# GSE81672_StressResistent_vs_Ctrl                          0.095424170                          0.14580608
# GSE81672_StressSusceptible_vs_Ctrl                        0.070465889                          0.17417176
# GSE84183_CMS_vs_Ctrl                                      0.255200419                          0.28516118
#                                       GSE116009_CUMS_vs_Control GSE132819_CSDS_vs_Ctrl GSE151807_CMS_vs_Ctrl
# GSE109315_StressResilient_Vs_Ctrl                0.0051408179             0.22137429            0.17328643
# GSE109315_StressSusceptible_Vs_Ctrl              0.0984006109             0.32121161            0.29989888
# GSE116009_CUMS_vs_Control                        1.0000000000             0.05707055           -0.02625782
# GSE132819_CSDS_vs_Ctrl                           0.0570705508             1.00000000            0.50950802
# GSE151807_CMS_vs_Ctrl                           -0.0262578247             0.50950802            1.00000000
# GSE56028_CMS_vs_Ctrl                             0.0591893210             0.33278465            0.19532472
# GSE59070_Stress8days_vs_Acute                    0.0010372019             0.32889175            0.39244423
# GSE59070_Stress13days_vs_Acute                   0.0499700996             0.38377007            0.48580557
# GSE81672_StressResistent_vs_Ctrl                 0.1803160449             0.30287598            0.33111125
# GSE81672_StressSusceptible_vs_Ctrl               0.0003110567             0.33813714            0.51098067
# GSE84183_CMS_vs_Ctrl                             0.0358286150             0.50015921            0.46356174
#                                 GSE56028_CMS_vs_Ctrl GSE59070_Stress8days_vs_Acute GSE59070_Stress13days_vs_Acute
# GSE109315_StressResilient_Vs_Ctrl             0.09810364                   0.193479354                      0.1295681
# GSE109315_StressSusceptible_Vs_Ctrl           0.11649547                   0.216821718                      0.1223603
# GSE116009_CUMS_vs_Control                     0.05918932                   0.001037202                      0.0499701
# GSE132819_CSDS_vs_Ctrl                        0.33278465                   0.328891751                      0.3837701
# GSE151807_CMS_vs_Ctrl                         0.19532472                   0.392444233                      0.4858056
# GSE56028_CMS_vs_Ctrl                          1.00000000                   0.187509944                      0.2065365
# GSE59070_Stress8days_vs_Acute                 0.18750994                   1.000000000                      0.7132741
# GSE59070_Stress13days_vs_Acute                0.20653652                   0.713274052                      1.0000000
# GSE81672_StressResistent_vs_Ctrl              0.16876515                   0.281650585                      0.3191741
# GSE81672_StressSusceptible_vs_Ctrl           -0.03534264                   0.189175477                      0.2284339
# GSE84183_CMS_vs_Ctrl                          0.55996659                   0.391838880                      0.4420614
#                                       GSE81672_StressResistent_vs_Ctrl GSE81672_StressSusceptible_vs_Ctrl
# GSE109315_StressResilient_Vs_Ctrl                         0.09542417                       0.0704658890
# GSE109315_StressSusceptible_Vs_Ctrl                       0.14580608                       0.1741717597
# GSE116009_CUMS_vs_Control                                 0.18031604                       0.0003110567
# GSE132819_CSDS_vs_Ctrl                                    0.30287598                       0.3381371429
# GSE151807_CMS_vs_Ctrl                                     0.33111125                       0.5109806699
# GSE56028_CMS_vs_Ctrl                                      0.16876515                      -0.0353426372
# GSE59070_Stress8days_vs_Acute                             0.28165059                       0.1891754770
# GSE59070_Stress13days_vs_Acute                            0.31917406                       0.2284339462
# GSE81672_StressResistent_vs_Ctrl                          1.00000000                       0.4450758481
# GSE81672_StressSusceptible_vs_Ctrl                        0.44507585                       1.0000000000
# GSE84183_CMS_vs_Ctrl                                      0.29184197                       0.0983400562
#                                     GSE84183_CMS_vs_Ctrl
# GSE109315_StressResilient_Vs_Ctrl             0.25520042
# GSE109315_StressSusceptible_Vs_Ctrl           0.28516118
# GSE116009_CUMS_vs_Control                     0.03582861
# GSE132819_CSDS_vs_Ctrl                        0.50015921
# GSE151807_CMS_vs_Ctrl                         0.46356174
# GSE56028_CMS_vs_Ctrl                          0.55996659
# GSE59070_Stress8days_vs_Acute                 0.39183888
# GSE59070_Stress13days_vs_Acute                0.44206140
# GSE81672_StressResistent_vs_Ctrl              0.29184197
# GSE81672_StressSusceptible_vs_Ctrl            0.09834006
# GSE84183_CMS_vs_Ctrl                          1.00000000

pdf("Heatmap_CorMatrix_HC_StressDatasets_GenesFDR05inAll.pdf", height=10, width=10)
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[which(MetaAnalysis_FoldChanges$x%in%TopGenesInAll11Datasets),-1]), use="pairwise.complete.obs"))
dev.off()

CorMatrixTopGenes<-cor(as.matrix(MetaAnalysis_FoldChanges[which(MetaAnalysis_FoldChanges$x%in%TopGenesInAll11Datasets),-1]), use="pairwise.complete.obs")

write.csv(CorMatrixTopGenes, "CorMatrix_HC_StressDatasets_GenesFDR05inAll.csv")


#################################

#I haven't updated this section yet (following debugging the SE calculations on 7/29):

#Messing around with heterogeneity statistics:

#Let's compare the results for a handful of genes that visually appear to have:

#extremely heterogeneous results (but still "sig" in meta-analysis): "S100a8","Gbx2", "Nkx2-1"
#homogenous results: "Fezf2", "Nmb", "Chrm2", "Npy", "Tet1"
#Results that are mostly heterogeneous/biased due to the CMS dataset: "Hectd2os", "Rps4l", "Vim"

# A data point that has a large value for Cook’s Distance indicates that it strongly influences the fitted values. A general rule of thumb is that any point with a Cook’s Distance over 4/n (where n is the total number of data points) is considered to be an outlier.
# https://www.statology.org/how-to-identify-influential-data-points-using-cooks-distance/

#Do we use the number of studies as n or the number of samples represented by those studies?
4/8
#[1] 0.5
4/6
#[1] 0.6666667
4/4
#[1] 1

# In general, large values of DFBETAS indicate observations that are influential in estimating a given parameter. Belsley, Kuh, and Welsch recommend 2 as a general cutoff value to indicate influential observations and 2/sqrt{n} as a size-adjusted cutoff.
# https://www.sfu.ca/sasdoc/sashtml/stat/chap55/sect38.htm

#Do we use the number of studies as n or the number of samples represented by those studies?
2/sqrt(8)
#[1] 0.7071068
2/sqrt(6)
#[1] 0.8164966
2/sqrt(4)
#[1] 1

# colnames(MetaAnalysis_FoldChanges_ForMeta )
# [1] "x"                            "StressResilient_Vs_Control"   "StressSusceptible_Vs_Control"
# [4] "CUMS_vs_Control"              "CMS_vs_Ctrl"                  "Stress8days_vs_Acute"        
# [7] "Stress13days_vs_Acute"        "StressResistent_vs_Ctrl"      "StressSusceptible_vs_Ctrl"  

library(MAd)

RunMetaAnalysisDiagnostics<-function(GeneSymbol){

print(GeneSymbol)  
  
effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
var<-as.numeric(MetAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])

FAILED <- FALSE
tryCatch(TempMeta<-rma(effect, var), error = function(e) {FAILED <<- TRUE})

if(FAILED){print("Meta-Analysis Failed")}else{
  TempMeta<-rma(effect, var)
  print(TempMeta)
  print("# of Comparisons:")
  print(8-sum(is.na(effect)))
  print("Cook's Distance:")
  print(data.frame(Study=colnames(MetaAnalysis_FoldChanges_ForMeta )[-1][is.na(effect)==FALSE],Distance=cooks.distance(TempMeta)))
  print("DF Betas:")
  print(data.frame(Study=colnames(MetaAnalysis_FoldChanges_ForMeta )[-1][is.na(effect)==FALSE], DfBetas=dfbetas(TempMeta)))
  rm(TempMeta)
}
rm(effect, var)
}


#Let's start with the pretty, homogenous results: "Fezf2", "Nmb", "Chrm2", "Npy", "Tet1"

RunMetaAnalysisDiagnostics("Fezf2")
# [1] "Fezf2"
# 
# Random-Effects Model (k = 8; tau^2 estimator: REML)
# 
# tau^2 (estimated amount of total heterogeneity): 0.0000 (SE = 0.0028)
# tau (square root of estimated tau^2 value):      0.0015
# I^2 (total heterogeneity / total variability):   0.03%
# H^2 (total variability / sampling variability):  1.00
# 
# Test for Heterogeneity: 
#   Q(df = 7) = 13.9445, p-val = 0.0522
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb    ci.ub     
# -0.1783  0.0241  -7.4134  <.0001  -0.2254  -0.1312  ***
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# [1] "# of Comparisons:"
# [1] 8
# [1] "Cook's Distance:"
# Study     Distance
# 1   StressResilient_Vs_Control 0.3080629254
# 2 StressSusceptible_Vs_Control 0.0174481871
# 3              CUMS_vs_Control 0.1542229719
# 4                  CMS_vs_Ctrl 0.0265373399
# 5         Stress8days_vs_Acute 0.2551501240
# 6        Stress13days_vs_Acute 0.2743558897
# 7      StressResistent_vs_Ctrl 0.0001236311
# 8    StressSusceptible_vs_Ctrl 0.0011287369
# [1] "DF Betas:"
# Study     intrcpt
# 1   StressResilient_Vs_Control -0.50413466
# 2 StressSusceptible_Vs_Control -0.06010801
# 3              CUMS_vs_Control  0.39311042
# 4                  CMS_vs_Ctrl  0.05476901
# 5         Stress8days_vs_Acute  0.17695102
# 6        Stress13days_vs_Acute  0.18381205
# 7      StressResistent_vs_Ctrl -0.01112600
# 8    StressSusceptible_vs_Ctrl  0.03363047


#Interestingly, even though this gene was really pretty, it *almost* shows significant heterogeneity.
#None of the studies we fit outlier status though

RunMetaAnalysisDiagnostics("Nmb")

# [1] "Nmb"
# 
# Random-Effects Model (k = 8; tau^2 estimator: REML)
# 
# tau^2 (estimated amount of total heterogeneity): 0 (SE = 0.0029)
# tau (square root of estimated tau^2 value):      0
# I^2 (total heterogeneity / total variability):   0.00%
# H^2 (total variability / sampling variability):  1.00
# 
# Test for Heterogeneity: 
#   Q(df = 7) = 2.6215, p-val = 0.9177
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb    ci.ub     
# -0.2189  0.0300  -7.2922  <.0001  -0.2777  -0.1601  ***
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# [1] "# of Comparisons:"
# [1] 8
# [1] "Cook's Distance:"
# Study    Distance
# 1   StressResilient_Vs_Control 0.000305984
# 2 StressSusceptible_Vs_Control 0.001138605
# 3              CUMS_vs_Control 0.031879143
# 4                  CMS_vs_Ctrl 0.044626264
# 5         Stress8days_vs_Acute 0.114780281
# 6        Stress13days_vs_Acute 0.010895421
# 7      StressResistent_vs_Ctrl 0.001190420
# 8    StressSusceptible_vs_Ctrl 0.007075822
# [1] "DF Betas:"
# Study     intrcpt
# 1   StressResilient_Vs_Control  0.01749240
# 2 StressSusceptible_Vs_Control -0.03374322
# 3              CUMS_vs_Control  0.17854731
# 4                  CMS_vs_Ctrl -0.21124929
# 5         Stress8days_vs_Acute  0.33879239
# 6        Stress13days_vs_Acute -0.10438114
# 7      StressResistent_vs_Ctrl -0.03450247
# 8    StressSusceptible_vs_Ctrl -0.08411791

RunMetaAnalysisDiagnostics("Chrm2")
# [1] "Chrm2"
# 
# Random-Effects Model (k = 7; tau^2 estimator: REML)
# 
# tau^2 (estimated amount of total heterogeneity): 0 (SE = 0.0140)
# tau (square root of estimated tau^2 value):      0
# I^2 (total heterogeneity / total variability):   0.00%
# H^2 (total variability / sampling variability):  1.00
# 
# Test for Heterogeneity: 
#   Q(df = 6) = 0.9835, p-val = 0.9862
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb    ci.ub     
# -0.2276  0.0614  -3.7067  0.0002  -0.3479  -0.1072  ***
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# [1] "# of Comparisons:"
# [1] 7
# [1] "Cook's Distance:"
# Study     Distance
# 1   StressResilient_Vs_Control 0.1916263657
# 2 StressSusceptible_Vs_Control 0.0042694361
# 3              CUMS_vs_Control 0.0079353172
# 5         Stress8days_vs_Acute 0.0490396479
# 6        Stress13days_vs_Acute 0.0007203939
# 7      StressResistent_vs_Ctrl 0.0025403690
# 8    StressSusceptible_vs_Ctrl 0.0047191532
# [1] "DF Betas:"
# Study     intrcpt
# 1   StressResilient_Vs_Control -0.43775149
# 2 StressSusceptible_Vs_Control  0.06534092
# 3              CUMS_vs_Control  0.08908040
# 5         Stress8days_vs_Acute  0.22144897
# 6        Stress13days_vs_Acute -0.02684015
# 7      StressResistent_vs_Ctrl -0.05040207
# 8    StressSusceptible_vs_Ctrl  0.06869609



#extremely heterogeneous results (but still "sig" in meta-analysis): "S100a8","Gbx2", "Nkx2-1"

RunMetaAnalysisDiagnostics("S100a8")
# [1] "S100a8"
# 
# Random-Effects Model (k = 8; tau^2 estimator: REML)
# 
# tau^2 (estimated amount of total heterogeneity): 0.0068 (SE = 0.0698)
# tau (square root of estimated tau^2 value):      0.0822
# I^2 (total heterogeneity / total variability):   3.87%
# H^2 (total variability / sampling variability):  1.04
# 
# Test for Heterogeneity: 
#   Q(df = 7) = 5.0905, p-val = 0.6489
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub     
# 1.0924  0.1402  7.7914  <.0001  0.8176  1.3672  ***
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# [1] "# of Comparisons:"
# [1] 8
# [1] "Cook's Distance:"
# Study    Distance
# 1   StressResilient_Vs_Control 0.002697505
# 2 StressSusceptible_Vs_Control 0.019525785
# 3              CUMS_vs_Control 0.048117234
# 4                  CMS_vs_Ctrl 1.160755929
# 5         Stress8days_vs_Acute 0.259622689
# 6        Stress13days_vs_Acute 0.080254443
# 7      StressResistent_vs_Ctrl 0.013280482
# 8    StressSusceptible_vs_Ctrl 0.024261336
# [1] "DF Betas:"
# Study     intrcpt
# 1   StressResilient_Vs_Control  0.05090108
# 2 StressSusceptible_Vs_Control  0.13608696
# 3              CUMS_vs_Control -0.22379451
# 4                  CMS_vs_Ctrl  1.12578581
# 5         Stress8days_vs_Acute -0.52335958
# 6        Stress13days_vs_Acute -0.26782301
# 7      StressResistent_vs_Ctrl -0.11391299
# 8    StressSusceptible_vs_Ctrl  0.15123787

#This one isn't significantly heterogeneous - looking at it, it is probably because the error bars are huge for the results.
#But the CMS study is deemed highly influential

RunMetaAnalysisDiagnostics("Gbx2")

# [1] "Gbx2"
# 
# Random-Effects Model (k = 8; tau^2 estimator: REML)
# 
# tau^2 (estimated amount of total heterogeneity): 0 (SE = 0.0607)
# tau (square root of estimated tau^2 value):      0
# I^2 (total heterogeneity / total variability):   0.00%
# H^2 (total variability / sampling variability):  1.00
# 
# Test for Heterogeneity: 
#   Q(df = 7) = 5.4714, p-val = 0.6026
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub     
# 0.4370  0.0604  7.2385  <.0001  0.3187  0.5554  ***
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# [1] "# of Comparisons:"
# [1] 8
# [1] "Cook's Distance:"
# Study     Distance
# 1   StressResilient_Vs_Control 0.2319234263
# 2 StressSusceptible_Vs_Control 0.1683312492
# 3              CUMS_vs_Control 0.0702309617
# 4                  CMS_vs_Ctrl 8.4538473167
# 5         Stress8days_vs_Acute 0.0658208328
# 6        Stress13days_vs_Acute 0.0065925693
# 7      StressResistent_vs_Ctrl 0.0005295941
# 8    StressSusceptible_vs_Ctrl 0.0009859968
# [1] "DF Betas:"
# Study     intrcpt
# 1   StressResilient_Vs_Control -0.28419455
# 2 StressSusceptible_Vs_Control -0.24477079
# 3              CUMS_vs_Control -0.19421137
# 4                  CMS_vs_Ctrl -2.90755005
# 5         Stress8days_vs_Acute  0.25655571
# 6        Stress13days_vs_Acute  0.08119464
# 7      StressResistent_vs_Ctrl -0.02301291
# 8    StressSusceptible_vs_Ctrl -0.03140059

#Ditto.

#Results that are mostly heterogeneous/biased due to the CMS dataset: "Hectd2os", "Rps4l", "Vim"

RunMetaAnalysisDiagnostics("Hectd2os")
[1] "Hectd2os"

# Random-Effects Model (k = 6; tau^2 estimator: REML)
# 
# tau^2 (estimated amount of total heterogeneity): 0.0000 (SE = 0.1299)
# tau (square root of estimated tau^2 value):      0.0022
# I^2 (total heterogeneity / total variability):   0.00%
# H^2 (total variability / sampling variability):  1.00
# 
# Test for Heterogeneity: 
#   Q(df = 5) = 5.4868, p-val = 0.3594
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb    ci.ub     
# -0.1611  0.0192  -8.3693  <.0001  -0.1988  -0.1234  ***
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# [1] "# of Comparisons:"
# [1] 6
# [1] "Cook's Distance:"
# Study     Distance
# 1   StressResilient_Vs_Control 1.321724e-04
# 2 StressSusceptible_Vs_Control 2.904133e-04
# 3              CUMS_vs_Control 5.308566e-04
# 4                  CMS_vs_Ctrl 1.024478e+02
# 7      StressResistent_vs_Ctrl 1.599074e-04
# 8    StressSusceptible_vs_Ctrl 1.293545e+00
# [1] "DF Betas:"
# Study      intrcpt
# 1   StressResilient_Vs_Control   0.01157052
# 2 StressSusceptible_Vs_Control   0.01704483
# 3              CUMS_vs_Control   0.02319696
# 4                  CMS_vs_Ctrl -10.18003552
# 7      StressResistent_vs_Ctrl   0.01264402
# 8    StressSusceptible_vs_Ctrl  -0.26620818

#Interesting- so this metric cannot pick up on the fact that all of the studies have super variable means with large error bars with the exception of the CMS study.
#... but the CMS study has such a large cooks distance, it is represented with an exponent (!)

RunMetaAnalysisDiagnostics("Rps4l")
# [1] "Rps4l"
# 
# Random-Effects Model (k = 8; tau^2 estimator: REML)
# 
# tau^2 (estimated amount of total heterogeneity): 0 (SE = 0.0041)
# tau (square root of estimated tau^2 value):      0
# I^2 (total heterogeneity / total variability):   0.00%
# H^2 (total variability / sampling variability):  1.00
# 
# Test for Heterogeneity: 
#   Q(df = 7) = 2.8271, p-val = 0.9005
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb    ci.ub     
# -0.1146  0.0147  -7.8197  <.0001  -0.1433  -0.0859  ***
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# [1] "# of Comparisons:"
# [1] 8
# [1] "Cook's Distance:"
# Study     Distance
# 1   StressResilient_Vs_Control 1.527351e-03
# 2 StressSusceptible_Vs_Control           NA
# 3              CUMS_vs_Control 1.981083e-03
# 4                  CMS_vs_Ctrl 1.693069e+01
# 5         Stress8days_vs_Acute 4.227787e-03
# 6        Stress13days_vs_Acute 6.750058e-05
# 7      StressResistent_vs_Ctrl 6.619277e-03
# 8    StressSusceptible_vs_Ctrl 9.679011e-03
# [1] "DF Betas:"
# Study      intrcpt
# 1   StressResilient_Vs_Control  0.039081333
# 2 StressSusceptible_Vs_Control           NA
# 3              CUMS_vs_Control  0.044509356
# 4                  CMS_vs_Ctrl -4.114692165
# 5         Stress8days_vs_Acute  0.065021437
# 6        Stress13days_vs_Acute  0.008215874
# 7      StressResistent_vs_Ctrl  0.081358940
# 8    StressSusceptible_vs_Ctrl  0.098381967

#Ditto.

RunMetaAnalysisDiagnostics("Vim")
# [1] "Vim"
# 
# Random-Effects Model (k = 8; tau^2 estimator: REML)
# 
# tau^2 (estimated amount of total heterogeneity): 0 (SE = 0.0066)
# tau (square root of estimated tau^2 value):      0
# I^2 (total heterogeneity / total variability):   0.00%
# H^2 (total variability / sampling variability):  1.00
# 
# Test for Heterogeneity: 
#   Q(df = 7) = 8.7731, p-val = 0.2694
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub     
# 0.0817  0.0187  4.3572  <.0001  0.0449  0.1184  ***
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# [1] "# of Comparisons:"
# [1] 8
# [1] "Cook's Distance:"
# Study     Distance
# 1   StressResilient_Vs_Control 4.440457e-03
# 2 StressSusceptible_Vs_Control 3.147907e-02
# 3              CUMS_vs_Control 2.001137e-03
# 4                  CMS_vs_Ctrl 8.396264e-02
# 5         Stress8days_vs_Acute 5.204104e-03
# 6        Stress13days_vs_Acute 4.190352e-02
# 7      StressResistent_vs_Ctrl 1.628028e-03
# 8    StressSusceptible_vs_Ctrl 9.301728e-05
# [1] "DF Betas:"
# Study      intrcpt
# 1   StressResilient_Vs_Control -0.066636752
# 2 StressSusceptible_Vs_Control -0.177423418
# 3              CUMS_vs_Control -0.044734064
# 4                  CMS_vs_Ctrl  0.065292204
# 5         Stress8days_vs_Acute  0.072139474
# 6        Stress13days_vs_Acute  0.204703498
# 7      StressResistent_vs_Ctrl  0.040348828
# 8    StressSusceptible_vs_Ctrl -0.009554527

#Whereas for this gene, neither the heterogeneity or distance metrics pick up on the fact that the CMS study is driving the results (because the Log2FC falls in the middle of the dataset)


#How about funnel plots? Those are meant for detecting publication bias (which doesn't really play out in transcriptional profiling studies the same way as in other kinds of meta-analyses), but maybe might show that the CMS study has weirdly low SE?

GeneSymbol="Rps4l"
effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
var<-as.numeric(MetAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])

plot(sqrt(var)~effect)

GeneSymbol="Vim"
effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
var<-as.numeric(MetAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])

plot(sqrt(var)~effect)

GeneSymbol="Hectd2os"
effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
var<-as.numeric(MetAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])

plot(sqrt(var)~effect)

#This doesn't really answer the question, since the CMS effect is almost always in the middle of these effects, which is how these genes ended up in the top results. 

#I think what I actually need is something like robust regression - i.e., an indication that the meta-analysis produces completely different results when the CMS dataset isn't included, but excluding other datasets doesn't have the same impact.

robust.rma.mv(TempMeta)
#hmm... looks like it only accepts rma.mv objects - so I need to rewrite the meta-analysis formula
#... also is pretty computationally expensive, and could lead to overfitting. First, let's assess the magnitude of the problem:


DistanceMetrics<-matrix(NA, length(MetaAnalysis_FoldChanges_ForMeta$x), 3)
CooksDistanceForStudies<-matrix(NA, length(MetaAnalysis_FoldChanges_ForMeta$x), 8)
#And then run a loop that run's a meta-analysis on the differential expression results (columns 2-7) for each gene (row):
for(i in c(1:length(MetaAnalysis_FoldChanges_ForMeta$x))){
  
  effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-1])
  var<-as.numeric(MetAnalysis_SV_ForMeta[i,-1])
  
  skip_to_next <- FALSE
  tryCatch(TempMeta<-rma(effect, var), error = function(e) {skip_to_next <<- TRUE})
  
  if(skip_to_next){}else{
    TempMeta<-rma(effect, var)
    DistanceMetrics[i, 1]<-TempMeta$QEp  
    DistanceMetrics[i, 2]<-TempMeta$I2 
    CooksDistance<-cooks.distance(TempMeta)
    DistanceMetrics[i, 3]<-max(CooksDistance) 
    CooksDistanceForStudies[i, as.numeric(names(cooks.distance(TempMeta)))]<-CooksDistance
    rm(TempMeta, CooksDistance)
  }
  rm(effect, var)
}

#That code took around 15 min to run.

colnames(DistanceMetrics)<-c("QEp", "I2", "MaxCooksDistance")
rownames(DistanceMetrics)<-MetaAnalysis_FoldChanges_ForMeta$x

colnames(CooksDistanceForStudies)<-colnames(MetaAnalysis_FoldChanges_ForMeta)[-1]
rownames(CooksDistanceForStudies)<-MetaAnalysis_FoldChanges_ForMeta$x

setwd("~/Documents/Teaching/Grinnell Interns 2022")
write.csv(DistanceMetrics, "DistanceMetrics.csv")
write.csv(CooksDistanceForStudies, "CooksDistanceForStudies.csv")

sum(DistanceMetrics[,1]<0.05, na.rm=TRUE)
#[1] 8160
sum(is.na(DistanceMetrics[,1])==FALSE)
#[1] 24414
#So approximately 1/4 of the results are statistically heterogeneous

hist(DistanceMetrics[,3], breaks=100)
#Wow, I can't visualize it because one of the genes has a cook's distance of 8000, LOL
hist(DistanceMetrics[DistanceMetrics[,3]<40,3], breaks=100)
#The vast majority are less than 10
hist(DistanceMetrics[DistanceMetrics[,3]<10,3], breaks=100)
#But at least half look like they are above 1, which is a very liberal cut-off.
summary(DistanceMetrics[,3])
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#     0.00     0.38     0.65    32.34     1.05 82020.11      474 
#Ok, it's more like 1/4 of the genes, but still a lot.

apply(CooksDistanceForStudies, 2, function(y) median(y, na.rm=TRUE))
# StressResilient_Vs_Control StressSusceptible_Vs_Control     CUMS_vs_Control         CMS_vs_Ctrl 
# 0.20187291                   0.05939162                   0.03014850                   0.22547694 
# Stress8days_vs_Acute        Stress13days_vs_Acute      StressResistent_vs_Ctrl    StressSusceptible_vs_Ctrl 
# 0.04559244                   0.03793348                   0.05946675                   0.04626252 

#Interestingly, with that metrix, the first dataset looks as bad as the CMS dataset.
boxplot(CooksDistanceForStudies)
#And almost all of the studies have an extreme cooks distance for at least a few genes 
boxplot(CooksDistanceForStudies, ylim=c(0,5), ylab="Cook's Distance per gene")
#Yeah, the CMS dataset is the worst, but the first stress resilient vs ctrl is up there too.

pdf("Boxplot_CooksDistance.pdf", width=10, height=5)
boxplot(CooksDistanceForStudies, ylim=c(0,5), ylab="Cook's Distance per gene")
dev.off()

apply(CooksDistanceForStudies, 2, function(y) quantile(y, probs=c(0.75), na.rm=TRUE))
# StressResilient_Vs_Control StressSusceptible_Vs_Control     CUMS_vs_Control      CMS_vs_Ctrl 
# 0.5750926                    0.1903901                    0.1001810                    0.7360000 
# Stress8days_vs_Acute        Stress13days_vs_Acute      StressResistent_vs_Ctrl    StressSusceptible_vs_Ctrl 
# 0.1473577                    0.1321165                    0.1715968                    0.1393737 

apply(CooksDistanceForStudies, 2, function(y) quantile(y, probs=c(0.90), na.rm=TRUE))
# StressResilient_Vs_Control StressSusceptible_Vs_Control   CUMS_vs_Control             CMS_vs_Ctrl 
# 1.0494830                    0.4322668                    0.2397078                    1.5842185 
# Stress8days_vs_Acute        Stress13days_vs_Acute      StressResistent_vs_Ctrl    StressSusceptible_vs_Ctrl 
# 0.3474975                    0.3559948                    0.4005313                    0.3155570 


#################################
setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01")

HRLR<-read.csv("F0_Meta_F2_4Models_withRnor6v88Coord_JustGenes.csv", header=TRUE, stringsAsFactors = FALSE)
str(HRLR)
colnames(HRLR)
HRLR_DF<-HRLR[,c(3:4, 6,11,16,21,29,31,35,36,38)]

HRLR_vs_metaOutput<-join(HRLR_DF, metaOutputFDR_DF, by="Symbol", type="inner")

plot(HRLR_vs_metaOutput$estimate~HRLR_vs_metaOutput$ChronicStress_estimate)
#also pretty underwhelming
plot(HRLR_vs_metaOutput$Coef.Lineage_AsFactorbLR~HRLR_vs_metaOutput$ChronicStress_estimate)
#ditto...

