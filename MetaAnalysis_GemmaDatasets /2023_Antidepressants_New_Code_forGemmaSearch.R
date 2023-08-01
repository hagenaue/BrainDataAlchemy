#Code for MetaData Extraction for Gemma Datasets
#Topic: Antidepressant effects in the hippocampus
#Megan Hagenauer, June 23, 2023

#####################

#General notes:

#This meta-analysis topic was originally pursued by Erin Hernandez (summer 2022), but after exclusions there were very few datasets left
#I'm going to try re-running it this summer (2023) with different search terms with the hope of finding more applicable datasets
#I also plan to fix the agilent datasets (which need additional QC and DE analysis - they look like they were log2 transformed twice) and that should add to the n

#####################

#To do this work, we will be using an R wrapper for Gemma's restful API.
#Here's the github site and documentation for it:
https://github.com/PavlidisLab/gemma.R
https://pavlidislab.github.io/gemma.R/articles/gemma.R.html


#####################

#1. Installing the Gemma API package:

if(!requireNamespace("devtools", quietly=T)){
  install.packages("devtools")
}

devtools:: install_github("PavlidisLab/gemma.R", force=T)

#####################

#2. Create a folder for output on your computer and then set the working directory to that own folder.
#I find that the easiest way to do this is via the GUI drop down menu (Session-> Set Working Directory-> Choose Directory), but then I save the code that appears in the Console so that I can easily navigate back.

#Here is an example of the code for setting a working directory to my own folder:
setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Teaching/BrainDataAlchemy/Grinnell Interns 2022/2022Cohort_GoogleDrive/GoogleDriveDownload_20230508/ErinHernandez_Antidepressants_Hippocampus/Methods_MH_2023")

#If you want to double-check which working directory you are in:
getwd()

#To list all the files in your working directory:
list.files()

#########################

#3. Searching for datasets that have metadata that includes particular terms:

#About the function:
help("gemma.R ::search_datasets")
??gemma.R ::search_datasets


#Some idiosyncracies:

#The boolean operator * can be used to capture all variations of a word that starts in a particular way
#e.g., hippocamp* can capture hippocampus, hippocampal, hippocampi

#but some other common search syntax doesn't work normally:

#using two words in a search captures results that include both terms:
#e.g., "hippocamp* antidepress*" would capture all results that include both hippocamp* and antidepress*
result <- gemma.R ::search_datasets("hippocamp* antidepress*", taxon = 'mouse',limit = 100) 

#but, strangely enough, I get slightly different results using "hippocamp* AND antidepress*"
result <- gemma.R ::search_datasets("hippocamp* AND antidepress*", taxon = 'mouse',limit = 100) 
#I don't know why this produces different results
#reviewing them, the version that doesn't include the AND seems more correct.

#using OR can make it so that you can get results that fit one set of terms or a different set of terms, e.g.,
result2 <- gemma.R ::search_datasets("hippocamp* antidepress*", taxon = 'mouse',limit = 100) 
result3 <- gemma.R ::search_datasets("dentate gyrus antidepress*", taxon = 'mouse',limit = 100) 

resultVenn<-gemma.R ::search_datasets("(dentate gyrus antidepress*) OR (hippocamp* antidepress*)", taxon = 'mouse',limit = 100) 
#17 obs
sum(result2$experiment.ShortName%in%result3$experiment.ShortName)
#[1] 9
sum(result3$experiment.ShortName%in%result2$experiment.ShortName)
#[1] 9
(16-9)+9+(10-9)
#[1] 17

SearchTerms<-"(dentate gyrus antidepress*) OR (hippocamp* antidepress*)"

resultVenn<-gemma.R ::search_datasets(SearchTerms, taxon = 'mouse',limit = 100) 

#...but for some reason you can't reverse the syntax (I suspect because the AND argument doesn't work as expected) 
SearchTerms<-"(dentate gyrus OR hippocamp*) AND (antidepress*)"
resultVenn<-gemma.R ::search_datasets(SearchTerms, taxon = 'mouse',limit = 100) 
#I don't know why this doesn't produce the same results


#########################

#Here's an example using some syntax that works:

#Search terms from Erin's reading, HDRF conferences, or from the FDA website (https://www.fda.gov/consumers/free-publications-women/depression-medicines#Tricyclic-Tetracyclic)

SearchTerms<-"
(dentate gyrus antidepress*) OR 
(hippocamp* antidepress*) OR 
(Ammon's Horn antidepress*) OR 
(dentate gyrus SSRI) OR
(hippocamp* SSRI) OR
(Ammon's Horn SSRI) OR
(dentate gyrus selective serotonin reuptake inhibitor*) OR 
(hippocamp* selective serotonin reuptake inhibitor*) OR 
(Ammon's Horn selective serotonin reuptake inhibitor*) OR 
(dentate gyrus fluoxetine) OR
(hippocamp* fluoxetine) OR
(Ammon's Horn fluoxetine) OR
(dentate gyrus sertraline) OR
(hippocamp* sertraline) OR
(Ammon's Horn sertraline) OR
(dentate gyrus paroxetine) OR
(hippocamp* paroxetine) OR
(Ammon's Horn paroxetine) OR
(dentate gyrus citalopram) OR
(hippocamp* citalopram) OR
(Ammon's Horn citalopram) OR
(dentate gyrus escitalopram) OR
(hippocamp* escitalopram) OR
(Ammon's Horn escitalopram) OR
(dentate gyrus fluvoxamine) OR
(hippocamp* fluvoxamine) OR
(Ammon's Horn fluvoxamine) OR
(dentate gyrus vilazodone) OR
(hippocamp* vilazodone) OR
(Ammon's Horn vilazodone) OR
(dentate gyrus vortioxetine) OR
(hippocamp* vortioxetine) OR
(Ammon's Horn vortioxetine) OR
(dentate gyrus tricyclic) OR
(hippocamp* tricyclic) OR
(Ammon's Horn tricyclic) OR
(dentate gyrus TCA*) OR
(hippocamp* TCA*) OR
(Ammon's Horn TCA*) OR
(dentate gyrus amitriptyline) OR
(hippocamp* amitriptyline) OR
(Ammon's Horn amitriptyline) OR
(dentate gyrus imipramine) OR
(hippocamp* imipramine) OR
(Ammon's Horn imipramine) OR
(dentate gyrus amoxapine) OR
(hippocamp* amoxapine) OR
(Ammon's Horn amoxapine) OR
(dentate gyrus desipramine) OR
(hippocamp* desipramine) OR
(Ammon's Horn desipramine) OR
(dentate gyrus nortriptyline) OR
(hippocamp* nortriptyline) OR
(Ammon's Horn nortriptyline) OR
(dentate gyrus clomipramine) OR
(hippocamp* clomipramine) OR
(Ammon's Horn clomipramine) OR
(dentate gyrus trimipramine) OR
(hippocamp* trimipramine) OR
(Ammon's Horn trimipramine) OR
(dentate gyrus protriptyline) OR
(hippocamp* protriptyline) OR
(Ammon's Horn protriptyline) OR
(dentate gyrus doxepin) OR
(hippocamp* doxepin) OR
(Ammon's Horn doxepin) OR
(dentate gyrus tetracyclic) OR
(hippocamp* tetracyclic) OR
(Ammon's Horn tetracyclic) OR
(dentate gyrus maprotiline) OR
(hippocamp* maprotiline) OR
(Ammon's Horn maprotiline) OR
(dentate gyrus trazadone) OR
(hippocamp* trazadone) OR
(Ammon's Horn trazadone) OR
(dentate gyrus nefazodone) OR
(hippocamp* nefazodone) OR
(Ammon's Horn nefazodone) OR
(dentate gyrus mirtazapine) OR
(hippocamp* mirtazapine) OR
(Ammon's Horn mirtazapine) OR
(dentate gyrus MAOI*) OR
(hippocamp* MAOI*) OR
(Ammon's Horn MAOI*) OR
(dentate gyrus monoamine oxidase inhibitor*) OR
(hippocamp* monoamine oxidase inhibitor*) OR
(Ammon's Horn monoamine oxidase inhibitor*) OR
(dentate gyrus phenelzine) OR
(hippocamp* phenelzine) OR
(Ammon's Horn phenelzine) OR
(dentate gyrus nialamide) OR
(hippocamp* nialamide) OR
(Ammon's Horn nialamide) OR
(dentate gyrus isocarboxazid) OR
(hippocamp* isocarboxazid) OR
(Ammon's Horn isocarboxazid) OR
(dentate gyrus hydracarbazine) OR
(hippocamp* hydracarbazine) OR
(Ammon's Horn hydracarbazine) OR
(dentate gyrus tranylcypromine) OR
(hippocamp* tranylcypromine) OR
(Ammon's Horn tranylcypromine) OR
(dentate gyrus selegiline) OR
(hippocamp* selegiline) OR
(Ammon's Horn selegiline) OR
(dentate gyrus SNRI*) OR
(hippocamp* SNRI*) OR
(Ammon's Horn SNRI*) OR
(dentate gyrus serotonin norepinephrine reuptake inhibitor*) OR
(hippocamp* serotonin norepinephrine reuptake inhibitor*) OR
(Ammon's Horn serotonin norepinephrine reuptake inhibitor*) OR
(dentate gyrus venlafaxine) OR
(hippocamp* venlafaxine) OR
(Ammon's Horn venlafaxine) OR
(dentate gyrus desvenlafaxine) OR
(hippocamp* desvenlafaxine) OR
(Ammon's Horn desvenlafaxine) OR
(dentate gyrus duloxetine) OR
(hippocamp* duloxetine) OR
(Ammon's Horn duloxetine) OR
(dentate gyrus levomilnacipran) OR
(hippocamp* levomilnacipran) OR
(Ammon's Horn levomilnacipran) OR
(dentate gyrus norepinephrine dopamine reuptake inhibitor*) OR
(hippocamp* norepinephrine dopamine reuptake inhibitor*) OR
(Ammon's Horn norepinephrine dopamine reuptake inhibitor*) OR
(dentate gyrus bupropion) OR
(hippocamp* bupropion) OR
(Ammon's Horn bupropion) OR
(dentate gyrus duloxetine) OR
(hippocamp* duloxetine) OR
(Ammon's Horn duloxetine) OR
(dentate gyrus levomilnacipran) OR
(hippocamp* levomilnacipran) OR
(Ammon's Horn levomilnacipran) OR
(dentate gyrus ketamine) OR
(hippocamp* ketamine) OR
(Ammon's Horn ketamine) OR
(dentate gyrus esketamine) OR
(hippocamp* esketamine) OR
(Ammon's Horn esketamine) OR
(dentate gyrus tianeptine) OR
(hippocamp* tianeptine) OR
(Ammon's Horn tianeptine) OR
(dentate gyrus brexanolone) OR
(hippocamp* brexanolone) OR
(Ammon's Horn brexanolone)"

results_Antidepressants_Mice<-gemma.R ::search_datasets(SearchTerms, taxon = 'mouse',limit = 100) 
#27 observations
results_Antidepressants_Rats<-gemma.R ::search_datasets(SearchTerms, taxon = 'rat',limit = 100) 
#11 observations

results_Antidepressants_Combined<-rbind.data.frame(results_Antidepressants_Mice, results_Antidepressants_Rats) 
#38 observations
results_Antidepressants_Combined_Unique<-unique(results_Antidepressants_Combined)
#38 observations

write.csv(results_Antidepressants_Combined_Unique, "results_Antidepressants_Combined_Unique.csv")

###################

#What if there are more than 100 results?

#From Gemma API vignette:
# Note that a single call of these functions will only return 20 results by default and a 100 results maximum, controlled by the limit argument. 
# In order to get all available results, offset argument should be used to make multiple calls.

SearchTerms<-"(dentate gyrus stress*) OR 
(hippocamp* stress*) OR 
(Ammon's Horn stress*)"

result <- gemma.R ::search_datasets(SearchTerms, taxon = 'mouse',limit = 1) 
print(attributes(result)$totalElements)
#[1] 111

results_Stress_Hippocampus_Mice<-gemma.R ::search_datasets(SearchTerms, taxon = 'mouse',limit = 100) 
results_Stress_Hippocampus_Mice2<-gemma.R ::search_datasets(SearchTerms, taxon = 'mouse',limit = 100, offset=100) 

results_Stress_Hippocampus_Mice_Combined<-rbind.data.frame(results_Stress_Hippocampus_Mice, results_Stress_Hippocampus_Mice2)

str(results_Stress_Hippocampus_Mice_Combined)
#Classes ‘data.table’ and 'data.frame':	111 obs. of  23 variables:

result <- gemma.R ::search_datasets(SearchTerms, taxon = 'mouse',limit = 1) 
print(attributes(result)$totalElements)
#[1]  111

SearchIterations<-floor(attributes(result)$totalElements/100)

TempResult<-gemma.R ::search_datasets(SearchTerms, taxon = 'mouse',limit = 100) 

for(i in c(1:SearchIterations)){
  TempResult<-rbind.data.frame(TempResult,(gemma.R ::search_datasets(SearchTerms, taxon = 'mouse',limit = 100, offset=i*100)))
}

str(TempResult)
#Classes ‘data.table’ and 'data.frame':	110 obs. of  23 variables:

result <- gemma.R ::search_datasets(SearchTerms, taxon = 'rat',limit = 1) 
print(attributes(result)$totalElements)
#[1] 36


#######################

#Snooping through those results:

#https://pavlidislab.github.io/Gemma/geeq.html
#We need to be careful with the external datasets - some of them are log2 transformed twice, and almost all of them lack batch information.

#How many are missing Raw Data?

table(results_Antidepressants_Combined_Unique$geeq.rawData)
# -1  1 
# 13 23 

#How many are troubled?
table(results_Antidepressants_Combined_Unique$experiment.Troubled)
#FALSE 
#38 


results_Antidepressants_Combined_Unique$experiment.ShortName

#This pulls up the metadata for the differential expression analyses for a dataset:
#If we are going to use this output, the factor of interest for the meta-analysis (i.e. an antidepressant) and a control group should be included in the analysis
gemma.R::get_dataset_differential_expression_analyses ("GSE27532")
GSE27532_Design<-gemma.R::get_dataset_differential_expression_analyses ("GSE27532")
#Checking to see whether an antidepressant is included as a factor in the experiment:
GSE27532_Design$experimental.factorValue
#[1] "desipramine"      "high swim stress-induced analgesia"             "desipramine_high swim stress-induced analgesia"
GSE27532_Design$baseline.factorValue
#[1] "reference substance role"          "low swim stress-induced analgesia" NA    
#Yep

#This pulls up the sample metadata for any particular dataset:
GSE27532_SampleMetaData<-gemma.R::get_dataset_samples("GSE27532")
str(GSE27532_SampleMetaData)
GSE27532_SampleMetaData$sample.Name
#The variables used for the differential expression analysis for all samples:
GSE27532_SampleMetaData$sample.FactorValues
#More details for all of the samples:
GSE27532_SampleMetaData$sample.Characteristics
  
#This pulls up the log2 expression for the dataset:
gemma.R::get_dataset_expression("GSE27532")
#Note that these values should probably range between -1 and 12 for RNA-Seq and 4-12 for microarray. 
#If you are only seeing values in a narrow range (e.g., between 2-3.58) that suggests there may be a problem, esp. if the data was external and no raw data was available
#I suspect this is happening because the external microarray data is often only available in log2 format and then accidentally log2-transformed again within Gemma's pipeline



#This pulls out the actual differential expression results for a dataset:
GSE27532_DEResults<-gemma.R::get_differential_expression_values ("GSE27532")
#Note that this object may contain different sets of results
str(GSE27532_DEResults)
# List of 3

# $ 546193:Classes ‘data.table’ and 'data.frame':	25659 obs. of  10 variables:
# ..$ Probe                : chr [1:25659] "ILMN_2616509" "ILMN_3120335" "ILMN_1225966" "ILMN_3132050" ...
# ..$ NCBIid               : chr [1:25659] "11864" "14402" "67903" "19417" ...
# ..$ GeneSymbol           : chr [1:25659] "Arnt2" "Gabrb3" "Gipc1" "Rasgrf1" ...
# ..$ GeneName             : chr [1:25659] "aryl hydrocarbon receptor nuclear translocator 2" "gamma-aminobutyric acid (GABA) A receptor, subunit beta 3" "GIPC PDZ domain containing family, member 1" "RAS protein-specific guanine nucleotide-releasing factor 1" ...
# ..$ pvalue               : num [1:25659] 1.87e-01 4.43e-01 7.41e-05 3.89e-01 8.04e-01 ...
# ..$ corrected_pvalue     : num [1:25659] 0.7263 0.8778 0.0322 0.8532 0.9686 ...
# ..$ rank                 : num [1:25659] 0.2569 0.5045 0.0023 0.4558 0.8299 ...
# ..$ contrast_67519_log2fc: num [1:25659] -0.1765 -0.0512 0.2352 0.059 0.0154 ...
# ..$ contrast_67519_tstat : num [1:25659] -1.605 -0.523 4.344 0.444 0.289 ...
# ..$ contrast_67519_pvalue: num [1:25659] 0.1264 0.6078 0.000416 0.6623 0.7756 ...
# ..- attr(*, ".internal.selfref")=<externalptr> 
#   ..- attr(*, "call")= chr "https://gemma.msl.ubc.ca/rest/v2/resultSets/546193"

# $ 546194:Classes ‘data.table’ and 'data.frame':	25656 obs. of  10 variables:
#   ..$ Probe                : chr [1:25656] "ILMN_2710604" "ILMN_2712697" "ILMN_1231518" "ILMN_2740764" ...
# ..$ NCBIid               : chr [1:25656] "13169" "23970" "213012" "56613" ...
# ..$ GeneSymbol           : chr [1:25656] "Dbnl" "Pacsin2" "Abhd10" "Rps6ka4" ...
# ..$ GeneName             : chr [1:25656] "drebrin-like" "protein kinase C and casein kinase substrate in neurons 2" "abhydrolase domain containing 10" "ribosomal protein S6 kinase, polypeptide 4" ...
# ..$ pvalue               : num [1:25656] 3.08e-04 6.71e-01 1.98e-01 6.94e-01 2.72e-05 ...
# ..$ corrected_pvalue     : num [1:25656] 0.005671 0.8613 0.4851 0.874 0.000923 ...
# ..$ rank                 : num [1:25656] 0.0542 0.7786 0.4086 0.7934 0.0294 ...
# ..$ contrast_67520_log2fc: num [1:25656] -0.1712 -0.00822 0.01658 0.00458 -0.192 ...
# ..$ contrast_67520_tstat : num [1:25656] -2.35 -0.1726 0.337 0.0655 -3.47 ...
# ..$ contrast_67520_pvalue: num [1:25656] 0.03079 0.865 0.7401 0.9485 0.00283 ...
# ..- attr(*, ".internal.selfref")=<externalptr> 
#   ..- attr(*, "call")= chr "https://gemma.msl.ubc.ca/rest/v2/resultSets/546194"

# $ 546195:Classes ‘data.table’ and 'data.frame':	25658 obs. of  10 variables:
#   ..$ Probe                      : chr [1:25658] "ILMN_2642800" "ILMN_3132361" "ILMN_2641278" "ILMN_1248696" ...
# ..$ NCBIid                     : chr [1:25658] "66042" "27057" "170731" "100604" ...
# ..$ GeneSymbol                 : chr [1:25658] "Sostdc1" "Ncoa4" "Mfn2" "Lrrc8c" ...
# ..$ GeneName                   : chr [1:25658] "sclerostin domain containing 1" "nuclear receptor coactivator 4" "mitofusin 2" "leucine rich repeat containing 8 family, member C" ...
# ..$ pvalue                     : num [1:25658] 0.0276 0.8884 0.5781 0.2474 0.8486 ...
# ..$ corrected_pvalue           : num [1:25658] 0.64 0.991 0.962 0.892 0.99 ...
# ..$ rank                       : num [1:25658] 0.043 0.896 0.601 0.277 0.857 ...
# ..$ contrast_67519_67520_log2fc: num [1:25658] 0.4492 0.0242 -0.0465 -0.0949 0.0115 ...
# ..$ contrast_67519_67520_tstat : num [1:25658] 2.403 0.142 -0.567 -1.197 0.194 ...
# ..$ contrast_67519_67520_pvalue: num [1:25658] 0.0276 0.8884 0.5781 0.2474 0.8486 ...
# ..- attr(*, ".internal.selfref")=<externalptr> 
#   ..- attr(*, "call")= chr "https://gemma.msl.ubc.ca/rest/v2/resultSets/546195"

#The results have result.ids that match the result.ids in GSE27532_Design
#Antidepressant treatment is result id 546193
GSE27532_DEResults_desipramine<-GSE27532_DEResults[[1]]
colnames(GSE27532_DEResults_desipramine)
# [1] "Probe"                 "NCBIid"                "GeneSymbol"            "GeneName"              "pvalue"                "corrected_pvalue"     
# [7] "rank"                  "contrast_67519_log2fc" "contrast_67519_tstat"  "contrast_67519_pvalue"
