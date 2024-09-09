#Sample code for MetaData extraction for Gemma datasets, part2: Filtering
#Megan Hagenauer, Jun 27 2024

############

#This code has not been cleaned up, fully annotated, functionalized, and made easily generalizable

###############

library(gemma.R)
library(dplyr)

###############

#After we've decided on our search terms:
    
MyQueryTerms<-"SSRI fluoxetine sertraline paroxetine citalopram escitalopram fluvoxamine vilazodone vortioxetine amitriptyline imipramine amoxapine desipramine nortriptyline clomipramine trimipramine protriptyline doxepin maprotiline trazadone nefazodone mirtazapine phenelzine nialamide isocarboxazid hydracarbazine tranylcypromine selegiline venlafaxine desvenlafaxine duloxetine levomilnacipran bupropion duloxetine levomilnacipran ketamine esketamine tianeptine brexanolone"

MyQueryTerms 

result_MyQueryTerms_NoFilter<- gemma.R ::get_datasets(query=MyQueryTerms) %>% 
  gemma.R:::get_all_pages() 

result_MyQueryTerms_NoFilter
#Classes ‘data.table’ and 'data.frame':	194 obs. of  23 variables:


#What are some characteristics that we might want to filter by?

result_MyQueryTerms_NoFilter$taxon.name

table(result_MyQueryTerms_NoFilter$taxon.name)
# human mouse   rat 
# 83    91    20
#We can only easily use the mouse and rat datasets

result_MyQueryTerms_NoFilter$experiment.troubled

table(result_MyQueryTerms_NoFilter$experiment.troubled)
# FALSE 
# 194 

#We only want the non-troubled datasets (so FALSE) - looks good!

table(result_MyQueryTerms_NoFilter$experiment.rawData)
# -1   1 
# 47 144 
#The experiments marked "-1" don't have raw data available
#So Gemma needed to import their summarized gene expression values from an external source
#That source may or may not have used standard normalization
#We will need to tread carefully with those datasets

#Filtering by taxa (mouse, rat) within the query itself:
result_MyQueryTerms_RatsMice<- gemma.R ::get_datasets(query=MyQueryTerms, taxa = c("mouse", "rat")) %>% 
  gemma.R:::get_all_pages() 

str(result_MyQueryTerms_RatsMice)
#Classes ‘data.table’ and 'data.frame':	111 obs. of  23 variables:

#Filtering by quality can probably be done within the search, but here's the subsetting version:

#result_MyQueryTerms_RatsMice[rows,columns]

#Examples of grabbing information from the first column:
result_MyQueryTerms_RatsMice$experiment.shortName
result_MyQueryTerms_RatsMice[,1]

#Example of grabbing info from the first row
result_MyQueryTerms_RatsMice[1,]
  
result_MyQueryTerms_Filtered<-result_MyQueryTerms_RatsMice[result_MyQueryTerms_RatsMice$experiment.troubled==FALSE,]

str(result_MyQueryTerms_Filtered)
#Classes ‘data.table’ and 'data.frame':	111 obs. of  23 variables:

#Here's how that same code would be written in the style of the tidyverse:
result_MyQueryTerms_Filtered<-result_MyQueryTerms_RatsMice %>% filter(experiment.troubled==FALSE)

str(result_MyQueryTerms_Filtered)
#Classes ‘data.table’ and 'data.frame':	111 obs. of  23 variables:

#Let's write that out for our records:
#Comma separate variable (.csv) file format
#CSV files are easily read into commonly used spreadsheet programs like Excel or GoogleSheets
write.csv(result_MyQueryTerms_Filtered, "result_MyQueryTerms_Filtered.csv")

#Where did my file go?
getwd()

###################

#Let's find out what brain regions we have:

MyResults<-result_MyQueryTerms_Filtered

#Note: some datasets have data from more than one "organism part"
#That means that if we survey all of the organism parts in our datasets
#we may end up with a vector of organism parts that is of unknown length 
#and potentially longer than our number of datasets

#Creating an empty vector to save our results 
#We can just add on to this using concatenate
#That way we don't need to prespecify length
OrganismPartAnnotations_All<-vector(mode="character")

#Looping over all of the datasets:

#How many datasets do we have?
length(MyResults$experiment.shortName)
nrow(MyResults)

for(i in c(1:nrow(MyResults))){
  
  ExperimentName<-MyResults$experiment.shortName[i]
  
  ExperimentAnnotations<-get_dataset_annotations(dataset=ExperimentName)
  rm(ExperimentName)
  
  #Making sure that there is "organism part" annotation for the dataset before attempting to extract it:
  
  #Do we have annotation for organism part? Here is a true/false vector
  #ExperimentAnnotations$class.name=="organism part"
  #This is the column of terms for the annotation
  #ExperimentAnnotations$term.name 
  #Example of subsetting for organism part:
  #ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"]
  
  #How many annotations do we have?
  length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"])
  
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"])>0){
    
  OrganismPartAnnotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"]
  rm(ExperimentAnnotations)
  
  OrganismPartAnnotations_All<-c(OrganismPartAnnotations_All, OrganismPartAnnotations)
  rm(OrganismPartAnnotations)
  
  #if there isn't organism part annotation, we just tell the loop to move on to the next dataset entry:
  
  }else{
    rm(ExperimentAnnotations)
  }
}

#This is what our result looks like:
OrganismPartAnnotations_All

#How many datasets have each organism part?  Let's make a summary table!

table(OrganismPartAnnotations_All)

#That is probably worth saving:

write.csv(table(OrganismPartAnnotations_All), "Table_OrganismPartAnnotations.csv")

#I was planning to filter by developmental stage annotations as well
#But snooping at a lot of the datasets, they seem to lack annotation for that.

########################

#If we then want to narrow our datasets down to a particular brain region:

#E.g., Hippocampal formation:
#I'm going to use the URI ontology term for the hippocampus so that I catch all of datasets annotated with the child terms as well:
result_MyQueryTerms_RatsMice_Hippocampus <- gemma.R ::get_datasets(query=MyQueryTerms, filter = 'allCharacteristics.valueUri in (http://purl.obolibrary.org/obo/UBERON_0002421)', taxa = c("mouse", "rat")) %>% 
  gemma.R:::get_all_pages() 

str(result_MyQueryTerms_RatsMice_Hippocampus)
#Classes ‘data.table’ and 'data.frame':	31 obs. of  23 variables:

#... then filter down again to high quality data:

result_MyQueryTerms_RatsMice_Hippocampus_Filtered<-result_MyQueryTerms_RatsMice_Hippocampus[result_MyQueryTerms_RatsMice_Hippocampus$experiment.troubled==FALSE,]

str(result_MyQueryTerms_RatsMice_Hippocampus_Filtered)

#Classes ‘data.table’ and 'data.frame':	31 obs. of  23 variables:

#How many of those are from questionable "external" datasets without raw data available:
table(result_MyQueryTerms_RatsMice_Hippocampus_Filtered$experiment.rawData)
#-1  1 
#11 20
#So 20 results that are more likely to be usable
#and 11 results that may or may not be usable depending on how the original authors of the study normalized their data

#We could write out these results to dig through them more easily in a spreadsheet program like Excel or Google Sheets

write.csv(result_MyQueryTerms_RatsMice_Hippocampus_Filtered, "result_MyQueryTerms_RatsMice_Hippocampus_Filtered.csv")

#We could then snoop again to see what the distribution is for organism part in our results:

MyResults<-result_MyQueryTerms_RatsMice_Hippocampus_Filtered

#I'm reusing an entire chunk of code from earlier - guess I should have functionalized it!

OrganismPartAnnotations_All<-vector(mode="character")

for(i in c(1:length(MyResults$experiment.shortName))){
  ExperimentName<-MyResults$experiment.shortName[i]
  
  ExperimentAnnotations<-get_dataset_annotations(dataset=ExperimentName)
  rm(ExperimentName)
  
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"])>0){
    
    OrganismPartAnnotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"]
    rm(ExperimentAnnotations)
    
    OrganismPartAnnotations_All<-c(OrganismPartAnnotations_All, OrganismPartAnnotations)
    rm(OrganismPartAnnotations)
    
  }else{
    rm(ExperimentAnnotations)
  }
}

table(OrganismPartAnnotations_All)
#25 are annotated with Ammon's Horn, 5 are dentate gyrus
#Some of these organism part annotations... aren't organism parts. I wonder if we should do this filtering by hand.

write.csv(OrganismPartAnnotations_All, "OrganismPartAnnotations_AfterFilteringToRegion.csv")

