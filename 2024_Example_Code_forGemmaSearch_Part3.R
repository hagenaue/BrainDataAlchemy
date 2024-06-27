#Sample code for MetaData extraction for Gemma datasets, part3: Annotating
#Megan Hagenauer, Jun 27 2024

############

#This code has not been cleaned up, fully annotated, functionalized, and made easily generalizable

###############

library(gemma.R)
library(dplyr)

###############

#First - just as a reminder, this is my toy Gemma search from earlier (part 1 & part 2 of the Gemma search code):
    
MyQueryTerms<-"SSRI fluoxetine sertraline paroxetine citalopram escitalopram fluvoxamine vilazodone vortioxetine amitriptyline imipramine amoxapine desipramine nortriptyline clomipramine trimipramine protriptyline doxepin maprotiline trazadone nefazodone mirtazapine phenelzine nialamide isocarboxazid hydracarbazine tranylcypromine selegiline venlafaxine desvenlafaxine duloxetine levomilnacipran bupropion duloxetine levomilnacipran ketamine esketamine tianeptine brexanolone"

#I searched for hippocampal datasets from rats and mice that fit my query terms:
result_MyQueryTerms_RatsMice_Hippocampus <- gemma.R ::get_datasets(query=MyQueryTerms, filter = 'allCharacteristics.valueUri in (http://purl.obolibrary.org/obo/UBERON_0002421)', taxa = c("mouse", "rat")) %>% 
  gemma.R:::get_all_pages() 

str(result_MyQueryTerms_RatsMice_Hippocampus)
#Classes ‘data.table’ and 'data.frame':	31 obs. of  23 variables:

#... and then filtered down to high quality data:

result_MyQueryTerms_RatsMice_Hippocampus_Filtered<-result_MyQueryTerms_RatsMice_Hippocampus[result_MyQueryTerms_RatsMice_Hippocampus$experiment.troubled==FALSE,]

str(result_MyQueryTerms_RatsMice_Hippocampus_Filtered)

######################################

#What if we wanted to add additional annotation to this basic output from the getDatasets function?
#e.g., info about the organism parts, developmental stage annotation, experimental factors used in the analyses
#This could make triaging (inclusion/exclusion decisions) so much easier!

#I wrote some code to do this.
#I should probably functionalize it...
#... but I haven't yet, so I've just renamed my results as something generic 
# to make it easier to functionalize later...
MyResults<-result_MyQueryTerms_RatsMice_Hippocampus_Filtered

#Let's make some empty vectors that are the same length as the columns in our results
#These empty vectors will be used to store our annotations while we loop through the rows of datasets:
OrganismParts<-vector(mode="character", length=nrow(MyResults))
CellTypes<-vector(mode="character", length=nrow(MyResults))
DevelopmentalStages<-vector(mode="character", length=nrow(MyResults))
Treatments<-vector(mode="character", length=nrow(MyResults))
Diseases<-vector(mode="character", length=nrow(MyResults))
DiseaseModels<-vector(mode="character", length=nrow(MyResults))
Genotypes<-vector(mode="character", length=nrow(MyResults))
Strains<-vector(mode="character", length=nrow(MyResults))
Sex<-vector(mode="character", length=nrow(MyResults))

#I'm going to loop over all of the rows (row number =i) in my results (i.e., dataset metadata)
#And collect all of this annotation information
#And then format it in a way so that it can be added into my simple dataframe of results
#And then outputted and read easily in a spreadsheet program like excel

for(i in c(1:nrow(MyResults))){
  
  #Pulling out the name for the dataset in a row (row number=i):
  ExperimentName<-MyResults$experiment.shortName[i]
  
  #Accessing the annotations for the dataset:
  ExperimentAnnotations<-get_dataset_annotations(dataset=ExperimentName)
  #The number and type of annotations for the datasets is quite variable
  
  rm(ExperimentName)
  
  #Determining whether there is any annotation for organism part:
  
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"])>0){
    
    #If there is organism part annotation, I'm grabbing it:
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"]
    
    #And then collapsing that vector of annotations into a single string 
    #that can be easily stashed in a single cell in a data.frame (or Excel spreadsheet) 
    #This will eventually become part of the the row for that dataset in the results
    # e.g., "annotation 1; annotation 2; annotation 3"
    OrganismParts[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  #If there isn't any annotation for organism part, we move on to the next type of annotation:
  }else{}
  
  #Now grabbing the annotation for cell type in a similar manner: 
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="cell type"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="cell type"]
    
    CellTypes[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  #Now grabbing the annotation for developmental stage in a similar manner:
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="developmental stage"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="developmental stage"]
    
    DevelopmentalStages[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  #Now grabbing the annotation for treatment in a similar manner:
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="treatment"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="treatment"]
    
    Treatments[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  #Now grabbing the annotation for disease in a similar manner:
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="disease"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="disease"]
    
    Diseases[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  #Now grabbing the annotation for disease model in a similar manner:
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="Disease model"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="Disease model"]
    
    DiseaseModels[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  #Now grabbing the annotation for genotype in a similar manner:
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="genotype"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="genotype"]
    
    Genotypes[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  #Now grabbing the annotation for strain in a similar manner:
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="strain"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="strain"]
    
    Strains[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  #Now grabbing the annotation for biological sex in a similar manner:
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="biological sex"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="biological sex"]
    
    Sex[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  rm(ExperimentAnnotations)
}

#Adding all of those vectors of annotation to my data.frame of results:
MyResults_Annotated<-cbind.data.frame(MyResults, 
                                OrganismParts,
                                CellTypes,
                                DevelopmentalStages,
                                Treatments,
                                Diseases,
                                DiseaseModels,
                                Genotypes,
                                Strains,
                                Sex)

#very pretty, very useful

#Let's add some empty columns for taking inclusion/exclusion notes too

Excluded<-vector(mode="character", length=nrow(MyResults))
WhyExcluded<-vector(mode="character", length=nrow(MyResults))

MyResults_Annotated<-cbind.data.frame(MyResults_Annotated, Excluded, WhyExcluded)

#And then write out the results so that we can snoop through them in a spreadsheet program like Excel:
write.csv(MyResults_Annotated, "MyResults_Annotated.csv")

####################

#Sanity check:
#I'm curious as to how many of the results I get using just the antidepressant terms and taxa are completely missing organism part info:

result_MyQueryTerms_RatsMice<- gemma.R ::get_datasets(query=MyQueryTerms, taxa = c("mouse", "rat")) %>% 
  gemma.R:::get_all_pages() 

str(result_MyQueryTerms_RatsMice)
#Classes ‘data.table’ and 'data.frame':	111 obs. of  23 variables:

#... then filter down again to high quality data:

result_MyQueryTerms_RatsMice_Filtered<-result_MyQueryTerms_RatsMice[result_MyQueryTerms_RatsMice$experiment.troubled==FALSE,]

str(result_MyQueryTerms_RatsMice_Filtered)
#Classes ‘data.table’ and 'data.frame':	111 obs. of  23 variables:

MyResults<-result_MyQueryTerms_RatsMice_Filtered
#... and re-ran the code above. 
#Definitely time to functionalize...

write.csv(MyResults_Annotated, "MyResults_Annotated_AllOrganismParts.csv")

#almost all studies were annotated with either an organism part or cell type
#there are a handful of studies that aren't
#Those studies may be worth looking at more carefully if there aren't many datasets for a research question/region
#Sometimes they're missing that info because there isn't a publication or metadata though
#So probably not worth digging through for any research question that already has many datasets
