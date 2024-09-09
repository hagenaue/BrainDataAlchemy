#Sample code for MetaData extraction for Gemma datasets
#Megan Hagenauer, Apr 23 2024

#####################

#General notes: The Gemma Database

#Gemma is a database of curated and re-analyzed gene expression studies
#It currently includes >18,000 re-analyzed gene expression studies (mostly from the Gene Expression Omnibus repository)
#The Gemma database can be accessed via a website point-and-click GUI (Graphical User Interface) here:
#https://gemma.msl.ubc.ca/home.html

#For more details and background, this is a recent paper overviewing the Gemma database:
#https://pubmed.ncbi.nlm.nih.gov/33599246/

#The datasets within Gemma can be searched and browsed via a point-and-click GUI on the website:
#https://gemma.msl.ubc.ca/browse/#/

# We are going to use R coding to conduct formal, replicable searches of Gemma's datasets using their API
  #An API is an "application programming interface"
  #It allows developers to integrate data, services and capabilities from other applications
  
# To do this, we will use the R wrapper that Gemma has made for their restful API.
  #Here's the github site and documentation for it:
# https://github.com/PavlidisLab/gemma.R
# https://pavlidislab.github.io/gemma.R/articles/gemma.R.html

  #Here is the reference manual on Bioconductor:
 #https://bioconductor.org/packages/release/bioc/manuals/gemma.R/man/gemma.R.pdf


  #Add PRISMA diagram code???
  
####################
  
  
#1. Installing necessary code packages (this only needs to be done once)
  
#Installing the Gemma API package from Bioconductor:
  #this installation code requires us to use another package to install it, "BiocManager":
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gemma.R", force=T)


#Another version of the installation code: 
#This code is to install the developers version of gemma.R package from Github
#This version will be more up-to-date (i.e., bugs fixed)
#...but may also no longer completely match the reference manual/vignettes provided on Bioconductor 

#this version of the installation code requires us to use another package to install it, "Devtools"
if(!requireNamespace("devtools", quietly=T)){
  install.packages("devtools")
}

devtools:: install_github("PavlidisLab/gemma.R", force=T)


#Installing other code packages that we will use:

# plyr (a code package for working with databases):
install.packages("plyr")

# dplyr (a code package for working with databases):

# The easiest way to get dplyr is to install the whole tidyverse:
install.packages("tidyverse")

# Alternatively, install just dplyr:
install.packages("dplyr")

####################

#2. Loading the necessary code packages:

library(gemma.R)
library(plyr)
library(dplyr)


#########################

#5. Basic query: Searching for datasets that have metadata that includes a particular term:

#About the function:
help("gemma.R ::get_datasets")
??gemma.R ::get_datasets

# The first argument in this function is query
# Inputting basic text to this argument means that a search will be performed of the entire Gemma record for that exact text
# i.e., it will search the dataset title, description, annotations, and all of the other fields

#For example, this code pulls up a data.frame with all of the datasets that include "hippocampus" in their Gemma record somewhere 
AllHippocampusDatasets<- gemma.R ::get_datasets(query ="hippocampus")

#We can quickly overview our table of results using the "structure" function:
str(AllHippocampusDatasets)
#Classes ‘data.table’ and 'data.frame':	20 obs. of  23 variables:

#You can also peruse the results in Rstudio by going to the Environment tab in the righthand corner and clicking on the object

#Note that these searches *default to being limited to 20 results*.
#To get the full set of results you have to add some code:

AllHippocampusDatasets<- gemma.R ::get_datasets(query ="hippocampus") %>% 
  gemma.R:::get_all_pages()
#If broad, the search may take a little bit of time (e.g., 1 minute)

str(AllHippocampusDatasets)
#Classes ‘data.table’ and 'data.frame':	1570 obs. of  23 variables:

#Much better!

#Note: This search would not capture any variation of the term hippocampus, only exact matches
#So it would miss records that include hippocampal formation, hippocampal tissue, hippocampi, etc

#To search for an exact phrase we use quotations "" around the phrase
#... which is a little awkward because the search is already in quotations
# so we add a backslash before the inner quotation marks... 

HippocampalFormationDatasets<- gemma.R ::get_datasets(query ="\"hippocampal formation\"") %>% 
  gemma.R:::get_all_pages()

str(HippocampalFormationDatasets)
#Classes ‘data.table’ and 'data.frame':	1260 obs. of  23 variables:

##########################

#6. More complex query: Searching for datasets that have metadata that includes a variety of terms:

#Ideally to capture all datasets with some version of hippocampus, we need a more sophisticated search.

#More detailed instructions for how to query
# https://gemma.msl.ubc.ca/resources/restapidocs/index.html
#Look up "QueryArg" on this page

# To create a more sophisticated search, we often use Boolean operators
# Boolean operators are used to express mathematical sets and database logic.
# They connect your search words together to either narrow or broaden your set of results
# This website also has a nice, brief introducton to Boolean operators:
# https://libguides.mit.edu/c.php?g=175963&p=1158594
# It includes some really nice visuals to illustrate AND, OR, and NOT on the website. 

#Here are some important points of interest:

# OR can be used to combine terms:

AllHippocampOR_Datasets<- gemma.R::get_datasets(query="hippocampus OR hippocampal OR hippocampi") %>% 
  gemma.R:::get_all_pages()

str(AllHippocampOR_Datasets)
#Classes ‘data.table’ and 'data.frame':	1721 obs. of  23 variables:

#So we captured more datasets with this search than just with "hippocampus"

#Note: These are all unique datasets:
length(unique(AllHippocampOR_Datasets$experiment.shortName))

#Just listing a variety of terms seems to function similarly to OR:
AllHippocampOR_Datasets<- gemma.R::get_datasets(query="hippocampus hippocampal hippocampi") %>% 
  gemma.R:::get_all_pages()

str(AllHippocampOR_Datasets)
#Classes ‘data.table’ and 'data.frame':	1721 obs. of  23 variables:


# AND is used to specify only Gemma records that include *both* terms
# e.g., if you search for "dentate gyrus" you get any record that includes both dentate or gyrus
# but if you want records that include dentate gyrus together, you need to combine them with and:

AllHippocampusDG_Datasets<- gemma.R::get_datasets(query="(dentate AND gyrus)") %>% 
  gemma.R:::get_all_pages()

str(AllHippocampusDG_Datasets)
#Classes ‘data.table’ and 'data.frame':	206 obs. of  23 variables:

#Example of combining multiple terms:
AllHippocampusDG_Datasets<- gemma.R::get_datasets(query="(hippocampus) OR (dentate AND gyrus)") %>% 
  gemma.R:::get_all_pages()

str(AllHippocampusDG_Datasets)
#Classes ‘data.table’ and 'data.frame':	1588 obs. of  23 variables:

#NOT is used to exclude records that include particular terms:
AllHippocampusNotDG_Datasets<- gemma.R::get_datasets(query="(hippocampus) NOT (dentate AND gyrus)") %>% 
  gemma.R:::get_all_pages()

str(AllHippocampusNotDG_Datasets)
#Classes ‘data.table’ and 'data.frame':	1549 obs. of  23 variables:

# * can supposedly be used to capture all variations of a word that starts in a particular way
#e.g., hippocamp* should capture hippocampus, hippocampal, hippocampi

AllHippocamp_Datasets<- gemma.R::get_datasets(query="hippocamp*") %>% 
  gemma.R:::get_all_pages()

str(AllHippocamp_Datasets)
#Classes ‘data.table’ and 'data.frame':	1721 obs. of  23 variables:
#The same number of records pulled when using "hippocampus hippocampal hippocampi"

#Note: In earlier versions of the Gemma.R package, "*" wasn't working properly, but Gemma fixed the bug.

#For contrasting with our earlier exact phrase code for hippocampal formation
#note that "hippocampal AND formation" do not produce the same search results:
HippocampalAndFormationDatasets<- gemma.R ::get_datasets(query ="hippocampal AND formation") %>% 
  gemma.R:::get_all_pages()

str(HippocampalAndFormationDatasets)
#Classes ‘data.table’ and 'data.frame':	1289 obs. of  23 variables:


##########################

#7. Leveling up: Leveraging ontologies instead of just search terms

#The searches described above suffer from two problems:

  #A. They require you to brainstorm every possible variation on "hippocampus" that could be used in the title, description, or metadata for a dataset
    #e.g., Ammon's horn, dentate gyrus, CA1, etc

  #B. The datasets that you discover may include the term hippocampus in their Gemma record but not actually reflect an experiment studying hippocampal tissue
    #e.g., In the abstract, the authors might describe how they are following up on previous investigations in the hippocampus by studying the amygdala...

#To get around these problems, we can make use of the fact that the datasets in Gemma are annotated using structured terms
#This annotation references formal, existing ontologies 
# Ontologies: "a set of concepts and categories in a subject area or domain that shows their properties and the relations between them."

#You can get a glimpse of the ontological terms used to annotate the datasets using the point-and-click dataset search GUI on the Gemma website:
# https://gemma.msl.ubc.ca/browse/#/
  #To see the ontological terms: after running a search, at the bottom of the screen is "Dataset download code"
  #Under this download code tab, the option "Gemma.R" provides the R code for the search.
  
#Many of the ontological terms are hierarchical
  #e.g., dentate gyrus is part of the hippocampus
  
#Within an ontological hierarchy:
  #the broader, umbrella term is called the "parent" (in the previous example, the hippocampus is the parent)
  #the term representing the more specific subset is called the "child" (in the previous example, the dentate gyrus is the child)
  
# Why hierarchical ontology is useful:
# When requesting Gemma records with particular annotation, a parent term can often pull up all records with the child terms

#This is an easy way to browse the hierarchies for the ontological terms:
# https://www.ebi.ac.uk/ols4
# Gemma considers the children of a term to be the terms listed under "is a"/"subclass of" and "part of"

# Here are several categories of terms that we regularly reference

#Tissues (Organism Parts):
# https://www.ebi.ac.uk/ols4/ontologies/uberon

#Experimental Factors:
# https://www.ebi.ac.uk/efo/
  
#Drugs - note: Gemma does *not* infer children from parent terms within this ontological framework:
# https://www.ebi.ac.uk/ols4/ontologies/chebi

###############################

#8. Example of searching for datasets using particular ontological annotations instead of just search terms:

#We can search for annotations related to hippocampus using the search_annotations function
#e.g., 
HippocampusAnnotations<-gemma.R ::search_annotations("hippocampus")
#Better:
HippocampusAnnotations<-gemma.R ::search_annotations("hippocamp*")
#for whatever reason, in this search the * does seem to work...

#To view the results you can click on the HippocampusAnnotations object in the Environment tab 

#That search provides ontological annotations that include terms that start with hippocamp*
#if we want to determine which of these annotations are actually used in the Gemma database, we can run this code to count their usage:

annots<-HippocampusAnnotations

annot_counts <- annots$value.URI %>% sapply(\(x){
  attributes(get_datasets(uris = x))$totalElements
})
annots$counts = annot_counts

#Then we can filter down to just the ontological terms that are actually used in the Gemma database:
annots_wRecords<-annots %>% filter(counts>0)

#You'll note that one of the frequently used terms is an "organism part" hippocampal formation
#That implies that this annotation is used to characterize the actual tissue used in the experiment
#It has the annotation ontological URI: http://purl.obolibrary.org/obo/UBERON_0002421

#If we use "UBERON_0002421" (hippocampal formation) in our search for datasets
# we will get not just results for hippocampal formation, but results matching its child terms as well

#If we wanted to see the hierarchy of ontological terms encompassed by "hippocampal formation" and/or encompassing "hippocampal formation"
#We can go to EMBL-EBI Ontology Lookup Service:
# https://www.ebi.ac.uk/ols4
#And put in the search box UBERON_0002421

#EMBL-EBI will tell us all of the various terms/concepts encompassed by "hippocampal formation" in the "Related from" section under "part of"
#e.g., layer of hippocampus, subiculum, fornix, Ammon's horn, etc

#Whereas to see if "hippocampal formation" can be broadened to include other useful nearby areas, we can click above it in the hierarchy ("cerebral cortex") to see our options

#We can also double check that those child terms for "hippocampal formation" are actually used within Gemma using this function:
get_child_terms("http://purl.obolibrary.org/obo/UBERON_0002421")

#We could search for all datasets annotated with this ontological term (or any of its child terms) using:

AllHippocampalDatasets_byURI<- gemma.R::get_datasets(uris="http://purl.obolibrary.org/obo/UBERON_0002421") %>% 
  gemma.R:::get_all_pages()

str(AllHippocampalDatasets_byURI)
#Classes ‘data.table’ and 'data.frame':	1254 obs. of  23 variables:

#The code can also be written like this:
AllHippocampalDatasets_byURI <- get_datasets(filter = 'allCharacteristics.valueUri in (http://purl.obolibrary.org/obo/UBERON_0002421)') %>% 
  gemma.R:::get_all_pages()

str(AllHippocampalDatasets_byURI)
#Classes ‘data.table’ and 'data.frame':	1254 obs. of  23 variables:

#Here is an example of datasets that are annotated with a term ("ammon's horn") that is a child to the broader hippocampal formation term

AmmonAnnotations<-gemma.R ::search_annotations("ammon's horn")

AllAmmonDatasets_byURI <- get_datasets(filter = 'allCharacteristics.valueUri in (http://purl.obolibrary.org/obo/UBERON_0001954)') %>% 
  gemma.R:::get_all_pages()

str(AllAmmonDatasets_byURI)
#Classes ‘data.table’ and 'data.frame':	1085 obs. of  23 variables:

#All of those datasets were also captured by the hippocampal formation term:

sum(AllAmmonDatasets_byURI$experiment.shortName %in% AllHippocampalDatasets_byURI$experiment.shortName)
#[1] 1085


#What if we want to conduct searches using a combination of ontological terms?

#If we wanted both hippocampus and amygdala:

AmygdalaAnnotations<-gemma.R ::search_annotations("amygdala")

#URI for amygdala as a biosource:
http://purl.obolibrary.org/obo/UBERON_0001876

AllHippocampalOrAmygdalaDatasets_byURI<- gemma.R::get_datasets(uris=("http://purl.obolibrary.org/obo/UBERON_0002421, http://purl.obolibrary.org/obo/UBERON_0001876")) %>% 
  gemma.R:::get_all_pages()

str(AllHippocampalOrAmygdalaDatasets_byURI)
#Classes ‘data.table’ and 'data.frame':	1356 obs. of  23 variables:

#The code can also be written like this:
AllHippocampalOrAmygdalaDatasets_byURI <- get_datasets(filter = 'allCharacteristics.valueUri in (http://purl.obolibrary.org/obo/UBERON_0002421, http://purl.obolibrary.org/obo/UBERON_0001876)') %>% 
  gemma.R:::get_all_pages()

str(AllHippocampalOrAmygdalaDatasets_byURI)
#Classes ‘data.table’ and 'data.frame':	1307 obs. of  23 variables:


########################################

#9. Another example: When searching for datasets using particular ontological annotations instead of search terms doesn't work quite as well...

#Example: searching for antidepressant studies

#When looking for records about the effects of antidepressant drugs, a query might give us any dataset that mentions "antidepressant-like effects" of an intervention in the abstract
#...Even if the intervention discussed is actually a genetic knock-out or behavioral intervention.
#In contrast, searching for records annotated with the ontological terms for specific antidepressants will only pull up records that use antidepressants in the experiment.

#We want to figure out all of the relevant annotations that might accompany our datasets of interest

AntidepressantAnnotations<-gemma.R ::search_annotations("antidepressant")

#This ontological term is listed under the category "treatment": http://purl.obolibrary.org/obo/CHEBI_35469

#When I look it up on the EMBL-EBI Ontology Lookup Service:
# https://www.ebi.ac.uk/ols4
#I get:
#  https://www.ebi.ac.uk/ols4/ontologies/chebi/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FCHEBI_35469
  #It tells me that it is a child term to "psychotropic drug" and lists a whole bunch of different antidepressants as "has role"
  #This is a much larger list than we generated from the literature during previous summers - nice!

  AllAntidepressantDatasets<-  gemma.R::get_datasets(filter="allCharacteristics.valueUri = http://purl.obolibrary.org/obo/CHEBI_35469") %>% 
    gemma.R:::get_all_pages()
  #but only one dataset has this annotation!
  
  #Likewise, if we go back and count the usage of the antidepressant annotations in Gemma:
  annots<-AntidepressantAnnotations
  
  annot_counts <- annots$value.URI %>% sapply(\(x){
    attributes(get_datasets(uris = x))$totalElements
  })
  annots$counts = annot_counts
  #... and then click on annots in the global environment to view it
  #There are almost no matches in the Gemma database!
  
  #And when we double-check in Gemma what "antidepressants" (CHEBI_35469) has as its children, we find...
  get_child_terms("http://purl.obolibrary.org/obo/CHEBI_35469")
  #character(0)
  #There are no children!
  
  #It turns ou that this is because of two reasons:
  #A. Gemma considers the children of a term to be the terms listed under "is a"/"subclass of" and "part of" - the specific antidepressants listed with CHEBI_35469 only "have role"
  #B. Gemma doesn't use the hierarchical information for drugs/molecules in the CHEBI ontology
  #Darn.
  
#Instead, we have to look up the individual antidepressants
  # To figure out which antidepressants to search for, we dig through the drugs listed as "have role" in CHEBI_35469 
  # We should also reference the literature
 #To make this go faster, we can list more than one when searching for annotations 

  #version of the search with all the antidepressants used by Erin Hernandez in her meta-analysis:
  #note: this code isn't happy with line breaks instead of spaces between drug names
  AntidepressantAnnotations_byName<-gemma.R ::search_annotations("SSRI fluoxetine sertraline paroxetine citalopram escitalopram fluvoxamine vilazodone vortioxetine amitriptyline imipramine amoxapine desipramine nortriptyline clomipramine trimipramine protriptyline doxepin maprotiline trazadone nefazodone mirtazapine phenelzine nialamide isocarboxazid hydracarbazine tranylcypromine selegiline venlafaxine desvenlafaxine duloxetine levomilnacipran bupropion duloxetine levomilnacipran ketamine esketamine tianeptine brexanolone")
    
  #Let's recycle that code from earlier for figuring out which of the identified annotations actually have associated records:
  annots<-AntidepressantAnnotations_byName
  
  annot_counts <- annots$value.URI %>% sapply(\(x){
    attributes(get_datasets(uris = x))$totalElements
  })
  annots$counts = annot_counts
  
  #Then we can filter down to just the ontological terms that are actually used in the Gemma database:
  annots_wRecords<-annots %>% filter(counts>0)
  
  #Some of these rows aren't actually useful to us:
  
  annots_wRecords_Good<-annots_wRecords[annots_wRecords$value.name!="major depressive disorder",]
  
  AntidepressantsByNameURIDatasets<- gemma.R::get_datasets(uris=annots_wRecords_Good$value.URI) %>% 
    gemma.R:::get_all_pages()
  
  str(AntidepressantsByNameURIDatasets)
  #Classes ‘data.table’ and 'data.frame':	77 obs. of  23 variables:
  
  ###########
