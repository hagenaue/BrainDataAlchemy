if(!requireNamespace("devtools", quietly=T)){
  install.packages("devtools")
}

devtools:: install_github("PavlidisLab/gemma.R", force=T)

#To change: Create a folder for your output on your computer and then set the working directory to your own folder (this one is Annaka’s):
setwd("C:/Users/annak/Documents/gemma datasets")

#For Amrita: This code may potentially change – if there are more results now, you may need to add another row of code with new offsets to capture additional results
##had to offset (skip specified num of objects found from database) to avoid the
##100 obj limit ex.) lps2 is the first 100 mouse lps studies, lps3 is the 101-200
##lps mouse studies and lps4 is the 201-229 mouse lps studies
##none of the other search terms hit the limit so they only have 1 and 2 (rat, mouse)
# 

lps_hipp1 <- gemma.R :: search_datasets("lps hippocamp*", taxon="rat", limit=100)
##1 row
lps_hipp2 <- gemma.R :: search_datasets("lps hippocamp*", taxon="mouse", limit=100)
#18 rows

lps_hippv2_1 <- gemma.R :: search_datasets("lipopolysaccharide hippocamp*", taxon= "rat", limit= 100)
#1 obs

lps_hippv2_2 <- gemma.R :: search_datasets("lipopolysaccharide hippocamp*", taxon= "mouse", limit= 100)
#17 obs

combined3<-rbind.data.frame(lps_hippv2_1, lps_hippv2_2, lps_hipp1, lps_hipp2)
#37 obs

unique_combined3 <- unique(combined3)
#20 obs

write.table(x=unique_combined3, file= "unique_data_Hippocampus2_20230512.txt", sep= "\t")

#This produces the same output as before, but with 3 new datasets that are inapplicable.

##########################

# Some of Annaka's original code:

# lps1 <- gemma.R :: search_datasets("LPS", taxon= "rat", limit= 100)
# str(lps1)
# #14 obs
# 
# lps2 <- gemma.R :: search_datasets("LPS", taxon= "mouse", limit= 100)
# #100 obs
# 
# lps3 <- gemma.R :: search_datasets("LPS", taxon= "mouse", offset=100,  limit= 100)
# #100 obs
# 
# lps4 <- gemma.R :: search_datasets("LPS", taxon= "mouse", offset=200,  limit= 100)
# #48 obs
# 
# lps_bacteria1 <- gemma.R :: search_datasets("bacterial infection", taxon= "rat", limit= 100)
# ##1 row
# lps_bacteria2 <- gemma.R :: search_datasets("bacterial infection", taxon="mouse", limit= 100)
# ##74 row
# #3 more than when Annaka ran it
# 
# lpsv2_1 <- gemma.R :: search_datasets("lipopolysaccharide", taxon= "rat", limit= 100)
# str(lps1)
# #12 obs
# 
# lpsv2_2 <- gemma.R :: search_datasets("lipopolysaccharide", taxon= "mouse", limit= 100)
# #100 obs
# 
# lpsv2_3 <- gemma.R :: search_datasets("lipopolysaccharide", taxon= "mouse", offset=100,  limit= 100)
# #100 obs
# 
# lpsv2_4 <- gemma.R :: search_datasets("lipopolysaccharide", taxon= "mouse", offset=200,  limit= 100)
# #20 obs
# 
# #This produced slightly fewer results than LPS
# 
# #This code will need to be updated to list the objects that you actually have in your code:
# combined<-rbind.data.frame(lps1, lps2, lps3, lps4, lps_bacteria1, lps_bacteria2, lpsv2_1, lpsv2_2, lpsv2_3, lpsv2_4)
# #569 obs
# 
# unique_combined <- unique(combined) 
# #343 rows
# 
# write.table(x=unique_combined, file= "unique_data_AllTissue_20230512.txt", sep= "\t")
# 
# 
# #lps_hipp1 <- gemma.R :: search_datasets("lps hippocampus", taxon="rat", limit=100)
# ##1 row
# #lps_hipp2 <- gemma.R :: search_datasets("lps hippocampus", taxon="mouse", limit=100)
# ##17 rows

# lps_ammon1 <- gemma.R :: search_datasets("lps Ammon's horn", taxon="rat", limit=100)
# #write.table(x=lps_hipp1, file= "lpshippocampussearch result1.txt", sep= "\t")
# ##1 row
# lps_ammon2 <- gemma.R :: search_datasets("lps Ammon's horn", taxon="mouse", limit=100)
# #10 rows
# 
# lps_ammonv2_1 <- gemma.R :: search_datasets("lipopolysaccharide Ammon's horn", taxon= "rat", limit= 100)
# #1 obs
# 
# lps_ammonv2_2 <- gemma.R :: search_datasets("lipopolysaccharide Ammon's horn", taxon= "mouse", limit= 100)
# #11 obs
# 
# combined4<-rbind.data.frame(lps_ammon1, lps_ammon2, lps_ammonv2_1, lps_ammonv2_2)
# #23 obs
# 
# unique_combined4 <- unique(combined4)
#12 obs
