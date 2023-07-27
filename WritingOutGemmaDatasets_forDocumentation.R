#Writing out the Gemma Summarized Experiment data into Excel-friendly file format
#We should save a copy of this for our records for each dataset
#This is because we want our analysis to be reproducible
#But the next time someone runs our code to extract the dataset from Gemma, Gemma may have updated something
#So we need documentation of the current version of the data that we used for our analysis
#Another way to do this would be to save our Workspace (.Rdata)
#But that is less accessible to folks who don't work in R

str(rowData(SummarizedExperiment_Filtered[[1]]))

write.csv(rowData(SummarizedExperiment_Filtered[[1]]), "Example_Annotation.csv")

str(colData(SummarizedExperiment_Filtered[[1]]))

write.csv(colData(SummarizedExperiment_Filtered[[1]])[,-1], "Example_MetaData.csv")

write.csv(ExpressionData_Filtered, "ExpressionData_Filtered.csv")

