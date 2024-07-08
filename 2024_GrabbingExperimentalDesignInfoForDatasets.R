
#These ids can also be read in from a .csv file using read.csv()
ExperimentIDs<-c("GSE135306", "GSE81672", "GSE180465", "GSE179667", "GSE146358", "GSE128255", "GSE187418")

FactorInfo<-vector(mode="character", length=length(ExperimentIDs))
BaselineFactorInfo<-vector(mode="character", length=length(ExperimentIDs))
TreatmentFactorInfo<-vector(mode="character", length=length(ExperimentIDs))

FactorInfoDF<-cbind.data.frame(ExperimentIDs, FactorInfo, TreatmentFactorInfo, BaselineFactorInfo)

for(i in c(1:length(ExperimentIDs))){
  Design<-gemma.R::get_dataset_differential_expression_analyses(ExperimentIDs[i])
  #print(ExperimentIDs[i])
  #print(paste(Design$factor.category, collapse="; "))
  FactorInfoDF[i,2]<-paste(Design$factor.category, collapse="; ")
  
  ExperimentalFactorVector<-vector(mode="character", length=1)
  
  for(j in c(1:length(Design$experimental.factors))){
    ExperimentalFactorVector<-c(ExperimentalFactorVector,Design$experimental.factors[[j]]$summary)
  }
  
  BaselineFactorVector<-vector(mode="character", length=1)
  
  for(k in c(1:length(Design$baseline.factors))){
    BaselineFactorVector<-c(BaselineFactorVector,Design$baseline.factors[[k]]$summary)
  }
  
  FactorInfoDF[i,3]<-paste(ExperimentalFactorVector, collapse="; ")
  FactorInfoDF[i,4]<-paste(BaselineFactorVector, collapse="; ")
  
  rm(Design, ExperimentalFactorVector, BaselineFactorVector)
  
}
