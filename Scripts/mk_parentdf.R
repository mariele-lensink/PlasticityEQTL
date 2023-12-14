library(data.table)
library(qtl2)
library(openxlsx)
library(readxl)

geno<- read.csv("BayxSha_geno.csv")
geno$ID<-as.numeric(gsub('RIL_([0-9]+)','\\1',geno$ID))
geno$X<-NULL
pheno<-fread("BayxSha_SApheno.csv")
pheno$V1<-NULL
colnames(pheno)<-sub("^X","",colnames(pheno))
probedf_full <-as.data.table(read.xlsx("../Row Labels in Order.xlsx"))
probedf_full<-probedf_full[!is.na(probedf_full$AGI),]
phenosubset<-pheno[,]

subsetnames<-which(colnames(pheno) %in% probedf_full$`›`)
subsetdt<-as.data.frame(pheno[,c(subsetnames), with = F])

newnames<-data.table(sapply(colnames(subsetdt), function(x){
  loc<-which(probedf_full$`›` == x)
  probedf_full$AGI[loc]
}))

colnames(subsetdt)<-unlist(newnames)
subsetdt<-as.data.table(subsetdt)
subsetdt<-cbind(pheno$ID,subsetdt)
colnames(subsetdt)[1]<-"ID"
SAsubsetdt<-data.table(treatment = "SA",subsetdt)
################################SW##########################################
swpheno<-as.data.table(read.csv("../SWqtl/BayxSha_SWpheno.csv"))
subsetdt<-as.data.frame(swpheno[,c(subsetnames), with = F])
colnames(subsetdt)<-unlist(newnames)
subsetdt<-as.data.table(subsetdt)
subsetdt<-cbind(pheno$ID,subsetdt)
colnames(subsetdt)[1]<-"ID"
SWsubsetdt<-data.table(treatment = "SW",subsetdt)
################################SW#########################################
SA1<-t(as.data.table(read_excel("ALL_211_RILs/SA1.xls",sheet = 2)))
SA2<-t(as.data.table(read_excel("ALL_211_RILs/SA2.xls",sheet = 2)))
SW1<-t(as.data.table(read_excel("ALL_211_RILs/sw1.xls",sheet = 2)))
SW2<-t(as.data.table(read_excel("ALL_211_RILs/sw2.xls",sheet = 2)))
dts<-list(SA1,SA2,SW1,SW2)
parentdt<-lapply(dts, function(x){
  ids<-rownames(x)[-1]
  ids<-sapply(ids, function(x){
    x<-substr(x,1,3)
  })
  repnum<-as.character(x[2,1])
  trt<-as.character(x[2,2])
  dt<-as.data.frame(x[,-c(1:5)])
  colnames(dt)<-dt[1,]
  dt<-dt[-1,]
  dtfinal<-data.table(genotype=ids,treatment=trt,rep=repnum,dt)
  
})
pSA1<-parentdt[[1]]
pSA2<-parentdt[[2]]
pSW1<-parentdt[[3]]
pSW2<-parentdt[[4]]

parentdtfinal<-rbindlist(parentdt)
subsetnames<-which(colnames(parentdtfinal) %in% probedf_full$Affy.probe.set.ID)
dt<-as.data.frame(parentdtfinal[,c(subsetnames), with = F])
newnames<-data.table(sapply(colnames(dt), function(x){
  loc<-which(probedf_full$Affy.probe.set.ID == x)
  probedf_full$AGI[loc]
}))

colnames(dt)<-c(unlist(newnames))
parentdtfinal<-cbind(parentdtfinal[,1:3],dt)
write.csv(parentdtfinal,"parents_pheno_alltrts_AGI.csv")


