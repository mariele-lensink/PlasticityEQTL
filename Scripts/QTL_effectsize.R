library(data.table)
library(parallel)
#read in qtl peak table
qtlssave<-fread("~/Desktop/finalpeaks_withmarkernames_10cMwindow_LODinfo.txt")
#read in map of genotypes for the rils
geno<-fread("Data/QTLAinput/BayxSha_geno.csv")
geno$V1<-NULL

#read in conversion table of transcripts, remove rows where AGI name is NA
probedf_full <-as.data.table(read.xlsx("Row Labels in Order.xlsx"))
probedf_full<-probedf_full[!is.na(probedf_full$AGI),]
colnames(probedf_full)[1]<-"Affy.probe.set.ID"
#read in phenotype files, transcript abundances for each treatment group
sa<-fread("BayxSha_SApheno.csv")
sa$V1<-NULL
sa<-data.table(treatment="SA",sa)
sw<-fread("SWqtl/BayxSha_SWpheno.csv")
sw$V1<-NULL
sw<-data.table(treatment="SW",sw)
delta<-fread("delta/BayxSha_deltapheno.csv")
delta$V1<-NULL
delta<-data.table(treatment="delta",delta)
pheno_alltrts<-rbind(sa,sw,delta)
#write them all to a csv so i can save this step and skip the previous ones from now on
write.csv(pheno_alltrts,"rils_pheno_alltrts.csv")
#change the names of the transcripts in the phenotype file to AGI names
colnames(pheno_alltrts)<-sub("^X","",colnames(pheno_alltrts))
subsetnames<-which(colnames(pheno_alltrts) %in% probedf_full$Affy.probe.set.ID)
dt<-as.data.frame(pheno_alltrts[,c(subsetnames),with = F])
newnames<-data.table(sapply(colnames(dt), function(x){
  loc<-which(probedf_full$Affy.probe.set.ID == x)
  probedf_full$AGI[loc]
}))
colnames(dt)<-c(unlist(newnames))

pheno_all<-data.table(pheno_alltrts[,1:2],dt)
#save this phenotype file with agi names to skip previous steps
write.csv(pheno_all,"rils_pheno_alltrts_AGI.csv")

#Use the lm() function in R to fit a linear model to the data, 
#where the dependent variable is the phenotype of interest and the 
#independent variable is the genotype at the QTL. 
#You can extract the genotype data from the qtl_results data table.

#Once you have fit the linear model, you can extract the effect size 
#of the QTL by examining the coefficient of the qtl_genotypes variable 
#in the model. The coefficient represents the change in the phenotype for 
#each additional copy of the QTL allele.

txdt<-as.data.table(read.csv("rils_pheno_alltrts_AGI.csv"))
txdt$X<-NULL
txdt<-txdt[ID!=417]
#edit genotype table so transcript names match and are all AGI
#edit geno names to match the rest of the input files
test<-sub("(At.+)\\.[0-9]+", "\\1",colnames(geno))
test<-test[-1]
cc<-fread("../gmap_conversion.txt")

newnames<-c(unlist(lapply(test, function(x){
 if (grepl("At",x)==FALSE) {
   loc<-which(cc$originalmarker == x)
   cc$marker[loc]
 } else {
   x
 }
})))
colnames(geno)<-c("ID",newnames)
#qtls, geno, txdt

clust<-makeCluster(4)
clusterExport(clust,c("qtls","txdt","geno"))
clusterEvalQ(clust,library(data.table)) 
start <- Sys.time()
effectoutput<-parApply(cl = clust, qtls[1:20,],1, function(x){
 dt<-data.table(transcript=txdt[treatment == x[7],get(x[1])],genotype=geno[,get(x[13])])
 lmout<-(lm(transcript ~ genotype,dt, na.action = "na.exclude"))
 summary(lmout)
})
print( Sys.time() - start )
stopCluster(clust)

saveRDS(effectoutput,"lm_summary")



