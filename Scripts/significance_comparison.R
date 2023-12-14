library(data.table)
library(readxl)
library(car)
library(ggplot2)
library(dplyr)
library(parallel)
library(ggtern)

SA1df <-read_excel("ALL_211_RILs/SA1.xls",col_names = FALSE)
SA2df <- read_excel("ALL_211_RILs/SA2.xls",col_names = FALSE)
SW2df <- read_excel("ALL_211_RILs/sw2.xls",col_names = FALSE)
SW1df <- read_excel("ALL_211_RILs/sw1.xls",col_names = FALSE)
###have to do the same thing to every code, this is repetitie but oh well
SA1<- data.table(t(SA1df[-c(1:5),]))
colnames(SA1)<-unlist(SA1[1,])
SA1<-data.table(sapply(SA1,as.numeric))
SA1<-data.table(genotype = as.numeric(SA1df[1,-c(1)]),treatment = "SA",rep="1",SA1)

SA2<- data.table(t(SA2df[-c(1:5),]))
colnames(SA2)<-unlist(SA2[1,])
SA2<-data.table(sapply(SA2,as.numeric))
SA2<-data.table(genotype = as.numeric(SA2df[1,-c(1)]),treatment = "SA",rep="2",SA2)


SW1<- data.table(t(SW1df[-c(1:5),]))
colnames(SW1)<-unlist(SW1[1,])
SW1<-data.table(sapply(SW1,as.numeric))
SW1<-data.table(genotype = as.numeric(SW1df[1,-c(1)]),treatment = "SW",rep="1",SW1)

SW2<- data.table(t(SW2df[-c(1:5),]))
colnames(SW2)<-unlist(SW2[1,])
SW2<-data.table(sapply(SW2,as.numeric))
SW2<-data.table(genotype = as.numeric(SW2df[1,-c(1)]),treatment = "SW",rep="2",SW2)
#bind them all together into one big table
batches<-data.table(rbindlist(list(SA1,SA2,SW1,SW2)))
batchestest<-fread("rils_pheno_alltrts_AGI.csv")

#sigtest on rils
clust<-makeCluster(4)
clusterExport(clust,c("batches"))
clusterEvalQ(clust,library(data.table)) 
start <- Sys.time()
lm_aov_sumsqs_rils<-rbindlist(parApply(cl=clust,batches[,-c(1:3),with = FALSE], 2, function(x){
  #run and anova on each gene according to batch number
  geneaov<-aov(x ~ as.factor(treatment)+as.factor(rep)+as.factor(genotype), data = batches)
  #make summary table of anova stats. we are looking at the sum of squares 
  #(sum of the sq of the variation), distance between each ind value and the mean (line of best fit) and then summed
  sumstats<-summary(geneaov)
  #this data object is a list of lists, so this command pulls out the sum of sqs batch & residuals
  sumsq<-sumstats[[1]][[2]]
  #table of sum of squares for each gene and its proportion, factor/residual
  sumtable<-data.table(treatmentv=sumsq[1],replicatev=sumsq[2],genotypev=sumsq[3],residuals=sumsq[4],
                       tprop=(sumsq[1]/sumsq[3]),rprop=(sumsq[2]/sumsq[4]),gprop=(sumsq[3]/sumsq[4]))
  ##im pretty confident this is wrong but also was 
  return(sumtable)
}))
print( Sys.time() - start )
stopCluster(clust)

clust<-makeCluster(4)
clusterExport(clust,c("batches"))
clusterEvalQ(clust,library(data.table)) 
start <- Sys.time()
#sigtest on parents
lm_aov_sumsqs_parents<-rbindlist(parApply(cl=clust,parentdtfinal[,-c("rep",'treatment','genotype'), with = F], 2, function(x){
  #run and anova on each gene according to batch number (batches, parentdtfinal)
  lm_genes<-lm(data = parentdtfinal,x ~ as.factor(treatment)+as.factor(rep)+
                 as.factor(genotype)+(as.factor(rep)/as.factor(treatment))+
                 (as.factor(genotype)*as.factor(treatment)))
  #make summary table of anova stats. we are looking at the sum of squares 
  #(sum of the sq of the variation), distance between each ind value and the mean (line of best fit) and then summed
  lm_aov1<-Anova(lm_genes,type="III")
  sumtable<-data.table(ss.treatment=lm_aov1[[1]][[2]],ss.replicate=lm_aov1[[1]][[3]],ss.genotype=lm_aov1[[1]][[4]],ss.replicate_n_treatment=lm_aov1[[1]][[5]],ss.treatment_i_genotype=lm_aov1[[1]][[6]],residualss=lm_aov1[[1]][[7]])
  #return sum of sq tables
  return(sumtable)
}))
print( Sys.time() - start )
stopCluster(clust)

write.csv(lm_aov_sumsqs,file = "Data/sigtests_output/lm_aov_sumsqs.rils.csv")
write.csv(lm_aov_sumsqs_parents,file = "Data/sigtests_output/lm_aov_sumsqs_parents.csv")
lm_aov_sumsqs_rils<-fread('sigtests_output/lm_aov_sumsqs.rils.csv')
lm_aov_sumsqs_rils$V1<-NULL
lm_aov_sumsqs_parents<-fread('sigtests_output/lm_aov_sumsqs_parents.csv')
lm_aov_sumsqs_parents$V1<-NULL

#Getting proportions
rilprops<-data.table()
lm_aov_sumsqs_rils$totalss<-lm_aov_sumsqs_rils$ss.treatment+lm_aov_sumsqs_rils$ss.genotype+lm_aov_sumsqs_rils$ss.treatment_i_genotype
rilprops$treatment<-lm_aov_sumsqs_rils$ss.treatment/lm_aov_sumsqs_rils$totalss
#rilprops$replicate<-lm_aov_sumsqs_rils$ss.replicate/lm_aov_sumsqs_rils$totalss
rilprops$genotype<-lm_aov_sumsqs_rils$ss.genotype/lm_aov_sumsqs_rils$totalss
#rilprops$rep.nested.treatment<-lm_aov_sumsqs_rils$ss.replicate_n_treatment/lm_aov_sumsqs_rils$totalss
rilprops$genotype.treatment.interaction<-lm_aov_sumsqs_rils$ss.treatment_i_genotype/lm_aov_sumsqs_rils$totalss
#rilprops$residuals<-lm_aov_sumsqs_rils$residualss/lm_aov_sumsqs_rils$totalss
rilprops$ID<-rilprops$ID<-1:nrow(rilprops)

parentprops<-data.table()
lm_aov_sumsqs_parents$totalss<-lm_aov_sumsqs_parents$ss.treatment+lm_aov_sumsqs_parents$ss.genotype+lm_aov_sumsqs_parents$ss.treatment_i_genotype
parentprops$treatment<-lm_aov_sumsqs_parents$ss.treatment/lm_aov_sumsqs_parents$totalss
#parentprops$replicate<-lm_aov_sumsqs_parents$ss.replicate/lm_aov_sumsqs_parents$totalss
parentprops$genotype<-lm_aov_sumsqs_parents$ss.genotype/lm_aov_sumsqs_parents$totalss
#parentprops$rep.nested.treatment<-lm_aov_sumsqs_parents$ss.replicate_n_treatment/lm_aov_sumsqs_parents$totalss
parentprops$genotype.treatment.interaction<-lm_aov_sumsqs_parents$ss.treatment_i_genotype/lm_aov_sumsqs_parents$totalss
#parentprops$residuals<-lm_aov_sumsqs_parents$residualss/lm_aov_sumsqs_parents$totalss
parentprops$ID<-1:nrow(parentprops)

qtls<-fread("finalpeaks_10cMwindow_june5.txt")
qtls<-unique(qtls[,.(name,type)])
result <- data.table(dcast(qtls, name ~ type, fun.aggregate = function(x) as.integer(length(x) > 0)))
cisgenes<-result[cis=='1'&trans=='0',name]
transgenes<-result[cis=='0'&trans=='1',name]
bothgenes<-result[cis=='1'&trans=='1',name]
pgenesdt<-fread("parents_pheno_alltrts_AGI.csv")
pgenes<-colnames(pgenesdt)[5:length(pgenesdt)]
rgenesdt<-fread("rils_pheno_alltrts_AGI.csv")
rgenes<-colnames(rgenes)[5:length(rgenes)]

parentprops
parentp<-ggtern(data=parentprops[,1:3],aes(treatment,genotype,genotype.treatment.interaction))+
  labs(x="",y="",z="",xarrow="Treatment",yarrow="Genotype",zarrow="Interaction")+
  stat_density_tern(aes(fill=..level..),geom = 'polygon',bins = 8,bdl = 0.08)+
  theme_light()+
  theme_arrowdefault()+
  theme(text = element_text(size=9),legend.position = 'none')

ggtern(data=rilprops[,1:3],aes(treatment,genotype,genotype.treatment.interaction))+
  labs(z="interaction")+
  labs(x="",y="",z="",xarrow="Treatment",yarrow="Genotype",zarrow="Interaction")+
  stat_density_tern(aes(fill=..level..),geom = 'polygon',bins = 100000,bdl = 0.00001)+
  theme_light()+
  theme_arrowdefault()+
  theme(text = element_text(size=9),legend.position = 'none')
library(patchwork)
parentp+rilp
dim(pheno)




