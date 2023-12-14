library(data.table)
library(Hmisc)
library(readxl)
library(gdata)
library(ggplot2)
library(ggpubr)
#uploading data
probedf_full <-data.frame(read_xlsx("Row Labels in Order.xlsx"))
parents<-fread("/Users/marielelensink/Documents-local/SalicylicAcid/ALL_211_RILs/Parental_averages.csv", colClasses = c("character","numeric","numeric","numeric","numeric"))
# make the delta value
parents$Bay0.delta<-((parents$Bay0.SA-parents$Bay0.SW)/((parents$Bay0.SA+parents$Bay0.SW)/2))
parents$Sha.delta<-((parents$Sha.SA-parents$Sha.SW)/((parents$Sha.SA+parents$Sha.SW)/2))
#save parent names to add
parentnames<-colnames(parents)
#transpose parent data table to match rils
parents<-t(parents)
#making the transcript names the column names
colnames(parents)<-parents[1,]
#getting rid of transcript names in first row
parents<-parents[-1,]
parents<-data.table(parents)
#add parent names back
parents<- data.table(ID = c(parentnames[-1]), parents)
################################################################################
#read in the expression data for the rils
SAexp<-fread("SAqtl/BayxSha_SApheno.csv", colClasses = "numeric")
SAexp$V1<-NULL
SWexp<-fread("SWqtl/BayxSha_SWpheno.csv", colClasses = "numeric")
SWexp$V1<-NULL
#make names of transcripts identical to parents
colnames(SAexp)<- sub("^X", "",colnames(SAexp))
colnames(SWexp)<- sub("^X", "",colnames(SWexp))
#some of the transcripts (64 of them) do not have tair names in the conversiontable
#we remove them but save them just in case
SAexpsubset_noconversion<-SAexp[,1:65]
SAexpsubset<-SAexp[,c(1,66:22810)]
SWexpsubset_noconversion<-SWexp[,1:65]
SWexpsubset<-SWexp[,c(1,66:22810)]
parentssubset_noconversion<-parents[,1:65]
parentssubset<-parents[,c(1,66:22810)]
#remove ril 417 because it has no genotype data and not used in qtl analysis
which(!ids %in% genotypes$ID)]
SAexpsubset<-SAexpsubset[ID != 417]
SWexpsubset<-SWexpsubset[ID != 417]
##############################################################################
#log10+1 transformation to all of the rils
SA_log10_plus1<-data.table(cbind(SAexpsubset[,1],log10(SAexpsubset[,2:22746]+1)))
SW_log10_plus1<-data.table(cbind(SWexpsubset[,1],log10(SWexpsubset[,2:22746]+1)))
#deltat<-(SAt-SWt)/((SAt+SWt)/2)

####################################################################################
#get stats
SA_stats<-data.table(t(apply(SA_log10_plus1[,2:22746],1,range)))
colnames(SA_stats)<-c("min","max")
SA_stats$ID<-"SA_RIL"
SW_stats<-data.table(t(apply(SW_log10_plus1[,2:22746],1,range)))
colnames(SW_stats)<-c("min","max")
SW_stats$ID<-"SW_RIL"
#delta_stats<-apply(deltat,1,range)
############
parentssubset<-data.table(cbind(parentssubset[,1],apply(parentssubset[,2:22746],2,as.numeric)))
p_stats<-data.table(cbind(parentssubset[,1],t(apply(parentssubset[,2:22746],1,range))))
colnames(p_stats)<-c("ID","min","max")
#remove delta for now
pstats<-p_stats[1:4,]
pstats$ID<-c("SA_parental","SA_parental","SW_parental","SW_parental")
####################################################################################
ranges<-rbind(pstats,SA_stats,SW_stats)
rangemelt<-melt(ranges,id.vars = 'ID')
ggplot(rangemelt,aes(x = as.factor(ID),y = value,color = variable)) +geom_jitter()+ggtitle("log10+1") 

comparedt<-data.table(t(rbind(parentssubset[3,2:1001],SA_log10_plus1[4,2:1001])))
colnames(comparedt)<-c('parent_SA','SA')
ggplot(comparedt,aes(parent_SA,SA))+geom_point()+
  ggtitle("comparing ril line (log10 transformation) to parent line")+
  stat_smooth(method = 'lm',formula = y~x,geom = 'smooth')+
  stat_cor(hjust = -4)

  

#ranges<-data.table(cbindX(prange,SArange,SWrange))
#ggplot(melt(ranges),aes(variable,value))+geom_boxplot()+geom_jitter()+ggtitle("log10")









#swap names of probe ids in peak table for transcript names (so they match the gmap)
transcriptnames<-sapply(colnames(SAexpsubset), function(x){
  loc<-which(probedf_full$Affy.probe.set.ID == x)
  probedf_full$AGI[loc]
})

ids<-ids[-208]
colnames(SAexpsubset)<-transcriptnames
SAexpsubset$ID<-ids
colnames(SWexpsubset)<-transcriptnames
SWexpsubset$ID<-ids
colnames(dexpsubset)<-transcriptnames
dexpsubset$ID<-ids


range(as.numeric(SAexpsubset[2,]))
range(log10(as.numeric(SAexpsubset[3,])))


###############REDO FOUND THE DATA#############################
SAexp<-fread("SAqtl/BayxSha_SApheno.csv", colClasses = "numeric")
SAexp$V1<-NULL
SWexp<-fread("SWqtl/BayxSha_SWpheno.csv", colClasses = "numeric")
SWexp$V1<-NULL

SAparents1<-as.data.table(read_excel("ALL_211_RILs/SA1.xls",sheet = 2),)
SAparents1<-SAparents1[,lapply(SAparents1[5:22814,-1],as.numeric)]
SAparents$Bay0_1<-(SAparents1$`Bay-0...6`+SAparents1$`Bay-0...7`)/2

###########delta averages
d<-read.csv("delta/BayxSha_deltapheno.csv")
d<-as.data.table(d)
d$X<-NULL
d$ID<-NULL
delta_avg_rils<-data.frame(apply(d,2,mean))
write.csv(delta_avg_rils,"delta_average_rils.csv")
