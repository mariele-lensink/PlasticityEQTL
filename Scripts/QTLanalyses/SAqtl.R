#load packages
library(qtl2)
library(ggplot2)
library(readxl)
library(Rcpp)
library(data.table)
library(abind)

###############FORMATTING INPUT FILES ####################################
SA1df <- read_excel("Data/originaldata/ALL_211_RILs/SA1.xls")
SA2df <- read_excel("Data/originaldata/ALL_211_RILs/SA2.xls")

###make phenotype tables for each salicylic acid exp###
#SA1
#remove extraneous infor from imported file
SA1<- data.frame(t(SA1df[-c(1:4),-1]))
#remove index and make numeric matrix (can do matrix math on it this way)
SA1<-data.table(sapply(SA1,as.numeric))
colnames(SA1)[]<- ""
SA1<-as.matrix(SA1)
#SA2
#deleting extraneous info from imported data
SA2<- data.frame(t(SA2df[-c(1:4),-1]))
#remove index and make numeric matrix (can do matrix math on it this way)
SA2<-data.table(sapply(SA2,as.numeric))
colnames(SA2)[]<- ""
SA2<-as.matrix(SA2)

###get average value between salicylic acid reps###
SAmean<-(SA1+SA2)/2
#turn into data table (easier to work with)
SApheno<-data.table(SAmean)

####format phenotype table for qtl analysis###
#getting the list of RILs to add to the mean SA pheno table
SA1temp<-data.frame(SA1df)
RILlist<-colnames(SA1df)
RILlist<-RILlist[-1]
#get list of transcripts from og data set to ad to mean SA pheno table
gene<-c(SA1temp[-c(1:4),1])
#assign transcript names to SA pheno table
colnames(SApheno)<-gene
#add an ID column with ril numbers
SApheno<-data.frame(ID = c(as.numeric(RILlist)),SApheno)

###had to find which ril is not in carolines files and remove it in mine (RIL 417), not in genotype file either###
SApheno[146,1]
#remove the row for RIL 417
SApheno<-SApheno[-146,]

###reorder the salicylic acid file by ril ###
SApheno.sortID<-SApheno[order(SApheno$ID),]
#make this a CSV
write.csv(SApheno.sortID,file = "BayxSha_SApheno.csv")
####################################################################################################
saphenoedit<-fread("SAqtl/BayxSha_SApheno.csv")
saphenoedit<-saphenoedit[,-1]
write.csv(saphenoedit,"SAqtl/BayxSha_SApheno.csv")
SApheno<-fread("Data/QTLAinput/BayxSha_SApheno.csv")
SApheno$V1<-NULL
####################################################################################################


###edit genotype file to match salicylic acid phenotype file###
##this can be reused for all of the other analysis
#get rid of "ril_" text on each ril number, easier to have everything simply numeric
geno<- read.csv("Data/QTLAinput/BayxSha_geno.csv")
geno$X<-NULL

#double check to make sure pheno and geno are the same
SApheno$ID[!(SApheno$ID %in% geno$ID)]
geno$ID[!(geno$ID %in% SApheno$ID)]
#make genotype file a csv
###in other qtl analyses, you dont have to edit this now###
write.csv(geno,file = "BayxSha_geno.csv")

#########################################################################
#########QTL ANALYSIS########################

setwd("/Users/marielelensink/Documents-local/CisTransPlasticity2/Data/QTLAinput/")
bayxsha_SA<-read_cross2("BayxSha_SA.yaml",quiet = FALSE)

#calculate conditional genotype probabilities given the marker data at each putative QTL position
##insert pseudomarkers into the genetic map
map_SA<-insert_pseudomarkers(bayxsha_SA$gmap,step=0)
#calculate genotype probabilities
pr_SA<- calc_genoprob(bayxsha_SA,map_SA,error_prob = 0.002)


####run the genome scan####

#scan1() takes as input the genotype probabilities, a matrix of phenotypes, 
#optional additive and interactive covariates, and the special X chromosome covariates.
out_SA <- scan1(pr_SA, bayxsha_SA$pheno)
#output is LOD scores.
#LOD stands for "logarithm of the odds." 
# positions X phenotypes

#find peaks
peaks_SA<-find_peaks(out_SA, map_SA, threshold=2, drop=1.5)
png("../../Figures/SA_freq_tx_per_qtl.png")
hist(table(peaks_SA$lodcolumn),
     main = "Frequency of transcripts per QTL in\nSalicylic Acid Treatment Group",
     xlab = "#of transcripts per peak",
     ylab = "frequency of QTLs controlling varying #'s of transcripts")
dev.off()
fwrite(peaks_SA,"../QTLAoutput/SA_peaks_june5.txt")

