#load packages
library(qtl2)
library(readxl)
library(Rcpp)
library(data.table)
library(ggplot2)


###############FORMATTING INPUT FILES ####################################
SW1df <- read_excel("../originaldata/ALL_211_RILs/sw1.xls")
SW2df <- read_excel("../originaldata/ALL_211_RILs/sw1.xls")

###make phenotype tables for each salicylic acid exp###
#SW1
#remove extraneous infor from imported file
SW1<- data.frame(t(SW1df[-c(1:4),-1]))
#remove index and make numeric matrix (can do matrix math on it this way)
SW1<-data.table(sapply(SW1,as.numeric))
colnames(SW1)[]<- ""
SW1<-as.matrix(SW1)
#SW2
#deleting extraneous info from imported data
SW2<- data.frame(t(SW2df[-c(1:4),-1]))
#remove index and make numeric matrix (can do matrix math on it this way)
SW2<-data.table(sapply(SW2,as.numeric))
colnames(SW2)[]<- ""
SW2<-as.matrix(SW2)

###get average value between salicylic acid reps###
SWmean<-(SW1+SW2)/2
#turn into data table (easier to work with)
SWpheno<-data.table(SWmean)

####format phenotype table for qtl analysis###
#getting the list of RILs to add to the mean SA pheno table
SW1temp<-data.frame(SW1df)
RILlist<-colnames(SW1df)
RILlist<-RILlist[-1]
#get list of transcripts from og data set to ad to mean SA pheno table
gene<-c(SW1temp[-c(1:4),1])
#assign transcript names to SA pheno table
colnames(SWpheno)<-gene
#add an ID column with ril numbers
SWpheno<-data.frame(ID = c(as.numeric(RILlist)),SWpheno)

###reorder the salicylic acid file by ril ###
SWpheno.sortID<-SWpheno[order(SWpheno$ID),]
#make this a CSV
write.csv(SWpheno.sortID,file = "SWqtl/BayxSha_SWpheno.csv")

########################################################################
swphenoedit<-fread("BayxSha_SWpheno.csv")
swphenoedit<-swphenoedit[,-1]
write.csv(swphenoedit,file = "SWqtl/BayxSha_SWpheno.csv")
#########################################################################
#########QTL ANALYSIS########################
setwd("/Users/marielelensink/Documents-local/CisTransPlasticity2/Data/QTLAinput/")
bayxsha_SW<-read_cross2("BayxSha_SW.yaml",quiet = FALSE)

#calculate conditional genotype probabilities given the marker data at each putative QTL position
##insert pseudomarkers into the genetic map
map_SW<-insert_pseudomarkers(bayxsha_SW$gmap,step=0)
#calculate genotype probabilities
pr_SW<- calc_genoprob(bayxsha_SW,map_SW,error_prob = 0.002)
#analyzing marker density
#determine grid of pseudomarkers
grid_SW <- calc_grid(bayxsha_SW$gmap, step=0)

####run the genome scan####

#scan1() takes as input the genotype probabilities, a matrix of phenotypes, 
#optional additive and interactive covariates, and the special X chromosome covariates.
out_SW <- scan1(pr_SW, bayxsha_SW$pheno)
#output is LOD scores.
#LOD stands for "logarithm of the odds." 

#find peaks
peaks_SW<-find_peaks(out_SW, map_SW, threshold=2, drop=1.5)
png("../../Figures/SW_freq_tx_per_qtl.png")
hist(table(peaks_SW$lodcolumn),
     main = "Frequency of transcripts per QTL in\nSilwet Treatment Group",
     xlab = "#of transcripts per peak",
     ylab = "frequency of QTLs controlling varying #'s of transcripts")
dev.off()
fwrite(peaks_SW,"../QTLAoutput/SW_peaks_june5.txt")
