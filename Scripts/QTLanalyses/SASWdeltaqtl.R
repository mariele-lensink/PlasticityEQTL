#load packages
library(qtl2)
library(readxl)
library(Rcpp)
library(data.table)


#delta = (SA-SW)/((SA+SW)/2)
###############CREATING DELTA PHENO TABLE FROM SA AND SW PHENO TABLES####################################
#read in CSVs
setwd("/Users/marielelensink/Documents-local/CisTransPlasticity2/Data/QTLAinput/")
SAdf<-fread("BayxSha_SWpheno.csv")
SWdf<-fread("BayxSha_SApheno.csv")

rils<-SAdf[,2]
sam<-as.matrix(SAdf)
swm<-as.matrix(SWdf)
sam<-sam[,-c(1:2)]
swm<-swm[,-c(1:2)]
#calculate difference divided by mean
delta <- as.data.table((sam-swm)/((sam+swm)/2))
delta<-cbind(rils,delta)
#delta- no normalization
diff<- as.data.table(sam-swm)
diff<-cbind(rils,diff)
#make this a CSV
write.csv(delta,file = "BayxSha_deltapheno.csv")
write.csv(diff,file = "BayxSha_diffpheno.csv")
#########################################################################
#########QTL ANALYSIS########################
bayxsha_delta<-read_cross2("BayxSha_delta.yaml",quiet = FALSE)
bayxsha_diff<-read_cross2("BayxSha_diff.yaml",quiet = FALSE)
#calculate conditional genotype probabilities given the marker data at each putative QTL position
##insert pseudomarkers into the genetic map
map_delta<-insert_pseudomarkers(bayxsha_delta$gmap,step=0)
map_diff<-insert_pseudomarkers(bayxsha_diff$gmap,step=0)
#calculate genotype probabilities
pr_delta<- calc_genoprob(bayxsha_delta,map_delta,error_prob = 0.002)
pr_diff<- calc_genoprob(bayxsha_diff,map_diff,error_prob = 0.002)
#analyzing marker density
#determine grid of pseudomarkers
grid_delta <- calc_grid(bayxsha_delta$gmap, step=0)
grid_diff <- calc_grid(bayxsha_diff$gmap, step=0)

####run the genome scan####

#scan1() takes as input the genotype probabilities, a matrix of phenotypes, 
#optional additive and interactive covariates, and the special X chromosome covariates.
out_delta <- scan1(pr_delta, bayxsha_delta$pheno)
out_diff <- scan1(pr_diff, bayxsha_diff$pheno)
#output is LOD scores.
#LOD stands for "logarithm of the odds." 
# positions X phenotypes

#find peaks
peaks_delta<-find_peaks(out_delta, map_delta, threshold=2, drop=1.5)
peaks_diff<-find_peaks(out_diff, map_diff, threshold=2, drop=1.5)

fwrite(peaks_delta,"../QTLAoutput/deltapeaks_june5.txt")
fwrite(peaks_delta,"../QTLAoutput/diffpeaks_june5.txt")

png("../../Figures/delta_freq_tx_per_qtl.png")
hist(table(peaks_delta$lodcolumn),
     main = "Frequency of transcripts per QTL in\nDelta Treatment Group",
     xlab = "#of transcripts per peak",
     ylab = "frequency of QTLs controlling varying #'s of transcripts")
dev.off()

png("../../Figures/diff_freq_tx_per_qtl.png")
hist(table(peaks_diff$lodcolumn),
     main = "Frequency of transcripts per QTL in\n Difference Treatment Group",
     xlab = "#of transcripts per peak",
     ylab = "frequency of QTLs controlling varying #'s of transcripts")
dev.off()
