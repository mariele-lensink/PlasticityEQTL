parents<-fread("Data/parents_pheno_alltrts_AGI.csv")
parents$V1<-NULL
parents$rep<-NULL
rils<-fread("Data/rils_pheno_alltrts_AGI.csv")
rils$V1<-NULL
colnames(rils)[2]<-"genotype"

delta_all<-rbind(rils[treatment=='delta'],parents[treatment=='delta'])
delta_all[,1:7]

parents[,1:5]
p_avg<-parents[,lapply(.SD,mean),by=.(genotype,treatment)]
dsha<-(p_avg[2,-c(1:2)]-p_avg[4,-c(1:2)])/((p_avg[2,-c(1:2)]+p_avg[4,-c(1:2)])/2)
delta_sha<-data.table(genotype= "Sha",treatment = "delta", dsha)

dbay<-(p_avg[1,-c(1:2)]-p_avg[3,-c(1:2)])/((p_avg[1,-c(1:2)]+p_avg[3,-c(1:2)])/2)
delta_bay<-data.table(genotype="bay",treatment="delta",dbay)

delta_all<-rbind(delta_sha,delta_bay,delta_all)
fwrite(delta_all,"rils.parents_delta.txt")
