#AOP1 AT4G03070
#AOP3 AT4G03050
#AOP2 AT4G03060

aopdtril<-pheno[treatment!='delta',.(treatment,ID,At4g03070,At4g03050,At4g03060)]
#this is the closest marker to the aop genes
#also the marker upstream and downstream are always the same
#so im pretty confident that the genes in between are the same genotype
rilgenos<-geno[,.(ID,At4g03210)]
colnames(rilgenos)<-c("ID","localgenotype")
aoprils<-merge(aopdtril,rilgenos,by='ID')
aoprils<-aoprils[localgenotype!='-1',.(ID,localgenotype,At4g03070,At4g03050,At4g03060)]
aoprilslong<-melt(aoprils,id.vars = c('ID','localgenotype'))


ggplot(aoprilslong,aes(x=variable,y=value,fill=localgenotype,color = localgenotype))+
  geom_jitter(position = position_dodge(width=0.2))



dt<-fread("parents_pheno_alltrts_AGI.csv")
aopdtparent<-dt[treatment!='delta',.(genotype,At4g03070,At4g03050,At4g03060)]
aopdtparent$V1<-NULL
aopdtparent[,1:4]
aopparent<-melt(aopdtparent,id.vars = "genotype")
ggplot(aopparent,aes(x=variable,y=value,fill=genotype,color = genotype))+
  geom_jitter(position = position_dodge(width=0.2))
