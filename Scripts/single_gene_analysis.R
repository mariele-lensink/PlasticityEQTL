library(ggplot2)
library(data.table)
library(gridExtra)
library(tidyr)
library(pbapply)

dt<-fread("finalpeaks_10cMwindow_trtgrps_rsq.txt")
alldt<-fread("finalpeaks_10cMwindow_trtgrps_rsq.txt")
highs<-dt[,.SD[rsq==max(rsq)],by = .(trtgroups,type)]

lows<-dt[,min(rsq),by = .(trtgroups,type)]
lowsdt<-dt[which(rsq %in% lows$V1),]
#interesting genes to look at
pr1<-dt[name == "At2g14610",]
pr5<-dt[name =="At1g75040"]
pad4<-dt[name == "At3g52430"]

#get transcript abundance data
pheno<-fread("rils_pheno_alltrts_AGI.csv")
pheno$V1<-NULL
pheno[1:5,1:5]
#read in map of genotypes for the rils
geno<-fread("BayxSha_geno.csv")
geno$V1<-NULL
test<-sub("(At.+)\\.[0-9]+", "\\1",colnames(geno))
test<-test[-1]
cc<-fread("gmap_conversion.txt")
newnames<-c(unlist(lapply(test, function(x){
  if (grepl("At",x)==FALSE) {
    loc<-which(cc$originalmarker == x)
    cc$marker[loc]
  } else {
    x
  }
})))
colnames(geno)<-c("ID",newnames)

#if transcript does not have genotype information, just get the genotype information of the nearest marker
gmap<-fread("gmap_edited_TAIR_allinfo.txt")
highsdt


#highsdt$closestmarkertotx
alldt$closestmarkertotx<-as.data.table(apply(alldt,1,FUN = function(x) {
  txchr<-as.numeric(substring(x[2],3,3))
  txpos<-as.numeric(substring(x[2],5,9))
  #make subset of data based on chromosome
  subsetgmap<-gmap[chr == txchr,]
  newloc<-subsetgmap$loc[which.min(abs(subsetgmap$loc-txpos))]
  subsetgmap[loc == newloc, marker]
}))



# Create a function to generate plots for each gene and group
create_gene_plot <- function(row) {
  #row<-highsdt[7,]
  gene <- as.character(row[2])
  trtment<-as.character(row[4])
  nearestmarker <-as.character(row[8])
  group <- as.character(row[7])
  regtype <- as.character(row[5])
  marker<- as.character(row[3])
  
  gene_tx <- pheno[, .(treatment = get(trtment), ID, transcript = get(gene))]
  gene_geno <- geno[, .(ID, genotype = get(nearestmarker))]
  genedt <- merge(gene_tx, gene_geno)
  genedt[genotype == -1] <- NA
  genedtmean <- na.omit(genedt[, mean(transcript), by = .(treatment, genotype)])
  colnames(genedtmean)[3] <- "transcription"
  
  saswmax<-max(genedtmean[treatment!="delta",transcription])
  saswmin<-min(genedtmean[treatment!="delta",transcription])
  deltamax<-max(genedtmean[treatment=="delta",transcription])
  deltamin<-min(genedtmean[treatment=="delta",transcription])
  scalenum<-(saswmax-saswmin)/(deltamax-deltamin)
  shiftnum<-saswmin-deltamin
  #function to scale the secondary axis
  scale_function<-function(x,scalevalue,shiftvalue){
    return ((x-shiftvalue)/scalevalue)
  }
  inv_scale_function<-function(x,scalevalue,shiftvalue){
    return ((x)*scalevalue+shiftvalue)
  }
  #change data from long to wide
  genedtmean_wide<-spread(genedtmean,treatment,transcription)
  
  p <- ggplot(genedtmean_wide, aes(x = genotype, y = SA))+
    geom_line(aes(group =1,color = "SA"))+
    geom_line(aes(y = SW, color = "SW",group = 1))+
    geom_line(aes(y = inv_scale_function(delta,scalenum,shiftnum),group = 1,color = "Delta"))+
    scale_y_continuous(,sec.axis = sec_axis(~scale_function(.,scalenum,shiftnum),name="Delta Scale"))+
    ggtitle(paste("Gene:", gene, "Group:", group,"\n", "Type:",regtype,"Marker:",marker))+
    theme(plot.title = element_text(size = 11, face = "bold"))+
    ylab("Transcript Abundance")
  
  
  return(p)
}

# Apply the create_gene_plot function to each row of highsdt and get a list of plots
plot_list <- apply(highsdt, 1, create_gene_plot)

# Save the grid of plots to a PDF file
pdf("single_gene_plots_fixed_final.pdf", width = 8, height = 14) # Adjust dimensions as needed
do.call("grid.arrange", c(plot_list, ncol = 2, top = "Gene Expression Plots"))
dev.off() # Close the PDF device




#for every individual, at a specific transcript, 
#gather the genotype information of the transcript(marker nearest to transcript if necesseary),
#the transcript abundance
#this makes the data set large so dont do the full thing
merge_pheno_geno<-function(dt){
  pbapply(dt,1,FUN = function(row){
    gene <- as.character(row[2])
    trtment<-as.character(row[4])
    nearestmarker <-as.character(row[8])
    group <- as.character(row[7])
    regtype <- as.character(row[5])
    marker<- as.character(row[3])
    gene_tx <- pheno[, .(treatment = trtment, ID, transcript = get(gene))]
    gene_geno <- geno[, .(ID, genotype = get(nearestmarker))]
    genedt <- merge(gene_tx, gene_geno)
    genedt[genotype == -1] <- NA
    genedt$txname<-rep(gene)
    genedt$closestmarkertotx<-rep(nearestmarker)
    genedt$trt<-rep(trtment)
    genedt$trtgroup<-rep(group)
    genedt$type<-rep(regtype)
    genedt$qtlmarker<-rep(marker)
    return(genedt)
  })
}



#get the 98th percentile of rsq for each treatment group to limit number of transcripts to look at
sathreshold<-quantile(alldt[trt=="SA",rsq],.99)
swthreshold<-quantile(alldt[trt=="SW",rsq],.99)
deltathreshold<-quantile(alldt[trt=="delta",rsq],.99)
sa99<-alldt[trt=="SA"&rsq>sathreshold,]
sw99<-alldt[trt=="SW"&rsq>swthreshold,]
delta99<-alldt[trt=="delta"&rsq>deltathreshold,]
dt99<-rbind(sa99,sw99,delta99)
dt99<-dt99[!duplicated(dt99),]

geneinfodttest<-rbindlist(merge_pheno_geno(dt99))
geneinfodt_nodup
highsgeneinfo<-rbindlist(merge_pheno_geno(highsdt))
highsgeneinfo<-highsgeneinfo[genotype == -1]<-NA
highsgeneinfo<-na.omit(highsgeneinfo)

#highplotlist<-lapply(highsgeneinfo, FUN = function(genedt){
genedt[genotype == -1] <- NA
genedtmean <- na.omit(genedt[, mean(transcript), by = .(treatment, genotype)])
colnames(genedtmean)[3] <- "transcription"

saswmax<-max(genedtmean[treatment!="delta",transcription])
saswmin<-min(genedtmean[treatment!="delta",transcription])
deltamax<-max(genedtmean[treatment=="delta",transcription])
deltamin<-min(genedtmean[treatment=="delta",transcription])
scalenum<-(saswmax-saswmin)/(deltamax-deltamin)
shiftnum<-saswmin-deltamin
#function to scale the secondary axis
scale_function<-function(x,scalevalue,shiftvalue){
  return ((x-shiftvalue)/scalevalue)
}
inv_scale_function<-function(x,scalevalue,shiftvalue){
  return ((x)*scalevalue+shiftvalue)
}
#change data from long to wide
genedtmean_wide<-spread(genedtmean,treatment,transcription)

p <- ggplot(genedtmean_wide, aes(x = genotype, y = SA))+
  geom_line(aes(group =1,color = "SA"))+
  geom_line(aes(y = SW, color = "SW",group = 1))+
  geom_line(aes(y = inv_scale_function(delta,scalenum,shiftnum),group = 1,color = "Delta"))+
  scale_y_continuous(,sec.axis = sec_axis(~scale_function(.,scalenum,shiftnum),name="Delta Scale"))+
  ggtitle(paste("Gene:", gene, "Group:", group,"\n", "Type:",regtype,"Marker:",marker))+
  theme(plot.title = element_text(size = 11, face = "bold"))+
  ylab("Transcript Abundance")


return(p)
})







#make environment the x axis
ggplot(geneinfodttest[treatment!="delta"&type=="cis"], aes(x=treatment,y=transcript,color=genotype,group = ID))+
  geom_line()+
  ggtitle("cis")
ggplot(geneinfodttest[treatment!="delta"&type=="trans"], aes(x=treatment,y=transcript,color=genotype,group = ID))+
  geom_line()+
  ggtitle("trans")


ggplot(highsgeneinfo[treatment!="delta"&type=="cis"], aes(x=treatment,y=transcript,color=genotype,group = ID))+
  geom_line()+
  ggtitle("cis")

