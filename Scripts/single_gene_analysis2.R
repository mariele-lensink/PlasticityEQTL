library(ggplot2)
library(data.table)
library(gridExtra)
library(tidyr)
library(pbapply)
library(ggh4x)
library(forcats)
library(ggforce)

dt<-fread("finalpeaks_10cMwindow_trtgrps_rsq.txt")
highs<-dt[,.SD[rsq==max(rsq)],by = .(trtgroups,type)]
colnames(highs)[4]<-'txname'

#get transcript abundance data
pheno<-fread("rils_pheno_alltrts_AGI.csv")
pheno$V1<-NULL
pheno[1:5,1:5]
#read in map of genotypes for the rils, make names compatible 
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
#read in gmap
gmap<-fread("gmap_edited_TAIR_allinfo.txt")
#if transcript does not have genotype information, 
#just get the genotype information of the nearest marker
highs$closestmarkertotx<-as.data.table(apply(highs,1,FUN = function(x) {
  txchr<-as.numeric(substring(x[4],3,3))
  txpos<-as.numeric(substring(x[4],5,9))
  #make subset of data based on chromosome
  subsetgmap<-gmap[chr == txchr,]
  newloc<-subsetgmap$loc[which.min(abs(subsetgmap$loc-txpos))]
  subsetgmap[loc == newloc, marker]
}))
#function to merge the phenotype and genotype information based on the qtl data
merge_pheno_geno<-function(qtldt){
  pbapply(qtldt,1,FUN = function(row){
    #row<-highs[1,]
    tx <- as.character(row[4])
    trtment<-row[6]
    nearestmarker <-as.character(row[8])
    group <- as.character(row[1])
    regtype <- as.character(row[2])
    marker<- as.character(row[5])
    gene_tx <- pheno[, .(ID,treatment, transcript = get(tx))]
    gene_geno <- geno[, .(ID, genotype = get(nearestmarker))]
    genedt <- merge(gene_tx, gene_geno)
    genedt$txname<-rep(tx)
    genedt$closestmarkertotx<-rep(nearestmarker)
    genedt$trtgroup<-rep(group)
    genedt$type<-rep(regtype)
    genedt$qtlmarker<-rep(marker)
    return(genedt)
  })
}
highsgeneinfo<-rbindlist(merge_pheno_geno(highs))

geneplots<-lapply(highsgeneinfolist,FUN = function(genetable){
  genetable[genotype==-1]<-NA
  genetable<-na.omit(genetable)
  genedtmean <- na.omit(genetable[, mean(transcript), by = .(treatment, genotype)])
  colnames(genedtmean)[3] <- "transcript_abundance"
  gpl<-ggplot(genedtmean,aes(x=genotype,y=transcript_abundance,color=treatment))+
  geom_line(aes(group=treatment))
  
  return(gpl)
  })
pdf("single_gene_plots_fixed_final.pdf", width = 8, height = 14) # Adjust dimensions as needed
do.call("grid.arrange", c(plot_list, ncol = 2, top = "Gene Expression Plots"))
dev.off() # Close the PDF device






#highsdtmean<-na.omit(highsgeneinfo[,mean(transcript),by = .(treatment,genotype,txname,type)])
colnames(highsdtmean)[5] <- "transcript_abundance"
#highsgeneinfo$id.trt<-paste(highsgeneinfo[,1],highsgeneinfo[,2])

highsgeneinfo[genotype==-1]<-NA
highsgeneinfo<-na.omit(highsgeneinfo)
##########################################################################################################
#environment single gene plots
##########################################################################################################
unique_combinations<-unique(highsgeneinfo[,.(trtgroup,txname,type,qtlmarker)])
highsnodelta<-highsgeneinfo[treatment!='delta']
highsnodelta<-highsnodelta[order(match(trtgroup,c("SA","SW","delta","SA SW","SA delta","SW delta","SA SW delta")))]
ggplot(highsnodelta, aes(x=treatment, y=transcript, color=genotype)) +
  geom_line(aes(group=interaction(ID,txname))) +
  facet_grid2(fct_relevel(trtgroup,"SA","SW","delta","SA SW","SA delta","SW delta","SA SW delta")~type,scales = 'free_y',independent = 'y')+
  theme_bw()+
  scale_color_manual(values = c("#56B4E9","#CC79A7"))+
  theme(legend.position = 'none',
        text = element_text(size = 9),
        axis.text = element_text(size = 9))+
  ylab("transcript abundace")
  # This line adds facets, in a grid organized by 'trtgroup' and 'type'
  # Replace with your actual main title
##########################################################################################################
#genotype single gene plots
##########################################################################################################
highs_avg<-highsgeneinfo[genotype!="-1",mean(transcript),by=c('treatment','genotype','type','trtgroup')]
highs_avg1<-highs_avg[,.(ID,treatment,genotype,txname,qtlmarker,trtgroup,type,avg_abundance)]
custom_colors <- c("SA" = "#377EB8", "delta" = "#E6550D", "SW" = "#4DAF4A")
ggplot(highs_avg1, aes(x=genotype, y=avg_abundance, color=treatment)) +
  geom_line(aes(group=treatment)) +
  facet_grid2(fct_relevel(trtgroup,"SA","SW","delta","SA SW","SA delta","SW delta","SA SW delta")~type,scales = 'free_y',independent = 'y')+
  theme_bw()+
  scale_x_discrete(labels=c('C'="Bay","L"="Sha"))+
  scale_color_manual(values = c(custom_colors))+
  theme(legend.position = 'none',
        text = element_text(size = 9),
        axis.text = element_text(size = 9))+
  ylab("transcript abundace")
##########################################################################################################
# Plot for non-delta treatments
p1 <- ggplot(subset(highs_avg1, treatment != "delta"), aes(x=genotype, y=avg_abundance, color=treatment)) +
  geom_line(aes(group=treatment)) +
  facet_grid2(fct_relevel(trtgroup,"SA","SW","delta","SA SW","SA delta","SW delta","SA SW delta") ~ type,scale= 'free') +
  theme_bw() +
  scale_x_discrete(labels=c('C'="Bay","L"="Sha")) +
  scale_color_manual(values = custom_colors) +
  theme(legend.position = 'none',
        text = element_text(size = 9),
        axis.text = element_text(size = 9)) +
  ylab("transcript abundance")

# Plot for delta treatment
p2 <- ggplot(subset(highs_avg1, treatment == "delta"), aes(x=genotype, y=avg_abundance, color=treatment)) +
  geom_line(aes(group=treatment)) +
  facet_grid2(fct_relevel(trtgroup,"SA","SW","delta","SA SW","SA delta","SW delta","SA SW delta") ~ type,scale = 'free') +
  theme_bw() +
  scale_x_discrete(labels=c('C'="Bay","L"="Sha")) +
  scale_color_manual(values = custom_colors) +
  theme(legend.position = 'none',
        text = element_text(size = 9),
        axis.text = element_text(size = 9)) +
  ylab("transcript abundance")

# Combine the plots
grid.arrange(p1,p2, ncol=2)















avg_data <- highsgeneinfo[treatment != 'delta',.(ID,transcript,treatment,txname,genotype,closestmarkertotx,trtgroup,type,qtlmarker),] 
avg_data<-avg_data[,mean(transcript),by=.(genotype,txname,trtgroup,type,treatment)]
colnames(avg_data)[6]<-"avg_transcript"  
#need to manually go through and clean out
unique_combinations<-unique(avg_data[,.(trtgroup,txname,type)])
# Open PDF file for saving plots
pdf("Average_Plots.pdf", width=4, height=4)
# Loop over each unique factor in the column (replace 'your_column' with the name of the column)
for(i in 1:nrow(unique_combinations)) {
  # Subset data for the specific factor
  sub_data <- avg_data[trtgroup == unique_combinations[i,"trtgroup"]&type== unique_combinations[i,"type"]&txname==unique_combinations[i,"txname"]]
  trt<-sub_data$trtgroup[1]
  cisortrans<-sub_data$type[1]
  transcriptname<-sub_data$txname[1]
  # Generate the ggplot for each subset
  p <- ggplot(sub_data, aes(x=treatment, y=avg_transcript, color=genotype)) +
    geom_point() +
    geom_line(aes(group=genotype))+
    ggtitle(paste("Average Plot for", transcriptname,"\nTreatment Groups: ",trt,"\nType:",cisortrans))
  print(p)
}
dev.off()









































saswmax<-max(highsdtmean[treatment!="delta",transcript_abundance])
saswmin<-min(highsdtmean[treatment!="delta",transcript_abundance])
deltamax<-max(highsdtmean[treatment=="delta",transcript_abundance])
deltamin<-min(highsdtmean[treatment=="delta",transcript_abundance])
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
highsdtmean_wide<-spread(highsdtmean,treatment,transcript_abundance)
highsdtmean


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
