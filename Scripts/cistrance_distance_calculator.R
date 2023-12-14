library(data.table)
library(ggplot2)
library(readxl)
library(stringr)

#probe id gene name conversion list from Dan
probedf_full<-data.table(read_xls("Data/originaldata/Row Labels in Order.xls"))
#get peak tables from qtl scripts
data.table(peaks_SA)
data.table(peaks_SW)
data.table(peaks_delta)
#reordering peaks by chromosome then position on chromosome (not sure if necessary but it makes me feel better)
SA<-peaks_SA[,2:7]
SA$trt<-"SA"
SW<-peaks_SW[,2:7]
SW$trt<-"SW"
delta<-peaks_delta[,2:7]
delta$trt<-"delta"
#combining the tables
peakdf<-as.data.table(rbind(SA,SW,delta))
#regex to remove x infront of probenames
peakdf$lodcolumn <- sub("^X", "",peakdf$lodcolumn)
#for now only keeping the probeids that exist in my conversion table, lost 157 peaks
peakdf_subset<-peakdf[peakdf$lodcolumn %in% probedf_full$`Affy probe set ID`]

#swap names of probe ids in peak table for transcript names (so they match the gmap)
newnames<-data.table(sapply(peakdf_subset$lodcolumn, function(x){
  loc<-which(probedf_full$`Affy probe set ID` == x)
  probedf_full$AGI[loc]
}))
colnames(newnames)<-"name"
peakdf<-cbind(newnames,peakdf_subset[,2:7])

######get marker information#########
gmap<-fread("Data/QTLAinput/BayxSha_gmap.csv")
#removing decimal and following digits
gmap$originalmarker<-gmap$marker
gmap$marker <- sub("\\.[0-9]+", "",gmap$marker)

#Get gaps between marker values
hist(gmap[chr==1,pos], breaks = 150)
hist(gmap[chr==2,pos], breaks = 150)
hist(gmap[chr==3,pos], breaks = 150)
hist(gmap[chr==4,pos], breaks = 150)
hist(gmap[chr==5,pos], breaks = 150)
#gonna revisit this but first just moving forward like normal

#####580 of the 83333 peaks do not have a tair name, simply removing for now###
wrongnamespeaks<-peakdf[!str_detect(peakdf$name,"At")]
peakdf<-peakdf[str_detect(peakdf$name,"At")]

###formatting peak table to make comparing to gmap easier
#adding the number of the chromosome where the transcript is located to compare to chromosome the peak is on
peakdf$chr_transcript<-as.numeric(substr(peakdf$name, 3,3))
#using the last part of the tair name as a pseudolocation
peakdf$loc<-as.numeric(substr(peakdf$name, 5,10))

#formatting gmap
#subset markers that dont have tair names, there are 38
neednewnames1<-gmap[!str_detect(gmap$marker,"At")]
gmap<-gmap[str_detect(gmap$marker,"At")]
###to clean up gmap#####
gmap[marker =="T1G11"]
neednewnames[38]
neednewnames[38,1]<- "At5g64810"
neednewnames
neednewnames1

gmap<-rbind(gmap,neednewnames)
############################### can carry on w/out fixed markers for now

#adding pseudolocation based on tair name for distance calc
gmap$loc<-as.numeric(substr(gmap$marker, 5,10))
#reorder gmap and pekas by tair name
gmap<-gmap[order(gmap$marker)]

gmap<-fread("gmap_edited_TAIR_allinfo.txt")
peakdf<-peakdf[order(peakdf$name)]
#separate out peaks that are on separate chromosomes from the transcript they control, they are automatically trans
transpeaks<-peakdf[chr!=chr_transcript]
transpeaks$newloc<-NA
transpeaks$distance<-100
calcpeaks<-data.table(peakdf[chr==chr_transcript])


#function that gets the closest marker to the peak marker and creates a new column for it, factor in chromosome!
calcpeaks$newloc<-as.data.table(apply(calcpeaks,1,FUN = function(x) {
  #make subset of data based on chromosome
  subsetgmap<-gmap[chr == x[2],]
  subsetgmap$loc[which.min(abs(subsetgmap$loc-as.numeric(x[9])))]
}))

#add a column with the distance of each peak to its marker
calcpeaks$distance<-as.data.table(apply(calcpeaks,1, FUN = function(x){
  subsetgmap<-gmap[chr==as.numeric(x[2])]
  abs(as.numeric(x[3])-as.numeric(subsetgmap[loc == as.numeric(x[10]),pos]))
}))
#add a column with the name of the marker linked to the transcript


moretrans<-calcpeaks[distance>=5]
cispeaks<-data.table(calcpeaks[distance<5])
cispeaks$type<-"cis"

transpeaks$newloc<-NA
transpeaks$distance<-NA
transpeaks<-rbind(transpeaks,moretrans)
transpeaks$type<-"trans"                   

peaksfinal<-as.data.table(rbind(transpeaks,cispeaks))

#aadding the marker name to each row 
peaksfinal$markername<-as.data.table(apply(peaksfinal,1, FUN = function(x){
  subsetgmap<-gmap[chr==as.numeric(x[2])]
  subsetgmap[pos ==as.numeric(x[3]),marker]
}))

fwrite(peaksfinal,"Data/finalpeaks_10cMwindow_june5.txt")



