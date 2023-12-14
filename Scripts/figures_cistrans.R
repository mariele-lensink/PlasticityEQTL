library(data.table)
library(ggplot2)
library(tidyr)
library(ggupset)
library(stringr)
library(forcats)
library(dplyr)
library(reshape2)
library(ggpmisc)
library(devtools)
library(plotly)
library(rcartocolor)

dt<-fread("finalpeaks_10cMwindow_june5.txt")

#calculating the proportion of transcripts regulated by cis vs trans qtl for each treatment group
x<-data.table(table(dt[,type,trt]))
colnames(x)[3]<-'txcount'
ggplot(x, aes(fill = type, x = trt, y= txcount))+
  geom_col(position = 'fill')+
  geom_text(aes(label= txcount),position = position_fill(vjust=0.5),size = 3)+
  scale_y_continuous()+
  ylab("Proportion of eQTL that are Cis vs Trans")+
  xlab("Treatment Group")+
  labs(fill="Regulatory Type")+
  theme_bw()+
  theme(legend.position="none",text = element_text(size = 9),axis.text = element_text(size=9))+
  scale_fill_grey(start = 0.5,end=0.8)
 
#get the sum stats and figures for how many qtl are linked to a given transcript
dt1<-dt[,.(qtlpertx=.N),by= .(name,trt)]
custom_colors <- c("SA" = "#377EB8", "delta" = "#E6550D", "SW" = "#4DAF4A")
ggplot(dt1,aes(x=qtlpertx,fill=factor(trt)))+
  geom_histogram(position = position_dodge2(padding = 0.05),binwidth = 1,)+
  scale_fill_manual(values=custom_colors) +
  scale_y_continuous(n.breaks = 6)+
  xlab('Number of QTL per Transcript')+
  ylab("Number of Transcripts")+
  labs(fill="Treatment Group")+
  theme_bw()+
  theme(legend.position = 'none',text = element_text(size = 9),axis.text = element_text(size=9))
#get summary stats for the above figure
summary(dt1[trt=='SA',qtlpertx])
###trying to make the same figure but with bar plot instead###
dtbar<-dt1[,.(Count=.N),by=.(trt,qtlpertx)]
ggplot(dtbar,aes(x=qtlpertx,y=Count,fill=trt))+
  geom_col(position=position_dodge())+
  scale_fill_manual(values=custom_colors)+
  xlab('Number of QTL per Transcript')+
  ylab("Number of Transcripts")+
  labs(fill="Treatment Group")+
  theme_bw()+
  theme(legend.position = 'none',text = element_text(size = 9),axis.text = element_text(size=9))


cisdt<-dt[type=='cis']
cisdt<-cisdt[, transcript.marker := paste0(name,".",markername)]
cisdt<-cisdt[,.(transcript.marker,trt)]
cisdt2<-cisdt[,.(alltrts=list(trt)),by = transcript.marker]
cisdt2$type<-'cis'
transdt<-dt[type =='trans',]
transdt<-transdt[, transcript.marker := paste0(name,".",markername)]
transdt<-transdt[,.(transcript.marker,trt)]
transdt2<-transdt[,.(alltrts=list(trt)),by = transcript.marker]
transdt2$type<-'trans'
alldt2<-rbind(cisdt2,transdt2)
#upset plot with transcript counts
ggplot(alldt2,aes(x=alltrts,fill = type))+
  geom_bar(position = 'dodge')+
  scale_x_upset()+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-.1,size = 3,position = position_dodge(width=1.1))+
  xlab("Treatment Group Membership(s)")+
  scale_x_upset(order_by = 'degree')+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  ylab('# of QTLs')+
  theme_bw()+
  theme(legend.position = 'none', text = element_text(size = 9),axis.text = element_text(size=9))+
  scale_fill_grey(start = 0.5,end=0.8)




effectsizes_dt<-readRDS("lm_summary")

rsqs<-unlist(lapply(effectsizes_dt,function(x){
 x[[8]]
}))

dt<-data.table(dt,rsq= rsqs)
dt<-dt[, transcript.marker := paste0(name,".",markername)]

#get the group of treatments that each marker has, when done append to 'dt' the original data table with all of the information
dt2<-dt[,.(name,markername,trt,type,rsq)]
dt2<-dt2[,.(transcript.marker,name,markername,trt,type,rsq)]
groupdt<-dt2[, .(group=list(trt)),by = .(type,transcript.marker)]
groupdt$trtgroups<-apply(groupdt,1,function(x){
  i<-toString(x[3])
  i<-gsub("^c","",i)
  gsub("[^[:alnum:][:space:]]","",i)
})
groupdt<-groupdt[,.(type,transcript.marker,trtgroups)]

dt3<-dt2 %>% inner_join(groupdt,by=c('transcript.marker'='transcript.marker', 'type'='type'))


#upset plot wih effect sizes
ggplot(dt2,aes(x=fct_infreq(trtgroups),y=rsq,factor = type))+
  geom_point(aes(color = type,alpha = 0.10),position=position_jitterdodge(0.2))+
  geom_violin()+
  xlab("Treatment Group Membership(s)")+
  ylab("Effect Size (rsq)")+
  theme_bw()+
  scale_color_grey(start = 0.5,end=0.8)+
  theme(text = element_text(size = 9),
        axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 45,vjust = 0.7),
        axis.title.x = element_text(margin = margin(b=1,t=0)),
        legend.position = 'none')
  
  
 


SASW<-dt2[trtgroups == "SA SW"]
SAdelta<-dt2[trtgroups =="SA delta"]
SWdelta<-dt2[trtgroups =="SW delta"]
dtall<-dt2[trtgroups =="SA SW delta"]
SASW_wide<-dcast(SASW,transcript.marker+type~trt,value.var = "rsq")
SAdelta_wide<-dcast(SAdelta,transcript.marker+type~trt,value.var = "rsq")
SWdelta_wide<-dcast(SWdelta,transcript.marker+type~trt,value.var = "rsq")
all_wide<-dcast(dtall,transcript.marker+type~trt,value.var = "rsq")

ggplot(SAdelta_wide,aes(delta,SA))+
  geom_point(aes(alpha=0.1))+
  geom_smooth(method = 'lm',formula = y~x)+
  facet_wrap(~type)+
  stat_poly_eq(size = 3)+
  ylim(0,1)+
  xlim(0,1)+
  theme_bw()+
  theme(text = element_text(size = 9),
        axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 0.9),
        legend.position = 'none',
        panel.spacing = unit(1, "lines"))
  


ggplot(all_wide,aes(x=SA,y=SW,z=delta,color=as.factor(type)))+
  theme_void()+
  axes_3D()+
  stat_3D()

plot_ly(all_wide,
        x=all_wide$SA,
        y=all_wide$SW,
        z=all_wide$delta,
        type="scatter3d",
        mode="markers",
        color=all_wide$type)
fwrite(dt2,"Data/finalpeaks_10cMwindow_trtgrps_rsq.txt")


  

