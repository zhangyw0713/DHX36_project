setwd("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/00_analysis/05_dstruct_analysis/localized")
read.table("dif_localized_newpara_5.tsv",sep = "\t",header = F)->df
colnames(df)<-c("peak","g4","g4sub","region")
length(unique(df[,1]))->difpeak
allpeak=4263
matrix(c((4263-difpeak)/4263*100,"non-DRRs",difpeak/4263*100,"DRRs"),nrow=2,ncol=2)->df1
data.frame(t(df1))->df2
colnames(df2)<-c("value","group")
bp<- ggplot(df2, aes(x="", y=as.numeric(as.character(value)), fill=group))+
  geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+ scale_fill_manual(values=c("#d8b365", "#5ab4ac"))+theme_minimal(base_size = 13)+ylab("Proportion of DRRs")
bp
ggsave("DRR_proportion_newpara_5.pdf",bp,width=5,height=5)


bp<- ggplot(region.df, aes(x="", y=as.numeric(as.character(freq)), fill=region))+
  geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+ scale_fill_brewer(palette="Dark2")+theme_minimal(base_size = 13)+ylab("Proportion of DHX36 binding regions")
bp
ggsave("DRR_region_proportion_newpara_5.pdf",bp,width=5,height=5)
#-----------------------bin shape--------------------------------------
library(tidyverse)
library(matrixStats)
library(Rmisc)

library(grid)
library(gridExtra)
setwd("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/00_analysis/05_dstruct_analysis/localized")
read.table("drr.shape_newpara_5.tsv",sep = "\t",header=F)->df

df[!duplicated(df[,1]),]->df1
df1[,c(1,3)]->ko.df
df1[,c(1,4)]->wt.df

separate(data = ko.df, col = V3, into = paste0('N', 1:101), sep = ",")->ko.df1
rownames(ko.df1)<-ko.df1[,1]
ko.df1[,-1]->ko.df2
matrix(nrow=nrow(ko.df2),ncol=5)->mm.ko
for (i in 1:nrow(ko.df2)) {

  length(str_split(ko.df[i,2],",")[[1]])-1->g4len
  seq(1,g4len,length.out=6)->bx
  ko.df2[i,1:g4len]->y
  x=c(1:g4len)
  yS <- binMeans(as.numeric(as.character(y)), x=x, bx=bx)
  yS ->mm.ko[i,]
}
apply(mm.ko, 2, mean,na.rm=T)->ko.drr.mean
apply(mm.ko, 2, sd,na.rm=T)


separate(data = wt.df, col = V4, into = paste0('N', 1:101), sep = ",")->wt.df1
rownames(wt.df1)<-wt.df1[,1]
wt.df1[,-1]->wt.df2
matrix(nrow=nrow(wt.df2),ncol=5)->mm.wt
for (i in 1:nrow(wt.df2)) {
  
  length(str_split(wt.df[i,2],",")[[1]])-1->g4len
  seq(1,g4len,length.out=6)->bx
  wt.df2[i,1:g4len]->y
  x=c(1:g4len)
  yS <- binMeans(as.numeric(as.character(y)), x=x, bx=bx)
  yS ->mm.wt[i,]
}
apply(mm.wt, 2, mean,na.rm=T)->wt.drr.mean


data.frame(wt.drr.mean) %>% mutate(group="wt",pos=1:5)->wt.out
data.frame(ko.drr.mean) %>% mutate(group="ko",pos=1:5)->ko.out
colnames(wt.out)[1]<-"value"
colnames(ko.out)[1]<-"value"
rbind(wt.out,ko.out)->ofplot



matrix(nrow=5,ncol=4)->deltadf
for (i in 1:5) {
  deltadf[i,1]<-i
  t.test(mm.wt[,i],mm.ko[,i],paired = T)->aa
  aa$estimate[[1]]->deltadf[i,2]
  aa$conf.int[1]->deltadf[i,3]
  aa$conf.int[2]->deltadf[i,4]
}
data.frame(deltadf)->deltadf

my_theme <- theme_bw()+
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_blank())
y_only_theme <- my_theme+
  theme(axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank())

ggplot(ofplot,aes(x = pos, y = value, colour = factor(group)))+
  geom_line(size = 1)+
  scale_y_continuous(limits = c(0.15,0.55))+y_only_theme->toppanel

ggplot(data=ofplot,aes(x=pos,y=as.numeric(value),fill=group)) +geom_bar(stat = "identity",position = position_dodge(),alpha=.8,width = 0.8)+theme_bw(base_size = 13)+xlab("Bins")+ylab("Average reactivity")+theme(legend.position="top",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())->toppanel
ggplot(data = deltadf, aes(x = X1))+
  geom_bar(aes(y = -1*X2), fill = "grey", position = position_dodge(), stat='identity',width = 0.8)+
  geom_errorbar(aes(ymin = -1*X3, ymax = -1*X4),width=.2, position=position_dodge(.9),alpha = 0.3, colour = "black")+theme_bw(base_size=13)+ylim(-0.3,0)+theme(legend.position="none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.title.x = element_blank())+ylab("Delta reactivity")->btpanel


pdf(file = 'drr.bin5.pdf', width = 5, height = 6)
grid.arrange(toppanel,btpanel,nrow=2,heights=c(3,1.5))
dev.off()







ggplot(data = deltadf, aes(x = X1))+
  geom_bar(aes(y = -1*X2), fill = "grey", position = position_dodge(), stat='identity',width = 0.8)+
  geom_ribbon(aes(ymin = -1*X3, ymax = -1*X4), alpha = 0.3, colour = NA)+theme_bw(base_size=13)+ylim(-0.3,0)+theme(legend.position="none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.title.x = element_blank())+ylab("Delta reactivity")->btpanel

geom_errorbar(aes(ymin = value - ci, ymax = value + ci),
              width=.2, position=position_dodge(.9))
grid.arrange(toppanel,btpanel,nrow=2)
 








#---------------------normalized regional enrichment-------------------------------


x1 <- matrix(c(407, 404, 1664, 568), ncol = 2)
chisq.test(x1)$p.value
x2 <- matrix(c(407, 404, 985, 376), ncol = 2)
chisq.test(x2)$p.value
x3 <- matrix(c(1664, 568, 985, 376), ncol = 2)
chisq.test(x3)


setwd("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/00_analysis/05_dstruct_analysis/localized")
matrix(ncol=1,nrow=3)->m1
985/1361->m1[1,1]
1664/2232->m1[2,1]
407/811->m1[3,1]

colnames(m1)<-"proportion"
rownames(m1)<-c("CDS","3UTR","5UTR")

data.frame(m1)->df
df$region<-rownames(df)


factor(df$region,levels = c("5UTR","CDS","3UTR"))->df$region

ggplot(df, aes(x = region, y=proportion,fill = region)) +
  geom_bar(stat="identity",width = 0.8)+ scale_fill_brewer(palette="Dark2") +
  theme_bw()+theme(legend.position="top")+annotate("segment", x=c(1,2,1),xend=c(2,3,3), y= c(0.6,0.65,0.7), yend=c(0.6,0.65,0.7))+xlab("Regions")+ylab("Proportion of DRR-containing sites")+theme(axis.title.x = element_text(size=13,face = "bold"),axis.title.y = element_text(size=13,face = "bold"),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13))+  annotate("text",x=1.5,y=0.62,label="3.06e-25")+  annotate("text",x=2.5,y=0.67,label="0.161")+  annotate("text",x=2,y=0.72,label="5.73e-37")+ylim(0,0.9)->p1
p1
ggsave(p1,filename = "DRRproportionpeak_newpara_5.pdf",width = 4,height = 5)

