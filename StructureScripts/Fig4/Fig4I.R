#-------------------------------------------5UTR-----------------------------------------------
setwd("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/00_analysis/06_difgene_structureCompare/utr5_structure")
library(ggplot2)
read.table("avg.down_wt.tsv",sep = "\t",header = F)->wt.down
read.table("avg.down_ko.tsv",sep = "\t",header = F)->ko.down
read.table("avg.up_wt.tsv",sep = "\t",header = F)->wt.up
read.table("avg.up_ko.tsv",sep = "\t",header = F)->ko.up
data.frame(wt.down,"wt","down")->d1
data.frame(ko.down,"ko","down")->d2
data.frame(wt.up,"wt","up")->d3
data.frame(ko.up,"ko","up")->d4
colnames(d1)<-c("pos","value","group","dif")
colnames(d2)<-c("pos","value","group","dif")
colnames(d3)<-c("pos","value","group","dif")
colnames(d4)<-c("pos","value","group","dif")
mean(d1[,2],na.rm = T)->d1[,5]
mean(d2[,2],na.rm = T)->d2[,5]
mean(d3[,2],na.rm = T)->d3[,5]
mean(d4[,2],na.rm = T)->d4[,5]

cbind(d1,d2)->downbd
cbind(d3,d4)->upbd
matrix(nrow = 2,ncol=4)->sm

t.test(downbd[,2],downbd[,7],paired=T)->dt
dt$estimate[[1]]->sm[1,2]
"down"->sm[1,1]
dt$conf.int[1]->sm[1,3]
dt$conf.int[2]->sm[1,4]

t.test(upbd[,2],upbd[,7],paired=T)->ut
ut$estimate[[1]]->sm[2,2]
"up"->sm[2,1]
ut$conf.int[1]->sm[2,3]
ut$conf.int[2]->sm[2,4]
data.frame(sm)->smdf
smdf->utr5df
rbind(d1,d2,d3,d4)->df
factor(df$group,levels=c("ko","wt"))->df$group

df[which(df$dif=="down"),"ymin"]<-0
df[which(df$dif=="down"),"ymax"]<-0.2
df[which(df$dif=="up"),"ymin"]<-0
df[which(df$dif=="up"),"ymax"]<-0.2


ggplot(data=df, aes(x=pos, y=value,group=group))+  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +geom_line(aes(color=group))+ geom_point(aes(color=group))+theme_classic()+xlab("Position")+ylab("Average reactivity")+theme(axis.title.x = element_text(size=13,face = "bold"),axis.title.y = element_text(size=13,face = "bold"),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13),legend.position="top")+facet_grid(cols = vars(dif),scale="free")+ geom_hline(aes(yintercept=V5,group=dif,color=group),alpha=0.5)+scale_color_manual(values=c("#66c2a5","#fc8d62"))->p1
p1
ggsave("utr5.pdf",p1,width=12,height=3.8)





#-------------------------------------------CDS-----------------------------------------------
setwd("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/00_analysis/06_difgene_structureCompare/cds_structure")
library(ggplot2)
read.table("avg.down_wt.tsv",sep = "\t",header = F)->wt.down
read.table("avg.down_ko.tsv",sep = "\t",header = F)->ko.down
read.table("avg.up_wt.tsv",sep = "\t",header = F)->wt.up
read.table("avg.up_ko.tsv",sep = "\t",header = F)->ko.up
data.frame(wt.down,"wt","down")->d1
data.frame(ko.down,"ko","down")->d2
data.frame(wt.up,"wt","up")->d3
data.frame(ko.up,"ko","up")->d4
colnames(d1)<-c("pos","value","group","dif")
colnames(d2)<-c("pos","value","group","dif")
colnames(d3)<-c("pos","value","group","dif")
colnames(d4)<-c("pos","value","group","dif")
mean(d1[,2],na.rm = T)->d1[,5]
mean(d2[,2],na.rm = T)->d2[,5]
mean(d3[,2],na.rm = T)->d3[,5]
mean(d4[,2],na.rm = T)->d4[,5]

cbind(d1,d2)->downbd
cbind(d3,d4)->upbd
matrix(nrow = 2,ncol=4)->sm

t.test(downbd[,2],downbd[,7],paired=T)->dt
dt$estimate[[1]]->sm[1,2]
"down"->sm[1,1]
dt$conf.int[1]->sm[1,3]
dt$conf.int[2]->sm[1,4]

t.test(upbd[,2],upbd[,7],paired=T)->ut
ut$estimate[[1]]->sm[2,2]
"up"->sm[2,1]
ut$conf.int[1]->sm[2,3]
ut$conf.int[2]->sm[2,4]
data.frame(sm)->smdf
smdf->cdsdf

rbind(d1,d2,d3,d4)->df
factor(df$group,levels=c("ko","wt"))->df$group


df[which(df$dif=="down"),"ymin"]<-0.05
df[which(df$dif=="down"),"ymax"]<-0.25
df[which(df$dif=="up"),"ymin"]<-0.05
df[which(df$dif=="up"),"ymax"]<-0.3





ggplot(data=df, aes(x=pos, y=value,group=group))+  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +geom_line(aes(color=group))+ geom_point(aes(color=group))+theme_classic()+xlab("Position")+ylab("Average reactivity")+theme(axis.title.x = element_text(size=13,face = "bold"),axis.title.y = element_text(size=13,face = "bold"),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13),legend.position="top")+facet_grid(cols = vars(dif),scale="free")+ geom_hline(aes(yintercept=V5,group=dif,color=group),alpha=0.5)+scale_color_manual(values=c("#66c2a5","#fc8d62"))->p1
p1
ggsave("cds.pdf",p1,width=12,height=3.8)

#------------------------------------------------------3UTR------------------------------------------------
setwd("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/00_analysis/06_difgene_structureCompare/utr3_structure")
library(ggplot2)
read.table("avg.down_wt.tsv",sep = "\t",header = F)->wt.down
read.table("avg.down_ko.tsv",sep = "\t",header = F)->ko.down
read.table("avg.up_wt.tsv",sep = "\t",header = F)->wt.up
read.table("avg.up_ko.tsv",sep = "\t",header = F)->ko.up
data.frame(wt.down,"wt","down")->d1
data.frame(ko.down,"ko","down")->d2
data.frame(wt.up,"wt","up")->d3
data.frame(ko.up,"ko","up")->d4
colnames(d1)<-c("pos","value","group","dif")
colnames(d2)<-c("pos","value","group","dif")
colnames(d3)<-c("pos","value","group","dif")
colnames(d4)<-c("pos","value","group","dif")
mean(d1[,2],na.rm = T)->d1[,5]
mean(d2[,2],na.rm = T)->d2[,5]
mean(d3[,2],na.rm = T)->d3[,5]
mean(d4[,2],na.rm = T)->d4[,5]

cbind(d1,d2)->downbd
cbind(d3,d4)->upbd
matrix(nrow = 2,ncol=4)->sm

t.test(downbd[,2],downbd[,7],paired=T)->dt
dt$estimate[[1]]->sm[1,2]
"down"->sm[1,1]
dt$conf.int[1]->sm[1,3]
dt$conf.int[2]->sm[1,4]

t.test(upbd[,2],upbd[,7],paired=T)->ut
ut$estimate[[1]]->sm[2,2]
"up"->sm[2,1]
ut$conf.int[1]->sm[2,3]
ut$conf.int[2]->sm[2,4]
data.frame(sm)->smdf

rbind(d1,d2,d3,d4)->df
factor(df$group,levels=c("ko","wt"))->df$group

df[which(df$dif=="down"),"ymin"]<-0.1
df[which(df$dif=="down"),"ymax"]<-0.45
df[which(df$dif=="up"),"ymin"]<-0.1
df[which(df$dif=="up"),"ymax"]<-0.45



ggplot(data=df, aes(x=pos, y=value,group=group))+  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +geom_line(aes(color=group))+ geom_point(aes(color=group))+theme_classic()+xlab("Position")+ylab("Average reactivity")+theme(axis.title.x = element_text(size=13,face = "bold"),axis.title.y = element_text(size=13,face = "bold"),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13),legend.position="top")+facet_grid(cols = vars(dif),scale="free")+ geom_hline(aes(yintercept=V5,group=dif,color=group),alpha=0.5)+scale_color_manual(values=c("#66c2a5","#fc8d62"))->p1
p1
ggsave("avg_react_uniq1.pdf",p1,width=12,height=3.8)

