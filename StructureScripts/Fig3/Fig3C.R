setwd("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/localized_combined/01_why3utr")
#tp
read.table("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/localized_combined/tp_uniq/bdsite2score_ko.bed",sep = "\t",header = F)->tpko
tpko[,c(1,2,3,4)]->tpko1
for (i in 1:nrow(tpko)) {
  tpko[i,7]->score
  str_split(score,",")->dd
  as.numeric(as.character(dd[[1]]))->dt
  mean(dt,na.rm=T)->mdt
  mdt->tpko1[i,5]
}

read.table("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/localized_combined/tp_uniq/bdsite2score_wt.bed",sep = "\t",header = F)->tpwt
tpwt[,c(1,2,3,4)]->tpwt1
for (i in 1:nrow(tpwt)) {
  tpwt[i,7]->score
  str_split(score,",")->dd
  as.numeric(as.character(dd[[1]]))->dt
  mean(dt,na.rm=T)->mdt
  mdt->tpwt1[i,5]
}
cbind(tpko1,tpwt1,"3UTR")->tpdf


#fp
read.table("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/localized_combined/fp_uniq/bdsite2score_ko.bed",sep = "\t",header = F)->fpko
fpko[,c(1,2,3,4)]->fpko1
for (i in 1:nrow(fpko)) {
  fpko[i,7]->score
  str_split(score,",")->dd
  as.numeric(as.character(dd[[1]]))->dt
  mean(dt,na.rm=T)->mdt
  mdt->fpko1[i,5]
}

read.table("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/localized_combined/fp_uniq/bdsite2score_wt.bed",sep = "\t",header = F)->fpwt
fpwt[,c(1,2,3,4)]->fpwt1
for (i in 1:nrow(fpwt)) {
  fpwt[i,7]->score
  str_split(score,",")->dd
  as.numeric(as.character(dd[[1]]))->dt
  mean(dt,na.rm=T)->mdt
  mdt->fpwt1[i,5]
}
cbind(fpko1,fpwt1,"5UTR")->fpdf


#cds
read.table("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/localized_combined/cds_uniq/bdsite2score_ko.bed",sep = "\t",header = F)->cdsko
cdsko[,c(1,2,3,4)]->cdsko1
for (i in 1:nrow(cdsko)) {
  cdsko[i,7]->score
  str_split(score,",")->dd
  as.numeric(as.character(dd[[1]]))->dt
  mean(dt,na.rm=T)->mdt
  mdt->cdsko1[i,5]
}

read.table("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/localized_combined/cds_uniq/bdsite2score_wt.bed",sep = "\t",header = F)->cdswt
cdswt[,c(1,2,3,4)]->cdswt1
for (i in 1:nrow(cdswt)) {
  cdswt[i,7]->score
  str_split(score,",")->dd
  as.numeric(as.character(dd[[1]]))->dt
  mean(dt,na.rm=T)->mdt
  mdt->cdswt1[i,5]
}
cbind(cdsko1,cdswt1,"CDS")->cdsdf

fpdf[,c(1,2,3,4,5,10,11)]->fpdf1
tpdf[,c(1,2,3,4,5,10,11)]->tpdf1
cdsdf[,c(1,2,3,4,5,10,11)]->cdsdf1
colnames(fpdf1)<-c("peak","transcript","start","end","ko","wt","region")
colnames(tpdf1)<-c("peak","transcript","start","end","ko","wt","region")
colnames(cdsdf1)<-c("peak","transcript","start","end","ko","wt","region")

read.csv("/lustre/home/zhangyw/data/bowtie2index/hg38_transcript_gencodeV33/3UTR_ENST_composition.csv",header=T)->tpcomp
read.csv("/lustre/home/zhangyw/data/bowtie2index/hg38_transcript_gencodeV33/5UTR_ENST_composition.csv",header=T)->fpcomp
read.csv("/lustre/home/zhangyw/data/bowtie2index/hg38_transcript_gencodeV33/CDS_ENST_composition.csv",header=T)->cdscomp


merge(fpdf1,fpcomp,by="transcript")->fpdf2
merge(tpdf1,tpcomp,by="transcript")->tpdf2
merge(cdsdf1,cdscomp,by="transcript")->cdsdf2
fpdf2$delta=fpdf2$ko-fpdf2$wt
tpdf2$delta=tpdf2$ko-tpdf2$wt
cdsdf2$delta=cdsdf2$ko-cdsdf2$wt
cor.test(fpdf2$delta,fpdf2$GC_content)
cor.test(fpdf2$delta,fpdf2$length,method = "pearson")
cor.test(tpdf2$delta,tpdf2$GC_content)
cor.test(tpdf2$delta,tpdf2$length)
cor.test(cdsdf2$delta,cdsdf2$GC_content)
cor.test(cdsdf2$delta,cdsdf2$length)
log(fpdf2$length)->fpdf2$loglen
log(tpdf2$length)->tpdf2$loglen
log(cdsdf2$length)->cdsdf2$loglen
lm(delta~GC_content+loglen,data=fpdf2)
lm(delta~GC_content+loglen,data=tpdf2)
lm(delta~GC_content+loglen,data=cdsdf2)




rbind(fpdf2,tpdf2,cdsdf2)->ddf2

ddf2$delta=ddf2$ko-ddf2$wt
cor.test(ddf2$delta,ddf2$GC_content)
cor.test(ddf2$delta,ddf2$length)

log(ddf2$length)->ddf2$loglen
lm(delta~GC_content+loglen,data=ddf2)->m1
barplot(m1$coefficients[2:3],ylim=c(0,0.3))


factor(ddf2$region,levels=c("5UTR","CDS","3UTR"))->ddf2$region
ggplot(ddf2,aes(x=region,y=as.numeric(as.character(GC_content))))+geom_boxplot(aes(fill = region),outlier.shape = NA,notch = TRUE,alpha=0.8,width=0.1)+scale_fill_brewer(palette="Dark2")+theme_bw()+xlab("Transcript regions")

ggplot(ddf2,aes(x=region,y=as.numeric(as.character(G_content))))+geom_boxplot(aes(fill = region),outlier.shape = NA,notch = TRUE,alpha=0.8,width=0.1)+scale_fill_brewer(palette="Dark2")+theme_bw()+xlab("Transcript regions")





read.table("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/localized_combined/01_why3utr/localize_model.tsv",sep = "\t",header=T)->df
library(pheatmap)
library(RColorBrewer)


range <- max(abs(df))
min(df)->minr

pheatmap(df,cluster_cols = F,cluster_rows = F,border_color="white" ,breaks = seq(-0.05, 0.25, length.out = 100),color = colorRampPalette(c("#8DA0CB", "white", "#F09596"))(100))










read.csv("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/localized_combined/01_why3utr/bdsite2seq_composition.csv",header=T)->local



str_split_fixed(local[,1],"_ENST",n=2)[,1]->local$peak
merge(local,ddf2,by="peak")->p1
merge(ddf2,local,by="peak")->p2
log(p2$length.x)->p2$loglen

pairwise.t.test(p2$GC_content.y,p2$region)

p2$region<-factor(p2$region,levels=(c("5UTR","CDS","3UTR")))
ggplot(p2,aes(region,GC_content.y, fill=region)) +
  geom_boxplot(outlier.shape = NA,notch = TRUE,width=0.8)  +
  theme(legend.position = "none")+theme_bw()+ylim(0.1,1)
ggsave("local_gc.pdf",width=6,height=6)

lm(delta~GC_content.y+GC_content.x+loglen,data=p2)->m1


p2[which(p2$region=="3UTR"),]->tpp2
p2[which(p2$region=="5UTR"),]->fpp2
p2[which(p2$region=="CDS"),]->cdsp2
log(tpp2$length.x)->tpp2$loglen
log(fpp2$length.x)->fpp2$loglen
log(cdsp2$length.x)->cdsp2$loglen
library(car)
lm(delta~GC_content.y+GC_content.x+loglen,data=tpp2)->m2
lm(delta~GC_content.y+GC_content.x+loglen,data=fpp2)->m3
lm(delta~GC_content.y+GC_content.x+loglen,data=cdsp2)->m4

vif(m2)
vif(m3)
vif(m4)
