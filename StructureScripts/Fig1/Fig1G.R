library(pqsfinder)
setwd("/lustre/zhangyw/myProject/Dhx36/PAR_CLIP/03_parse/newexon/regional_g4_count")

read.table("gqs.tsv",sep="\t",header = F)->gqs
matrix(nrow = nrow(gqs),ncol=11)->df1
for (i in 1:nrow(gqs)) {
  gqs[i,1]->seq
  pqs<-pqsfinder(DNAString(seq),strand="+",run_max_len = 3L,loop_min_len = 1L,loop_max_len = 21L,min_score = 1L)
  elementMetadata(pqs)->metadata
  max(score(pqs))->ss
  as.character(seq)->df1[i,1]
  ss->df1[i,2]
  which(score(pqs)==max(score(pqs)))->tt
  as.character(pqs[[tt[1]]])->frag
  metadata[[3]][tt[1]]->tract
  metadata[[9]][tt[1]]->l1
  metadata[[10]][tt[1]]->l2
  metadata[[11]][tt[1]]->l3
  as.character(gqs[i,2])->df1[i,3]
  as.character(gqs[i,3])->df1[i,4]
  as.character(gqs[i,5])->df1[i,5]
  max(l1,l2,l3)->maxl
  frag->df1[i,6]
  tract->df1[i,7]
  l1->df1[i,8]
  l2->df1[i,9]
  l3->df1[i,10]
  maxl->df1[i,11]
}
write.table(df1,"pqs_regional_score.tsv",sep="\t",row.names = F,col.names = F,quote = F)

read.table("pqs_regional_score.tsv",sep="\t",header=F)->df
head(df)
pairwise.t.test(df$V3,df$V2)$p.value
library(ggplot2)

factor(df$V2,levels=c("5UTR","CDS","3UTR"))->df$V2
ggplot(df,aes(x=V2,y=V3,))+geom_violin(aes(fill = V2),alpha=0.8, trim=T)+geom_boxplot(aes(fill = V2),outlier.shape = NA,notch = TRUE,alpha=0.8,width=0.1)+annotate("segment", x=c(1,2,1),xend=c(2,3,3), y= c(78,83,88), yend=c(78,83,88))+ scale_fill_brewer(palette="Dark2")+theme_bw()+annotate("text", x=1.5, y=80, label="P=1.43e-81")+annotate("text", x=2, y=90, label="P=1.51e-12")+annotate("text", x=2.5, y=85, label="P=4.11e-56")+ylab("GQS score")+theme(axis.title.x = element_text(size=13,face = "bold"),axis.title.y = element_text(size=13,face = "bold"),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13))+ ylim(0,92)+coord_cartesian(ylim = c(0,92))+xlab("Transcript regions")
