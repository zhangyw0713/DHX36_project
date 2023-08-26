read.table("/lustre/home/zhangyw/data/Homo_info/gencode.v33.pcg.list",sep = "\t",header = F)->pcg
colnames(pcg)<-c("enst","gene","name")
pcg[,2:3]->pcg
unique(pcg)->pcg
read.table("/lustre/zhangyw/myProject/Dhx36/new_RNAseq/00_deseq2/DHX36_wt_vs_ko_raw.tsv",sep="\t",header=T)->df1
read.table("/lustre/zhangyw/myProject/Dhx36/new_RNAseq/chrRNAseq/03_dif/chrRNAseq_wt_vs_ko_raw.tsv",sep="\t",header=T)->df
merge(df1,df,by="row.names")->upc
upc[,c(1,3,9)]->upcd
colnames(upcd)<-c("gene","total","chr")
2^(-1*upcd$total)->upcd$totalfc
2^(-1*upcd$chr)->upcd$chrfc
upcd$fcfc=upcd$totalfc/upcd$chrfc
upcd->allcd
upcd[which(upcd$fcfc>1.5),]->testfc
testfc[,1]->post_up_gene
#str_split_fixed(testfc[,1],"\\.",n=2)[,1]->post_up_gene
length(post_up_gene)
upcd[which(upcd$fcfc<1/1.5),]->testfc
#str_split_fixed(testfc[,1],"\\.",n=2)[,1]->post_down_gene
testfc[,1]->post_down_gene
length(post_down_gene)

testfc[which(testfc$gene %in% post_down_gene),]->downcd


merge(allcd,pcg,by="gene")->allcd.pcg

read.table("/lustre/zhangyw/myProject/Dhx36/PAR_CLIP/03_parse/newexon/finalfile.pcgtrans.bed",sep = "\t",header = F)->bd
bd[which(bd$V14=="3UTR"),]->bdtp
allcd.pcg[which(allcd.pcg$gene %in% bdtp[,11] ),]->bd.df
#allcd.pcg[which(allcd.pcg$gene %in% bd[,11]),]->bd.df

read.table("/lustre/zhangyw/myProject/Dhx36/m6a/yth/03_merge/YTHDF1/finalfile.trans.bed",sep="\t",header=F)->df1
df1[which(df1$V14=="1 3UTR"),]->df1
bd.df[which(bd.df$gene %in% df1[,11]),]->bd.df1
bd.df[-which(bd.df$gene %in% df1[,11]),]->bd.nodf1

data.frame(bd.df1,"YTHDF1(+)")[,6:8]->df2
data.frame(bd.nodf1,"YTHDF1(-)")[,6:8]->df3
t.test(bd.df1$fcfc,bd.nodf1$fcfc)
colnames(df3)<-c("fcfc","name","group")
colnames(df2)<-c("fcfc","name","group")
rbind(df2,df3)->newdf
newdf[which(newdf$fcfc<4),]->ss


ggplot(newdf,aes(x=group,y=as.numeric(as.character(fcfc)),))+geom_violin(aes(fill = group),alpha=0.8, trim=T)+geom_boxplot(aes(fill = group),outlier.shape = NA,notch = TRUE,alpha=0.8,width=0.1)+theme_bw()+ coord_cartesian(ylim = c(0,4))+ylim(0,4)+theme(legend.position = "none",axis.title.x = element_text(size=13,face = "bold"),axis.title.y = element_text(size=13,face = "bold"),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13))+ylab("pFC(KO/WT)")+ scale_fill_manual(values=c("#F29696", "#879AE8"))

p1

ggplot(newdf,aes(x=group,y=as.numeric(as.character(fcfc)),))+geom_boxplot(aes(fill = group),outlier.shape = NA,notch = TRUE,alpha=1)+ coord_cartesian(ylim = c(0,2.5))+theme_bw()+ylim(0,2.5)+theme(legend.position = "none",axis.title.x = element_text(size=13,face = "bold"),axis.title.y = element_text(size=13,face = "bold"),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13))+ylab("pFC(KO/WT)")+ scale_fill_manual(values=c("#F29696", "#879AE8"))


newdf[which(newdf$fcfc<4),]->newdf1
t.test(newdf$fcfc~newdf$group)
boxplot(newdf$fcfc~newdf$group,outline=F)

ggplot(newdf,aes(x=fcfc,color=group))+stat_ecdf()+xlim(0,3)+theme_bw()+theme(axis.title.x = element_text(size=13,face = "bold"),axis.title.y = element_text(size=13,face = "bold"),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13))+xlab("pFC(KO/WT)")+ylab("Accumulative proportion")+annotate("text",x=1.5,y=0.5,label="DF1(+) vs DF1(-):0.029")+scale_color_brewer(palette="Dark2")->p1
p1


setwd("/lustre/zhangyw/myProject/Dhx36/new_RNAseq/chrRNAseq/03_dif/new/bd_feature_25percent")
ggsave("3UTR_YTHDF1_fold_change_vs_YTHDF1_unbound1.pdf",p1,width=6.5,height=4.5)














read.table("/lustre/zhangyw/myProject/Dhx36/m6a/miclip_293T/finalfile.trans.bed",sep="\t",header=F)->df1
df1[which(df1$V14=="1 3UTR"),]->m6a
intersect(m6a[,11],bdtp[,11])->m6a.bd
allcd.pcg[which(allcd.pcg$gene %in% m6a.bd),]->m6a.pcg

read.table("/lustre/zhangyw/myProject/Dhx36/m6a/yth/03_merge/YTHDF1/finalfile.trans.bed",sep="\t",header=F)->df1
df1[which(df1$V14=="1 3UTR"),]->df1
intersect(df1[,11],bdtp[,11])->df1.bd
allcd.pcg[which(allcd.pcg$gene %in% df1.bd),]->df1.pcg




data.frame(allcd.pcg,"Overall")[,6:8]->df1
data.frame(df1.pcg,"YTHDF1")[,6:8]->df2
data.frame(m6a.pcg,"m6A_modify")[,6:8]->df3


t.test(df1$fcfc,df2$fcfc)$p.value
t.test(df1$fcfc,df3$fcfc)$p.value
t.test(df2$fcfc,df3$fcfc)

colnames(df1)<-c("fcfc","name","group")
colnames(df2)<-c("fcfc","name","group")
colnames(df3)<-c("fcfc","name","group")


rbind(df1,df2,df3)->newdf



ggplot(newdf,aes(x=1/fcfc,color=group))+stat_ecdf()+xlim(0,4)+theme_bw()+theme(axis.title.x = element_text(size=13,face = "bold"),axis.title.y = element_text(size=13,face = "bold"),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13))+xlab("Fold change of totalRNA / Fold change of chrRNA (WT/KO)")+ylab("Accumulative proportion")+annotate("text",x=1.5,y=0.5,label="Overall vs DF1:1.26e-40")+annotate("text",x=1.5,y=0.4,label="Overall vs m6A_modify:7.21e-38")+annotate("text",x=1.5,y=0.3,label="DF1 vs m6A_modify:0.712")+scale_color_brewer(palette="Dark2")->p1
p1

ggsave("3UTR_YTHDF1orM6a_fold_change_vs_overall1.pdf",p1,width=6.5,height=4.5)
