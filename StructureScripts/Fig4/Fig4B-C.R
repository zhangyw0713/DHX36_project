
#-------------------------------------------figure B--------------------------------------------------
read.table("/lustre/zhangyw/myProject/Dhx36/new_RNAseq/00_deseq2/DHX36_wt_vs_ko_raw.tsv",sep="\t",header=T)->df1
-1*df1$log2FoldChange->df1$log2Fc_ko_wt
df1[which(df1$padj<0.05 & df1$log2Fc_ko_wt>=log2(1.5)),]->up
df1[which(df1$padj<0.05 & df1$log2Fc_ko_wt<=log2(1/1.5)),]->down

setwd("/lustre/zhangyw/myProject/Dhx36/new_RNAseq/chrRNAseq/03_dif")
read.table("chrRNAseq_wt_vs_ko_raw.tsv",sep="\t",header=T)->df
-1*df$log2FoldChange->df$log2Fc_ko_wt

merge(df1,df,by="row.names")->upc
upc[,c(1,8,15)]->upcd
colnames(upcd)<-c("gene","total","chr")
2^upcd$total->upcd$totalfc
2^upcd$chr->upcd$chrfc
upcd$fcfc=upcd$totalfc/upcd$chrfc
upcd[which(upcd$fcfc>1.5),]->testfc
testfc[,1]->post_up_gene
#str_split_fixed(testfc[,1],"\\.",n=2)[,1]->post_up_gene
length(post_up_gene)
upcd[which(upcd$fcfc<1/1.5),]->testfc
#str_split_fixed(testfc[,1],"\\.",n=2)[,1]->post_down_gene
testfc[,1]->post_down_gene
length(post_down_gene)

ggplot(data=upcd, aes(x=total, y=chr)) +
  geom_point(data=subset(upcd,upcd$fcfc <= 1.5 & upcd$fcfc >=1/1.5),color="gray",alpha=0.6) +
  geom_point(data=subset(upcd,upcd$fcfc > 1.5),color="#fc8d62",alpha=0.9) +
  geom_point(data=subset(upcd,upcd$fcfc < 1/1.5),color="#94C47D",alpha=0.9) +
  xlab("FC in whole cell (KO/WT)") + ylab("FC in chromatin fraction (KO/WT)") +   #定义X轴和Y轴的名称
  theme_classic(base_size=15)+   #去除灰色背景并设置字体大小
  theme(panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+ #去除背景格线
  theme(plot.title = element_text(size=15,hjust = 0.5),legend.position='none')+xlim(-4,4)+ylim(-4,4)+
  annotate(geom="text", x=2.5, y=-3, label="Post-transcriptionally upregulated genes",color="#fc8d62")+annotate(geom="text", x=-2, y=3, label="Post-transcriptionally downregulated genes",color="#94C47D")
 ggsave("fc_choose_ko_wt.pdf",width=8.3,height=8)
 
 
#--------------------------------------------Figure C--------------------------------------------- 
setwd("/lustre/zhangyw/myProject/Dhx36/new_RNAseq/chrRNAseq/03_dif")
read.table("/lustre/zhangyw/myProject/Dhx36/PAR_CLIP/03_parse/newexon/finalfile.trans.bed",sep = "\t",header=F)->bd
bd[,11]->bd.gene
allcd.pcg[which(allcd.pcg$gene %in% bd.gene),]->bd.fc
t.test(bd.fc$fcfc,allcd.pcg$fcfc)$p.value
allcd.pcg[,c(1,4,5,6,7)]->allcd.pcgnew
write.table(allcd.pcgnew,"fcfc_ko_wt.tsv",sep = "\t",row.names = F,col.names = F,quote = F)

data.frame(allcd.pcg,"Overall")[,6:8]->df1
data.frame(bd.fc,"DHX36-bound")[,6:8]->df2
colnames(df1)<-c("fcfc","name","group")
colnames(df2)<-c("fcfc","name","group")
rbind(df1,df2)->newdf

ggplot(newdf,aes(x=fcfc,color=group))+stat_ecdf()+xlim(0,2)+theme_bw()+theme(axis.title.x = element_text(size=13,face = "bold"),axis.title.y = element_text(size=13,face = "bold"),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13))+xlab("pFC (KO/WT)")+ylab("Accumulative proportion")+annotate("text",x=1.4,y=0.25,label="Overall vs DHX36-binding:2.24e-26")+scale_color_brewer(palette="Dark2")->p1
p1
ggsave("overall_dhx36-binding_fold_change.pdf",p1,width=6.5,height=4.5)