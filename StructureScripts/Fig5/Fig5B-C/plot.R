setwd("/lustre/zhangyw/myProject/Dhx36/PAR_CLIP/03_parse/newexon/m6a_enrich")
#df1
read.table("random_df1_count2.tsv",sep = "\t",header=F)->down1
as.numeric(as.character(down1$V1))->down1$V1
mean(down1[,2]/2296)->m1
sd(down1[,2]/2296)->m1sd
matrix(c(m1,"Simulated",m1sd,552/2296,"True 3UTR",0),ncol = 3,nrow = 2,byrow = T)->m2
data.frame(m2)->df
colnames(df)<-c("value","label","sd")

as.numeric(as.character(df$value))->df$value
as.numeric(as.character(df$sd))->df$sd
df$label<-factor(df$label,levels=c("True 3UTR","Simulated"))
ggplot(df) +
  geom_bar( aes(x=label, y=value), stat="identity", width=0.8,fill="indianred", alpha=0.9,color="black") +
  geom_errorbar( aes(x=label, ymin=value-sd, ymax=value+sd), width=0.2, colour="grey20", alpha=0.9, size=0.8)+theme_bw()+theme_classic()+ylab("Proportion of peaks with m6A reader")+xlab("Group")+theme(axis.title.x = element_text(size=13,face = "bold"),axis.title.y = element_text(size=13,face = "bold"),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13))+coord_cartesian(ylim=c(0.15,0.265))+annotate("text", x=1.5, y=0.265, label="***")+annotate("segment", x=1,xend=2, y=0.26, yend=0.26)->p1
ggsave("true2simulated.3utr.df1.pdf",p1,width = 4.5,height = 5.5)

#m6a
read.table("random_m6a_count_newmiCLIP.tsv",sep = "\t",header=F)->down1
as.numeric(as.character(down1$V2))->down1$V2
mean(down1[,2]/2286)->m1
sd(down1[,2]/2286)->m1sd
matrix(c(m1,"Simulated",m1sd,438/2286,"True 3UTR",0),ncol = 3,nrow = 2,byrow = T)->m2
data.frame(m2)->df
colnames(df)<-c("value","label","sd")

as.numeric(as.character(df$value))->df$value
as.numeric(as.character(df$sd))->df$sd
df$label<-factor(df$label,levels=c("True 3UTR","Simulated"))
ggplot(df) +
  geom_bar( aes(x=label, y=value), stat="identity", width=0.8,fill="#2ca25f", alpha=0.9,color="black") +
  geom_errorbar( aes(x=label, ymin=value-sd, ymax=value+sd), width=0.2, colour="grey20", alpha=0.9, size=0.8)+theme_bw()+theme_classic()+ylab("Proportion of peaks with m6A modification")+xlab("Group")+theme(axis.title.x = element_text(size=13,face = "bold"),axis.title.y = element_text(size=13,face = "bold"),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13))+coord_cartesian(ylim=c(0.1,0.21))+annotate("text", x=1.5, y=0.21, label="***")+annotate("segment", x=1,xend=2, y=0.21, yend=0.21)->p1
p1
ggsave("true2simulated.3utr.m6a.miclip.pdf",p1,width = 4.5,height = 5.5)
