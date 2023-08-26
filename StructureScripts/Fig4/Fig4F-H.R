setwd("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/00_analysis/01_regional_compare/change2fcfc")
read.table("react2fc2_ko_wt.tsv",header=F,sep="\t")->df
df$V5/df$V3->df$reactfc
df$V5/df$V3->df$ginifc
cor.test(log2(df$reactfc),log2(df$V7),method = "pearson")




cor.test(log2(df$V4),log2(df$V5),method="pearson")
ggplot(df, aes(x=log2(df$reactfc), y=log2(df$V7))) + stat_density_2d(geom = 'polygon', aes(alpha = ..level..),fill="#00AFBB", bins = 4) +theme_bw(base_size = 13)+geom_smooth(method="lm",color="black")+annotate("text", x = -0.3, y = 1, label = "R=-0.16, p<2.2e-16")+xlim(-1.5,0.5)+xlab("Log2(Fold change average reactivity)(KO/WT)")+ylab("pFC(KO/WT)")->p1
p1
ggsave("overall.cor_ko_wt.pdf",p1,width=7.5,height=6)



read.table("fp.react2fc2_ko_wt.tsv",header=F,sep="\t")->df

df[which(df$V4>0 & df$V5>0),]->df
df$V5/df$V3->df$reactfc
df$V5/df$V3->df$ginifc
cor.test(log2(df$reactfc),log2(df$V7),method = "pearson")


ggplot(df, aes(x=log2(df$reactfc), y=log2(df$V7))) + stat_density_2d(geom = 'polygon', aes(alpha = ..level..),fill="#1b9e77", bins = 4) +theme_bw(base_size = 13)+geom_smooth(method="lm",color="black")+annotate("text", x = -1, y = 1, label = "R=-0.05, p=0.001")+xlim(-1.5,1.3)+xlab("Log2(Fold change average reactivity)(KO/WT)")+ylab("pFC(KO/WT)")->p1
p1
ggsave("5utr.cor_ko_wt.pdf",p1,width=7.5,height=6)




read.table("tp.react2fc2_ko_wt.tsv",header=F,sep="\t")->df

df[which(df$V4>0 & df$V5>0),]->df
df$V5/df$V3->df$reactfc
df$V5/df$V3->df$ginifc
cor.test(log2(df$reactfc),log2(df$V7),method = "pearson")


ggplot(df, aes(x=log2(df$reactfc), y=log2(df$V7))) + stat_density_2d(geom = 'polygon', aes(alpha = ..level..),fill="#7570b3", bins = 4) +theme_bw(base_size = 13)+geom_smooth(method="lm",color="black")+annotate("text", x = -1, y = 1, label = "R=-0.18, p=2.75e-39")+xlim(-1.5,0.3)+xlab("Log2(Fold change average reactivity)(KO/WT)")+ylab("pFC(KO/WT)")->p1
p1
ggsave("3utr.cor_ko_wt.pdf",p1,width=7.5,height=6)





read.table("cds.react2fc2_ko_wt.tsv",header=F,sep="\t")->df

df[which(df$V4>0 & df$V5>0),]->df
df$V5/df$V3->df$reactfc
df$V5/df$V3->df$ginifc
cor.test(log2(df$reactfc),log2(df$V7),method = "pearson")


ggplot(df, aes(x=log2(df$reactfc), y=log2(df$V7))) + stat_density_2d(geom = 'polygon', aes(alpha = ..level..),fill="#d95f02", bins = 4) +theme_bw(base_size = 13)+geom_smooth(method="lm",color="black")+annotate("text", x = -1, y = 1, label = "R=-0.07, p=1.35e-06")+xlim(-1.5,1.2)+xlab("Log2(Fold change average reactivity)(KO/WT)")+ylab("pFC(KO/WT)")->p1
p1
ggsave("cds.cor_ko_wt.pdf",p1,width=7.5,height=6)



