setwd("RBNS/00_analysis")
#--------------------------------------------------K-------------------------------------------------
read.csv("K/tables/enrichments/RHAU_enrichment_R.6mers.w_all_Zscores.txt",sep = "\t",header = T)->df

df[-1,c(1,seq(2,ncol(df)-1,2))]->df1
library(dplyr)

colnames(df1)<-c("motif","5nm","20nm","80nm","320nm","1300nm")

df1[1:10,] %>%
  pivot_longer(cols = -motif, 
               names_to = "concentration", 
               values_to = "Rvalue",
               names_transform = list(time = as.integer),
               values_drop_na = TRUE) ->df2

df2->dfk



#-------------------------------------------Li-----------------------------------
setwd("/lustre/zhangyw/forCollaborator/forKit/RBNS/00_analysis")
read.csv("Li/tables/enrichments/RHAU_enrichment_R.6mers.w_all_Zscores.txt",sep = "\t",header = T)->df
df[-1,c(1,seq(2,ncol(df)-1,2))]->df1
library(dplyr)

colnames(df1)<-c("motif","5nm","20nm","80nm","320nm","1300nm")

df1[1:10,] %>%
  pivot_longer(cols = -motif, 
               names_to = "concentration", 
               values_to = "Rvalue",
               names_transform = list(time = as.integer),
               values_drop_na = TRUE) ->df2

df2->dfl


#----------------------------------------combine------------------------------------------
"K+"->dfk$ion
"Li+"->dfl$ion
rbind(dfk,dfl)->c1
factor(c1$concentration,levels=c("5nm","20nm","80nm","320nm","1300nm"))->c1$concentration
factor(c1$motif,levels=unique(c1$motif))->c1$motif
ggplot(data=c1,mapping=aes(x=concentration,y=as.numeric(as.character(Rvalue)),group=motif))+
  geom_line(aes(color=motif,lty=ion))+
  geom_point(aes(color=motif,shape=ion))+theme_bw()+theme(axis.title.x = element_text(size=13,face = "bold"),axis.title.y = element_text(size=13,face = "bold"),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13))+ylab("R value")+ylim(c(0.5,22))->p1
p1

#------------------------------------G4 proportion------------------------------------------

read.table("g4_compare.tsv",sep="\t",header=F)->df
df[which(df[,2]=="K"),]->k
df[which(df[,2]=="Li"),]->li

library(ggplot2)
data.frame(df)->df1
log2(df$V2)->df$fold

ggplot(df,aes(x=V1,y=fold, fill=V1)) +
geom_boxplot(width=0.6) +
geom_point()+ 
geom_line(aes(group=V3,color=V3)) +
theme_bw(base_size = 13)+ylim(-1,5.5)+ylab("GQS enrichment (pulldown/input, log2)")
