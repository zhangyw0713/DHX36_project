setwd("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/00_getDRRdistribution")
read.table("/lustre/zhangyw/myProject/Dhx36/PAR_CLIP/03_parse/newexon/finalfile.pcgtrans.uniq.bed",sep="\t",header=F)->pcg
dir("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/dstruct_output_5/")->sfile
str_split_fixed(sfile,"_",2)->id
id[,1]->newid

read.table("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/localized_combined/binding.site2tPosition_full_final3.bed",sep = "\t",header = F)->pos
pcg[which(pcg$V4 %in% pos$V4),]->pcg1



setwd("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/00_getDRRdistribution")
library(ggplot2)
library(dplyr)
### Read in file
m6a.dist <- read.delim ("metagene.tsv", header = F)

utr5.m6a.dist <- m6a.dist[m6a.dist$V4 < 1, 4]
cds.m6a.dist <- m6a.dist [m6a.dist$V4 < 2 & m6a.dist$V4 >= 1, 4]-1
utr3.m6a.dist <- m6a.dist[m6a.dist$V4 >= 2, 4]-2

data.frame(utr5.m6a.dist,"utr5",1)->df1
data.frame(cds.m6a.dist,"cds",1)->df2
data.frame(utr3.m6a.dist,"utr3",1)->df3

colnames(df1)<-c("prop","region","weight")
colnames(df2)<-c("prop","region","weight")
colnames(df3)<-c("prop","region","weight")
rbind(df1,df2,df3)->ddf

ggplot(m6a.dist) + geom_density(aes(x = V4)) + xlim(0, 3) + 
  theme_bw() + geom_vline(xintercept = 1:2, col = "grey")

plot(density(m6a.dist$V4))
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")

ggplot() + geom_density(data = ddf, aes(prop, weight = weight, colour = factor(region), y = after_stat(count * mean(density) / mean(count) / n_distinct(group))))+ facet_grid(cols  = vars(region))+theme_classic() ->p1
p1
ggsave("metagene.pdf",width=10,height=5)


#---------------------------------------------DRRs gain/loss structure-------------------------------------
matrix(c((57476-42390)/57476*100,"KO_loss",42390/57476*100,"KO_gain"),nrow=2,ncol=2)->df1
data.frame(t(df1))->df2
colnames(df2)<-c("value","group")
bp<- ggplot(df2, aes(x="", y=as.numeric(as.character(value)), fill=group))+
  geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+ scale_fill_manual(values=c("#d8b365", "#5ab4ac"))+theme_minimal(base_size = 13)+ylab("DRR structure change")
bp
ggsave("DRRchange.pdf",bp,width=5,height=5)


data <- data.frame(
  name=c("KO_loss","KO_gain") ,  
  value=c(0.2625*100,100*(1-0.2625))
)

# Barplot
ggplot(data, aes(x=name, y=value)) + 
  geom_bar(stat = "identity",width = 0.7)+theme_classic()+ scale_y_continuous(expand = c(0, 0),limits = c(0,80))

ggsave("gainloss.proportion.pdf")
#--------------------------------------------length corrected---------------------------------------------------




## A simple histogram
hist(m6a.dist$V4, breaks = 200, col = "black")
m6a.dist->dist
colnames(dist)<-c("transcript","start","end","rel_location","utr5_size","cds_size","utr3_size","region")
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")

## Rescale regions
## Determine scale factor
utr5.SF <- median(dist$utr5_size, na.rm = T)/median(dist$cds_size, na.rm = T)
utr3.SF <- median(dist$utr3_size, na.rm = T)/median(dist$cds_size, na.rm = T)

# assign the regions to new dataframes
utr5.dist <- dist[dist$rel_location < 1, ]
cds.dist <- dist [dist$rel_location < 2 & dist$rel_location >= 1, ]
utr3.dist <- dist[dist$rel_location >= 2, ]

# rescale 5'UTR and 3'UTR
library("scales")
utr5.dist$rel_location <- rescale(utr5.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
utr3.dist$rel_location <- rescale(utr3.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))

# Combine the regions in a new dataframe and plot
all.regions <- c(utr5.dist$rel_location, cds.dist$rel_location, utr3.dist$rel_location)
hist.data <- hist(all.regions, breaks = 200, col = "black", xlim = c(0,3)) # plot and save to variable
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")

## Alternate representations of the metagene
## A line plot
plot(hist.data$breaks[1:length(hist.data$breaks)-1], hist.data$counts, type = 'l')
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")

## A smooth density plot
plot(density(all.regions))
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")

## An absolute distance plot
hist(dist$utr3_st, xlim = c(-500,500), breaks = 2000, col = "black")
abline(v=0, col = "red")



#------------------------------------------calculate mean reactivity-----------------------------------------


setwd("/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/00_getDRRdistribution")
read.table("metagene.onlyregional.tsv",sep = "\t",header=F)->dist
read.table("/lustre/zhangyw/myProject/Dhx36/PAR_CLIP/03_parse/newexon/peak_distribution/relative_position_full.new.tsv",sep = "\t",header=F)->bdsite



ggplot(tp) + geom_density(aes(x = V4)) + xlim(0, 3) + 
  theme_bw() + geom_vline(xintercept = 1:2, col = "grey")+facet_grid(rows = vars(V9))



mutate(dist,weight="1")->newdist
ggplot(newdist,aes(x=as.numeric(as.character(V4)),y=as.numeric(as.character(weight))))+
  stat_density2d(aes(fill=..density..), geom="bar") +
  scale_fill_gradient(low="blue", high="green")


dist[,c(4,9)]->newdist

pheatmap(newdist)


newdist[which(newdist$V9=="Only_3utr"),]->tp

ggplot(tp,aes(x = V4, y = weight)) +
  geom_bin2d(
    binwidth = c(0.01, 0.01)
  ) +
  scale_fill_gradient(low = "red", high = "yellow") 

newdist[which(newdist$V9=="Only_5utr"),]->fp

ggplot(fp,aes(x = V4, y = weight)) +
  geom_bin2d(
    binwidth = c(0.01, 0.01)
  ) +
  scale_fill_gradient(low = "red", high = "yellow") 


library(aplot)


newdist[which(newdist$V9=="Only_cds"),]->cds
bdsite[which(bdsite$V3 %in% cds$V1),]->ss
ss$V4+1->ss$V4
ggplot(ss) + geom_density(aes(x = V4)) + xlim(0, 3) + 
  theme_classic()+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + geom_vline(xintercept = 1:2, col = "grey")->cdsp1
ggplot(cds,aes(x = V4, y = weight)) +
  geom_bin2d(
    binwidth = c(0.01, 0.01)
  ) +
  scale_fill_gradient(low = "#6980E2", high = "#F8766D")+ 
  theme_bw() ->cdsp2

insert_bottom(p1,p,height=1)->pcds
ggsave("onlyCDS.DRR.distribution.pdf",pcds,height=3,width = 8)


newdist[which(newdist$V9=="Only_3utr"),]->tp
bdsite[which(bdsite$V3 %in% tp$V1),]->ss
ss$V4+2->ss$V4
ggplot(ss) + geom_density(aes(x = V4)) + xlim(0, 3) + 
  theme_classic()+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + geom_vline(xintercept = 1:2, col = "grey")->tpp1
ggplot(tp,aes(x = V4, y = weight)) +
  geom_bin2d(
    binwidth = c(0.01, 0.01)
  ) +
  scale_fill_gradient(low = "#6980E2", high = "#F8766D")+ 
  theme_bw() ->tpp2

insert_bottom(p1,p,height=1)->ptp
ggsave("only3UTR.DRR.distribution.pdf",ptp,height=3,width = 8)



newdist[which(newdist$V9=="Only_5utr"),]->fp
bdsite[which(bdsite$V3 %in% fp$V1),]->ss
ggplot(ss) + geom_density(aes(x = V4)) + xlim(0, 3) + 
  theme_classic()+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + geom_vline(xintercept = 1:2, col = "grey")->fpp1
ggplot(fp,aes(x = V4, y = weight)) +
  geom_bin2d(
    binwidth = c(0.01, 0.01)
  ) +
  scale_fill_gradient(low = "#6980E2", high = "#F8766D")+ 
  theme_bw() ->fpp2

insert_bottom(p1,p,height=1)->pfp
ggsave("only5UTR.DRR.distribution.pdf",pfp,height=3,width = 8)


fpp1 %>% insert_bottom(fpp2,height=1) %>% insert_bottom(cdsp1,height=1) %>% insert_bottom(cdsp2,height=1) %>% insert_bottom(tpp1,height=1) %>% insert_bottom(tpp2,height=1)




#------------------------------------------------------------------normalized---------------------------------------------------------------------

matrix(c((1635)/811,"5UTR",20059/1361,"CDS",34921/2232,"3UTR"),nrow=3,ncol=2)->df1

matrix(ncol=1,nrow=3)->m1

colnames(m1)<-"proportion"
rownames(m1)<-c("CDS","3UTR","5UTR")

data.frame(m1)->df
df$region<-rownames(df)


factor(df$region,levels = c("5UTR","CDS","3UTR"))->df$region

ggplot(df, aes(x = region, y=proportion,fill = region)) +
  geom_bar(stat="identity",width = 0.8)+ scale_fill_brewer(palette="Dark2") +
  theme_bw()+theme(legend.position="top")+xlab("Regions")+ylab("Normalized density")+theme(axis.title.x = element_text(size=13,face = "bold"),axis.title.y = element_text(size=13,face = "bold"),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13))->p1
p1
ggsave(p1,filename = "normalize_drr_density.pdf",width = 4,height = 5)










