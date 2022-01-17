library(pqsfinder)
setwd("D:\\research\\DHX36\\pgsfinder")

read.table("mouse_gqs_200.txt",sep="\t",header = F)->gqs
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
write.table(df1,"mouse_pqs_score_extend_3L_1L_200.tsv",sep="\t",row.names = F,col.names = F,quote = F)





read.table("pqs_score_extend.tsv",sep = "\t",header = F)->df1
df1[which(df1$V7>4),]->df2
max(as.numeric(as.character(df1[,11])))
min(df1[,11])
quantile(df1[,11])
pqsfinder(DNAString("CTGGGCAAAGGAAATGACAAGGGGACGGGGTCT"),strand="+",run_max_len = 3L,loop_min_len = 1L,loop_max_len = 21L,min_score = 1L)
pqsfinder(DNAString("AGGCGACGGTGGGGAAGATGGCGTACCAGAGCTTGCGGCTGGAGTACCTGCAGATCCCACCGGTCAGCCGCGCCTACACCACTGCCTGCGTCCTCACCACC"),run_max_len = 3L,min_score = 1L)->pqs
pqs
elementMetadata(pqs)
max(score(pqs))->ss
which(score(pqs)==max(score(pqs)))->tt
as.character(pqs[[tt[1]]])->frag
metadata[[3]][tt[1]]->tract
metadata[[9]][tt[1]]->l1
metadata[[10]][tt[1]]->l2
metadata[[11]][tt[1]]->l3
as.character(gqs[i,6])->df1[i,3]
as.character(gqs[i,2])->df1[i,4]
as.character(gqs[i,5])->df1[i,5]
max(l1,l2,l3)->maxl
frag->df1[i,6]
tract->df1[i,7]
l1->df1[i,8]
l2->df1[i,9]
l3->df1[i,10]
maxl->df1[i,11]