for f in Dhx36_E335A_1 Dhx36_E335A_2; do
perl /lustre/home/zhangyw/bin/ctk/fastq_filter.pl -v -if sanger -f mean:0-24:20 -of fastq ../$f.fastq - | gzip -c > $f.filter.fastq.gz
done

for f in Dhx36_E335A_1 Dhx36_E335A_2; do
cutadapt -f fastq -n 1 --quality-cutoff 5 -m 20 -a TAATATCGTATGCCGTCTTCTGCTTG -o $f.filter.trim.fastq.gz ../$f.fastq > $f.cutadpt.log
done

for f in Dhx36_E335A_1 Dhx36_E335A_2; do
perl /lustre/home/zhangyw/bin/ctk/fastq2collapse.pl ../00_trim/$f.filter.trim.fastq.gz - | gzip -c > $f.filter.trim.c.fastq.gz
done

for f in Dhx36_E335A_1 Dhx36_E335A_2; do
bwa aln -t 4 -n 0.06 -q 20 /lustre/home/zhangyw/data/bwaIndex/hg38/hg38_gencodeV33.fa ../01_collapse/$f.filter.trim.c.fastq.gz > $f.sai
bwa samse /lustre/home/zhangyw/data/bwaIndex/hg38/hg38_gencodeV33.fa $f.sai ../01_collapse/$f.filter.trim.c.fastq.gz | gzip -c > $f.sam.gz
done

for f in Dhx36_E335A_1 Dhx36_E335A_2; do
perl /lustre/home/zhangyw/bin/ctk/parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file $f.mutation.txt ../02_align/$f.sam.gz - | gzip -c > $f.tag.bed.gz
done

for f in Dhx36_E335A_1 Dhx36_E335A_2; do
python /lustre/home/zhangyw/bin/ctk/joinWrapper.py $f.mutation.txt $f.tag.uniq.bed 4 4 N $f.tag.uniq.mutation.txt
done

for f in Dhx36_E335A_1 Dhx36_E335A_2; do
perl /lustre/home/zhangyw/bin/ctk/tag2collapse.pl -v -big -weight --weight-in-name --keep-max-score --keep-tag-name $f.tag.bed $f.tag.uniq.bed
done

perl /lustre/home/zhangyw/bin/ctk/bed2rgb.pl -v -col "188,0,0" Dhx36_E335A_1.tag.uniq.bed Dhx36_E335A_1.tag.uniq.rgb.bed
perl /lustre/home/zhangyw/bin/ctk/bed2rgb.pl -v -col "188,0,0" Dhx36_E335A_2.tag.uniq.bed Dhx36_E335A_2.tag.uniq.rgb.bed

perl /lustre/home/zhangyw/bin/ctk/bed2annotation.pl -dbkey hg38 -ss -big -region -v -summary Dhx36_E335A_2.tag.uniq.annot.summary.txt Dhx36_E335A_2.tag.uniq.rgb.bed Dhx36_E335A_2.tag.uniq.annot.txt
perl /lustre/home/zhangyw/bin/ctk/bed2annotation.pl -dbkey hg38 -ss -big -region -v -summary Dhx36_E335A_1.tag.uniq.annot.summary.txt Dhx36_E335A_1.tag.uniq.rgb.bed Dhx36_E335A_1.tag.uniq.annot.txt


#CIMS rep1
perl /lustre/home/zhangyw/bin/ctk/getMutationType.pl -t t2c --summary Dhx36_E335A_1.tag.uniq.mutation.stat.txt Dhx36_E335A_1.tag.uniq.mutation.txt Dhx36_E335A_1.tag.uniq.t2c.bed
perl /lustre/home/zhangyw/bin/ctk/CIMS.pl -big -n 10 -p -outp Dhx36_E335A_1.tag.uniq.t2c.pos.stat.txt -v Dhx36_E335A_1.tag.uniq.rgb.bed Dhx36_E335A_1.tag.uniq.t2c.bed Dhx36_E335A_1.tag.uniq.t2c.CIMS.txt
awk '{if($9<=0.05) {print $0}}' Dhx36_E335A_1.tag.uniq.t2c.CIMS.txt | sort -k 9,9n -k 8,8nr -k 7,7n > Dhx36_E335A_1.tag.uniq.t2c.CIMS.s13.txt
cut -f 1-6 Dhx36_E335A_1.tag.uniq.t2c.CIMS.s13.txt > Dhx36_E335A_1.tag.uniq.t2c.CIMS.s13.bed

#CIMS rep2
perl /lustre/home/zhangyw/bin/ctk/getMutationType.pl -t t2c --summary Dhx36_E335A_2.tag.uniq.mutation.stat.txt Dhx36_E335A_2.tag.uniq.mutation.txt Dhx36_E335A_2.tag.uniq.t2c.bed
perl /lustre/home/zhangyw/bin/ctk/CIMS.pl -big -n 10 -p -outp Dhx36_E335A_2.tag.uniq.t2c.pos.stat.txt -v Dhx36_E335A_2.tag.uniq.rgb.bed Dhx36_E335A_2.tag.uniq.t2c.bed Dhx36_E335A_2.tag.uniq.t2c.CIMS.txt
awk '{if($9<=0.05) {print $0}}' Dhx36_E335A_2.tag.uniq.t2c.CIMS.txt | sort -k 9,9n -k 8,8nr -k 7,7n > Dhx36_E335A_2.tag.uniq.t2c.CIMS.s13.txt
cut -f 1-6 Dhx36_E335A_2.tag.uniq.t2c.CIMS.s13.txt > Dhx36_E335A_2.tag.uniq.t2c.CIMS.s13.bed


#get overlapped crosslink sites
bedtools intersect -a /lustre/zhangyw/myProject/Dhx36/PAR_CLIP/03_parse/rep2/Dhx36_E335A_2.tag.uniq.t2c.CIMS.s13.bed -b /lustre/zhangyw/myProject/Dhx36/PAR_CLIP/03_parse/rep1/Dhx36_E335A_1.tag.uniq.t2c.CIMS.s13.bed -wo >Dhx36_E355A_repOverlap_t2c_CIMS.s13.bed


#annotation
bedtools intersect -a Dhx36_E355A_repOverlap_t2c_CIMS.s13.bed -b /lustre/home/zhangyw/data/Homo_info/region/sorted_3_utr.bed -s -wao |perl -lane 'print if $F[-1] != 0' | awk '{print $0,"3UTR"}'> annotated_result_3utr_2.bed
bedtools intersect -a Dhx36_E355A_repOverlap_t2c_CIMS.s13.bed -b /lustre/home/zhangyw/data/Homo_info/region/cds_sorted_hg38.bed -s -wao |perl -lane 'print if $F[-1] != 0' | awk '{print $0,"CDS"}'> annotated_result_cds_2.bed
bedtools intersect -a Dhx36_E355A_repOverlap_t2c_CIMS.s13.bed -b /lustre/home/zhangyw/data/Homo_info/region/sorted_5_utr.bed -s -wao |perl -lane 'print if $F[-1] != 0' | awk '{print $0,"5UTR"}'> annotated_result_5utr_2.bed
cat annotated_result_5utr_2.bed  annotated_result_3utr_2.bed annotated_result_cds_2.bed > finalfile2.bed


