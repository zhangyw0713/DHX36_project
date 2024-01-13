#Script for structure-seq data analysis
#Note: the commands for WT samples were shown here. When dealing with KO samples, please revise the filename. 

#Quality control
fastqc -t 6 file.fq

#Adapter removal
cutadapt -q 30 -m 20 -n 1 -a NNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o WT-NAI-rep1_1_trimmed.fq WT-NAI-rep1_L3_1.fq 2>wt_nai_rep1.log
cutadapt -q 30 -m 20 -n 1 -a NNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o WT-NAI-rep2_1_trimmed.fq WT-NAI-rep2_L2_1.fq 2>wt_nai_rep2.log
cutadapt -q 30 -m 20 -n 1 -a NNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o WT-DMSO-rep1_1_trimmed.fq WT-DMSO-rep1_L3_1.fq 2>wt_dmso_rep1.log
cutadapt -q 30 -m 20 -n 1 -a NNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o WT-DMSO-rep2_1_trimmed.fq WT-DMSO-rep2_L2_1.fq 2>wt_dmso_rep2.log
cutadapt -q 30 -m 20 -n 1 -a NNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o KO-NAI-rep1_1_trimmed.fq KO-NAI-rep1_L3_1.fq 2>ko_nai_rep1.log
cutadapt -q 30 -m 20 -n 1 -a NNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o KO-NAI-rep2_1_trimmed.fq KO-NAI-rep2_L2_1.fq 2>ko_nai_rep2.log
cutadapt -q 30 -m 20 -n 1 -a NNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o KO-DMSO-rep1_1_trimmed.fq KO-DMSO-rep1_L3_1.fq 2>ko_dmso_rep1.log
cutadapt -q 30 -m 20 -n 1 -a NNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o KO-DMSO-rep2_1_trimmed.fq KO-DMSO-rep2_L2_1.fq 2>ko_dmso_rep2.log

#mapping to transcriptome reference
bowtie2 -a -p 12 -x transcriptome_reference -q WT-DMSO-rep1_1_trimmed.fq >WT_DMSO_rep1.sam 2> WT_DMSO_rep1.log
bowtie2 -a -p 12 -x transcriptome_reference -q WT-DMSO-rep2_1_trimmed.fq >WT_DMSO_rep2.sam 2> WT_DMSO_rep2.log
bowtie2 -a -p 12 -x transcriptome_reference -q WT-NAI-rep1_1_trimmed.fq >WT_NAI_rep1.sam 2> WT_NAI_rep1.log
bowtie2 -a -p 12 -x transcriptome_reference -q WT-NAI-rep2_1_trimmed.fq >WT_NAI_rep2.sam 2> WT_NAI_rep2.log
bowtie2 -a -p 12 -x transcriptome_reference -q KO-DMSO-rep1_1_trimmed.fq >KO_DMSO_rep1.sam 2> KO_DMSO_rep1.log
bowtie2 -a -p 12 -x transcriptome_reference -q KO-DMSO-rep2_1_trimmed.fq >KO_DMSO_rep2.sam 2> KO_DMSO_rep2.log
bowtie2 -a -p 12 -x transcriptome_reference -q KO-NAI-rep1_1_trimmed.fq >KO_NAI_rep1.sam 2> KO_NAI_rep1.log
bowtie2 -a -p 12 -x transcriptome_reference -q KO-NAI-rep2_1_trimmed.fq >KO_NAI_rep2.sam 2> KO_NAI_rep2.log

#filter sam
python /lustre/home/zhangyw/bin/StructureFold2-master/sam_filter.py

#count RT stop
python /lustre/home/zhangyw/bin/StructureFold2-master/sam_to_rtsc.py transcriptome_reference.fasta

#RT stop correlation
python /lustre/home/zhangyw/bin/StructureFold2-master/rtsc_correlation.py  ko_dmso_rep1_filtered.rtsc ko_dmso_rep2_filtered.rtsc ko_nai_rep1_filtered.rtsc ko_nai_rep2_filtered.rtsc wt_dmso_rep1_filtered.rtsc wt_dmso_rep2_filtered.rtsc wt_nai_rep1_filtered.rtsc wt_nai_rep2_filtered.rtsc 

#Calculate RT coverage of each transcript
python /lustre/home/zhangyw/bin/StructureFold2-master/rtsc_coverage.py gencodev33_transcriptome.fa -f  ko_nai_rep1_filtered.rtsc ko_nai_rep2_filtered.rtsc wt_nai_rep1_filtered.rtsc wt_nai_rep2_filtered.rtsc -ol -bases AGCT -ot 1 -on NAI.overlap_1.0.lite.txt


#Merge replicates
python /lustre/home/zhangyw/bin/StructureFold2-master/rtsc_combine.py wt_dmso_rep1_filtered.rtsc wt_dmso_rep2_filtered.rtsc -name wt_dmso_combined.rtsc
python /lustre/home/zhangyw/bin/StructureFold2-master/rtsc_combine.py wt_nai_rep1_filtered.rtsc wt_nai_rep2_filtered.rtsc -name wt_nai_combined.rtsc
python /lustre/home/zhangyw/bin/StructureFold2-master/rtsc_combine.py ko_dmso_rep1_filtered.rtsc ko_dmso_rep2_filtered.rtsc -name ko_dmso_combined.rtsc
python /lustre/home/zhangyw/bin/StructureFold2-master/rtsc_combine.py ko_nai_rep1_filtered.rtsc ko_nai_rep2_filtered.rtsc -name ko_nai_combined.rtsc


#Calculate SHAPE reactivity. A library size correction factor should be added to line62 of rtsc_to_react.py.
python /lustre/home/zhangyw/bin/StructureFold2-master/rtsc_to_react.py wt_dmso_combined.rtsc wt_nai_combined.rtsc gencodev33_transcriptome.fa -bases ATCG -name wt_combined.react
python /lustre/home/zhangyw/bin/StructureFold2-master/rtsc_to_react.py ko_dmso_combined.rtsc ko_nai_combined.rtsc gencodev33_transcriptome.fa -bases ATCG -scale wt_combined.scale -name ko_combined.react 


#Splice structure-seq data by transcript regions
python /lustre/home/zhangyw/bin/StructureFold2-master/splice_reacts_by_FASTA.py wt_combined.react 5UTR_ENST.fa CDS_ENST.fa 3UTR_ENST.fa &
python /lustre/home/zhangyw/bin/StructureFold2-master/splice_reacts_by_FASTA.py ko_combined.react 5UTR_ENST.fa CDS_ENST.fa 3UTR_ENST.fa &


#Calculate average reactivity and Gini index
python /lustre/home/zhangyw/bin/StructureFold2-master/react_statistics.py  -react wt_combined.react -restrict NAI.overlap_1.0.lite.txt -name statistic_wt_full.sig.react 
python /lustre/home/zhangyw/bin/StructureFold2-master/react_statistics.py  -react wt_combined_CDS.react -restrict NAI.overlap_1.0.lite.txt -name statistic_wt_CDS.sig.react 
python /lustre/home/zhangyw/bin/StructureFold2-master/react_statistics.py  -react wt_combined_tpUTR.react -restrict NAI.overlap_1.0.lite.txt -name statistic_wt_tp.sig.react 
python /lustre/home/zhangyw/bin/StructureFold2-master/react_statistics.py  -react wt_combined_fpUTR.react -restrict NAI.overlap_1.0.lite.txt -name statistic_wt_fp.sig.react 
python /lustre/home/zhangyw/bin/StructureFold2-master/react_statistics.py  -react ko_combined.react -restrict NAI.overlap_1.0.lite.txt -name statistic_ko_full.sig.react 
python /lustre/home/zhangyw/bin/StructureFold2-master/react_statistics.py  -react ko_combined_CDS.react -restrict NAI.overlap_1.0.lite.txt -name statistic_ko_CDS.sig.react 
python /lustre/home/zhangyw/bin/StructureFold2-master/react_statistics.py  -react ko_combined_tpUTR.react -restrict NAI.overlap_1.0.lite.txt -name statistic_ko_tp.sig.react 
python /lustre/home/zhangyw/bin/StructureFold2-master/react_statistics.py  -react ko_combined_fpUTR.react -restrict NAI.overlap_1.0.lite.txt -name statistic_ko_fp.sig.react 
