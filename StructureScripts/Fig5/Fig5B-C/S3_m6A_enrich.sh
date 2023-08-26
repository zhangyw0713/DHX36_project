for i in `ls /lustre/zhangyw/myProject/Dhx36/PAR_CLIP/03_parse/newexon/m6a_enrich/extend/|grep shuffled`
do

	bedtools intersect -a extend/$i -b /lustre/zhangyw/myProject/Dhx36/m6a/miclip_293T/combined_merged.bed -wo -s > m6a_overlap.bed
	num2=`cut -f 4 m6a_overlap.bed |sort|uniq|wc -l`
	out2=$i$'\t'$num2>>random_m6a_count.tsv
	echo $out2 >>random_m6a_count_newmiCLIP.tsv
done


