use List::Util qw(sum);
$bed="cut -f 1-6 ../m6A.clipsite1.bed >m6a.site && bedtools intersect -a m6a.site -b /lustre/home/zhangyw/data/Homo_info/region/exon.final.bed -wao -s >m6a.site2exon.bed";
system(`$bed`);

open(info,"/lustre/home/zhangyw/data/Homo_info/gtf2bed/gencode.v33.annotation.bed");
while(<info>)
{
	chomp;
	@token=split(/\t/,$_);
	$id2info{$token[3]}=$_;
}
close(info);

open(file,"m6a.site2exon.bed");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	$hash{$token[3]."\t".$token[9]}=$token[1]."\t".$token[2]."\t".$token[11]."\t".$token[9]."\t".$token[0];
}
close(file);
open(res,">m6a.site2tPosition.20.bed");
while(my ($k,$v)=each %hash)
{
	@token=split(/\t/,$k);
	$tid=$token[1];
	$peak_id=$token[0];
	$info=$id2info{$tid};
	@token1=split(/\t/,$v);
	$chr=$token1[4];
	$strand=$token1[2];
	@token2=split(/\t/,$info);
	$exonlength=$token2[10];
	@length=split(/\,/,$exonlength);
	$fulllength=sum(@length);
	$exonstart=$token2[11];
	@estart=split(/\,/,$exonstart);
	# for($)
	if($fulllength>40)
	{
		if($strand eq "+")
		{

			$gene_start=$token2[1];
			$gene_end=$token2[2];
			$POS_u=$token1[0];
			$bd_start=$token1[0]-20;
			$bd_end=$token1[0]+20;
			$exon=0;
			for($a=0;$a<@estart;$a++)
			{
				$exon_start=$estart[$a]+$gene_start;
				$exon_end=$exon_start+$length[$a];
				if($POS_u>=$exon_start and $POS_u<=$exon_end)
				{
					if($bd_start<$exon_start)
					{
						# $peakstart=$exon+1;
						$peakend=$fulllength+1;
					}else{
						$peakstart=$exon+$bd_start-$exon_start+1;
						$peakend=$peakstart+40;
					}
					
					# if($peakend > $fulllength)
					# {
						# $peakend=$fulllength;
						# $peakstart=$fulllength-100;
					# }
					
				}
				$exon=$exon+$length[$a];
			}
			if($peakend<=$fulllength)
			{
				$peakstart=$peakstart-1;
				$peakend=$peakend-1;
				print res $chr."\t".$token1[0]."\t".$token1[1]."\t".$k."\t"."\.\t".$strand."\t".$tid."\t".$peakstart."\t".$peakend."\n";
			}
		}else{
			$gene_start=$token2[1];
			$gene_end=$token2[2];
			$POS_u=$token1[0];
			$bd_start=$POS_u-20;
			$bd_end=$POS_u+20;
			$exon=0;
			for($a=0;$a<@estart;$a++)
			{
				$exon_start=$estart[$a]+$gene_start;
				$exon_end=$exon_start+$length[$a];
				if($POS_u>=$exon_start and $POS_u<=$exon_end)
				{
					if($bd_start<=$gene_start)
					{
						# $peakstart=$exon+1;
						# $peakend=$peakstart+100;
						# $peakend=$fulllength;
						$peakstart=0;
					}else{
						$peakend=$fulllength-($exon+$bd_start-$exon_start)+1;
						$peakstart=$peakend-40;
					}
					
					# if($peakstart <1)
					# {
						# $peakstart=1;
						# $peakend=101;
					# }
					
				}
				$exon=$exon+$length[$a];
			}
			# $peakstart=$fulllength-$peakend+1;
			# $peakend=$peakstart+100;
			# if($peakend>$fulllength)
			# {
				# $peakend=$fulllength;
			# }
			if($peakstart>=1)
			{
			$peakstart=$peakstart-1;
				$peakend=$peakend-1;
			print res $chr."\t".$token1[0]."\t".$token1[1]."\t".$k."\t"."\.\t".$strand."\t".$tid."\t".$peakstart."\t".$peakend."\n";
			}	
		}
	}
	
}
close(res);