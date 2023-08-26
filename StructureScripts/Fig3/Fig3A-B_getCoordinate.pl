use List::Util qw(sum);
$bed="cut -f 1-6 ../finalfile.pcgtrans.bed >binding2.site && bedtools intersect -a binding2.site -b /lustre/home/zhangyw/data/Homo_info/region/exon.final.bed -wao -s >binding.site2exon_full.bed";
system(`$bed`);

open(info,"/lustre/home/zhangyw/data/Homo_info/gtf2bed/gencode.v33.annotation.bed");
while(<info>)
{
	chomp;
	@token=split(/\t/,$_);
	$id2info{$token[3]}=$_;
}
close(info);

open(file,"binding.site2exon.bed");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	$hash{$token[3]."\t".$token[9]}=$token[1]."\t".$token[2]."\t".$token[11]."\t".$token[9]."\t".$token[0];
}
close(file);
open(res,">binding.site2tPosition.bed");
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
	# if($tid eq "ENST00000388968.8")
	# {
		# print $fulllength."\n";
	# }
	# for($)
	if($fulllength>100)
	{
		if($strand eq "+")
		{

			$gene_start=$token2[1];
			$gene_end=$token2[2];
			$POS_u=$token1[0];
			# print $POS_u."\n";
			$bd_start=$token1[0]-50;
			$bd_end=$token1[0]+50;
			$exon=0;
			for($a=0;$a<@estart;$a++)
			{
				$exon_start=$estart[$a]+$gene_start;
				$exon_end=$exon_start+$length[$a];
				if($POS_u>=$exon_start and $POS_u<=$exon_end)
				{
					if($bd_start<$gene_start)
					{
						$peakstart=$exon+1;
						$peakend=$peakstart+100;
						$distance=$POS_u-$exon_start;
						
						# print $distance."\n";
						
					}else{
						$peakstart=$exon+$bd_start-$exon_start+1;
						$peakend=$peakstart+100;
						$distance=50;
						# print  $chr."\t".$token1[0]."\t".$token1[1]."\t".$k."\t"."\.\t".$strand."\t".$tid."\t".$peakstart."\t".$peakend."\t$distance\n";
						# print $distance."\n";
					}
					
					if($peakend > $fulllength)
					{
						$distance=$distance+($peakend-$fulllength);
						# print  $chr."\t".$token1[0]."\t".$token1[1]."\t".$k."\t"."\.\t".$strand."\t".$tid."\t".$peakstart."\t".$peakend."\t$distance\n";
						$peakend=$fulllength;
						$peakstart=$fulllength-100;
						
						# print  $chr."\t".$token1[0]."\t".$token1[1]."\t".$k."\t"."\.\t".$strand."\t".$tid."\t".$peakstart."\t".$peakend."\t$distance\n"
					}
					
				}
				$exon=$exon+$length[$a];
				
			}
			print res $chr."\t".$token1[0]."\t".$token1[1]."\t".$k."\t"."\.\t".$strand."\t".$tid."\t".$peakstart."\t".$peakend."\t$distance\n";
		}else{
			$gene_start=$token2[1];
			$gene_end=$token2[2];
			$POS_u=$token1[0];
			$bd_start=$POS_u-50;
			$bd_end=$POS_u+50;
			$exon=0;
			for($a=0;$a<@estart;$a++)
			{
				$exon_start=$estart[$a]+$gene_start;
				$exon_end=$exon_start+$length[$a];
				if($POS_u>=$exon_start and $POS_u<=$exon_end)
				{
					if($bd_start<$gene_start)
					{
						$peakend=$fulllength-1;
						$peakstart=$peakend-100;
						$distance=50+($gene_start-$bd_start);

						# print  $chr."\t".$token1[0]."\t".$token1[1]."\t".$k."\t"."\.\t".$strand."\t".$tid."\t".$peakstart."\t".$peakend."\t$distance\n";

					}else{
						$peakend=$fulllength-($exon+$bd_start-$exon_start);
						$peakstart=$peakend-100;
						$distance=50;
						# print $peak_id."\t".$tid."\n";
					}
					
					if($peakstart <1)
					{
						$distance=$distance-(1-$peakstart);
						
						$peakstart=1;
						$peakend=101;

					}
					
				}
				$exon=$exon+$length[$a];
			}

			print res $chr."\t".$token1[0]."\t".$token1[1]."\t".$k."\t"."\.\t".$strand."\t".$tid."\t".$peakstart."\t".$peakend."\t$distance\n";
				
		}
	}
	
}
close(res);
 



