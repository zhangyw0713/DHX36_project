open(tp,"/lustre/home/zhangyw/data/Homo_info/region/sorted_3_utr.bed");
while(<tp>)
{
	chomp;
	@token=split(/\t/,$_);
	# if($token[5] eq "+")
	# {
		# $start{$token[3]}=$token[1];
		# $end{$token[3]}=$token[2];
	# }else{
		$p1{$token[3]}=$token[1];
		$p2{$token[3]}=$token[2];
	# }

}
close(tp);

open(peak,"/lustre/zhangyw/myProject/Dhx36/PAR_CLIP/03_parse/newexon/3UTR_finalfile.trans.bed");
while(<peak>)
{
	chomp;
	@token=split(/\t/,$_);
	$p2t{$token[3]}=$token[9];
}
close(peak);

for($i=21;$i<23;$i++)
{
	$input=$i."_shuffled.bed";
	$output="extend/".$input;
	open(file,"$input");
	open(res,">$output");
	while(<file>)
	{
		chomp;
		@token=split(/\t/,$_);
		$enst=$p2t{$token[3]};
		if($token[5] eq "+")
		{
			$start=$p1{$enst};
			$end=$p2{$enst};
			$pstart=$token[1]-50;
			$pend=$token[2]+50;
			if($pstart<$start)
			{
				$pstart=$start;
				$pend=$start+100;
			}
			if($pend>$end)
			{
				$pend=$end;
			}
		}else{
			$start=$p1{$enst};
			$end=$p2{$enst};
			$pstart=$token[1]-50;
			$pend=$token[2]+50;
			if($pend>$end)
			{
				$pend=$end;
				$pstart=$end-100;
			}
			if($pstart<$start)
			{
				$pstart=$start;
			}
		}
		print res $token[0]."\t".$pstart."\t".$pend."\t".$token[3]."\t".$token[4]."\t".$token[5]."\n";
	}
	close(file);
	close(res);
}