
open(shape,"/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/wt_combined_format.out");
while(<shape>)
{
	chomp;
	@token=split(/\t/,$_);
	$id=$token[0];
	$score=$token[1];
	$t2s{$id}=$score;
}
close(shape);



open(file,"m6a.site2seq.20.bed");
open(res,">m6asite2score_wt_combined.20.bed");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	$tid=$token[7];
	$start=$token[8];
	$end=$token[9];
	$null=0;
	if(exists($t2s{$tid}))
	{
		$score=$t2s{$tid};
		@ss=split(/\,/,$score);
		@seg=@ss[$start-1..$end-1];
		for($a=0;$a<@seg;$a++)
		{
			if($seg[$a] == 0)
			{
				$null++;
			}
		}
		# if($null<95)
		# {
			$segseq=join(",",@seg);
			@token1=split(/\=/,$token[3]);
			@token2=split(/\]/,$token1[1]);
			$ph=$token2[0];
			print res $token[3]."\t".$tid."\t".$start."\t".$end."\t".$ph."\t".$token[10]."\t".$segseq."\n";
		# }
	}

}
close(file);
close(res);





open(shape,"/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/ko_combined_format.out");
while(<shape>)
{
	chomp;
	@token=split(/\t/,$_);
	$id=$token[0];
	$score=$token[1];
	$t2s1{$id}=$score;
}
close(shape);



open(file,"m6a.site2seq.20.bed");
open(res,">m6asite2score_ko_combined.20.bed");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	$tid=$token[7];
	$start=$token[8];
	$end=$token[9];
	$null=0;
	if(exists($t2s1{$tid}))
	{
		$score=$t2s1{$tid};
		@ss=split(/\,/,$score);
		@seg=@ss[$start-1..$end-1];
		for($a=0;$a<@seg;$a++)
		{
			if($seg[$a] == 0)
			{
				$null++;
			}
		}
		# if($null<95)
		# {
			$segseq=join(",",@seg);
			@token1=split(/\=/,$token[3]);
			@token2=split(/\]/,$token1[1]);
			$ph=$token2[0];
			print res $token[3]."\t".$tid."\t".$start."\t".$end."\t".$ph."\t".$token[10]."\t".$segseq."\n";
		# }
	}

}
close(file);
close(res);
