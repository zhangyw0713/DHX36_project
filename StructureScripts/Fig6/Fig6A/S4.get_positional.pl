open(file,"m6asite2score_wt_combined.20.bed");
open(res,">positional_wt_base.tsv");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	$seq=$token[5];
	$icshape=$token[6];
	$icshape=~s/-1/NA/g;
	$id=$token[1]."_".$token[2]."_".$token[3];
	@base=split(//,$seq);
	@struc=split(/\,/,$icshape);

	for($i=0;$i<@base;$i++)
	{
		${"g4".$i}{$base[$i]}++;
		$pos=$i-20;
		print res $token[0]."\t".$pos."\t".$base[$i]."\t".$struc[$i]."\tall\n";
	}

}
close(file);
close(res);

open(file,"m6asite2score_ko_combined.20.bed");
open(res,">positional_ko_base.tsv");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	$seq=$token[5];
	$icshape=$token[6];
	$icshape=~s/-1/NA/g;
	$id=$token[1]."_".$token[2]."_".$token[3];
	@base=split(//,$seq);
	@struc=split(/\,/,$icshape);

	for($i=0;$i<@base;$i++)
	{
		${"g4".$i}{$base[$i]}++;
		$pos=$i-20;
		print res $token[0]."\t".$pos."\t".$base[$i]."\t".$struc[$i]."\tall\n";
	}

}
close(file);
close(res);