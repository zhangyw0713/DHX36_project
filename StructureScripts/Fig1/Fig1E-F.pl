system(`python G4classify.py DHX36BindingSite.fa`);
open(g4,"final_result");
while(<g4>)
{
	chomp;
	@token=split(/\t/,$_);
	$key=$token[1];
	if($token[4] ne "Others")
	{
	$g4{$key}=0;
	}
}
close(g4);

open(file,"../finalfile.pcgtrans.bed");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	$region=$token[13];
	$rr{$region}++;
	if(exists($g4{$token[3]}))
	{
		$g4{$region}++;
	}
}
close(file);

open(res,">regional_g4.tsv");
while(my ($k,$v)=each %rr)
{
	$all=$v;
	$g4=$g4{$k};
	$nong4=$all-$g4;
	print res $k."\t".$g4."\t"."g4"."\n";
	print res $k."\t".$nong4."\t"."nong4"."\n";
}
close(res);