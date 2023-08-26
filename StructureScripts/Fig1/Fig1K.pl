open(info,"finalfile.pcgtrans.bed");     #peak annotation file
while(<info>)
{
	chomp;
	@token=split(/\t/,$_);
	$p2r{$token[3]}=$token[13];
}
close(info);

open(motif,"enrichments/RHAU_top10");   #sequence of top 10 motifs
while(<motif>)
{
	chomp;
	@token=split(/\t/,$_);
	$motif{$token[0]}=0;
}
close(motif);

open(seq,"binding.site2seq.bed");    #peak to sequence file
while(<seq>)
{
	chomp;
	@token=split(/\t/,$_);
	$p2seq{$token[3]}=$token[11];
}
close(seq);
@motifs=keys(%motif);

open(res,">seq2motif.tsv");
while(my ($k,$v)=each %p2seq)
{
	$region=$p2r{$k};
	$seq=$v;
	$count=0;
	for($i=0;$i<10;$i++)
	{
		$mm=$motifs[$i];
		if(grep(/$mm/,$seq))
		{
			$count++;
		}
	}
	print res $k."\t".$region."\t".$count."\n";
}
close(res);