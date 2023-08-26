open(res,">avg_ko.tsv");
for($i=-20;$i<21;$i++)
{
	$a=0;
	$b=0;
	open(file,"positional_ko_base.tsv");
	while(<file>)
	{
		chomp;
		@token=split(/\t/,$_);
		if($token[4] eq "all")
		{
			if(($token[3] ne "NA") and ($token[1]==$i))
			{
				# print $token[1]."\n";
				$a++;
				$b=$b+$token[3];
			}
		}
	}
	close(file);
	$avg=$b/$a;
	print  res $i."\t".$avg."\n";
}
close(res);

open(res,">avg_wt.tsv");
for($i=-20;$i<21;$i++)
{
	$a=0;
	$b=0;
	open(file,"positional_wt_base.tsv");
	while(<file>)
	{
		chomp;
		@token=split(/\t/,$_);
		if($token[4] eq "all")
		{
			if(($token[3] ne "NA") and ($token[1]==$i))
			{
				# print $token[1]."\n";
				$a++;
				$b=$b+$token[3];
			}
		}
	}
	close(file);
	$avg=$b/$a;
	print  res $i."\t".$avg."\n";
}
close(res);
