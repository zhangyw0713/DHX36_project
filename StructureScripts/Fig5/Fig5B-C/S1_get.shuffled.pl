for($i=0;$i<100;$i++)
{
	# $i=22;
	open(file,"/lustre/zhangyw/myProject/Dhx36/PAR_CLIP/03_parse/newexon/3UTR_finalfile.trans.bed");
	$out=$i."_shuffled.bed";
	open(res,">$out");
	
	while(<file>)
	{
		chomp;
		@token=split(/\t/,$_);
		open(res1,">tmp.peak.bed");
		$r1=join("\t",@token[0..5]);
		print res1 $r1."\n";
		close(res1);
		$enst=$token[9];
		$cmd1="grep $enst /lustre/home/zhangyw/data/Homo_info/region/3_utr >tmp.3utr.bed";
		# print $cmd1."\n";
		system(`$cmd1`);
		$cmd2="bedtools shuffle -i tmp.peak.bed -g /lustre/home/zhangyw/data/Homo_info/genome.size -incl tmp.3utr.bed -noOverlapping >tmp.shuffled.bed";
		# print "start\n";
		system(`$cmd2`);
		# print "done\n";
		open(f,"tmp.shuffled.bed");
		while(<f>)
		{
			chomp;
			print res $_."\n";
		}
		close(f);
		# system(`cat $out tmp.shuffled.bed >>$out`);
	}
	close(res);
}