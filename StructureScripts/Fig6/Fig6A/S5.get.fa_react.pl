#get sequence and reactivity for structure folding 
system(`mkdir m20_wt`);
open(file,"m6asite2score_wt_combined.20.bed");
open(res1,">m20_wt/bdsite_ko_uniq.react");
open(res2,">m20_wt/bdsite_uniq.fa");
open(res3,">m20_wt/transcript_uniq.list");

while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	$id=$token[0];
	$react=$token[6];
	$react=~s/\,/\t/g;
	$seq=$token[5];
	$id1{$token[0]}=$id;
	$seq1{$token[0]}=">".$id."\n".$seq;
	$react1{$token[0]}=$id."\n".$react;
	# print res3 $id."\n";
	# print res2 ">".$id."\n".$seq."\n";
	
	# print res1 $id."\n".$react."\n";
}
close(file);
while(my ($k,$v)=each %id1)
{
	print res1 $react1{$k}."\n";
	print res2 $seq1{$k}."\n";
	print res3 $v."\n";

}


close(res1);
close(res2);
close(res3);


