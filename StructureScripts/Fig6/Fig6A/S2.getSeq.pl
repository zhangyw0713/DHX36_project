open(file,"/lustre/zhangyw/myProject/Dhx36/icshape/0000_pipe_Res/GTF/gencodev33_transcriptome.fa");
local $/=">";
while(<file>)
{
	chomp;
	@token=split(/\n/,$_);
	$id=$token[0];
	shift @token;
	# $seq=~s//\n/g;
	$newseq=join("",@token);
	@head=split(/\t/,$id);
	$tid=$head[0];
	# print $newseq."\n";
	$hash{$tid}=$newseq;
	
}
close(file);

open(file,"m6a.site2tPosition.20.bed");
local $/="\n";
open(res,">m6a.site2seq.20.bed");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	$tid=$token[7];
	$start=$token[8]-1;
	$end=$token[9]-1;
	# print $tid."\n";
	if(exists($hash{$tid}))
	{
		@seq=split(//,$hash{$tid});
		# print $tid."\t".$seq[1]."\n";
		print res $_."\t";
		for($a=$start;$a<$end+1;$a++)
		{
			print res $seq[$a];
		}
		print res "\n";
	}
}
close(res);
close(file);










# open(info,"/lustre/home/zhangyw/data/bowtie2index/hg38_transcript_gencodeV33/genelist.txt");
# while(<info>)
# {
	# chomp;
	# @token=split(/\t/,$_);
	# $num2tid{$token[0]}=$token[1];
# }
# close(info);

# open(seq,"/lustre/home/zhangyw/data/bowtie2index/hg38_transcript_gencodeV33/hg38_transcript_reform.fa
# ");
# while(<seq>)
# {
	# chomp;
	# @token=split(/\t/,$_);
	# $tid=$num2tid{$token[0]};
	# $hash{$tid}=$token[1];
# }
# close(seq);

# open(file,"binding.site2tPosition.bed");
# open(res,">binding.site2seq_forG4_noCom.bed");
# while(<file>)
# {
	# chomp;
	# @token=split(/\t/,$_);
	# $tid=$token[6];
	# @token1=split(//,$hash{$tid});
	# $strand=$token[5];
	# print res $_."\t";
	# if($strand eq "+")
	# {
		# for($i=$token[7]-1;$i<$token[8];$i++)
		# {
			# if($token1[$i] eq "T")
			# {
				# print res "U";
			# }else{
			# print res $token1[$i];
			# }
		# }

	# }else{
		# for($i=$token[8]-1;$i>$token[7]-2;$i=$i-1)
		# {
			# if($token1[$i] eq "T")
			# {
				# print res "U";
			# }else{
			# print res $token1[$i];
			# }
		# }
	# }
	# print res "\n";
# }
# close(file);
# close(res);
