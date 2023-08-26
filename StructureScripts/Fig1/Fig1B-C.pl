open(file,"/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/localized_combined/final_result");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	$hash{$token[1]."\t".$token[2]}=$token[0];
}
close(file);
open(fp,"/lustre/home/zhangyw/data/bowtie2index/hg38_transcript_gencodeV33/5UTR.fa");
local $/=">";
while(<fp>)
{
	chomp;
	@token=split(/\n/,$_);
	$head=$token[0];
	@token1=split(/\s/,$head);
	$id=$token1[1];
	
	shift @token;
	$seq=join("",@token);
	$fp{$id}=$seq;
	# print $id."\t".$seq."\n";
}
close(fp);
open(tp,"/lustre/home/zhangyw/data/bowtie2index/hg38_transcript_gencodeV33/3UTR.fa");
local $/=">";
while(<tp>)
{
	chomp;
	@token=split(/\n/,$_);
	$head=$token[0];
	@token1=split(/\s/,$head);
	$id=$token1[1];
	
	shift @token;
	$seq=join("",@token);
	$tp{$id}=$seq;
}
close(tp);
open(cds,"/lustre/home/zhangyw/data/bowtie2index/hg38_transcript_gencodeV33/CDS.fa");
local $/=">";
while(<cds>)
{
	chomp;
	@token=split(/\n/,$_);
	$head=$token[0];
	@token1=split(/\s/,$head);
	$id=$token1[1];
	
	shift @token;
	$seq=join("",@token);
	$cds{$id}=$seq;
}
close(cds);
open(file,"../finalfile.bed");
open(res,">relative_position.new.tsv");
local $/="\n";
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	@token1=split(/_/,$token[9]);
	$tid=$token1[0];
	$key=$token[3]."\t".$token1[0];
	$seq100=$hash{$key};
	@token2=split(//,$seq100);
	$query=join("",@token2[50..57]);
	$region=$token[12];
	if($query ne "")
	{
	if($region eq "5UTR")
	{
		$seq=$fp{$tid};
		@token3=split(/$query/,$seq);
		@token4=split(//,$token3[0]);
		$length=@token4;
		@token5=split(//,$seq);
		$full=@token5;
		$pro=($length)/$full;
		
		print res $query."\t".$key."\t".$pro."\t5UTR\n";
	}
		if($region eq "3UTR")
	{
		$seq=$tp{$tid};
		@token3=split(/$query/,$seq);
		@token4=split(//,$token3[0]);
		$length=@token4;
		@token5=split(//,$seq);
		$full=@token5;
		$pro=($length)/$full;
		
		print res $query."\t".$key."\t".$pro."\t3UTR\n";
	}
		if($region eq "CDS")
	{
		$seq=$cds{$tid};
		@token3=split(/$query/,$seq);
		@token4=split(//,$token3[0]);
		$length=@token4;
		@token5=split(//,$seq);
		$full=@token5;
		$pro=($length)/$full;
		
		print res $query."\t".$key."\t".$pro."\tCDS\n";
	}
	
	
	}
}
close(file);
close(res);
