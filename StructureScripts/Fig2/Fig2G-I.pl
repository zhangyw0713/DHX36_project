open(file,"/lustre/home/zhangyw/data/bowtie2index/mm10_transcript_gencodevm24/3UTR_ENST_composition.csv");
while(<file>)
{
	chomp;
	@token=split(/\,/,$_);
	$tp{$token[0]}=$token[8];
}
close(file);
open(file,"/lustre/home/zhangyw/data/bowtie2index/mm10_transcript_gencodevm24/CDS_ENST_composition.csv");
while(<file>)
{
	chomp;
	@token=split(/\,/,$_);
	$cds{$token[0]}=$token[8];
}
close(file);
open(file,"/lustre/home/zhangyw/data/bowtie2index/mm10_transcript_gencodevm24/5UTR_ENST_composition.csv");
while(<file>)
{
	chomp;
	@token=split(/\,/,$_);
	$fp{$token[0]}=$token[8];
}
close(file);

open(file,"/lustre/zhangyw/for503/xiaona/Dhx36_reanalysis/02_annotation/00_rep_merge/newexon/finalfile.single.pcgtrans.bed");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	$hash{$token[9]}=0;
}
close(file);

opendir(dir,"dstruct_newparameter_output_10_dif");
while(my $file=readdir(dir))
{
	@token1=split(/_/,$file);
	if(exists($hash{$token1[0]}))
	{
	open(file,"dstruct_newparameter_output_10_dif/$file");
	while(<file>)
	{
		chomp;
		@token=split(/\t/,$_);
		if(!grep(/start/,$_))
		{
			$d=-1*$token[6];
			if($token[6]>=0)
			{
				${$token1[0]}{$token[0]."\t".$token[1]}=$d."\t"."Gain";
			}else{
				${$token1[0]}{$token[0]."\t".$token[1]}=$d."\t"."Loss";
			}
		}
	}
	close(file);
	}
}
closedir(dir);





open(res,">metagene.tsv");
opendir(dir,"/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/dstruct_newparameter_output_10/");
while(my $file=readdir(dir))
{
	@token=split(/_/,$file);
	$tid=$token[0];
	# open(file,"")
	if(exists($tp{$tid}) and exists($cds{$tid}) and exists($fp{$tid}) and exists($hash{$tid}))
	{
		$fplen=$fp{$tid};
		$cdslen=$cds{$tid};
		$cdsend=$fplen+$cdslen;
		$tplen=$tp{$tid};
		$tpend=$fplen+$cdslen+$tplen;
		open(file,"/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/dstruct_newparameter_output_10/$file");
		while(<file>)
		{
			chomp;
			if(!grep(/start/,$_))
			{
				@token=split(/\t/,$_);
				$pos=($token[0]+$token[1])/2;
				$new=${$tid}{$token[0]."\t".$token[1]};
				if($pos<=$fplen)
				{
					$ratio=$pos/$fplen;
					print res $tid."\t".$token[0]."\t".$token[1]."\t".$ratio."\t".$fplen."\t".$cdslen."\t".$tplen."\t5UTR\t$new\n";
				}
				
				if($pos>$fplen and $pos<$cdsend)
				{
					$ratio=($pos-$fplen)/$cdslen+1;
					print res $tid."\t".$token[0]."\t".$token[1]."\t".$ratio."\t".$fplen."\t".$cdslen."\t".$tplen."\tCDS\t$new\n";
				}
				if($pos>$cdsend)
				{
					$ratio=($pos-$cdsend)/$tplen+2;
					print res $tid."\t".$token[0]."\t".$token[1]."\t".$ratio."\t".$fplen."\t".$cdslen."\t".$tplen."\t3UTR\t$new\n";
				}
			}
		}
		close(file);
	}
}
close(res);
closedir(dir);