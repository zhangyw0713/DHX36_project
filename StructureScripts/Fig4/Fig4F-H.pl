open(ab,"/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/abundance/maxAbundance.list");
while(<ab>)
{
	chomp;
	@token=split(/\,/,$_);
	$hash{$token[1]}=$token[2];
}
close(ab);

open(wt,"/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/statistic_wt_tp.sig.react.csv");
while(<wt>)
{
	chomp;
	@token=split(/\,/,$_);
	$wt{$token[0]}=$token[2]."\t".$token[4];
}
close(wt);

open(ko,"/lustre/zhangyw/myProject/Dhx36/new_NAIdata/000_rtsc2react_sfalpha/statistic_ko_tp.sig.react.csv");
while(<ko>)
{
	chomp;
	@token=split(/\,/,$_);
	$ko{$token[0]}=$token[2]."\t".$token[4];
}
close(ko);

open(fcfc,"/lustre/zhangyw/myProject/Dhx36/new_RNAseq/chrRNAseq/03_dif/new/fcfc_ko_wt.tsv");
open(res,">tp.react2fc2_ko_wt.tsv");
while(<fcfc>)
{
	chomp;
	@token=split(/\t/,$_);
	if(exists($hash{$token[0]}))
	{
		$enst=$hash{$token[0]};
		if(exists($wt{$enst}) and exists($ko{$enst}))
		{
			$wtavg=$wt{$enst};
			$koavg=$ko{$enst};
			# $delta=$wtavg-$koavg;
			# if($koavg!=0)
			# {
			# $sfc=$wtavg/$koavg;
			$fcfc=$token[3];
			print res $enst."\t".$token[0]."\t".$wtavg."\t".$koavg."\t".$fcfc."\n";
			# }
		}
	}
}
close(fcfc);
close(res);



