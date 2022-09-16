#!/bin/bash
## cT=$1
## method=$2

## True label
resolution=1000000
for cT in H1ESC HFF #GM12878 H1ESC HAP1 HFF IMR90
do
    for method in Bulk #higashi Bulk # BandNorm 3DVI schicluster
    do
	echo $cT
	echo $method
	if [ "$method" == "Bulk" ]; then
	    inPath="/results/$method"

	    ## KR normalization
	    cat $inPath/binPair/${cT}/${cT}_chr*.binPair | awk -v OFS="\t" -v halfBin=500000 '{print $1, $2 + halfBin, $3, $4 + halfBin, int($5)}' >$inPath/binPair/$cT/$method.$cT.allChrom.binPair
	    
	else
	    inPath="/results/Duan2020/"
	    ## KR normalization
	    if [ "$method" == "higashi" ];then
		
		cat $inPath/binPair/$cT/*${method}*_chr*.binPair | awk -v OFS="\t" -v halfBin=500000 'function ceil(x, y){y=int(x); return(x>y?y+1:y)} {print $1, $2 + halfBin, $3, $4 + halfBin, ceil($5/10)}' >$inPath/binPair/$cT/$method.$cT.allChrom.binPair
	    else
		cat $inPath/binPair/$cT/*${method}*_chr*.binPair | awk -v OFS="\t" -v halfBin=500000 'function ceil(x, y){y=int(x); return(x>y?y+1:y)} {print $1, $2 + halfBin, $3, $4 + halfBin, ceil($5)}' >$inPath/binPair/$cT/$method.$cT.allChrom.binPair
	    fi

	fi


	mkdir -p $inPath/KRnorm
	python3 KR_norm_mHiC.py -r $resolution -l "/p/keles/yezheng/volumeA/HiC_essentialData/hg19.chrom.sizes.short" -c "whole" -tr 10 -f "$inPath/binPair/$cT/$method.$cT.allChrom.binPair" -o "$inPath/KRnorm"


	##Fit-Hi-C
	fithicPath="/p/keles/yezheng/volumeA/mHiC/mHiC_scripts/FitHiC/fithic/"
	mkdir -p $inPath/FitHiC
	## generate fragment file
	/p/keles/yezheng/volumeA/Softwares/miniconda2/bin/python2 $fithicPath/utils/Create_FitHiC_Fragments.py --ChrSizeFile /p/keles/yezheng/volumeA/HiC_essentialData/hg19.chrom.sizes.short --OutFile $inPath/FitHiC/$cT.fragments --binsize $resolution

	rm -rf $inPath/FitHiC/$cT/$method/fithic.spline_pass2.significances.txt
	
	/p/keles/yezheng/volumeA/Softwares/miniconda2/bin/python2 $fithicPath/fithic/fithic-runner.py -f $inPath/FitHiC/$cT.fragments -i $inPath/binPair/$cT/$method.$cT.allChrom.binPair -o $inPath/FitHiC/$cT/$method -b 200  -U 200000000 -L $resolution  -t $inPath/KRnorm/$method.$cT.allChrom.binPair.KRnorm.bias 

	gunzip $inPath/FitHiC/$cT/$method/fithic.spline_pass2.significances.txt.gz
    done
done


