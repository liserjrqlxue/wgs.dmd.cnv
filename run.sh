#!/bin/bash
set -e

sampleID=$1
bam=$2
outdir=${3:=test/$sampleID}

inputDir=$(dirname $(dirname $bam))
sex=$inputDir/QC/sex.txt
depth=$(cut -f 1 $sex)
gender=$(cut -f 7 $sex)

mkdir -p $outdir

prefix=$outdir/$sampleID

region="chrX:31000000-33500000"


samtools depth -aa -r $region -o $prefix.DMD.depth.txt $bam

conda activate r


\time -v Rscript dmd.cnv.cal.R $sampleID $gender $depth $prefix.DMD.depth.txt $outdir
