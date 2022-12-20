#!/bin/bash
set -e
set -x

sampleID=$1
workdir=$2
cnv=${3:-$workdir/$sampleID/CNV/cnvnator/$sampleID.nator.step6}
outdir=${4:-$workdir/$sampleID/CNV/cnvnator}

export PATH=$(dirname $(readlink -e $0)):$PATH

mkdir -p $outdir

bam=$workdir/$sampleID/bam_chr/$sampleID.final.merge.bam
sex=$workdir/$sampleID/QC/sex.txt
depthX=$(cut -f 1 $sex)
gender=$(cut -f 7 $sex)

\time -v filterCNV  -depthX $depthX -gender "$gender" -id $sampleID -prefix $outdir/$sampleID -cnv $cnv -bam $bam -skip
