#!/bin/bash
set -e
set -x

sampleID=$1
workdir=$2
outdir=${3:-$workdir/$sampleID/CNV/cnvnator}

export PATH=$(dirname $(readlink -e $0)):$PATH

mkdir -p $outdir

cnv=$workdir/$sampleID/CNV/cnvnator/$sampleID.nator.step6
bam=$workdir/$sampleID/bam_chr/$sampleID.final.merge.bam
sex=$workdir/$sampleID/QC/sex.txt
depthX=$(cut -f 1 $sex)
gender=$(cut -f 7 $sex)

\time -v filterCNV -cnv $cnv -bam $bam -depthX $depthX -gender $gender -id $sampleID -prefix $outdir/$sampleID -skip
