package main

import (
	"flag"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"log"
	"os"
	"path/filepath"
)

// os
var (
	ex, _   = os.Executable()
	exPath  = filepath.Dir(ex)
	binPath = filepath.Join(exPath, "..", "..", "bin")
	etcPath = filepath.Join(exPath, "..", "..", "etc")
)

// flag
var (
	id = flag.String(
		"id",
		"",
		"sampleID",
	)
	bam = flag.String(
		"bam",
		"",
		"bam file",
	)
	gender = flag.String(
		"gender",
		"",
		"gender,{M|F}",
	)
	depthX = flag.Float64(
		"depthX",
		0,
		"average depth of chrX",
	)
	prefix = flag.String(
		"prefix",
		"",
		"output prefix",
	)
	// config
	skip = flag.Bool(
		"skip",
		false,
		"")
	filter = flag.Bool(
		"filter",
		false,
		"output filter exon cnv",
	)
	region = flag.String(
		"region",
		"chrX:31000000-33500000",
		"calculate region",
	)
	width = flag.Int(
		"width",
		50,
		"bin width",
	)
	thresholdCoverage = flag.Float64(
		"coverage",
		0.75,
		"coverage threshold",
	)
	thresholdPercent = flag.Float64(
		"percent",
		75.0,
		"percent threshold",
	)
	exons = flag.String(
		"exon",
		filepath.Join(etcPath, "exon.info.txt"),
		"exon info",
	)
	control = flag.String(
		"control",
		filepath.Join(etcPath, "control.txt"),
		"control to correct ratio",
	)
)

func main() {
	flag.Parse()
	if *bam == "" && *prefix == "" {
		flag.Usage()
		log.Fatalln("-bam or -prefix required!")
	}

	if *prefix == "" {
		*prefix = *bam
	}

	var qc = &QC{
		ID:       *id,
		gender:   *gender,
		depthX:   *depthX,
		binWidth: *width,
	}

	// prepare
	var (
		exonInfo     = loadExon(*exons)
		outputDepth  = *prefix + ".DMD.depth.txt"
		outputBin    = *prefix + ".DMD.Bin.txt"
		outputCNV    = *prefix + ".DMD.CNV.txt"
		outputMerge  = *prefix + ".DMD.merged.CNV.txt"
		outputFilter = *prefix + ".DMD.filtered.CNV.txt"
		outputQC     = *prefix + ".DMD.QC.txt"
		skipDepth    = false
		skipDNAcopy  = false
	)

	if *skip {
		if osUtil.FileExists(outputCNV) {
			skipDNAcopy = true
		}
		if osUtil.FileExists(outputDepth) {
			skipDepth = true
		}
	}

	var depthInfo = bam2depth(*bam, outputDepth, *region, *control, qc, skipDepth)
	var binInfo = depth2bin(depthInfo, qc)
	var cnvInfo = bin2cnv(binInfo, outputBin, outputCNV, skipDNAcopy)
	var mergeInfo = mergeCNV(cnvInfo, qc)

	annotateInfos(mergeInfo, binInfo, exonInfo, qc)
	writeInfos(mergeInfo, outputMerge)
	if *filter {
		var filterInfo = filterInfos(mergeInfo, *thresholdCoverage, *thresholdPercent)
		writeInfos(filterInfo, outputFilter)
	}

	writeQC(qc, outputQC)

}
