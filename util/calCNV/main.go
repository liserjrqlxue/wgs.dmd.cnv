package main

import (
	"flag"
	"log"
	"os"
	"path/filepath"

	"wgsDmdCnv/util"

	"github.com/liserjrqlxue/goUtil/osUtil"
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
		"",
	)
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
	segment = flag.String(
		"segment",
		filepath.Join(binPath, "cbs.R"),
		"segment R script",
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

	var qc = &util.QC{
		ID:       *id,
		Gender:   *gender,
		DepthX:   *depthX,
		BinWidth: *width,
	}

	// prepare
	var (
		exonInfo     = util.LoadExon(*exons)
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

	var depthInfo = util.Bam2depth(*bam, outputDepth, *region, *control, qc, skipDepth)
	var binInfo = util.Depth2bin(depthInfo, qc)
	var cnvInfo = util.Bin2cnv(binInfo, outputBin, outputCNV, *segment, skipDNAcopy)
	var mergeInfo = util.MergeCNV(cnvInfo, qc)

	util.AnnotateInfos(mergeInfo, binInfo, exonInfo, qc)
	util.WriteInfos(mergeInfo, outputMerge)
	if *filter {
		var filterInfo = util.FilterInfos(mergeInfo, *thresholdCoverage, *thresholdPercent)
		util.WriteInfos(filterInfo, outputFilter)
	}

	util.WriteQC(qc, outputQC)
}
