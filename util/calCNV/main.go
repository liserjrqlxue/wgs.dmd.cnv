package main

import (
	"flag"
	"log"
	"os"
	"path/filepath"
	"strconv"
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
	depth = flag.String(
		"depth",
		"",
		"depth file",
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
	filter = flag.Bool(
		"filter",
		false,
		"output filter exon cnv",
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
	if *depth == "" {
		flag.Usage()
		log.Fatalln("-cna required!")
	}

	if *prefix == "" {
		*prefix = *depth + ".bin" + strconv.Itoa(*width)
	}

	var qc = &QC{
		ID:       *id,
		depthX:   *depthX,
		binWidth: *width,
	}

	// prepare
	var exonInfo = loadExon(*exons)

	var depthInfo = loadDepth(*depth, *control, qc)
	var binInfo = depth2bin(depthInfo, qc)
	var cnvInfo = bin2cnv(binInfo, *prefix, qc)
	var mergeInfo = mergeCNV(cnvInfo)

	annotateInfos(mergeInfo, binInfo, exonInfo, qc)
	writeInfos(mergeInfo, *prefix+".merged.CNV.txt")
	if *filter {
		var filterInfo = filterInfos(mergeInfo, *thresholdCoverage, *thresholdPercent)
		writeInfos(filterInfo, *prefix+".filtered.CNV.txt")
	}

	writeQC(qc, *prefix+".QC.txt")

}
