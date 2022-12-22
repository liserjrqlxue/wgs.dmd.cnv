package main

import (
	"flag"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"log"
	"os"
	"path/filepath"
	"wgsDmdCnv/util"
)

// os
var (
	ex, _   = os.Executable()
	exPath  = filepath.Dir(ex)
	etcPath = filepath.Join(exPath, "..", "..", "etc")
)

// flag
var (
	id = flag.String(
		"id",
		"",
		"sampleID",
	)
	cnv = flag.String(
		"cnv",
		"",
		"cnv input",
	)
	cnvType = flag.String(
		"cnvType",
		"nator",
		"nator or lumpy",
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
	control = flag.String(
		"control",
		filepath.Join(etcPath, "control.txt"),
		"control to correct ratio",
	)
	percentThreshold = flag.Float64(
		"percent",
		60,
		"percent threshold",
	)
)

func main() {
	flag.Parse()
	if *cnv == "" {
		flag.Usage()
		log.Fatalln("-cnv required!")

	}
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
		outputDepth     = *prefix + ".DMD.depth.txt"
		outputCNV       = *prefix + ".DMD.CNV.txt"
		outputFilterCNV = *prefix + ".DMD.CNV.filtered.txt"
		skipDepth       = false
	)

	if *skip {
		if osUtil.FileExists(outputDepth) {
			skipDepth = true
		}
	}

	log.Println("load cnv")
	var (
		cnvInfo []*util.Info
		title   []string
	)

	switch *cnvType {
	case "nator":
		cnvInfo, title = util.LoadNatorStep6(*cnv)
	case "lumpy":
		cnvInfo, title = util.LoadLumpyStep6(*cnv)
	default:
		log.Fatalf("-cnvType [%s] invalid, only support [nator,lumpy]", *cnvType)
	}

	log.Println("load depth")
	var depthInfo = util.Bam2depth(*bam, outputDepth, *region, *control, qc, skipDepth)
	log.Println("annotate cnv")
	util.AnnotateInfos(cnvInfo, depthInfo, nil, qc)
	log.Println("filter cnv")
	var filterInfo = util.FilterInfos(cnvInfo, -1, *percentThreshold)

	log.Println("write cnv")
	title = append(title, "depth", "ratio", "fixRatio", "factor", "percent")
	util.WriteCNV(cnvInfo, title, outputCNV)
	util.WriteCNV(filterInfo, title, outputFilterCNV)
	log.Println("all done")
}
