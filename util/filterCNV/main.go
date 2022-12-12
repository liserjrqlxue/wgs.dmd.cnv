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
		outputDepth = *prefix + ".DMD.depth.txt"
		outputCNV   = *prefix + ".DMD.CNV.txt"
		skipDepth   = false
	)

	if *skip {
		if osUtil.FileExists(outputDepth) {
			skipDepth = true
		}
	}

	log.Println("load cnv")
	var cnvInfo, title = util.LoadNatorStep6(*cnv)
	log.Println("load depth")
	var depthInfo = util.Bam2depth(*bam, outputDepth, *region, *control, qc, skipDepth)
	util.AnnotateInfos(cnvInfo, depthInfo, nil, qc)

	title = append(title, "depth", "ratio", "fixRatio", "factor")
	util.WriteCNV(cnvInfo, title, outputCNV)
}
