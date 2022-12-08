package main

import (
	"flag"
	"github.com/liserjrqlxue/goUtil/fmtUtil"
	math2 "github.com/liserjrqlxue/goUtil/math"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/liserjrqlxue/goUtil/stringsUtil"
	"github.com/liserjrqlxue/goUtil/textUtil"
	"log"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

// os
var (
	ex, _   = os.Executable()
	exPath  = filepath.Dir(ex)
	etcPath = filepath.Join(exPath, "..", "..", "etc")
)

// flag
var (
	cna = flag.String(
		"cna",
		"",
		"cnv result",
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
	if *cna == "" {
		flag.Usage()
		log.Fatalln("-cna required!")
	}

	if *prefix == "" {
		*prefix = *cna
	}

	var (
		depthData      = textUtil.File2Slice(*depth, "\t")
		controlData, _ = textUtil.File2MapArray(*control, "\t", nil)
		cnaData, _     = textUtil.File2MapArray(*cna, "\t", nil)

		exonInfo  []*Exon
		depthInfo []*Info
		binInfo   []*Info
		cnvInfo   []*Info

		depths      []float64
		factors     []float64
		depthRatios []float64
		fixRatios   []float64
		bin         = *width
		n           int
	)

	if len(depthData) != len(controlData) {
		log.Fatalf("nrow error:[depth:%d]vs.[control:%d]\n", len(depthData), len(controlData))
	}

	for i, str := range depthData {
		var info = &Info{
			chr:     str[0],
			start:   stringsUtil.Atoi(str[1]),
			end:     stringsUtil.Atoi(str[1]),
			numMark: 1,
			depth:   stringsUtil.Atof(str[2]),
			factor:  stringsUtil.Atof(controlData[i]["mean"]),
		}
		info.depthRatio = info.depth / *depthX
		info.fixRatio = info.depthRatio / info.factor

		depths = append(depths, info.depth)
		factors = append(factors, info.factor)
		depthRatios = append(depthRatios, info.depthRatio)
		fixRatios = append(fixRatios, info.fixRatio)

		depthInfo = append(depthInfo, info)
	}

	n = len(depthInfo) / bin
	for i := 0; i < n; i++ {
		var info = &Info{
			chr:        depthInfo[i*bin].chr,
			start:      depthInfo[i*bin].start,
			end:        depthInfo[(i+1)*bin-1].end,
			numMark:    bin,
			depth:      math2.Mean(depths[i*bin : (i+1)*bin-1]),
			factor:     math2.Mean(factors[i*bin : (i+1)*bin-1]),
			depthRatio: math2.Mean(depthRatios[i*bin : (i+1)*bin-1]),
			fixRatio:   math2.Mean(fixRatios[i*bin : (i+1)*bin-1]),
		}
		binInfo = append(binInfo, info)
	}

	for _, str := range textUtil.File2Slice(*exons, "\t") {
		var exon = &Exon{
			chr:    str[0],
			start:  stringsUtil.Atoi(str[1]) + 1,
			end:    stringsUtil.Atoi(str[2]),
			lenght: stringsUtil.Atoi(str[3]),
			gene:   str[4],
			trans:  str[5],
			strand: str[6],
			exon:   str[7],
		}
		exonInfo = append(exonInfo, exon)
	}

	for _, datum := range cnaData {
		var info = &Info{
			ID:      datum["ID"],
			chr:     datum["chr"],
			start:   stringsUtil.Atoi(datum["loc.start"]),
			end:     stringsUtil.Atoi(datum["loc.end"]),
			numMark: stringsUtil.Atoi(datum["num.mark"]),
			segMean: stringsUtil.Atof(datum["seg.mean"]),
			ratio:   stringsUtil.Atof(datum["ratio"]),
		}
		cnvInfo = append(cnvInfo, info)
	}

	var mergeInfo []*Info
	var i2 = 0 // 下一个
	for i, info := range cnvInfo {
		if i < i2 {
			continue
		}
		if i == len(cnvInfo)-1 {
			mergeInfo = append(mergeInfo, info)
			continue
		}
		i2 = i + 1
		for isSameLevel(info.ratio, cnvInfo[i2].ratio, *gender) {
			info = mergeInfos(info, cnvInfo[i2], *gender)
			i2 = i2 + 1
			if i2 >= len(cnvInfo) {
				break
			}
		}
		mergeInfo = append(mergeInfo, info)
	}

	var output = osUtil.Create(*prefix + ".merge.txt")
	var filterOutput *os.File
	if *filter {
		filterOutput = osUtil.Create(*prefix + ".merge.filter.txt")
	}
	fmtUtil.FprintStringArray(output, infoTitle, "\t")
	if *filter {
		fmtUtil.FprintStringArray(filterOutput, infoTitle, "\t")
	}

	for i, info := range mergeInfo {
		var (
			exonLength    = 0
			exonCnvLength = 0
			coverages     []string

			cnvRatios      []float64
			cnvDepths      []float64
			cnvFactors     []float64
			cnvDepthRatios []float64
			cnvFixRatios   []float64
		)

		info.percent = math.Min(info.percent, 100)

		// depth + factor
		for _, info2 := range binInfo {
			if info2.start <= info.end && info2.end >= info.start {
				cnvDepths = append(cnvDepths, info2.depth)
				cnvFactors = append(cnvFactors, info2.factor)
				cnvRatios = append(cnvRatios, info2.ratio)
				cnvDepthRatios = append(cnvDepthRatios, info2.depthRatio)
				cnvFixRatios = append(cnvFixRatios, info2.fixRatio)
			}
		}
		info.depth = math2.Mean(cnvDepths)
		info.factor = math2.Mean(cnvFactors)
		info.ratio = math2.Mean(cnvRatios)
		info.depthRatio = math2.Mean(cnvDepthRatios)
		info.fixRatio = math2.Mean(cnvFixRatios)

		log.Printf("Compare:\t%d\tdepthRatio:[%f:%f/%f]\tfixRatio:[%f:%f/%f]", i, info.depthRatio, info.depth, *depthX, info.fixRatio, info.depthRatio, info.factor)

		// exon info
		for _, e := range exonInfo {
			if e.start <= info.end && e.end >= info.start {
				info.annos = append(info.annos, e.exon)
				var (
					hitStart = e.start
					hitEnd   = e.end
				)
				if hitStart < info.start {
					hitStart = info.start
				}
				if hitEnd > info.end {
					hitEnd = info.end
				}
				exonLength += e.end - e.start + 1
				exonCnvLength += hitEnd - hitStart + 1
				var coverage = float64(hitEnd-hitStart+1) / float64(e.end-e.start+1)
				info.coverages = append(info.coverages, coverage)
				coverages = append(coverages, strconv.FormatFloat(coverage, 'f', 2, 32))
			}
		}
		var n = len(info.annos)
		if n > 0 {
			if n == 1 {
				info.allAnno = info.annos[0]
			} else {
				info.allAnno = info.annos[n-1] + "_" + info.annos[0]
			}
			info.anno = strings.Join(info.annos, ",")
			info.coverage = strings.Join(coverages, ",")
			info.allCoverage = float64(exonCnvLength) / float64(exonLength)
		}

		fmtUtil.Fprintln(output, info.String())
		if *filter && info.allCoverage >= *thresholdCoverage && info.percent >= *thresholdPercent {
			fmtUtil.Fprintln(filterOutput, info.String())
		}
	}

	simpleUtil.DeferClose(output)
	if *filter {
		simpleUtil.CheckErr(filterOutput.Close())
	}

}

func mergeInfos(x, y *Info, gender string) *Info {
	var info = &Info{
		ID:      x.ID,
		chr:     x.chr,
		start:   x.start,
		end:     y.end,
		numMark: x.numMark + y.numMark,
		ratio:   weightMean(x.ratio, y.ratio, x.numMark, y.numMark),
		percent: weightMean(x.percent, y.percent, x.numMark, y.numMark),
	}
	if gender == "M" {
		info.percent = math.Abs(info.ratio-1) * 100
	} else {
		info.percent = math.Abs(info.ratio-1) * 200
	}
	info.segMean = math.Log2(info.ratio)
	return info
}

func weightMean(x, y float64, wx, wy int) float64 {
	return (float64(wx)*x + float64(wy)*y) / float64(wx+wy)
}

var (
	cnF0 = 0.15
	cnF1 = 0.75
	cnF2 = 1.25

	cnM0 = 0.3
	cnM1 = 1.5
)

func isSameLevel(x, y float64, gender string) bool {
	if gender == "M" {
		if CalLevelM(x) == CalLevelM(y) {
			return true
		}
	} else {
		if CalLevelF(x) == CalLevelF(y) {
			return true
		}
	}
	return false
}

func CalLevelF(ratio float64) int {
	if ratio <= cnF0 {
		return 0
	} else if ratio < cnF1 {
		return 1
	} else if ratio < cnF2 {
		return 2
	} else {
		return 3
	}
}

func CalLevelM(ratio float64) int {
	if ratio <= cnM0 {
		return 0
	} else if ratio < cnM1 {
		return 1
	} else {
		return 3
	}
}
