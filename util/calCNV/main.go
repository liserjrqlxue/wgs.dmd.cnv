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
	gender = flag.String(
		"gender",
		"",
		"gender,{M|F}",
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
	exon = flag.String(
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

	var bin = *width

	var exonInfo []*Exon
	for _, str := range textUtil.File2Slice(*exon, "\t") {
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

	var controlData, _ = textUtil.File2MapArray(*control, "\t", nil)
	var (
		controls []*Info
		factors  []float64
	)
	for _, datum := range controlData {
		var info = &Info{
			chr:     datum["chromosome"],
			start:   stringsUtil.Atoi(datum["pos"]),
			end:     stringsUtil.Atoi(datum["pos"]),
			numMark: 1,
			factor:  stringsUtil.Atof(datum["mean"]),
		}
		factors = append(factors, info.factor)
		controls = append(controls, info)
	}

	var (
		binControls []*Info
	)
	var n = len(controls) / bin
	for i := 0; i < n; i++ {
		var info = &Info{
			chr:     binControls[i*bin].chr,
			start:   binControls[i*bin].start,
			end:     binControls[(i+1)*bin-1].end,
			numMark: bin,
			factor:  math2.Mean(factors[i*bin : (i+1)*bin-1]),
		}
		binControls = append(binControls, info)
	}

	var cnaData, _ = textUtil.File2MapArray(*cna, "\t", nil)
	var cnvInfo []*Info
	for _, datum := range cnaData {
		var info = &Info{
			ID:      datum["ID"],
			chr:     datum["chr"],
			start:   stringsUtil.Atoi(datum["loc.start"]),
			end:     stringsUtil.Atoi(datum["loc.end"]),
			numMark: stringsUtil.Atoi(datum["num.mark"]),
			segMean: stringsUtil.Atof(datum["seg.mean"]),
			ratio:   stringsUtil.Atof(datum["ratio"]),
			percent: stringsUtil.Atof(datum["percent"]),
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
	defer simpleUtil.DeferClose(output)
	var filterOutput *os.File
	if *filter {
		filterOutput = osUtil.Create(*prefix + ".merge.filter.txt")

	}
	fmtUtil.FprintStringArray(output, []string{"ID", "chr", "start", "end", "numMark", "segMean", "ratio", "factor", "percent", "allAnno", "allCoverage", "anno", "coverage"}, "\t")
	for _, info := range mergeInfo {
		var (
			exonLength    = 0
			exonCnvLength = 0
			coverages     []string
		)

		// factor
		var cnvFactors []float64
		for _, binControl := range binControls {
			if binControl.start <= info.end && binControl.end >= info.start {
				cnvFactors = append(cnvFactors, binControl.factor)
			}
		}
		info.factor = math2.Mean(cnvFactors)

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
		if info.allCoverage >= *thresholdCoverage && info.percent >= *thresholdPercent {
			fmtUtil.Fprintln(filterOutput, info.String())
		}
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
