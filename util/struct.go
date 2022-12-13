package util

import (
	"fmt"
	math2 "github.com/liserjrqlxue/goUtil/math"
	"math"
	"strconv"
	"strings"
)

type QC struct {
	ID         string
	Gender     string
	BinWidth   int
	DepthX     float64
	ratioCV    float64
	fixCV      float64
	binRatioCV float64
	binFixCV   float64
}

var qcTitle = []string{
	"ID",
	"binWdith",
	"DepthX",
	"ratioCV",
	"fixCV",
	"binRatioCV",
	"binFixCV",
}

func (qc *QC) String() string {
	return fmt.Sprintf(
		"%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f",
		qc.ID,
		qc.BinWidth,
		qc.DepthX,
		qc.ratioCV,
		qc.fixCV,
		qc.binRatioCV,
		qc.binFixCV,
	)
}

type Info struct {
	ID          string
	chr         string
	start       int
	end         int
	numMark     int
	segMean     float64
	depth       float64
	ratio       float64
	depthRatio  float64
	fixRatio    float64
	factor      float64
	percent     float64
	allAnno     string
	allCoverage float64
	anno        string
	coverage    string

	annos     []string
	coverages []float64

	raw map[string]string
}

var infoTitle = []string{
	"ID",
	"chr",
	"start",
	"end",
	"numMark",
	"segMean",
	"depth",
	"ratio",
	"depthRatio",
	"fixRatio",
	"factor",
	"percent",
	"allAnno",
	"allCoverage",
	"anno",
	"coverage",
}

func (info *Info) String() string {
	return fmt.Sprintf(
		"%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%.4f\t%s\t%s",
		info.ID,
		info.chr,
		info.start,
		info.end,
		info.numMark,
		info.segMean,
		info.depth,
		info.ratio,
		info.depthRatio,
		info.fixRatio,
		info.factor,
		info.percent,
		info.allAnno,
		info.allCoverage,
		info.anno,
		info.coverage,
	)
}

// CalPercent calculate percent from fixRatio
func (info *Info) CalPercent(gender string) {
	if gender == "M" {
		info.percent = math.Abs(info.fixRatio-1) * 100
	} else {
		info.percent = math.Abs(info.fixRatio-1) * 200
	}
	info.percent = math.Min(info.percent, 100)
}

// CalRatio calculate mean depth, mean factor, depthRatio and fixRatio from raw depth and raw factor
func (info *Info) CalRatio(depthX float64, rawInfos []*Info) {
	var (
		cnvDepths  []float64
		cnvFactors []float64
	)
	// depth + factor
	for _, info2 := range rawInfos {
		if info2.start <= info.end && info2.end >= info.start {
			cnvDepths = append(cnvDepths, info2.depth)
			cnvFactors = append(cnvFactors, info2.factor)
			if info.end < info2.end {
				info.end = info2.end
			}
		}
	}
	info.depth = math2.Mean(cnvDepths)
	info.factor = math2.Mean(cnvFactors)

	info.depthRatio = info.depth / depthX
	info.fixRatio = info.depthRatio / info.factor
}

// AnnotateExons exon info -> anno and coverage
func (info *Info) AnnotateExons(exonInfos []*Exon) {
	var (
		exonLength    = 0
		exonCnvLength = 0
		coverages     []string
	)

	for _, e := range exonInfos {
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
}

type Exon struct {
	chr    string
	start  int
	end    int
	lenght int
	gene   string
	trans  string
	strand string
	exon   string
}
