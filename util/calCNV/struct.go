package main

import (
	"fmt"
)

type QC struct {
	ID         string
	binWidth   int
	depthX     float64
	ratioCV    float64
	fixCV      float64
	binRatioCV float64
	binFixCV   float64
}

var qcTitle = []string{
	"ID",
	"binWdith",
	"depthX",
	"ratioCV",
	"fixCV",
	"binRatioCV",
	"binFixCV",
}

func (qc *QC) String() string {
	return fmt.Sprintf(
		"%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f",
		qc.ID,
		qc.binWidth,
		qc.depthX,
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
