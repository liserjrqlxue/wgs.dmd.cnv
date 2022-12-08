package main

import "fmt"

type Info struct {
	ID          string
	chr         string
	start       int
	end         int
	numMark     int
	segMean     float64
	ratio       float64
	factor      float64
	percent     float64
	allAnno     string
	allCoverage float64
	anno        string
	coverage    string

	annos     []string
	coverages []float64
}

func (info *Info) String() string {
	return fmt.Sprintf(
		"%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%.4f\t%s\t%s",
		info.ID, info.chr, info.start, info.end, info.numMark,
		info.segMean, info.ratio, info.factor, info.percent,
		info.allAnno, info.allCoverage, info.anno, info.coverage,
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
