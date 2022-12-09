package main

import (
	"github.com/liserjrqlxue/goUtil/fmtUtil"
	math2 "github.com/liserjrqlxue/goUtil/math"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/liserjrqlxue/goUtil/stringsUtil"
	"github.com/liserjrqlxue/goUtil/textUtil"
	"log"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
)

// INPUT

// load exon info to exonInfos
func loadExon(file string) (exonInfos []*Exon) {
	for _, str := range textUtil.File2Slice(file, "\t") {
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
		exonInfos = append(exonInfos, exon)
	}

	return
}

// load depth file to depthInfos
func loadDepth(depth, control string, qc *QC) (depthInfos []*Info) {
	var (
		depthData      = textUtil.File2Slice(depth, "\t")
		controlData, _ = textUtil.File2MapArray(control, "\t", nil)
		xDepth         = qc.depthX
	)

	for i, str := range depthData {
		var info = &Info{
			chr:     str[0],
			start:   stringsUtil.Atoi(str[1]),
			end:     stringsUtil.Atoi(str[1]),
			numMark: 1,
			depth:   stringsUtil.Atof(str[2]),
			factor:  stringsUtil.Atof(controlData[i]["mean"]),
		}
		info.depthRatio = info.depth / xDepth
		info.fixRatio = info.depthRatio / info.factor

		depthInfos = append(depthInfos, info)
	}
	return
}

// load cna data to cnvInfos
func loadCNA(file string) (cnvInfos []*Info) {
	var cnaData, _ = textUtil.File2MapArray(file, "\t", nil)
	for _, datum := range cnaData {
		var info = &Info{
			chr:     datum["chrom"],
			start:   stringsUtil.Atoi(datum["loc.start"]),
			end:     stringsUtil.Atoi(datum["loc.end"]),
			numMark: stringsUtil.Atoi(datum["num.mark"]),
			segMean: stringsUtil.Atof(datum["seg.mean"]),
		}
		info.ratio = math.Pow(2, info.segMean)
		cnvInfos = append(cnvInfos, info)
	}

	return
}

// OUTPUT

// write []*Info to file
func writeInfos(infos []*Info, path string) {
	var file = osUtil.Create(path)
	fmtUtil.FprintStringArray(file, infoTitle, "\t")
	for _, info := range infos {
		fmtUtil.Fprintln(file, info.String())
	}
	simpleUtil.CheckErr(file.Close())
}

// write *QC to file
func writeQC(qc *QC, path string) {
	log.Println("write QC")
	var file = osUtil.Create(path)
	fmtUtil.FprintStringArray(file, qcTitle, "\t")
	fmtUtil.Fprintln(file, qc.String())
	simpleUtil.CheckErr(file.Close())
}

func bam2depth(bam, depth, region, control string, qc *QC, skip bool) (depthInfos []*Info) {
	if !skip {
		samtoolsDepth(region, depth, bam)
	}
	return loadDepth(depth, control, qc)
}

// samtools depth -aa -r $region -o $output $bam
func samtoolsDepth(region, output, bam string) {
	var cmd = exec.Command(
		"samtools",
		"depth",
		"-aa",
		"-r", region,
		"-o", output,
		bam,
	)
	log.Println(cmd.String())
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	simpleUtil.CheckErr(cmd.Run())
}

// depthInfos to binInfos, and update qc
func depth2bin(depthInfos []*Info, qc *QC) (binInfos []*Info) {
	var (
		binWidth = qc.binWidth
		ID       = qc.ID

		depths       []float64
		factors      []float64
		depthRatios  []float64
		fixRatios    []float64
		binRatios    []float64
		binFixRatios []float64
	)

	for _, info := range depthInfos {
		depths = append(depths, info.depth)
		factors = append(factors, info.factor)
		depthRatios = append(depthRatios, info.depthRatio)
		fixRatios = append(fixRatios, info.fixRatio)
	}

	var n = len(depthInfos) / binWidth
	for i := 0; i < n; i++ {
		var (
			s = i * binWidth
			e = s + binWidth - 1
		)
		var info = &Info{
			ID:         ID,
			chr:        depthInfos[s].chr,
			start:      depthInfos[s].start,
			end:        depthInfos[e].end,
			numMark:    binWidth,
			depth:      math2.Mean(depths[s:e]),
			factor:     math2.Mean(factors[s:e]),
			depthRatio: math2.Mean(depthRatios[s:e]),
			fixRatio:   math2.Mean(fixRatios[s:e]),
		}
		binInfos = append(binInfos, info)

		binRatios = append(binRatios, info.depthRatio)
		binFixRatios = append(binFixRatios, info.fixRatio)
	}

	qc.ratioCV = calCV(depthRatios)
	qc.fixCV = calCV(fixRatios)
	qc.binRatioCV = calCV(binRatios)
	qc.binFixCV = calCV(binFixRatios)

	return
}

func bin2cnv(infos []*Info, outputBin, outputCNV string, skip bool) (cnvInfos []*Info) {
	writeInfos(infos, outputBin)
	if !skip {
		runDNAcopy(outputBin, outputCNV)
	}
	return loadCNA(outputCNV)
}

func runDNAcopy(input, output string) {
	var DNAcopy = exec.Command(
		"Rscript",
		filepath.Join(binPath, "dmd.cnv.cal.R"),
		input, output,
	)
	log.Println(DNAcopy.String())
	DNAcopy.Stdout = os.Stdout
	DNAcopy.Stderr = os.Stderr
	simpleUtil.CheckErr(DNAcopy.Run())
}

func mergeCNV(cnvInfos []*Info) (mergeCnvInfos []*Info) {
	var i2 = 0 // 下一个
	for i, info := range cnvInfos {
		if i < i2 {
			continue
		}
		if i == len(cnvInfos)-1 {
			mergeCnvInfos = append(mergeCnvInfos, info)
			continue
		}
		i2 = i + 1
		for isSameLevel(info.ratio, cnvInfos[i2].ratio, *gender) {
			info = mergeInfos(info, cnvInfos[i2], *gender)
			i2 = i2 + 1
			if i2 >= len(cnvInfos) {
				break
			}
		}
		mergeCnvInfos = append(mergeCnvInfos, info)
	}
	return
}

func annotateInfos(infos, rawInfos []*Info, exonInfos []*Exon, qc *QC) {
	for _, info := range infos {
		annotateInfo(info, rawInfos, exonInfos, qc)
	}
}

func filterInfos(infos []*Info, coverage, percent float64) (filterInfos []*Info) {
	for _, info := range infos {
		if info.allCoverage >= coverage && info.percent >= percent {
			filterInfos = append(filterInfos, info)
		}
	}
	return
}

func annotateInfo(info *Info, rawInfos []*Info, exonInfos []*Exon, qc *QC) {
	info.ID = qc.ID
	info.percent = math.Min(info.percent, 100)

	// depth + factor + ratio
	reCalRatio(info, rawInfos, qc)

	// exon info -> anno and coverage
	annoteateExon(info, exonInfos)
}

func reCalRatio(info *Info, rawInfos []*Info, qc *QC) {
	var (
		xDepth     = qc.depthX
		cnvDepths  []float64
		cnvFactors []float64
	)
	// depth + factor
	for _, info2 := range rawInfos {
		if info2.start <= info.end && info2.end >= info.start {
			cnvDepths = append(cnvDepths, info2.depth)
			cnvFactors = append(cnvFactors, info2.factor)
		}
	}
	info.depth = math2.Mean(cnvDepths)
	info.factor = math2.Mean(cnvFactors)

	info.depthRatio = info.depth / xDepth
	info.fixRatio = info.depthRatio / info.factor
}

func annoteateExon(info *Info, exonInfos []*Exon) {
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

func calCV(x []float64) float64 {
	var mean, sd = math2.MeanStdDev(x)
	return sd / mean
}

// calculate ratio level

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
