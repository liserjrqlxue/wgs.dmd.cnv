package util

import (
	"bufio"
	"fmt"
	"github.com/liserjrqlxue/goUtil/fmtUtil"
	math2 "github.com/liserjrqlxue/goUtil/math"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/liserjrqlxue/goUtil/stringsUtil"
	"github.com/liserjrqlxue/goUtil/textUtil"
	"io"
	"log"
	"math"
	"os"
	"os/exec"
	"regexp"
	"strings"
)

// INPUT

// LoadExon load exon info to exonInfos
func LoadExon(file string) (exonInfos []*Exon) {
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
		xDepth         = qc.DepthX
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

// load wgs nator.step6 cnv result
var isEX = regexp.MustCompile(`DMD_EX`)

func LoadNatorStep6(path string) (cnvInfos []*Info, title []string) {
	var file = osUtil.Open(path)
	var buf = bufio.NewReader(file)
	var line, err = buf.ReadString('\n')
	simpleUtil.CheckErr(err)
	title = strings.Split(strings.TrimSuffix(line, "\n"), "\t")
	for err == nil {
		line, err = buf.ReadString('\n')
		if err != nil {
			break
		}
		var datum = make(map[string]string)
		for j, k := range strings.Split(strings.TrimSuffix(line, "\n"), "\t") {
			datum[title[j]] = k
		}
		if isEX.MatchString(datum["OMIM_EX"]) {
			var info = &Info{
				ID:    datum["Sample"],
				chr:   datum["Chr"],
				start: stringsUtil.Atoi(datum["Start"]),
				end:   stringsUtil.Atoi(datum["End"]),
				ratio: stringsUtil.Atof(datum["normalized_RD"]),
				raw:   datum,
			}
			cnvInfos = append(cnvInfos, info)
		}
	}
	if err != io.EOF {
		log.Fatalf("bufio.Reader.ReadString('\n') error:[%+v]", err)
	}
	simpleUtil.CheckErr(file.Close())
	return
}

func LoadLumpyStep6(path string) (cnvInfos []*Info, title []string) {
	var file = osUtil.Open(path)
	var buf = bufio.NewReader(file)
	var line, err = buf.ReadString('\n')
	simpleUtil.CheckErr(err)
	title = strings.Split(strings.TrimSuffix(line, "\n"), "\t")
	for err == nil {
		line, err = buf.ReadString('\n')
		if err != nil {
			break
		}
		var datum = make(map[string]string)
		for j, k := range strings.Split(strings.TrimSuffix(line, "\n"), "\t") {
			datum[title[j]] = k
		}
		if isEX.MatchString(datum["OMIM_exon"]) {
			var info = &Info{
				ID:    datum["Sample"],
				chr:   datum["CHROM"],
				start: stringsUtil.Atoi(datum["POS"]),
				end:   stringsUtil.Atoi(datum["END"]),
				ratio: stringsUtil.Atof(datum["RATIO"]),
				raw:   datum,
			}
			cnvInfos = append(cnvInfos, info)
		}
	}
	if err != io.EOF {
		log.Fatalf("bufio.Reader.ReadString('\n') error:[%+v]", err)
	}
	simpleUtil.CheckErr(file.Close())
	return
}

// OUTPUT

// WriteInfos write []*Info to file
func WriteInfos(infos []*Info, path string) {
	var file = osUtil.Create(path)
	fmtUtil.FprintStringArray(file, infoTitle, "\t")
	for _, info := range infos {
		fmtUtil.Fprintln(file, info.String())
	}
	simpleUtil.CheckErr(file.Close())
}

// WriteCNV write []*Info to file
func WriteCNV(infos []*Info, title []string, path string) {
	var file = osUtil.Create(path)
	fmtUtil.FprintStringArray(file, title, "\t")
	for _, info := range infos {
		info.raw["depth"] = float2str(info.depth)
		info.raw["ratio"] = float2str(info.ratio)
		info.raw["fixRatio"] = float2str(info.fixRatio)
		info.raw["factor"] = float2str(info.factor)
		info.raw["percent"] = float2str(info.percent)
		var line []string
		for _, s := range title {
			line = append(line, info.raw[s])
		}
		fmtUtil.FprintStringArray(file, line, "\t")
	}
	simpleUtil.CheckErr(file.Close())
}

// write *QC to file
func WriteQC(qc *QC, path string) {
	log.Println("write QC")
	var file = osUtil.Create(path)
	fmtUtil.FprintStringArray(file, qcTitle, "\t")
	fmtUtil.Fprintln(file, qc.String())
	simpleUtil.CheckErr(file.Close())
}

func Bam2depth(bam, depth, region, control string, qc *QC, skip bool) (depthInfos []*Info) {
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
func Depth2bin(depthInfos []*Info, qc *QC) (binInfos []*Info) {
	var (
		binWidth = qc.BinWidth
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

func Bin2cnv(infos []*Info, outputBin, outputCNV, rScript string, skip bool) (cnvInfos []*Info) {
	WriteInfos(infos, outputBin)
	if !skip {
		runDNAcopy(outputBin, outputCNV, rScript)
	}
	return loadCNA(outputCNV)
}

func runDNAcopy(input, output, rScript string) {
	var DNAcopy = exec.Command(
		"Rscript",
		rScript,
		input, output,
	)
	log.Println(DNAcopy.String())
	DNAcopy.Stdout = os.Stdout
	DNAcopy.Stderr = os.Stderr
	simpleUtil.CheckErr(DNAcopy.Run())
}

func MergeCNV(cnvInfos []*Info, qc *QC) (mergeCnvInfos []*Info) {
	var gender = qc.Gender
	var i2 = 0 // ?????????
	for i, info := range cnvInfos {
		if i < i2 {
			continue
		}
		if i == len(cnvInfos)-1 {
			mergeCnvInfos = append(mergeCnvInfos, info)
			continue
		}
		i2 = i + 1
		for isSameLevel(info.ratio, cnvInfos[i2].ratio, gender) {
			info = mergeInfos(info, cnvInfos[i2])
			i2 = i2 + 1
			if i2 >= len(cnvInfos) {
				break
			}
		}
		mergeCnvInfos = append(mergeCnvInfos, info)
	}
	return
}

func AnnotateInfos(infos, rawInfos []*Info, exonInfos []*Exon, qc *QC) {
	for _, info := range infos {
		info.UpdateInfo(rawInfos, exonInfos, qc)
	}
}

func FilterInfos(infos []*Info, coverage, percent float64) (filterInfos []*Info) {
	for _, info := range infos {
		if info.allCoverage >= coverage && info.percent >= percent {
			filterInfos = append(filterInfos, info)
		}
	}
	return
}

func mergeInfos(x, y *Info) *Info {
	var info = &Info{
		ID:      x.ID,
		chr:     x.chr,
		start:   x.start,
		end:     y.end,
		numMark: x.numMark + y.numMark,
		ratio:   weightMean(x.ratio, y.ratio, x.numMark, y.numMark),
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

func float2str(f float64) string {
	return fmt.Sprintf("%.4f", f)
}
