package main

import (
	"bufio"
	"flag"
	"github.com/liserjrqlxue/goUtil/fmtUtil"
	"github.com/liserjrqlxue/goUtil/math"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/liserjrqlxue/goUtil/stringsUtil"
	"github.com/liserjrqlxue/goUtil/textUtil"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
)

// flag
var (
	input = flag.String(
		"input",
		"",
		"input info",
	)
	depthDir = flag.String(
		"depthDir",
		"depth",
		"default depth path:-depthDir/sampleID.DMD.depth.txt",
	)
	output = flag.String(
		"output",
		"",
		"output path,default -input.control.txt",
	)
	region = flag.String(
		"region",
		"chrX:31000000-33500000",
		"calculate region",
	)
	depthXthreshold = flag.Float64(
		"depthXthreshold",
		20,
		"depth threshold of chrX",
	)
	skip = flag.Bool(
		"skip",
		false,
		"",
	)
)

func main() {
	flag.Parse()
	if *input == "" {
		flag.Usage()
		log.Fatalln("-input required!")
	}
	if *output == "" {
		*output = *input + ".control.txt"
	}

	var (
		info, _ = textUtil.File2MapArray(*input, "\t", nil)
		out     = osUtil.Create(*output)
		IDs     []string
		depthXs []float64
		files   []*os.File
		scans   []*bufio.Scanner
	)
	for _, m := range info {
		var sampleID = m["sampleID"]
		if sampleID == "" {
			log.Fatalln("[sampleID] required for -input")
		}
		var depthX = stringsUtil.Atof(m["depthX"])
		if depthX < *depthXthreshold {
			log.Printf(
				"skip [%s] for chrX depth:[%.4f]<%.4f",
				sampleID, depthX, *depthXthreshold,
			)
			continue
		}
		IDs = append(IDs, sampleID)
		depthXs = append(depthXs, stringsUtil.Atof(m["depthX"]))

		if m["depth"] == "" {
			m["depth"] = filepath.Join(*depthDir, sampleID+".DMD.depth.txt")
		}

		if !osUtil.FileExists(m["depth"]) && !osUtil.FileExists(m["bam"]) {
			log.Fatalf("either [depth] or [bam] should exists for %s\n", sampleID)
		}

		if !*skip || !osUtil.FileExists(m["depth"]) {
			samtoolsDepth(*region, m["depth"], m["bam"])
		}

		var file = osUtil.Open(m["depth"])
		var scan = bufio.NewScanner(file)
		files = append(files, file)
		scans = append(scans, scan)
	}

	fmtUtil.Fprintf(out, "chromosome\tpos\tmean\n")
	for scans[0].Scan() {
		var (
			line  []string
			ratio []float64
		)
		for i, scan := range scans {
			if i > 0 && !scan.Scan() {
				log.Fatalf("can not scan.Scan() for [%d:%s]\n", i, IDs[i])
			}
			line = strings.Split(scan.Text(), "\t")
			ratio = append(ratio, stringsUtil.Atof(line[2])/depthXs[i])
			simpleUtil.CheckErr(scan.Err())
		}
		fmtUtil.Fprintf(out, "%s\t%s\t%.4f\n", line[0], line[1], math.Mean(ratio))
	}

	simpleUtil.CheckErr(out.Close())
	for _, file := range files {
		simpleUtil.CheckErr(file.Close())
	}
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
