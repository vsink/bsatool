package main

import (
	"bufio"
	"math/rand"
	"os/exec"
	"runtime"
	"runtime/pprof"
	"time"

	"fmt"

	"log"

	"encoding/gob"

	"net/http"
	"os"

	"regexp"

	"strings"

	"compress/lzw"

	"sort"
	"strconv"

	"./amino"

	kingpin "gopkg.in/alecthomas/kingpin.v2"

	"html/template"

	"./codon"

	"github.com/360EntSecGroup-Skylar/excelize"
	"github.com/mitchellh/hashstructure"
	"github.com/pkg/browser"
	"github.com/xrash/smetrics"

	"path/filepath"
)

const (
	logo = `
 _______    ________     __  ___________  ______      ______    ___       
|   _  "\  /"       )   /""\("     _   ")/    " \    /    " \  |"  |      
(. |_)  :)(:   \___/   /    \)__/  \\__/// ____  \  // ____  \ ||  |      
|:     \/  \___  \    /' /\  \  \\_ /  /  /    ) :)/  /    ) :)|:  |      
(|  _  \\   __/  \\  //  __'  \ |.  | (: (____/ //(: (____/ //  \  |___   
|: |_)  :) /" \   :)/   /  \\  \\:  |  \        /  \        /  ( \_|:  \  
(_______/ (_______/(___/    \___)\__|   \"_____/    \"_____/    \_______) 
` +
		`BSATool - Bacterial Snp Annotation Tool ver.` + version + "\n" +
		`      Laboratory of Social and Epidemic Infections
 Scientific Centre for Family Health and Human Reproduction Problems
     	(c) V.Sinkov, Irkutsk, Russia, 2018                                   
                                                  
	`
	list   = "list"
	ncFlag = "NC"
	aaFlag = "AA"
	tLMN   = "LMN"   //locus:Mutation:NAME
	tPMLN  = "PMLN"  // position:Mutation:locus:NAME
	tPMN   = "PMN"   // position:Mutation:NAME
	tLSAAN = "LSAAN" //locus:shortAA:codon:shortAA:name
	tLLAAN = "LLAAN" // locus:longAA:codon:longAA:name
	tLCN   = "LCN"   //locus:codon:name
	// tPMNT  = "PMNL"  //position_Ref>Alt{Tab}Name;tag (position, mutation,name, tag )
	vcfExt = "*.vcf"
	// pbtmpl      = `{{counters . }}{{ bar . "⎜" "agtc"   "⬤" "TCAG" "⎜"}}{{percent . }}{{rtime . " ETA %s  "}}{{speed .}} `
	lmnRegExp  = `^(\w+)\W+(\d+)(\D)>(\D)\s+(.*)` // LMN Regular expression
	pmlnRegExp = `^(\d+)_(\D)>(\D)\W+(\w+)\W+(.*)`
	pmnRegExp  = `(\d+)_(\D)>(\D)\W+(.*)$`
	// pmntRegExp  = `^(\d+)_(\D)>(\D)(.*)\b(;|:)\b(.*)`
	lsaanRegExp = `^(\w+)\W+(\D{1})(\d+)(\D{1})\W+(.*)`
	llaanRegExp = `^(\w+)\W+(\D{3})(\d+)(\D{3})\W+(.*)`
	lcnRegExp   = `^(\w+)\W+codon(\d+)\W+(.*)`
	version     = "0.28.18a"
)

var (

	//
	//Database flags
	app       = kingpin.New(logo, "BSATool - Bacterial Snp Annotation Tool")
	appAuthor = app.Author("V.Sinkov")
	appVer    = app.Version(version)

	gbWeb          = kingpin.Flag("web", " Open results in web browser").Short('w').Bool()
	gbXLSX         = kingpin.Flag("xlsx", " Export to XLSX format").Short('x').String()
	gbVerbose      = kingpin.Flag("verbose", "Show additional information ").Short('v').Default("false").Bool()
	gbIndex        = kingpin.Flag("index", "Calculate Complex Index for amino acid changes").Default("false").Bool()
	gbPort         = kingpin.Flag("port", "Use your own localhost:port (default:8080)").Default("8080").String()
	gbNoSeq        = kingpin.Flag("noseq", "Don't show nucleotides").Default("false").Bool()
	gbDebug        = kingpin.Flag("debug", "Debug mode").Default("false").Bool()
	gbLog          = kingpin.Flag("log", "write log file").Default("false").Bool()
	gbExcludeGenes = kingpin.Flag("exclude-genes", "file with genes which should be excluded from mkseq").String()
	gbExcludeSnp   = kingpin.Flag("exclude-snp", "file with genes which should be excluded from mkseq").String()

	mkdb      = kingpin.Command("mkdb", "Create database")
	dbName    = mkdb.Flag("out", "Name of database").Short('o').Required().String()
	dbGenbank = mkdb.Flag("gb", "Name of genbank file").Short('i').Required().String()

	//Annotation flags

	annAction = kingpin.Command("annotate", "annotate vcf").Alias("annotation").Alias("ann").Alias("a")
	annDB     = annAction.Flag("db", "Name of database").Short('b').Required().String()
	annVCF    = annAction.Flag("vcf", "Input VCF file").Short('i').Required().String()
	// annWeb           = annAction.Flag("web", "").Bool()
	annMakeSeq       = annAction.Flag("mkseq", "NC or AA").Short('m').String()
	annMakeSeqRef    = annAction.Flag("ref", "Generate reference sequence").Short('r').Default("false").Bool()
	annWithFilenames = annAction.Flag("wfn", "Show filenames in list annotated VCF's").Short('n').Bool()
	annInDel         = annAction.Flag("indel", "indel detection").Bool()
	annBench         = annAction.Flag("annprof", "cpuprofile").String()
	annSeqLen        = annAction.Flag("len", "indel detection").Int()
	// annBench         = annAction.Flag("cpuprofile", "cpuprofile").String()

	// compute statistic options

	statAction = kingpin.Command("stat", "Calculates statistic tests").Alias("s")
	statDB     = statAction.Flag("db", "Database file").Short('b').Required().String()
	statTask   = statAction.Flag("action", "Type of action:share, snp,dnds, ginfo,bed, matrix,range, circos").Short('a').Required().String()
	// statWeb     = annAction.Flag("web", "").Short('w').Bool()
	statInFile  = statAction.Flag("in", "Input file").Short('i').String()
	statOutFile = statAction.Flag("out", "Output file").Short('o').String()
	statTypeOf  = statAction.Flag("type", "Type of matrix (binary, gc3, dnds, table, nc, locus. freq, jw").Short('t').String()
	statInRule  = statAction.Flag("rule", "Input rule file").Short('r').String()
	statAll     = statAction.Flag("all", "show all dNdS results").Default("false").Bool()
	statBench   = statAction.Flag("statprof", "cpuprofile").String()
	statMakeSeq = statAction.Flag("mkseq", "make seq from snp filelist").Bool()
	// statRefGenomeName   = statAction.Flag("ref", "set custom genome name").String()
	statCircosTypeOf    = statAction.Flag("typeof", "make circos file without IGR regions").String()
	statCircosBandColor = statAction.Flag("color", "set color to band in circos").String()
	statCircosGenome    = statAction.Flag("genome", "set color to band in circos").String()
	statShowAnnotation  = statAction.Flag("annotation", "show annotations to genes").Bool()
	statNbrOfNSP        = statAction.Flag("snp-number", "number of snp for MST matrix").Int()
	statNonRandomize    = statAction.Flag("nonrandomize", "Randomize  mst snp matrix").Default("false").Bool()
	statGroupFromFile   = statAction.Flag("group", "File with filenames and their groups").String()
	// statMkSeq   = statAction.Flag("mkseq", "").Bool()
	// statTH      = statAction.Flag("th", "").Int()

	infoAction = kingpin.Command("info", "Get information")
	infoLocus  = infoAction.Flag("locus", "locus name").String()
	infoDB     = infoAction.Flag("db", "Database file").Short('b').Required().String()
	infoRanges = infoAction.Flag("range", "start:end").String()
	infoCodons = infoAction.Flag("codons", "start:end").String()
	infoShowAs = infoAction.Flag("showas", "Show as:direct (from left to right),gene(direction as in gene) ").String()

	devAction = kingpin.Command("dev", "Developer mode.").Alias("debug")
	devDB     = devAction.Flag("db", "Database file").Short('b').Required().String()
	devTask   = devAction.Flag("action", "Action...").Short('a').Required().String()
	devPwd    = devAction.Flag("pwd", "Password").Short('p').Required().String()

	betaAction  = kingpin.Command("beta", "Beta test mode for testing new functions")
	betaTask    = betaAction.Flag("complex", "Action...").Short('a').Required().String()
	betaDB      = betaAction.Flag("db", "Database").Short('b').Required().String()
	betaInFile  = betaAction.Flag("in", "Input file").Short('i').String()
	betaOutFile = betaAction.Flag("out", "Output file").Short('o').String()
)

var (
	resStart, resEnd int
	genomeSeqSlice   []string // информация об генах, загруженная из базы
	// var genomeCoordinates []genomeCoordInfo
	allGenesVal []geneInfo
	rulesArr    []rulesInfo

	// var snpCacheFromChan []vcfInfoQuery //сюда помещаются результаты запросов через vcfQuery тип
	snpCacheMap = make(map[string][]snpInfo)
	listOfFiles []string

	gInfo genomeInfo

	// var geneData = make(map[string]geneInfo)
	geneCoordinates = make(map[string]gCoords)

	colorReset  = "\x1b[39;49;0m"
	colorRed    = "\x1b[31;1m"
	colorGreen  = "\x1b[32;1m"
	colorYellow = "\x1b[33;1m"
	colorDbBlue = "\x1b[35;1m"
	colorBlue   = "\x1b[36;1m"
	colorWhite  = "\x1b[37;1m"
)

type (
	gCoords struct {
		Start, End int
		Type       string
	}

	geneInfo struct {
		Locus, Name, Product, Direction, GeneID, ProteinID, Note, GOA, TypeOf string
		Start, End, NucCore                                                   int
		PDB, InterPro, ProSite                                                []string
	}

	snpInfo struct {
		/*
					APos абсолютная позиция в геноме
						PosInGene позиция в гене
			PosInCodonG позиция в буквы в кодоне (0-первая, 1-средняя, 2-последняя)

		*/
		APos, PosInGene, PosInCodonG, CodonNbrInG, GeneLen, Start, End, NucCore, TangIdxVal int
		RefCodon, AltCodon, RefAA, AltAA, Locus,
		Direction, NucInPos, Product, Name,
		RefAAShort, AltAAShort, Mutation, Tang, Alt, Note, ReportType, ProteinID, GeneID, GOA, TiTv, TypeOf, ComplexIndex string
		InterPro, PDB, ProSite []string
	}

	snpCheckInfo struct {
		Locus, PosInGene, CodonNbrInG, Ref, Alt, Name, TypeOf, APos, AASref, AASalt, AALref, AALalt, Raw, Tag string
	}

	seqInfo struct {
		Name, Seq     string //
		UsedPositions []int
		// Len       int
	}
	genomeInfo struct {
		Organism, Strain, Version string
		Start, End, NucCore       int
	}

	statInfo struct {
		Pos, Count, Perc        int
		FilesWith, FilesWithout string
	}
	checkSNP struct {
		FileName, FoundSNP, ExactWithRule string
	}

	// DnDsRes is....
	DnDsRes struct {
		N, S, PN, PS, DN, DS, DNDS float64
		ND, NS                     int
		Locus, Product             string
	}

	// JaroWinklerInfo is....
	JaroWinklerInfo struct {
		Locus, Product string
		JWDist         float64
	}

	// GC3Type is....
	GC3Type struct {
		GC3Alt, GC3Ref, GCalt, GCref float64
		Locus                        string
	}

	rangePosInfo struct {
		Start, End, Len, Doubles int
		Seq, Prod, Gname, Note   string
		GC                       float64
	}

	// g.Start, g.End, g.Locus, g.Product, g.Direction, g.TypeOf

	genomeMapInfo struct {
		Start, End             int
		Locus, Product, TypeOf string
	}

	rangeArray struct {
		Start, End int
	}

	rulesInfo struct {
		Name     string
		Variants []string
		Lenght   int
	}

	// ------------------------------------------------------------------------

	vcfInfoQuery struct {
		File    string
		SnpInfo []snpInfo
	}

	vcfQuery struct {
		OutChan chan vcfInfoQuery
		File    string
		Print   bool
	}

	// ------------------------------------------------------------------------

	allPositionsInGene struct {
		pos             int
		alt, ref, locus string
	}

	// ------------------------------------------------------------------------

	snpInfoQuery struct {
		OutChan chan snpInfo
		apos    int
		g       geneInfo
		alt     string
		index   bool
	}
)

// var featAllGenesVal []featInfo

// var allPos []allPositionsInGene

func (q *vcfQuery) request() {
	q.OutChan <- vcfInfoQuery{File: q.File, SnpInfo: parserVCF(q.File, false, allGenesVal)}

}
func (q *vcfQuery) getOutChan() <-chan vcfInfoQuery {
	return q.OutChan
}

func (q *snpInfoQuery) request() {
	// q.OutChan = make(chan snpInfo)
	q.OutChan <- getSnpInfo(q.apos, q.g, q.alt, q.index)

}

func (q *snpInfoQuery) getOutChan() <-chan snpInfo {
	return q.OutChan
}

// ------------------------------------------------------------------------

func init() {
	listOfFiles = getListofVCF()
}

func main() {

	// парсинг флагов
	defer os.Exit(0)
	// flag.Parse()
	kingpin.New("BSATool", "BSATool - Bacterial Snp Annotation Tool")
	kingpin.Version(version).Author("V.Sinkov")

	switch kingpin.Parse() {

	// bsatool mkdb/create/makedb -i /--genbank FILE --out /-o FILE

	case "mkdb":

		if _, err := os.Stat(*dbGenbank); os.IsNotExist(err) {
			fmt.Printf("The %v file is not exist!\n", *dbGenbank)
			os.Exit(3)
		}
		allGenesVal, genomeSeqSlice = geneBankFileParser(*dbGenbank, *gbVerbose)
		writeDB(*dbName, allGenesVal)

	case "annotate":

		allGenesVal = readDB(*annDB)
		var (
			// excludedGenes []string

			exGenes = make(map[int]int)
			exSNPs  = make(map[int]int)
		)
		// go run bsatool.go annotate --vcf list  --db test_core -w

		// cpuprofile := *annBench
		// if cpuprofile != "" {
		// 	f, err := os.Create(cpuprofile)
		// 	if err != nil {
		// 		log.Fatal(err)

		// 	}
		// 	pprof.StartCPUProfile(f)
		// 	defer pprof.StopCPUProfile()
		// }

		if *gbExcludeGenes != "" {
			if _, err := os.Stat(*gbExcludeGenes); os.IsNotExist(err) {
				fmt.Println("No file with excluded genes was found")
				os.Exit(3)

			} else {
				// exFileGenes = true
				exGenes = loadExcludeGenes(*gbExcludeGenes)
			}

			if *gbExcludeSnp != "" {
				if _, err := os.Stat(*gbExcludeSnp); os.IsNotExist(err) {
					fmt.Println("No file with excluded SNPs was found")
					os.Exit(3)

				} else {
					// exSNPpos = true
					// fmt.Println(exSNPpos)
					exSNPs = loadExcludeSNP(*gbExcludeSnp)

				}

			}

		}
		// fmt.Println(exGenes)

		if *annVCF == "list" || *annVCF == "*" || *annVCF == "all" {
			if *gbWeb == false && *annMakeSeq == "" {
				parserBulkVCF(*annWithFilenames)
				// go run bsatool.go annotate --vcf list --mkseq=NC  --db test_core
			} else if *gbWeb == true && *annMakeSeq == "NC" {
				createNCWebServer(*gbPort, exGenes, exSNPs)
			} else if *gbWeb == false && *annMakeSeq == "NC" {
				// seq := makeSeq(*annMakeSeq, *gbVerbose, *annMakeSeqRef)
				seq := makeSeq(*annMakeSeq, *gbVerbose, *annMakeSeqRef, exGenes, exSNPs)
				// for _, val := range seq {
				// 	fmt.Println(val.Seq)
				// }
				for i := 0; i < len(seq); i++ {
					fmt.Println(seq[i].Seq)
				}
				if *gbDebug == true {
					fmt.Println(seq[0].UsedPositions)
				}

			}
		} else {
			// go run bsatool.go annotate --vcf 161_RuU_m.vcf  --db test_core -w
			qSNP := &vcfQuery{File: *annVCF, OutChan: make(chan vcfInfoQuery)}
			go qSNP.request()
			// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
			snpRes := <-qSNP.OutChan

			if _, err := os.Stat(*annVCF); os.IsNotExist(err) {
				fmt.Printf("The %v file is not exist!\n", *annVCF)
				os.Exit(3)
			}
			if *gbWeb == false && *annMakeSeq == "" && *gbXLSX == "" {

				for i := 0; i < len(snpRes.SnpInfo); i++ {
					printTextResults(snpRes.SnpInfo[i], *gbVerbose)
				}

				// parserVCF(*annVCF, true, allGenesVal)

			} else if *gbWeb == true && *annMakeSeq == "" {
				// snps := parserVCF(*annVCF, false, allGenesVal)
				printWebResults(snpRes.SnpInfo, *gbPort)
			}
			if *gbXLSX != "" {
				// snps := parserVCF(*annVCF, false, allGenesVal)
				exportToExcel(snpRes.SnpInfo, *gbXLSX)

			}
		}
	case "stat":
		cpuprofile := *statBench
		if cpuprofile != "" {
			f, err := os.Create(cpuprofile)
			if err != nil {
				log.Fatal(err)

			}
			pprof.StartCPUProfile(f)
			defer pprof.StopCPUProfile()
		}
		allGenesVal = readDB(*statDB)

		switch *statTask {
		case "share":
			// go run bsatool.go   stat -b test_core  -a share  -v

			getShareSNP(*gbVerbose, *gbWeb, listOfFiles)

		case "snp":
			//  go run bsatool.go   stat -b test_core  -a snp  -v
			snpStat()
		case "dnds":

			//  go run bsatool.go   stat -b test_core  -a dnds  -i 161_RuU_m.vcf -v
			if _, err := os.Stat(*statInFile); os.IsNotExist(err) {
				fmt.Println("No input file found")
				os.Exit(3)
			}
			// qSNP := &vcfQuery{File: *statInFile, OutChan: make(chan vcfInfoQuery)}
			// go qSNP.request()
			// // snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
			// snpRes := <-qSNP.OutChan
			// snpCacheMap[snpRes.File] = snpRes.SnpInfo
			// calcDnDsVal(*statInFile, *gbVerbose)
			calcDnDs(*statInFile)
		case "ginfo":
			// go run bsatool.go   stat -b test_core  -a ginfo
			testGeneInfo(allGenesVal)
		case "circos":
			//  go run bsatool.go   stat -b test_core  -a circos --typeof=[cds,igr] --color=lblue --genome NC_000962
			toCircos(allGenesVal)
		case "bed":
			// go run bsatool.go   stat -b test_core  -a bed

			toBED()
		case "matrix":

			// go run bsatool.go stat -a matrix --db test_core -t binary (binary, gc3, dnds, table, nc, locus. freq, jw) -o test.csv

			if *statTypeOf != "" && *statOutFile != "" {
				makeMatrix(*statTypeOf, *statOutFile, *gbVerbose)
			} else if *statTypeOf != "" && *statOutFile == "" {
				fmt.Println("Please, use key -o (--out) to save data.")
			}
		case "range":
			// go run bsatool.go stat -a range --db test_core -i BWAIntersections.txt

			// for _, val := range *statOptions {
			// var opt []string
			// opt =strings.Split(*statOptions, ",")

			//go run bsatool.go stat -a range --db test_core -i BWAIntersections.txt

			if *statInFile != "" {
				file := *statInFile
				res := getRangeFromFile(file, *gbVerbose, *gbNoSeq)
				printSequenceRange(res, *gbWeb, *gbPort)
			}
		case "check":
			// go run bsatool.go stat -a check --db test_core -i drugs2.txt -r rule.txt -w

			if *statMakeSeq == true {
				makeSeqFromSNPListFile(*statInFile)

			} else {
				var useRule bool
				if *statInRule != "" {
					useRule = true
					rulesArr = checkRuleFromFile(*statInRule)
				}
				if *statInFile != "" {

					checkSNPfromFile(*statInFile, *gbVerbose, *gbWeb, useRule)
				}
			}

		case "rule":
			// go run bsatool.go stat -a check --db test_core -i drugs2.txt -r rule.txt -w

			if *statInRule != "" {
				rulesArr = checkRuleFromFile(*statInRule)
				// fmt.Println(rulesArr)
			}
		}

	case "dev":
		allGenesVal = readDB(*devDB)

		switch *devTask {
		case "dump":
			if *devPwd == "vsink" {
				fmt.Println(allGenesVal)
				fmt.Println(geneCoordinates)
			}
		}
		fmt.Println("DevMode")
	case "beta":

		switch *betaTask {

		case "complex":

			if *betaInFile != "" {
				calcComplexIdxFromFile(*betaInFile)
			}
		}
	case "info":
		allGenesVal = readDB(*infoDB)
		var (
			seq, direction string
		)
		if *infoLocus != "" {

			start, end := getGenePosByName(*infoLocus)

			if start != 0 && end != 0 {
				switch *infoShowAs {
				case "", "direct":
					seq = getNucFromGenome(start-1, end)
					direction = "in forward direction"
				case "gene":
					seq = getGeneSequence(*infoLocus)
					direction = "in gene direction"
				}

				prod := getProductByName(*infoLocus)
				// вывод последовательности кодонов заданной ключом --codons=Start:End
				if *infoCodons != "" {
					// seqCodonsMap := map[int]string{}

					seqByCodons := make([]string, (len(seq)/3)+1)
					codonNbr := 1
					var codonBuffer strings.Builder
					for _, c := range seq {
						codonBuffer.WriteString(string(c))
						if codonBuffer.Len() == 3 {
							// seqCodonsMap[codonNbr] = codonBuffer.String()
							seqByCodons[codonNbr] = codonBuffer.String()
							// fmt.Println(codonNbr, codonBuffer.String(), seqByCodons[codonNbr])
							codonBuffer.Reset()
							codonNbr = codonNbr + 1
						}
					}
					codonsCoords := strings.Split(*infoCodons, ":")
					startCodon, _ := strconv.Atoi(codonsCoords[0])
					endCodon, _ := strconv.Atoi(codonsCoords[1])
					for i := startCodon; i <= endCodon; i++ {
						if *gbVerbose == false {
							codonBuffer.WriteString(seqByCodons[i])
						} else {
							// codonBuffer.WriteString(seqByCodons[i])
							codonBuffer.WriteString(fmt.Sprintf("[%v]%v", i, seqByCodons[i]))
						}
						// fmt.Println(i, seqCodonsMap[i])

					}

					fmt.Printf(">%v %v [%v:%v codons %v]\n%v\n", *infoLocus, prod, startCodon, endCodon, direction, codonBuffer.String())

				} else {
					gc, gc1, gc2, gc3 := codon.GcCodonCalc(seq)

					fmt.Printf(">%v %v (%v-%v %v)\n%v\ngc:%v gc1:%v gc2:%v gc3:%v Len:%v\n", *infoLocus, prod, start, end, direction, seq, gc, gc1, gc2, gc3, len(seq))
					// snp := checkSNPnotation("Rv1326c	1038G>A	lineage1.1.3 Indo-Oceanic EAI6 RD239")
					// seqArr := strings.Split(seq, "")
					// posInGene, _ := strconv.Atoi(snp.PosInGene)
					// for i := 0; i < len(seqArr); i++ {
					// 	if *gbVerbose == true {

					// 		seqArr[i] = fmt.Sprintf("%v->%v,", i+1, seqArr[i])
					// 	}
					// 	if i+1 == posInGene {
					// 		seqArr[i] = fmt.Sprintf("[%v_%v>%v]", posInGene, strings.ToUpper(seqArr[i]), snp.Alt)
					// 	}
					// }
					// // seqArr[posInGene+1] = fmt.Sprintf("(%v_%v>%v)", posInGene, strings.ToUpper(seqArr[posInGene+1]), strings.ToUpper(snp.Alt))
					// altSeq := strings.Join(seqArr, "")
					// fmt.Println(strings.TrimRight(altSeq, ","))

					// fmt.Println(seqArr[posInGene+1], snp)
					// fmt.Println(seq, gc, gc1, gc2, gc3, prod)
				}
			} else {
				fmt.Println(colorRed, *infoLocus, " not found!")
			}
		} else if *infoRanges != "" {
			coords := strings.Split(*infoRanges, ":")

			start, _ := strconv.Atoi(coords[0])
			end, _ := strconv.Atoi(coords[1])
			locname, _ := getGeneNameByPos(start, end)
			if start != 0 && end != 0 {
				seq := getNucFromGenome(start-1, end)
				// gc, gc1, gc2, gc3 := codon.GcCodonCalc(seq)
				// prod := getProductByName(*infoLocus)
				// fmt.Printf("%v (%v-%v)\n%v\ngc:%v gc1:%v gc2:%v gc3:%v Len:%v\n", prod, start, end, seq, gc, gc1, gc2, gc3, len(seq))
				// fmt.Println(seq, gc, gc1, gc2, gc3, prod)
				fmt.Printf(">%v|%v|%v (%v:%v)\n%v\n", gInfo.Strain, gInfo.Version, locname, start, end, seq)

			}

		}

	}

}

func getSnpInfo(apos int, g geneInfo, alt string, flgTang bool) snpInfo {

	var (
		snp                         snpInfo  // структура SNP
		codonPositions              []string // срез для разбивки кодона побуквенно
		altCodonPositions           []string // срез для разбивки кодона побуквенно альтернативным нуклеотидом
		locReportType, typeOf, titv string
		geneLen                     int
		codon                       string // переменная для кодона
		altCodon                    string // аналогично для альтернативного кодона
		mut                         string
		tangIdx                     string
		tangIdxVal                  int
	)
	// var trouble int
	lStart := g.Start // переменная начала гена
	lEnd := g.End
	posInGene := ((apos - lStart) + 1)           // позиция снипа в гене
	codonNbrInG := ((posInGene - 1) / 3) + 1     // номер кодона=номеру аминокислоты в трансляции
	posInCodonG := (codonNbrInG * 3) - posInGene // позиция в буквы в кодоне (0-первая, 1-средняя, 2-последняя)
	CPosInGene := ((lEnd - apos) + 1)            // комплементарный ген. позиция в гене
	CCodonNbrInG := ((CPosInGene - 1) / 3) + 1   // комплементарный ген. номер кодона = номеру аминокислоты
	// финт, который делал в snpMiner2, сдвиг на 1 букву. взял оттуда

	if posInCodonG == 2 {
		posInCodonG = 0
	} else if posInCodonG == 0 {
		posInCodonG = 2
	}

	/*
		определяем границы кодона, в зависимости от положения нуклеотида
		0- три буквы справа
		1-одна справа, одна слева
		2-три буквы слева
	*/
	if lStart > lEnd {
		geneLen = (lStart - lEnd) + 1

	} else if lEnd > lStart {
		geneLen = (lEnd - lStart) + 1

	}
	if posInCodonG == 0 {
		codon = getNucFromGenome((posInGene+lStart)-2, ((posInGene+lStart)-1)+2)
	} else if posInCodonG == 1 {
		codon = getNucFromGenome((posInGene+lStart)-3, ((posInGene+lStart)-1)+1)
	} else if posInCodonG == 2 {
		codon = getNucFromGenome((posInGene+lStart)-4, ((posInGene + lStart) - 1))
	}

	/*
		ревертируем кодон для комплементарных генов

	*/
	nucG := getNucFromGenomePos((posInGene + lStart) - 1)
	typeOf = g.TypeOf

	if g.Direction == "r" {
		alt = getComplement(alt)
		nucG = getComplement(nucG)
		codon = getComplement(codon)
		posInGene = CPosInGene
		codonNbrInG = CCodonNbrInG

		if posInCodonG == 2 {
			posInCodonG = 0
		} else if posInCodonG == 0 {
			posInCodonG = 2
		}

	}
	//

	codonPositions = strings.Split(codon, "")
	codonPositions[posInCodonG] = strings.ToUpper(codonPositions[posInCodonG])
	codon = strings.Join(codonPositions, "")
	altCodonPositions = codonPositions
	altCodonPositions[posInCodonG] = alt
	altCodonPositions[posInCodonG] = strings.ToUpper(altCodonPositions[posInCodonG])
	altCodon = strings.Join(altCodonPositions, "")
	aaRef, aaRefShort := amino.Codon2AA(codon)
	aaAlt, aaAltShort := amino.Codon2AA(altCodon)
	if aaRefShort == aaAltShort {
		mut = "synonymous"
		tangIdx = "0000"
	} else if aaRefShort != aaAltShort && aaAltShort != "X" {
		mut = "missense"

		// tangIdx = strconv.FormatFloat(amino.GetTangInx(aaRefShort, aaAltShort), 'f', 2, 64)
		tangIdx, tangIdxVal = amino.GetComplexIndex(aaRefShort, aaAltShort, *gbVerbose)
	} else if aaRefShort != aaAltShort && aaAltShort == "X" {
		mut = "nonsense"
		tangIdx = "0000"
	}

	// if strings.ToUpper(alt) == strings.ToUpper(nucG) {
	// 	trouble = 1
	// }

	if flgTang == true {
		locReportType = "T1"
	} else if flgTang == false {
		locReportType = "T0"
	}
	// amino.GetTangInx(aaRefShort, aaAltShort)
	titv = checkTiTv(nucG, alt)
	// fmt.Println(lStart, lEnd, g.Direction, geneLen, g.Locus, g.Product)

	snp = snpInfo{APos: apos, PosInGene: posInGene, PosInCodonG: posInCodonG,
		RefCodon: codon, RefAA: aaRef, NucInPos: strings.ToUpper(nucG), Locus: g.Locus,
		Direction: g.Direction, Name: g.Name, Product: g.Product,
		Start: g.Start, End: g.End, CodonNbrInG: codonNbrInG, AltCodon: altCodon,
		AltAA: aaAlt, RefAAShort: aaRefShort, AltAAShort: aaAltShort,
		Mutation: mut, Tang: tangIdx, TangIdxVal: tangIdxVal, Alt: alt,
		Note: g.Note, ReportType: locReportType, ProteinID: g.ProteinID,
		GeneID: g.GeneID, GOA: g.GOA, GeneLen: geneLen, TiTv: titv, TypeOf: typeOf}

	return snp
}

func getNucFromGenome(start int, end int) string {
	/*
		возвращает последовательность нуклеотидов из генома:

	*/
	var (
		result string
		slice  []string
	)
	slice = genomeSeqSlice[start:end]
	result = strings.Join(slice, "")
	return result

}

// func getGeneSeq(locus string) {
// 	var (
// 		seq string
// 	)
// 	start, end := getGenePosByName(locus)
// 	dir := getDirectionByName(locus)
// 	if dir == "f" {
// 		seq = getNucFromGenome(start-1, end)
// 	} else if dir == "r" {
// 		seq = getNucFromGenome(start-1, end)

// 	}

// }

func getNucFromGenomePos(pos int) string {
	/*
		возвращает последовательность нуклеотидов из генома:

	*/
	var (
		result string
		slice  []string
	)
	slice = genomeSeqSlice[pos-1 : pos]
	result = strings.Join(slice, "")
	return result

}

func geneBankFileParser(file string, verbose bool) (g []geneInfo, genomeSplice []string) {
	//функция для считывания файла в формате генбанка и занесения строк в массив linesFromGB

	var (
		qenomeSeq    = regexp.MustCompile(`\s+(\d+)\s+(\w{6}.*)`)
		regDelSpaces = regexp.MustCompile(`(\s+)(.*)`)
		checkOrigin  = regexp.MustCompile(`(\d+\s\w{10}\s\w{10}\s\w{10})`)
		// var checkCDS = regexp.MustCompile(`(CDS\s+complement\W\d+\W+\d+\W)|(CDS\s+\d+\W+\d+)`)
		// var endOfCDS = regexp.MustCompile(`(gene\s+\d+\W+\d+)|(gene\s+complement\W\d+\W+\d+)`)
		startOrigin    = regexp.MustCompile(`(ORIGIN\s+)|(translation=)`)
		feature        = regexp.MustCompile(`^\s{5}(\w+)\s+complement\W(\d)+\W+(\d)+\W|^\s{5}(\w+)\s+(\d)+\W+(\d+)`)
		makeAnchors    = regexp.MustCompile(`^(/)`)
		genomeSource   = regexp.MustCompile(`organism=\W(.*)\W`)
		genomeStrain   = regexp.MustCompile(`strain=\W(.*)\W`)
		genomeVersion  = regexp.MustCompile(`^VERSION\s+(.*)`)
		genomeStartEnd = regexp.MustCompile(`source\W+(\d+)..(\d+)`)

		gStart, gEnd int
		gVersion     string

		//  source          1..4411532
		//                  /organism="Mycobacterium tuberculosis H37Rv"
		//                  /mol_type="genomic DNA"
		//                  /strain="H37Rv"
		//                  /db_xref="taxon:83332"
		// var replaceCDS = regexp.MustCompile(`(CDS\s+.*)`)

		resString []string

		extractedData, extractedDataFeat []string
		splitedGenome                    []string
		originBlock, CDSblock, firstCDS  int

		changedStr, organismName, organismStrain string
		noteBuffer                               strings.Builder
	)
	// green := color.New(color.FgGreen)
	fmt.Println(logo)
	f, err := os.Open(file) // открываем файл

	if err != nil {
		fmt.Println(err)

	}
	defer f.Close()
	extractedData = append(extractedData, "START_BLOCK") // начало блока записи
	extractedDataFeat = append(extractedDataFeat, "START_BLOCK")
	scanner := bufio.NewScanner(f) //  новый сканер

	for scanner.Scan() {

		scanTxt := scanner.Text()
		// размер генома
		for _, gStartEnd := range genomeStartEnd.FindAllStringSubmatch(scanTxt, -1) {
			gStart, _ = strconv.Atoi(gStartEnd[1])
			gEnd, _ = strconv.Atoi(gStartEnd[2])
		}
		//  организм
		for _, gName := range genomeSource.FindAllStringSubmatch(scanTxt, -1) {
			organismName = gName[1]
			fmt.Printf("Organism:%v\n", organismName)
		}

		// штамм
		for _, gStrain := range genomeStrain.FindAllStringSubmatch(scanTxt, -1) {
			organismStrain = gStrain[1]
			fmt.Printf("Strain:%v\n%v%v...\r", organismStrain, "Working on ", file)
		}
		for _, genomevesion := range genomeVersion.FindAllStringSubmatch(scanTxt, -1) {
			gVersion = genomevesion[1]
			fmt.Printf("Version:%v\n", gVersion)
		}

		for _, sOrigin := range startOrigin.FindAllStringSubmatch(scanTxt, -1) {
			if len(sOrigin[1]) != 0 || len(sOrigin[2]) != 0 {
				CDSblock = 0
			}
		}

		for _, rfeat := range feature.FindAllStringSubmatch(scanTxt, -1) {
			if rfeat[1] != "gene" && rfeat[4] != "gene" && strings.Contains(rfeat[0], "source") == false {
				CDSblock = 1

				if rfeat[1] != "" {
					scanTxt = strings.Replace(string(scanTxt), rfeat[1], fmt.Sprintf("*%v", rfeat[1]), -1)
				} else if rfeat[4] != "" {
					scanTxt = strings.Replace(string(scanTxt), rfeat[4], fmt.Sprintf("*%v", rfeat[4]), -1)

				}

			} else if rfeat[1] == "gene" || rfeat[4] == "gene" && strings.Contains(rfeat[0], "source") == false {

				CDSblock = 0

			}
		}

		for _, origninFound := range checkOrigin.FindAllStringSubmatch(scanTxt, -1) {
			if len(origninFound[1]) != 0 {
				originBlock = 1
				CDSblock = 0

			} else {
				originBlock = 0
			}
		}

		if CDSblock == 1 {
			if firstCDS == 0 {

				firstCDS = 1

			}
			changedStr = regDelSpaces.ReplaceAllString(string(scanTxt), "$2 ") // удаляем пробелы

			changedStr = strings.Replace(string(changedStr), "/note=", "!note=", -1) // меняем / на ! для дальнейшего парсинга

			changedStr = strings.Replace(string(changedStr), "/product=", "!product=", -1) // см выше.

			changedStr = strings.Replace(string(changedStr), "\"", "", -1)

			changedStr = makeAnchors.ReplaceAllString(string(changedStr), "!!")
			if strings.Index(changedStr, "!!") != 0 && strings.Index(changedStr, "*") != 0 {

				noteBuffer.WriteString(changedStr)

			} else {
				if noteBuffer.Len() != 0 {

					extractedData = append(extractedData, strings.TrimSpace(strings.TrimPrefix(strings.Replace(noteBuffer.String(), "!", "\n!", -1), "!")))
					noteBuffer.Reset()

				}

			}

			if strings.Index(changedStr, "!!") == 0 {

				extractedData = append(extractedData, strings.TrimSpace(strings.TrimPrefix(changedStr, "!!")))

			} else if strings.Index(changedStr, "*") == 0 {

				if firstCDS == 1 {

					extractedData = append(extractedData, "START_OF")
					firstCDS = 2
				} else if firstCDS == 2 {
					extractedData = append(extractedData, "END_OF")
					extractedData = append(extractedData, "START_OF")
				}
				extractedData = append(extractedData, strings.TrimPrefix(changedStr, "*"))

			}

		}

		if originBlock == 1 {

			for _, genomeMatch := range qenomeSeq.FindAllStringSubmatch(scanTxt, -1) {

				resString = append(resString, strings.Replace(genomeMatch[2], " ", "", -1))
			}
		}

	}
	extractedData = append(extractedData, "END_OF")
	extractedData = append(extractedData, "END_BLOCK")
	splitedGenome = strings.SplitAfter(strings.Join(resString, ""), "")

	var (
		lDir, lLoc, lName, lProd, lNote, gID, pID, lGOA, lType string
		lPDBArr, lInterProArr, lProSiteArr                     []string
		nucCore, igrCount                                      int
		cdsStEnd                                               = regexp.MustCompile(`^(\w+)\s+(\d+)\W+(\d+)|^(\w+)\s+complement\W(\d+)\W+(\d+)`)
		cdsLocus                                               = regexp.MustCompile(`locus_tag=(.*)`)
		cdsName                                                = regexp.MustCompile(`gene=(.*)`)
		cdsProd                                                = regexp.MustCompile(`product=(.*)`)
		cdsNote                                                = regexp.MustCompile(`note=(.*)`)
		cdsgID                                                 = regexp.MustCompile(`db_xref=GeneID:(\d+)`)
		cdsprotID                                              = regexp.MustCompile(`protein_id=(.*)`)
		cdsGOA                                                 = regexp.MustCompile(`db_xref=GOA:(.*)`)
		cdsPDB                                                 = regexp.MustCompile(`db_xref=PDB:(.*)`)
		cdsInterPro                                            = regexp.MustCompile(`db_xref=InterPro:(.*)`)
		cdsProSite                                             = regexp.MustCompile(`inference=protein motif:PROSITE:(.*)`)

		startOfBlock, endOfBlock, cdsStart, cdsEnd, lStart, lEnd int
		// var cdsOpen, cdsClosed int
		cdsCount       = map[string]int{}
		igensS, igensE []int

		leftGene []string

	// fmt.Println(extractedData)
	)
	for _, val := range extractedData {
		// fmt.Println(val)
		if strings.Index(val, "START_BLOCK") != -1 {
			startOfBlock = 1
			endOfBlock = 0
		} else if strings.Index(val, "END_BLOCK") != -1 {
			endOfBlock = 1
			startOfBlock = 0
		}
		if strings.Index(val, "START_OF") != -1 {
			cdsStart = 1
			cdsEnd = 0
			// cdsCount++
		} else if strings.Index(val, "END_OF") != -1 {
			cdsStart = 0
			cdsEnd = 1
		}

		if startOfBlock == 1 && endOfBlock == 0 {

			if cdsStart == 1 && cdsEnd == 0 {

				for _, cdsStEndMatch := range cdsStEnd.FindAllStringSubmatch(val, -1) {

					if cdsStEndMatch[2] != "" {

						lType = cdsStEndMatch[1]
						lStart, _ = strconv.Atoi(strings.TrimSpace(cdsStEndMatch[2]))
						lEnd, _ = strconv.Atoi(strings.TrimSpace(cdsStEndMatch[3]))

						lDir = "f"

						cdsCount[cdsStEndMatch[1]] = cdsCount[cdsStEndMatch[1]] + 1

					} else if cdsStEndMatch[5] != "" {

						lType = cdsStEndMatch[4]
						lStart, _ = strconv.Atoi(strings.TrimSpace(cdsStEndMatch[5]))
						lEnd, _ = strconv.Atoi(strings.TrimSpace(cdsStEndMatch[6]))

						lDir = "r"

						cdsCount[cdsStEndMatch[4]] = cdsCount[cdsStEndMatch[4]] + 1

					}

				}

				for _, cdsNameMatch := range cdsName.FindAllStringSubmatch(val, -1) {
					lName = strings.TrimSpace(strings.Replace(cdsNameMatch[1], " ", "", -1))

				}

				for _, cdsLocusMatch := range cdsLocus.FindAllStringSubmatch(val, -1) {
					lLoc = strings.TrimSpace(strings.Replace(cdsLocusMatch[1], " ", "", -1))

				}
				for _, cdsProdMatch := range cdsProd.FindAllStringSubmatch(val, -1) {
					lProd = strings.TrimSpace(cdsProdMatch[1])

				}
				for _, cdsNoteMatch := range cdsNote.FindAllStringSubmatch(val, -1) {
					lNote = strings.TrimSpace(cdsNoteMatch[1])

				}
				for _, cdsgIDMatch := range cdsgID.FindAllStringSubmatch(val, -1) {
					gID = strings.TrimSpace(strings.Replace(cdsgIDMatch[1], " ", "", -1))

				}
				for _, cdsgprotIDMatch := range cdsprotID.FindAllStringSubmatch(val, -1) {
					pID = strings.TrimSpace(strings.Replace(cdsgprotIDMatch[1], " ", "", -1))

				}
				for _, cdsGOAMatch := range cdsGOA.FindAllStringSubmatch(val, -1) {
					lGOA = strings.TrimSpace(strings.Replace(cdsGOAMatch[1], " ", "", -1))

				}
				for _, cdsPDBMatch := range cdsPDB.FindAllStringSubmatch(val, -1) {
					lPDBArr = append(lPDBArr, strings.TrimSpace(strings.Replace(cdsPDBMatch[1], " ", "", -1)))

				}
				for _, cdsInterProMatch := range cdsInterPro.FindAllStringSubmatch(val, -1) {
					lInterProArr = append(lInterProArr, strings.TrimSpace(strings.Replace(cdsInterProMatch[1], " ", "", -1)))

				}
				for _, cdsProSiteMatch := range cdsProSite.FindAllStringSubmatch(val, -1) {
					lProSiteArr = append(lProSiteArr, strings.TrimSpace(strings.Replace(cdsProSiteMatch[1], " ", "", -1)))

				}
				// fmt.Println(lName)
				// fmt.Println(lProSiteArr)
				// if *flgDev == true {
				// 	fmt.Println(val)
				// }

				if gID != "" && lGOA == "" {
					nucCore = 1

				} else if gID == "" && lGOA != "" {
					nucCore = 1

				} else if gID == "" && lGOA == "" {
					nucCore = 2
				}
				// fmt.Println(nucCore)

				// if lName == "" {
				// 	lName = lLoc
				// }

			} else if cdsStart == 0 && cdsEnd == 1 {

				if lLoc == "" && lName == "" {
					lLoc = fmt.Sprintf("%v_%v_%v", lType, lStart, lEnd)
					// fmt.Println(lLoc)
					// } else if lLoc == "" && lName == "" {
					// 	lLoc = fmt.Sprintf("%v_%v_%v", lType, lStart, lEnd)

				} else if lLoc == "" && lName != "" {
					lLoc = lName
				}

				// if lType == "ncRNA" {
				// 	fmt.Println(lName, lStart, lEnd, lType)
				// }
				// if lLoc == "" && lName == "" {
				// 	lLoc = fmt.Sprintf("%v_%v_%v", lType, lStart, lEnd)

				// }

				if lProd == "" && lNote != "" {
					lProd = lNote
				} else if lProd == "" && lNote == "" {
					lProd = lType
				}

				g := geneInfo{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir,
					Product: lProd, Name: lName, GeneID: gID, ProteinID: pID, Note: lNote, GOA: lGOA, TypeOf: lType, PDB: lPDBArr, InterPro: lInterProArr, ProSite: lProSiteArr}

				allGenesVal = append(allGenesVal, g)

				// geneData[strings.ToUpper(lLoc)] = g
				if lType == "CDS" {
					// fmt.Println(lLoc, lStart, lEnd)
					geneCoordinates[lLoc] = gCoords{Start: lStart, End: lEnd, Type: "CDS"}

					// if lLoc == "Rv1300" {
					// 	fmt.Println(g, geneCoordinates["RV1300"])
					// }
				}

				if *gbVerbose == true {

					fmt.Printf("l:%v s:%v e:%v d:%v p:%v gId:%v pId:%v n:%v GOA:%v T:%v InterPro:%v PDB:%v ProSite:%v\n", lLoc, lStart, lEnd, lDir, lProd, gID, pID, lNote, lGOA, lType, lInterProArr, lPDBArr, lProSiteArr)

				}
				lName, lStart, lEnd, lLoc, lDir, lProd, gID, pID, lGOA, lNote, lType = "", 0, 0, "", "", "", "", "", "", "", ""
				lPDBArr = nil
				lInterProArr = nil
				lProSiteArr = nil

			}
		}

	}

	for key, val := range cdsCount {
		fmt.Printf("\nFound: %v %v", val, key)
	}
	if nucCore == 0 {
		fmt.Println("\nGene Info: Genbank Id Nomenclature", nucCore)
	} else if nucCore == 1 {
		fmt.Println("\nGene Info: high-quality Gene Ontology (GO) annotations", nucCore)
	} else if nucCore == 2 {
		fmt.Println("\nGene Info: unknown", nucCore)
	}

	// sort.Slice(genomeCoordinates, func(i, j int) bool {
	// 				return genomeCoordinates[i].Start < genomeCoordinates[j].Start
	// 			})

	gInfo = genomeInfo{NucCore: nucCore, Organism: organismName, Start: gStart, End: gEnd, Strain: organismStrain, Version: gVersion}
	// go process("Анализ файла закончен!")

	sort.Slice(allGenesVal, func(i, j int) bool {
		return allGenesVal[i].Start < allGenesVal[j].Start
	})

	for i, g := range allGenesVal {

		igensS = append(igensS, g.Start)
		igensE = append(igensE, g.End)
		leftGene = append(leftGene, g.Locus)

		if i >= 1 && i < len(allGenesVal)-1 {

			igrCount++
			checkEndOfCDS := (igensE[i-1] + 1) - g.Start

			if checkEndOfCDS < 0 {
				igr := geneInfo{Start: igensE[i-1] + 1, End: igensS[i] - 1, Locus: fmt.Sprintf("IGR_%v_%v", igensE[i-1]+1, igensS[i]-1), Direction: "f", Product: "Intergenic region (IGR)", Name: "IGR", TypeOf: "IGR"}
				allGenesVal = append(allGenesVal, igr)

				// geneData[strings.ToUpper(fmt.Sprintf("IGR_%v_%v", igensE[i-1]+1, igensS[i]-1))] = igr
				geneCoordinates[fmt.Sprintf("IGR_%v_%v", igensE[i-1]+1, igensS[i]-1)] = gCoords{Start: igensE[i-1] + 1, End: igensS[i] - 1, Type: "IGR"}
			}

		}

	}
	fmt.Println("Found IGR regions:", igrCount)
	// fmt.Println(genomeCoordinates)

	sort.Slice(allGenesVal, func(i, j int) bool {
		return allGenesVal[i].Start < allGenesVal[j].Start
	})

	return allGenesVal, splitedGenome

}

func writeDB(file string, gene []geneInfo) {

	gobFile, err := os.Create(file)
	if err != nil {
		log.Println(err.Error())
	}
	defer gobFile.Close()
	compressedGobFile := lzw.NewWriter(gobFile, lzw.LSB, 8)
	defer compressedGobFile.Close()
	gobParser := gob.NewEncoder(compressedGobFile)
	gobParser.Encode(&gene)
	gobParser.Encode(&genomeSeqSlice)
	gobParser.Encode(&gInfo)
	gobParser.Encode(&geneCoordinates)

	// gobParser.Encode(&featAllGenesVal)
	// gobParser.Encode(&genomeCoordinates)

	fmt.Println(file, " was created successfully.")

}

func readDB(file string) []geneInfo {
	var gene []geneInfo
	gobFile, err := os.Open(file)
	if err != nil {
		log.Println(err.Error())
	}
	defer gobFile.Close()
	compressedGobFile := lzw.NewReader(gobFile, lzw.LSB, 8)
	defer compressedGobFile.Close()
	gobParser := gob.NewDecoder(compressedGobFile)
	gobParser.Decode(&gene)
	gobParser.Decode(&genomeSeqSlice)
	gobParser.Decode(&gInfo)
	gobParser.Decode(&geneCoordinates)

	// gobParser.Decode(&featAllGenesVal)
	// gobParser.Decode(&genomeCoordinates)
	return gene

}

// func codonReverse(codon string) string {
// 	var (
// 		lCodonSplit []string
// 		result      string
// 	)
// 	lCodonSplit = strings.Split(codon, "")
// 	// fmt.Println(lCodonSplit)
// 	for _, n := range lCodonSplit {

// 		switch n {
// 		case "a":
// 			result = "t" + result
// 		case "t":
// 			result = "a" + result
// 		case "c":
// 			result = "g" + result
// 		case "g":
// 			result = "c" + result
// 		case "A":
// 			result = "T" + result
// 		case "T":
// 			result = "A" + result
// 		case "C":
// 			result = "G" + result
// 		case "G":
// 			result = "C" + result
// 		}

// 	}
// 	return result

// }

//

func parserVCF(f string, print bool, genes []geneInfo) []snpInfo {
	var (
		vcf = regexp.MustCompile(`^\S+\s+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)`)
		// var indel = regexp.MustCompile(`^^\S+\W+(\d+)\W+(\w+)\s+(\w+).*(INDEL).*DP=(\d+)`)
		validateVCF = regexp.MustCompile(`(##fileformat)=VCF`)
		vcfValid    bool

		//
		snpFromVCF []snpInfo
	)
	cpuprofile := *annBench
	if cpuprofile != "" {
		f, err := os.Create(cpuprofile)
		if err != nil {
			log.Fatal(err)

		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	file, err := os.Open(f)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		for _, vcfvalid := range validateVCF.FindAllStringSubmatch(scanner.Text(), -1) {

			if vcfvalid[0] == "##fileformat=VCF" {
				vcfValid = true
			}

		}
		if vcfValid == false {
			fmt.Printf("\n%v is not VCF file!!! Check it!\n", file.Name())
			break
		}
		if *annInDel == true {

		}

		for _, match := range vcf.FindAllStringSubmatch(scanner.Text(), -1) {

			if vcfValid == true {
				apos, _ := strconv.Atoi(match[1])
				ref := match[2]
				alt := match[3]

				// for _, g := range genes {
				for z := 0; z < len(genes); z++ {
					// g := genes[z]
					// lStart := genes[z].Start
					// lEnd := genes[z].End

					if apos >= genes[z].Start && apos <= genes[z].End {

						qSnpInfo := &snpInfoQuery{OutChan: make(chan snpInfo), apos: apos, g: genes[z], alt: alt, index: *gbIndex}
						go qSnpInfo.request()
						snp := <-qSnpInfo.OutChan

						// go qSNP.request()
						// // snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
						// snpRes := <-qSNP.OutChan
						// snpCacheMap[snpRes.File] = snpRes.SnpInfo

						// snp := getSnpInfo(apos, g, alt, *gbIndex)
						snp.InterPro = genes[z].InterPro
						snp.PDB = genes[z].PDB
						snp.ProSite = genes[z].ProSite

						// snp.FileName = f
						// fmt.Println(g.PDB)
						// br := testing.Benchmark(snp)
						// fmt.Println(br)

						if len(ref) == 1 && len(alt) == 1 {
							snpFromVCF = append(snpFromVCF, snp)
							if print == true {
								// printResults(snp)
								printTextResults(snp, *gbVerbose)
								// fmt.Println(snp)
							}
						}

					}

				}
			}
		}

	}

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}

	return snpFromVCF

}

func about() {
	fmt.Println("\n", logo)
}

// func process(str string) {
// 	fmt.Printf("%v\r", str)
// }

func parserBulkVCF(withFilenames bool) {

	files := &listOfFiles
	for i, file := range *files {

		snpsChan := make(chan []snpInfo)

		go func() {
			snpsChan <- makeSnps(file)
		}()
		snps := <-snpsChan
		if *gbVerbose == true {
			fmt.Printf(" File %v (%v from %v) was succefull annotated. \r", file, i+1, len(*files))
		}
		if withFilenames == true {
			fmt.Printf("\n\n%v:\n\n", file)

		}
		for _, val := range snps {
			printTextResults(val, *gbVerbose)
		}

	}

}

func makeSeq(typeof string, verbose bool, ref bool, exGenes map[int]int, exSNPs map[int]int) []seqInfo {

	var (
		AllPos, SelectedPos []int
		ResSeq              []seqInfo
		// passSNP = make(map[string]int)
		uniqSNP  = make(map[int]int)
		nbrOfSNP = 2000
		// posCount             = make(map[int]int)
	)
	// files := getListofVCF()

	// queryChan := make(chan vcfInfoQuery)

	if *annSeqLen != 0 {
		nbrOfSNP = *annSeqLen
	}

	files := &listOfFiles
	for i, file := range *files {
		// создаем запрос в виде типа vcfQuery, который передается через канал на выход <-qSNP.OutChan
		qSNP := &vcfQuery{File: file, OutChan: make(chan vcfInfoQuery), Print: verbose}
		go qSNP.request()
		// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
		snpRes := <-qSNP.OutChan
		snpCacheMap[snpRes.File] = snpRes.SnpInfo
		if verbose == true {
			fmt.Printf("Reading files:  %v from %v \r", i+1, len(*files))
		}

	}

	for _, snps := range snpCacheMap {

		for _, val := range snps {
			if len(exGenes) != 0 {

				for key, value := range exGenes {
					if val.APos >= key && val.APos <= value {
						// fmt.Println(val.Locus, val.Start, val.End, val.Product)
						// passSNP["genes"] = passSNP["genes"] + 1
						uniqSNP[val.APos] = 2
						continue
					} else if exSNPs[val.APos] == 1 {
						// passSNP["snp"] = passSNP["snp"] + 1
						uniqSNP[val.APos] = 2
						continue

						// fmt.Println(val.APos)

					} else {
						if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
							uniqSNP[val.APos] = 1

						}
					}
				}
			} else {
				if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
					uniqSNP[val.APos] = 1

				}
				// AllPosUnsort = append(AllPosUnsort, val.APos)
			}

			// if exGenes[val.Locus] != 1 && exSNPs[string(val.APos)] != 1 {
			// 	AllPosUnsort = append(AllPosUnsort, val.APos)
			// }
		}
	}

	// fmt.Println(uniqSNP)

	for key, value := range uniqSNP {

		if value == 1 {

			AllPos = append(AllPos, key)
			// } else if value == 2 {
			// 	fmt.Println(key)

		}

	}
	// go process("Working...          ")
	// AllPos = unique(AllPosUnsort)

	sort.Ints(AllPos)

	rand.Seed(time.Now().UnixNano())

	if nbrOfSNP > len(AllPos) {
		nbrOfSNP = len(AllPos) - 1
	}

	for i := 1; i <= nbrOfSNP; i++ {
		rnd := rand.Intn(len(AllPos)-i) + i
		// fmt.Println(AllPos[rnd])
		SelectedPos = append(SelectedPos, AllPos[rnd])
	}
	// fmt.Println(AllPos[0])

	// rand.Intn(max - min) + min

	sort.Ints(SelectedPos)

	if ref == true {
		switch typeof {
		case ncFlag:
			var refBuffer strings.Builder
			refBuffer.WriteString(fmt.Sprintf(">%v\n", "REFERENCE"))
			for _, allpos := range SelectedPos {
				refBuffer.WriteString(getNucFromGenomePos(allpos))
			}
			ResSeq = append(ResSeq, seqInfo{Name: "reference", Seq: refBuffer.String(), UsedPositions: SelectedPos})

		}
	}

	// fmt.Println(AllPos)
	for fname, snps := range snpCacheMap {

		// if verbose == true {
		// 	fmt.Printf("Generating sequences: Working on  %v from %v \r", i+1, len(snps.File))
		// }
		pos := make(map[int]string)
		var buffer strings.Builder

		buffer.WriteString(fmt.Sprintf(">%v\n", strings.ToUpper(fname)))

		switch typeof {
		case ncFlag:
			for _, val := range snps {
				pos[val.APos] = val.Alt
				// if exGenes[val.Locus] == 1 {
				// 	fmt.Println(val.Locus)
				// }

			}
			for _, allpos := range SelectedPos {
				// posCount[allpos] = posCount[allpos] + 1
				if pos[allpos] != "" {

					buffer.WriteString(pos[allpos])
				} else {

					buffer.WriteString(getNucFromGenomePos(allpos))
				}

			}

			ResSeq = append(ResSeq, seqInfo{Name: fname, Seq: buffer.String(), UsedPositions: SelectedPos})
			// fmt.Println(len(buffer.String()))

		case aaFlag:

			fmt.Println("Reserved")

		}
		// if *gbDebug == true {
		// 	fmt.Printf("%v\t:\nThere was passed %v SNPs from exclude gene file\n And %v SNPs from exclude snp file\n", fname, passSNP["genes"], passSNP["snp"])
		// }
	}

	return ResSeq
}

func unique(list []int) []int {
	uniqueSet := make(map[int]bool, len(list))
	for _, x := range list {
		uniqueSet[x] = true
	}
	result := make([]int, 0, len(uniqueSet))
	for x := range uniqueSet {
		result = append(result, x)
	}
	return result
}

func getInterGen(pos int) {
	var (
		igensS []int
		igensE []int
	)
	for i, g := range allGenesVal {

		igensS = append(igensS, g.Start)
		igensE = append(igensE, g.End)
		// igenStart, _ := strconv.Atoi(g.End)
		if i >= 1 && i < len(allGenesVal)-1 {
			// 	// igenEnd, _ := strconv.Atoi(g.Start)
			fmt.Printf("igen:s%v e:%v\n", igensS[i], igensE[i-1])
		}
		fmt.Printf("l:%v s:%v e:%v d:%v %v\n", g.Locus, g.Start, g.End, g.Direction, i)

	}
}

func printTextResults(snps snpInfo, verbose bool) {

	const fullAnnotations = " {{if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T0\"))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T1\"))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Tang}}({{.TangIdxVal}})\t{{.Product}}\n" +
		"{{else if (and (ne .TypeOf \"CDS\") (eq .ReportType \"T0\"))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.APos}}{{.NucInPos}}>{{.Alt}}\t-\t-\t-\t-\t{{.Product}}\n" +
		"{{else if (and (ne .TypeOf \"CDS\") (eq .ReportType \"T1\"))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.APos}}{{.NucInPos}}>{{.Alt}}\t-\t-\t-\t-\t-\t{{.Product}}\n" +
		"{{end}}"

	const cdsAnnotations = "{{if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T0\"))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T1\"))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Tang}}({{.TangIdxVal}})\t{{.Product}}\n" +
		"{{end}}"

	t := template.New("report")
	if verbose == true {
		t, _ = t.Parse(fullAnnotations)
		t.Execute(os.Stdout, snps)

	} else {
		t, _ = t.Parse(cdsAnnotations)
		t.Execute(os.Stdout, snps)
	}

}

func exportToExcel(snps []snpInfo, file string) {

	xlsx := excelize.NewFile()
	// Create a new sheet.
	index := xlsx.NewSheet("SNPs")
	// Set value of a cell.
	var locuses []snpInfo

	// xlsx.SetCellValue("Sheet1", "B2", 100)
	// Set active sheet of the workbook.
	xlsx.SetActiveSheet(index)
	// Save xlsx file by the given path.

	for _, locus := range snps {

		if locus.Locus != "" {
			locuses = append(locuses, locus)

		}
	}

	for i, snp := range locuses {
		// fmt.Println(snp.Locus)
		// fmt.Println(snp.APos, i)

		if i > 0 {
			xlsx.SetCellValue("SNPs", fmt.Sprintf("A%v", i+1), snp.Locus)
			xlsx.SetCellValue("SNPs", fmt.Sprintf("B%v", i+1), snp.Name)
			xlsx.SetCellValue("SNPs", fmt.Sprintf("C%v", i+1), snp.APos)
			xlsx.SetCellValue("SNPs", fmt.Sprintf("D%v", i+1), fmt.Sprintf("%v%v>%v", snp.PosInGene, snp.NucInPos, snp.Alt))
			xlsx.SetCellValue("SNPs", fmt.Sprintf("E%v", i+1), fmt.Sprintf("%v/%v", snp.RefCodon, snp.AltCodon))
			xlsx.SetCellValue("SNPs", fmt.Sprintf("F%v", i+1), fmt.Sprintf("%v%v%v", snp.RefAA, snp.CodonNbrInG, snp.AltAA))
			xlsx.SetCellValue("SNPs", fmt.Sprintf("G%v", i+1), snp.Mutation)
			xlsx.SetCellValue("SNPs", fmt.Sprintf("H%v", i+1), snp.Tang)
			xlsx.SetCellValue("SNPs", fmt.Sprintf("I%v", i+1), snp.TangIdxVal)
			xlsx.SetCellValue("SNPs", fmt.Sprintf("J%v", i+1), snp.Product)

			// fmt.Println(xlsx.GetCellValue("Sheet2", fmt.Sprintf("A%v", i)))
			// fmt.Println(fmt.Sprintf("A%v - %v", i, snp.Locus))
		} else {
			xlsx.SetCellValue("SNPs", "A1", "Locus")
			xlsx.SetCellValue("SNPs", "B1", "Gene")
			xlsx.SetCellValue("SNPs", "C1", "Position")
			xlsx.SetCellValue("SNPs", "D1", "Mutation")
			xlsx.SetCellValue("SNPs", "E1", "Codons")
			xlsx.SetCellValue("SNPs", "F1", "AA")
			xlsx.SetCellValue("SNPs", "G1", "Type")
			xlsx.SetCellValue("SNPs", "H1", "Complex Index")
			xlsx.SetCellValue("SNPs", "I1", "Index Sum")
			xlsx.SetCellValue("SNPs", "J1", "Product")

		}

	}
	err := xlsx.SaveAs(file)
	if err != nil {
		fmt.Println(err)
	} else {
		fmt.Printf("The %v file was created successfully\n", file)
	}
}

func printWebResults(snps []snpInfo, port string) {

	var htmlTitle = `   <!DOCTYPE html>
				<html>

			<table width="100%" cellspacing="0" cellpadding="4" border="1">
			<tbody>
			<tr>
			
			{{if eq .NucCore 1}}
			
				<td>Locus</td><td>Gene</td><td>Pos.</td><td>Mutation</td><td>Codons</td>
				<td>AA</td><td>Type</td><td>Product</td><td><p title="Gene Ontology Annotation (GOA) Database">GOA</p></td><td><p title="The Universal Protein Resource (UniProt) is a comprehensive resource for protein sequence and annotation data">UniProt</p></td><td><p title="InterPro provides functional analysis of proteins by classifying them into families and predicting domains and important sites">InterPro</p></td><td><p title="Protein Data Bank (PDB)">PDB</p></td><td>ProSite</td>
			{{else if eq .NucCore 0}}
				<td>Locus</td><td>Gene</td><td>Pos.</td><td>Mutation</td><td>Codons</td>
				<td>AA</td><td>Type</td><td>Product</td><td>GeneID</td><td>UniProt</td>
			{{else if  eq .NucCore 2}}
				<td>Locus</td><td>Gene</td><td>Pos.</td><td>Mutation</td><td>Codons</td>
				<td>AA</td><td>Type</td><td>Product</td><td>-</td><td>UniProt</td>
			{{end}}
			
			</tr>
			
`
	var htmlTemplate = `
			
			{{/* 1 */}}
				{{range $element := .}}
				{{/* 2 */}}
			
			{{if .GeneID}}
				<tr>
					
				<td><p title="{{.Note}}">{{.Locus}}</p></td><td>{{.Name}}</td><td>{{.APos}}</td><td>{{.PosInGene}}{{.NucInPos}}>{{.Alt}}</td>
				<td>{{.RefCodon}}/{{.AltCodon}}</td><td><p title="{{.RefAA}}{{.CodonNbrInG}}{{.AltAA}}">{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}</p></td>
				<td><p title="Complex Index: {{.Tang}}({{.TangIdxVal}})">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}"target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ncbi.nlm.nih.gov/gene/{{.GeneID}}={{.GeneID}}" target="_blank">{{.GeneID}}</a>
				</td><td><a href="http://www.uniprot.org/uniprot/?query={{.ProteinID}}&sort=score">{{.ProteinID}}</td>
				</tr>
					{{/* 2 */}}
			{{else}}
				<tr>
				
					{{/* 3 */}}
			{{if eq .TypeOf "CDS"}}			
				<td><p title="{{.Note}}">{{.Locus}}</p></td><td>{{.Name}}</td><td>{{.APos}}</td><td>{{.PosInGene}}{{.NucInPos}}>{{.Alt}}</td>
				<td>{{.RefCodon}}/{{.AltCodon}}</td><td><p title="{{.RefAA}}{{.CodonNbrInG}}{{.AltAA}}">{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}</p></td>			
					{{/* 4 */}}
			{{if eq .Mutation "missense"}}
					{{if eq .TangIdxVal 4 }}
				<td bgcolor="#CECEF6"><u><p title="Complex Index: {{.Tang}}({{.TangIdxVal}})">{{.Mutation}}</p></u></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}" target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ebi.ac.uk/QuickGO/GProtein?ac={{.GOA}}" target="_blank">{{.GOA}}</a>
					{{else}}
	<td bgcolor="#CECEF6"><p title="Complex Index: {{.Tang}}({{.TangIdxVal}})">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}" target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ebi.ac.uk/QuickGO/GProtein?ac={{.GOA}}" target="_blank">{{.GOA}}</a>
					{{end}}
					{{/* 4 */}}
			{{else if eq .Mutation "nonsense"}}
				<td bgcolor="#F78181"><p title="Complex Index: {{.Tang}}({{.TangIdxVal}})">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}" target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ebi.ac.uk/QuickGO/GProtein?ac={{.GOA}}" target="_blank">{{.GOA}}</a>
					{{/* 4 */}}
			{{else}}
				<td bgcolor="#ddffcc"><p title="Complex Index: {{.Tang}}({{.TangIdxVal}})">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}" target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ebi.ac.uk/QuickGO/GProtein?ac={{.GOA}}" target="_blank">{{.GOA}}</a>
					{{/* 4 */}}
			{{end}}
				</td><td><a href="http://www.uniprot.org/uniprot/?query={{.ProteinID}}&sort=score" target="_blank">{{.ProteinID}}</td>
				<td>
				{{ range $value := .InterPro }}
   				<a href="http://www.ebi.ac.uk/interpro/entry/{{$value}}" target="_blank">{{$value}}</a>
				{{ end }}
				</td>
				<td>
				{{ range $value := .PDB }}
   				<a href="https://www.rcsb.org/structure/{{$value}}" target="_blank">{{$value}}</a>
				{{ end }}			
				</td>
				<td>
				{{ range $value := .ProSite }}
   				<a href="https://prosite.expasy.org/{{$value}}" target="_blank">{{$value}}</a>
				{{ end }}			
				</td>
				</tr>
				{{/* 3 */}}
			{{end}}
				{{/* 2 */}}
			{{end}}
				{{/* 1 */}}		
			{{end}}
			
			</tbody>
			</table>
			<table width="100%" cellspacing="0" cellpadding="4" border="0">
			<tr>
			<td><a href="http://bsatool.ru" target="_blank">Created by BSATool (Bacterial Snp Annotation Tool)</a></td>
			</tr>
			</table>
			`

	t := template.New("t")
	tt := template.New("t1")
	// t1 := template.New("t1")
	t, err := t.Parse(htmlTemplate)
	t1, err := tt.Parse(htmlTitle)

	// t1 := template.New("t1")
	// t2, err := t1.Parse(htmlTitle)

	if err != nil {
		panic(err)
	}

	browser.OpenURL(fmt.Sprintf("localhost:%v", port))

	http.HandleFunc("/", func(w http.ResponseWriter, r *http.Request) {
		// err = t.Execute(w, &gInfo)

		err = t1.Execute(w, &gInfo)
		if err != nil {
			panic(err)
		}
		err = t.Execute(w, snps)
		if err != nil {
			panic(err)
		}
		go func() {
			defer os.Exit(0)
		}()
	})

	locPort := fmt.Sprintf(":%v", port)
	http.ListenAndServe(locPort, nil)
}

func checkSNPnotation(inString string) snpCheckInfo {
	var (
		snpCheck snpCheckInfo
	)
	rLMN := regexp.MustCompile(lmnRegExp) // L:1 PiG:2 REF:3 ALT:4 NAME:5
	// LocusMutationName(LMN)
	rPMLN := regexp.MustCompile(pmlnRegExp) // APOS:1 REF:2 ALT:3 L:4 NAME:5
	//PositionMutationLocusName (PMLN)
	rPMN := regexp.MustCompile(pmnRegExp)     // APOS:1 REF:2 ALT:3 NAME:4
	rLSAAN := regexp.MustCompile(lsaanRegExp) // L:1 AA_REF:2 POS:3 AA_ALT:4 NAME:5
	// LocusShortAminoAcidName (LSAAN)
	rLLAAN := regexp.MustCompile(llaanRegExp) // L:1 LAA_REF:2 PiG:3 LAA_ALT:4 NAME:5
	// LocusLongAminoAcidName (LLAAN)
	rLCN := regexp.MustCompile(lcnRegExp)
	// rPMNT := regexp.MustCompile(pmntRegExp) //Pos(1)_Ref(2)>Alt(3){Tab}NAME(4):|;(5)Tag(6)

	for _, matchLMN := range rLMN.FindAllStringSubmatch(inString, -1) {
		snpCheck = snpCheckInfo{Locus: matchLMN[1], PosInGene: matchLMN[2], Ref: matchLMN[3], Alt: matchLMN[4], Name: matchLMN[5], TypeOf: tLMN, Raw: matchLMN[0]}

		// fmt.Printf("%v\n", snpCheck)
		// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
	}
	for _, matchPMLN := range rPMLN.FindAllStringSubmatch(inString, -1) {
		snpCheck = snpCheckInfo{APos: matchPMLN[1], Ref: matchPMLN[2], Alt: matchPMLN[3], Locus: matchPMLN[4], Name: matchPMLN[5], TypeOf: tPMLN, Raw: matchPMLN[0]}

		// fmt.Printf("%v\n", snpCheck)
		// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
	}
	for _, matchPMN := range rPMN.FindAllStringSubmatch(inString, -1) {

		snpCheck = snpCheckInfo{APos: matchPMN[1], Ref: matchPMN[2], Alt: matchPMN[3], Name: matchPMN[4], TypeOf: tPMN, Raw: matchPMN[0]}

		// fmt.Printf("%v\n", snpCheck)
		// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
	}

	// for _, matchPMNT := range rPMNT.FindAllStringSubmatch(scanner.Text(), -1) {

	// 	snpCheck = snpCheckInfo{APos: matchPMNT[1], Ref: matchPMNT[2], Alt: matchPMNT[3], Name: matchPMNT[4], TypeOf: tPMNT, Raw: matchPMNT[0], Tag: matchPMNT[6]}

	// 	parsedSNP = append(parsedSNP, snpCheck)

	// 	// fmt.Println("!!!")
	// 	// fmt.Printf("%v\n", snpCheck.Tag)
	// 	// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
	// }

	for _, matchLSAAN := range rLSAAN.FindAllStringSubmatch(inString, -1) {
		snpCheck = snpCheckInfo{Locus: matchLSAAN[1], AASref: matchLSAAN[2], CodonNbrInG: matchLSAAN[3], AASalt: matchLSAAN[4], Name: matchLSAAN[5], TypeOf: tLSAAN, Raw: matchLSAAN[0]}

		// fmt.Printf("%v\n", snpCheck)
		// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
	}
	for _, matchLLAAN := range rLLAAN.FindAllStringSubmatch(inString, -1) {
		snpCheck = snpCheckInfo{Locus: matchLLAAN[1], AALref: matchLLAAN[2], CodonNbrInG: matchLLAAN[3], AALalt: matchLLAAN[4], Name: matchLLAAN[5], TypeOf: tLLAAN, Raw: matchLLAAN[0]}

		// fmt.Printf("%v\n", snpCheck)
		// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
	}
	for _, matchLCN := range rLCN.FindAllStringSubmatch(inString, -1) {
		snpCheck = snpCheckInfo{Locus: matchLCN[1], CodonNbrInG: matchLCN[2], Name: matchLCN[3], TypeOf: tLCN, Raw: matchLCN[0]}

		// fmt.Printf("%v\n", snpCheck)
		// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
	}

	return snpCheck
}

func readSNPFromFile(f string) []snpCheckInfo {
	var (
		parsedSNP    []snpCheckInfo
		snpCheckChan = make(chan snpCheckInfo)
	)
	file, err := os.Open(f)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		if !strings.HasPrefix(scanner.Text(), "#") {

			go func() {

				snpCheckChan <- checkSNPnotation(scanner.Text())
			}()
			snpCheck := <-snpCheckChan
			parsedSNP = append(parsedSNP, snpCheck)

		}
	}

	return parsedSNP
}

func getShareSNP(verbose bool, web bool, files []string) {
	// var countSNPs = 1
	var (
		pos   = map[uint64]int{}
		alt   = map[uint64]snpInfo{}
		share []uint64
		// var g geneInfo
		snpToWeb     []snpInfo
		snpToConsole []snpInfo
	// var snpCacheFromChan []vcfInfoQuery
	)
	upperLimit := len(files)
	// bar := pb.StartNew(len(files))
	// bar := pb.ProgressBarTemplate(pbtmpl).Start(len(files))
	// bar.SetWidth(90)
	// tmp := 0
	for i, file := range files {
		// tmp++
		// fmt.Printf("processed %v files\r", cnt+1)
		if verbose == true {
			fmt.Printf("Pass: %v from %v \r", i+1, len(files))
		}

		qSNP := &vcfQuery{File: file, OutChan: make(chan vcfInfoQuery), Print: verbose}
		go qSNP.request()
		// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
		snpRes := <-qSNP.OutChan
		snpCacheMap[snpRes.File] = snpRes.SnpInfo
	}

	for _, snps := range snpCacheMap {
		// fmt.Println(tmp, "2163790")
		for _, val := range snps {
			hash := getHashSNP(val)
			pos[hash] = pos[hash] + 1
			// fmt.Println(pos[val.APos+val.Start], val.Locus, upperLimit)
			// fmt.Println(val.APos, file, pos[val.APos])

			alt[hash] = val

			// if hash == "efe3a457c5b080a8f7abff7c37883febc1d722a9" {
			// 	fmt.Println(pos[hash], hash, val)
			// }
			// if pos[val.APos] > upperLimit {
			// 	fmt.Println(pos[val.APos], val.APos)
			// }
			// 	g = geneInfo{Locus: val.Locus, Start: val.Start, End: val.End, Name: val.Name, Product: val.Product, Direction: val.Direction, GeneI } else if lPos > upperLimit {
			// 	fmt.Println(i, alt[i], lPoD: val.GeneID, ProteinID: val.ProteinID}
			// }
			// hash = ""

		}

		for i, lPos := range pos {
			if lPos == upperLimit {
				share = append(share, i)

				//s)
				// } else if lPos > upperLimit {
				// 	fmt.Println(lPos, i, alt[i])
			}
			// fmt.Println(lPos, i)
			// if i == 5683364590609038770 {
			// fmt.Println(i, alt[i], lPos)
			// }
			// fmt.Println(i, alt[i], lPos)
		}

		// sort.Slice(share, func(i, j int) bool {
		// 	return share[i] < share[j]
		// })

	}
	// fmt.Println(share)
	for _, sharePos := range share {
		// countSNPs++
		// fmt.Println(countSNPs)
		// fmt.Println(alt[sharePos])

		if web == true {
			snpToWeb = append(snpToWeb, alt[sharePos])

			// printResults(alt[sharePos])

		} else {
			// printResults(alt[sharePos])
			snpToConsole = append(snpToConsole, alt[sharePos])
		}
	}

	// fmt.Println()
	// if *gbDebug == true {
	// 	fmt.Printf("f:%v snp: %v\n%v\n", len(files), countSNPs, files)
	// }

	//--------------------------------------------
	if web == true && len(snpToWeb) != 0 {
		sort.Slice(snpToWeb, func(i, j int) bool {
			return snpToWeb[i].Start < snpToWeb[j].Start
		})
		printWebResults(snpToWeb, *gbPort)
	} else if web == false && len(snpToConsole) != 0 {

		sort.Slice(snpToConsole, func(i, j int) bool {
			return snpToConsole[i].Start < snpToConsole[j].Start
		})
		for _, res := range snpToConsole {
			// printResults(res)

			printTextResults(res, verbose)
		}
	}
	//---------------------------------------------

}

func getSNPFromSnpInfo(snps []snpInfo, c chan int) {
	for _, val := range snps {
		c <- val.APos
	}
}

func snpStat() {
	// var countSNPs = 1
	var (
		pos                  = make(map[int]int)
		alt                  = make(map[int]snpInfo)
		f                    = make(map[int][]string)
		positions, posUnsort []int
	)
	// files := getListofVCF()
	files := &listOfFiles

	upperLimit := len(*files)

	for i, file := range *files {

		fmt.Printf("Reading: %v (%v from %v)%v \r", file, i+1, len(*files), strings.Repeat(" ", 60))
		// }
		// bar.Increment()
		// snps := parserVCF(file, false, allGenesVal)
		qSNP := &vcfQuery{File: file, OutChan: make(chan vcfInfoQuery)}
		go qSNP.request()
		// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
		snpRes := <-qSNP.OutChan
		snpCacheMap[snpRes.File] = snpRes.SnpInfo
	}

	for fname, snps := range snpCacheMap {
		for _, val := range snps {
			if pos[val.APos] <= upperLimit {
				pos[val.APos] = pos[val.APos] + 1        //count
				alt[val.APos] = val                      // pos
				f[val.APos] = append(f[val.APos], fname) //files
				posUnsort = append(posUnsort, val.APos)  //array of positions
			}

		}
		positions = unique(posUnsort)
		sort.Ints(positions)

	}
	var stat []statInfo

	for _, p := range positions {
		perc := (pos[p] * 100) / upperLimit
		filesNotInList := compareSlices(*files, f[p])
		stat = append(stat, statInfo{Pos: p, Count: pos[p], Perc: perc, FilesWith: strings.Join(f[p], ",\n"), FilesWithout: strings.Join(filesNotInList, ",\n")})

		// }
	}
	printWebStat(stat, *gbPort)
	// fmt.Println(stat)
}

func calcJaroWinklerDist(file string, print bool) []JaroWinklerInfo {

	var (
		altPositions = make(map[string][]allPositionsInGene)
		validData    []string
		altSequences []string
		jwRes        []JaroWinklerInfo
		// var dndsArray []DnDsRes

		jarwinkl float64
	)
	// snps := parserVCF(file, false, allGenesVal)
	if len(snpCacheMap) == 0 {
		qSNP := &vcfQuery{File: file, OutChan: make(chan vcfInfoQuery)}
		go qSNP.request()
		// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
		snpRes := <-qSNP.OutChan
		snpCacheMap[snpRes.File] = snpRes.SnpInfo
	}

	for fname, snps := range snpCacheMap {
		if fname == file {
			for _, val := range snps {
				// fmt.Println(val.Locus, val.PosInGene, val.Alt)
				if val.TypeOf == "CDS" {
					altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, alt: val.Alt, ref: val.NucInPos, locus: val.Locus})
				}
			}
		}
	}

	// выбор генов, в которых обнаружены мутации
	for key, val := range altPositions {
		// fmt.Println(key, len(val))
		if len(val) > 1 {
			validData = append(validData, key)
		}
	}

	sort.Strings(validData)

	for _, val := range validData {

		// start, end := getGenePosByName(val)
		altS := makeAltString(val, altPositions[val])
		altSequences = append(altSequences, altS)
		refS := getGeneSequence(val)
		//  jarwinkl, _ = textdistance.JaroDistance(refS, altS)
		jarwinkl = smetrics.JaroWinkler(refS, altS, 0.7, 4)
		jwRes = append(jwRes, JaroWinklerInfo{Locus: val, JWDist: jarwinkl})
		// fmt.Println(val, "\n", refS)
		// prod := getProductByName(val)
		// 	N, S, PN, PS, DN, DS, DNDS float64
		// ND, NS                     int

		// if dnds.ND != 0 && dnds.NS != 0 {
		// 	i++
		if print == true {
			fmt.Printf("L:%v JW:%v\n", val, jarwinkl)
		}
		// 		// 	fmt.Printf("\nl:%v\tdn/ds:%.2f (%v:%v)\nREF:%v\nALT:%v\n", val, dnds.DNDS, start, end, refS, altS)
		// 		// 	// fmt.Printf("Calculated %v dn/ds from %v\r", i, len(validData))
		// 		// } else if print == true && *gbNoSeq == true {
		// 		fmt.Printf("\nl:%v\tdn/ds:%.2f (%v:%v)\t%v", val, dnds.DNDS, start, end, prod)
		// 	}

		// 	dndsArray = append(dndsArray, DnDsRes{Locus: val, Product: prod, N: dnds.N, S: dnds.S, ND: dnds.ND, NS: dnds.NS, PN: dnds.PN, PS: dnds.PS, DN: dnds.DN, DS: dnds.DS, DNDS: dnds.DNDS})
		// }

	}
	// for i := 1; i < len(altSequences); i++ {
	// 	// for j := 1; j < len(altSequences); j++ {
	// 	// fmt.Println(altSequences[i], altSequences[j])
	// 	jarwinkl, _ = textdistance.JaroDistance(altSequences[i], altSequences[i-1])
	// 	if print == true {
	// 		fmt.Println(jarwinkl)
	// 		jwArray = append(jwArray, jarwinkl)
	// 	}
	// 	// }
	// }
	// fmt.Println(jwRes)
	return jwRes

}

func calcGC3Val(snps []snpInfo) []GC3Type {
	// var wg sync.WaitGroup
	var (
		altPositions = make(map[string][]allPositionsInGene)
		// var validData []string
		gcArray []GC3Type
		// wg.Add(1)

		allLocusUnsort, allLocuses []string
	)

	for _, val := range snps {

		if val.TypeOf == "CDS" {
			// fmt.Println(file, val.Locus)
			altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, alt: val.Alt, ref: val.NucInPos, locus: val.Locus})
			// fmt.Println(val.PosInGene)
			allLocusUnsort = append(allLocusUnsort, val.Locus)
		}

	}

	allLocuses = removeStringDuplicates(allLocusUnsort)

	for _, loc := range allLocuses {

		// start, end := getGenePosByName(loc)
		// fmt.Println(start, end, loc)
		// fmt.Println(start, end, altPositions[	loc], geneCoordinates[strings.ToUpper(loc)])
		// refS := getNucFromGenome(start-1, end)
		refS := getGeneSequence(loc)
		// fmt.Println(refS, end-start)
		altS := makeAltString(loc, altPositions[loc])

		// fmt.Println(altPositions[loc], len(refS), altS, loc)
		// if len(altS) == 0 {
		// 	altS = refS
		// 	fmt.Println("!!!!!")
		// }

		gcAlt, _, _, gc3Alt := codon.GcCodonCalc(altS)
		gcRef, _, _, gc3Ref := codon.GcCodonCalc(refS)

		gcArray = append(gcArray, GC3Type{Locus: loc, GC3Alt: gc3Alt, GC3Ref: gc3Ref, GCalt: gcAlt, GCref: gcRef})
		// fmt.Println(gcArray)
	}
	// wg.Wait()

	return gcArray

}

func compareSlices(slice1 []string, slice2 []string) []string {
	diffStr := []string{}
	m := map[string]int{}

	for _, s1Val := range slice1 {
		m[s1Val] = 1
	}
	for _, s2Val := range slice2 {
		m[s2Val] = m[s2Val] + 1
	}

	for mKey, mVal := range m {
		if mVal == 1 {
			diffStr = append(diffStr, mKey)
		}
	}

	return diffStr
}

func printWebStat(stat []statInfo, port string) {
	var htmlTemplate = `   <!DOCTYPE html>
				<html>
				<head>
				<style>
				.col {
				word-wrap: break-word; /* Перенос слов */
				}
				</style>
			</head>
			<table width="100%" cellspacing="0" cellpadding="4" border="1">
			<tbody>
			<tr>			
			<td>Pos</td><td>Count</td><td>Percent</td><td>+</td><td>-</td>
			</tr>
			<tr>
			{{range $element := .}}
			{{if eq .Perc 100}}				
			<td>{{.Pos}}</td><td>{{.Count}}</td><td bgcolor="#ffe6e6">{{.Perc}}%</td><td>ALL FILES</td><td>NONE</td>
			{{else if eq .Perc 5}}				
			<td>{{.Pos}}</td><td>{{.Count}}</td><td bgcolor="#ddffcc">{{.Perc}}%</td><td>{{.FilesWith}}</td><td>{{.FilesWithout}}</td>
			{{else}}
			<td>{{.Pos}}</td><td>{{.Count}}</td><td>{{.Perc}}%</td><td>{{.FilesWith}}</td><td>{{.FilesWithout}}</td>
			{{end}}

			</tr>	
			{{end}}		
			</tbody>
			</table>
			</table>
			<table width="100%" cellspacing="0" cellpadding="4" border="0">
			<tr>
			<td><a href="http://bsatool.ru" target="_blank">Created by BSATool (Bacterial Snp Annotation Tool)</a></td>
			</tr>
			</table>
			
			
`
	t := template.New("t")
	// t1 := template.New("t1")
	t, err := t.Parse(htmlTemplate)

	// t1 := template.New("t1")
	// t2, err := t1.Parse(htmlTitle)

	if err != nil {
		panic(err)
	}

	browser.OpenURL(fmt.Sprintf("localhost:%v", port))

	http.HandleFunc("/", func(w http.ResponseWriter, r *http.Request) {
		// err = t.Execute(w, &gInfo)

		err = t.Execute(w, stat)
		if err != nil {
			panic(err)
		}
		go func() {
			defer os.Exit(0)
		}()
	})
	// if *flgPort != 0 {
	// 	locPort := fmt.Sprintf(":%v", *flgPort)
	// 	http.ListenAndServe(locPort, nil)
	// } else {
	// 	http.ListenAndServe(":8080", nil)
	// }
	locPort := fmt.Sprintf(":%v", port)
	http.ListenAndServe(locPort, nil)
}

func createNCWebServer(port string, exGenes map[int]int, exSNPs map[int]int) {
	/*

	 */

	seq := makeSeq(ncFlag, *gbVerbose, *annMakeSeqRef, exGenes, exSNPs)
	var htmlTemplate = `
<!DOCTYPE html>
<html>
<head>
<style>
.col {
word-wrap: break-word; /* Перенос слов */
}
</style>
</head>
<body>
<div class="col">
{{range $element := .}}

<p>{{.Seq}}</p>
		{{end}}
</div>
</body>


`

	t := template.New("t")
	t, err := t.Parse(htmlTemplate)
	if err != nil {
		panic(err)
	}

	http.HandleFunc("/", func(w http.ResponseWriter, r *http.Request) {

		err = t.Execute(w, seq)
		if err != nil {
			panic(err)
		}

		go func() {
			defer os.Exit(0)
		}()
	})

	// if *flgPort != 0 {
	// 	locPort := fmt.Sprintf(":%v", *flgPort)
	// 	http.ListenAndServe(locPort, nil)
	// } else {
	// 	http.ListenAndServe(":8080", nil)
	// }
	// resp, err := http.Get(fmt.Sprintf("http://localhost:%v", port))
	// if err != nil {
	// 	fmt.Println("Failed:", err)

	// }
	// resp.Body.Close()
	// browser.OpenURL(fmt.Sprintf("localhost:%v", port))

	locPort := fmt.Sprintf(":%v", port)

	http.ListenAndServe(locPort, nil)

	// fmt.Println("Для выхода из программы, нажмите Ctl+C")

}

func getListofVCF() []string {
	/*
		возвращает список VCF файлов в папке в виде массива
	*/
	files, err := filepath.Glob(vcfExt)

	if err != nil {
		log.Fatal(err)
	}
	sort.Strings(files)
	return files
}

func checkSNPfromFile(f string, verbose bool, web bool, useRule bool) {
	/*
	   tLMN    = "LMN"
	   tPMLN   = "PMLN"
	   tLSAAN  = "LSAAN"
	   tLLAAN  = "LLAAN"
	*/
	mapOfVCF := make(map[string][]snpInfo)
	mapofSNP := make(map[string][]string)
	// var snpArr []checkSNP
	// files := getListofVCF()
	files := &listOfFiles

	locSNPcheck := readSNPFromFile(f)
	var chkSNP []checkSNP

	// var buffer strings.Builder
	// bar := pb.ProgressBarTemplate(pbtmpl).Start(len(files))
	// bar.SetWidth(90)
	for i, file := range *files {
		if verbose == true {
			fmt.Printf("Reading %v (%v:%v)%v\r", file, i+1, len(*files), strings.Repeat(" ", 60))
		}
		// snps := parserVCF(file, false, allGenesVal)

		// bar.Increment()
		// i++
		// fmt.Printf("processed %v files\r", i)
		qSNP := &vcfQuery{File: file, OutChan: make(chan vcfInfoQuery)}
		go qSNP.request()
		// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
		snpRes := <-qSNP.OutChan
		snpCacheMap[snpRes.File] = snpRes.SnpInfo

	}

	for fname, snps := range snpCacheMap {
		mapOfVCF[fname] = snps
	}
	// bar.Finish()
	// fmt.Println()

	for file, snp := range mapOfVCF {
		if verbose == true {
			fmt.Printf("Working on %v                     \r", file)
		}
		// buffer.WriteString(fmt.Sprintf("%v ", file))

		for _, val := range locSNPcheck {

			lGpoS, _ := strconv.Atoi(val.PosInGene)
			CodonNbrInG, _ := strconv.Atoi(val.CodonNbrInG)
			lAPos, _ := strconv.Atoi(val.APos)
			// fmt.Println(val.TypeOf)
			for _, snpFromFile := range snp {

				switch val.TypeOf {
				case tLMN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && strings.ToUpper(val.Alt) == strings.ToUpper(snpFromFile.Alt) && lGpoS == snpFromFile.PosInGene {
						// fmt.Println(val.Locus, "\t", val.PosInGene, "\t", snpFromFile.PosInGene, "\t", snpFromFile.Locus)
						// chkSNP = checkSNP{FileName: file, FoundSNP: fmt.Sprintf("%v[%v:%v%v>%v]\t", val.Name, val.Locus, val.PosInGene, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt))}
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:%v%v>%v]", val.Name, val.Locus, val.PosInGene, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))
						// buffer.WriteString(fmt.Sprintf("%v_%v:%v%v>%v\t", val.Name, val.Locus, val.PosInGene, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))
						//  fmt.Println(strings.ToUpper(buffer.String()))
					}
				case tPMLN:

					if lAPos == snpFromFile.APos && strings.ToUpper(val.Alt) == strings.ToUpper(snpFromFile.Alt) {
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:%v_%v>%v]", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))
						// chkSNP = checkSNP{FileName: file, FoundSNP: fmt.Sprintf("%v[%v:%v_%v>%v]\t", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt))}
						// snpArr = append(snpArr, checkSNP{FileName: file, FoundSNP: fmt.Sprintf("%v[%v:%v_%v>%v]\t", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt))})
						// buffer.WriteString(fmt.Sprintf("%v_%v:%v_%v>%v\t", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))
					}
				case tLSAAN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && CodonNbrInG == snpFromFile.CodonNbrInG &&

						val.AASalt == snpFromFile.AltAAShort {
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:%v%v%v]", val.Name, val.Locus, val.AASref, CodonNbrInG, val.AASalt))
						// buffer.WriteString(fmt.Sprintf("%v_%v:%v%v%v\t", val.Name, val.Locus, val.AASref, CodonNbrInG, val.AASalt))
					}
				case tLLAAN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && CodonNbrInG == snpFromFile.CodonNbrInG &&
						val.AALalt == snpFromFile.AltAA {
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:%v%v%v]", val.Name, val.Locus, val.AALref, CodonNbrInG, val.AALalt))
						// buffer.WriteString(fmt.Sprintf("%v_%v:%v%v%v\t", val.Name, val.Locus, val.AALref, CodonNbrInG, val.AALalt))

					}
				case tLCN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && CodonNbrInG == snpFromFile.CodonNbrInG {
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:codon%v]", val.Name, val.Locus, CodonNbrInG))
						// buffer.WriteString(fmt.Sprintf("%v_%v:codon%v\t", val.Name, val.Locus, CodonNbrInG))
					}
				case tPMN:
					// fmt.Println(lAPos == snpFromFile.APos)
					if lAPos == snpFromFile.APos && strings.ToUpper(val.Alt) == strings.ToUpper(snpFromFile.Alt) {
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:%v_%v>%v]", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))
						// chkSNP = checkSNP{FileName: file, FoundSNP: fmt.Sprintf("%v[%v:%v_%v>%v]\t", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt))}
						// snpArr = append(snpArr, checkSNP{FileName: file, FoundSNP: fmt.Sprintf("%v[%v:%v_%v>%v]\t", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt))})
						// buffer.WriteString(fmt.Sprintf("%v_%v:%v_%v>%v\t", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))

					}
				}
			}

			// fmt.Println(mapofSNP)

			// if len(mapofSNP) == 0 {
			// 	fmt.Println("Не было обнаружено ни одного совпадения в файле ", file, " между позициями в исследуемых файлах и списком СНИПов. Это может быть следствием того, что либо в ваших данных отсутствуют данные позиции, либо результатом использования различных референсных геномов при картировании и при анализе.  Обратите внимание на то, что имена локусов также могут быть различными и приведите ваш список сравнения соответственно используемому референсному геному. ")
			// 	// }
			// }
		}

	}

	if useRule == true {

		found := make(map[string][]string)

		for key, val := range mapofSNP {
			for _, tag := range rulesArr {

				for i := 0; i <= len(tag.Variants)-1; i++ {
					if strings.Contains(strings.ToUpper(strings.Join(val, ",")), tag.Variants[i]) == true {

						found[key] = appendIfMissing(found[key], strings.ToUpper(tag.Variants[i]))

					}
				}

			}

		}

		for key, val := range mapofSNP {
			chkSNP = append(chkSNP, checkSNP{FileName: key, FoundSNP: strings.Join(val, ","), ExactWithRule: strings.Join(found[key], ",")})
		}

	} else {
		for key, val := range mapofSNP {
			chkSNP = append(chkSNP, checkSNP{FileName: key, FoundSNP: strings.Join(val, ",")})
		}
	}

	if web == true {

		printSNPfromFile(chkSNP, *gbPort)
	} else {
		for _, key := range chkSNP {
			fmt.Printf("%v\t%v\t%v\n", key.FileName, key.FoundSNP, key.ExactWithRule)
		}
	}

}

func appendIfMissing(slice []string, val string) []string {
	sort.Slice(slice, func(i, j int) bool {
		return slice[i] < slice[j]
	})

	for _, ele := range slice {
		if ele == val {
			return slice
		}
	}

	return append(slice, val)
}

func printSNPfromFile(stat []checkSNP, port string) {

	var htmlTemplate = `   <!DOCTYPE html>
			<html>
			<head>
			<meta charset="utf-8">			
			</head>		
			<table width="100%" cellspacing="0" cellpadding="4" border="1">
			<body>
			<tr>		
		 
			{{range $element := .}}	
				{{ $length := len .ExactWithRule}}{{ if eq $length 0 }}
			   				<td>{{.FileName}}</td><td>{{.FoundSNP}}</td>		
					{{else}} 
							<td>{{.FileName}}</td><td>{{.FoundSNP}}</td><td>{{.ExactWithRule}}</td>
				{{end}}
				
			</tr>
			{{end}}	
			</body>
			</table>
		
			<table width="100%" cellspacing="0" cellpadding="4" border="0">
			<tr>
			<td><a href="http://bsatool.ru" target="_blank">Created by BSATool (Bacterial Snp Annotation Tool)</a></td>
			</tr>
			</table>			
			`

	t := template.New("t")
	// t1 := template.New("t1")
	t, err := t.Parse(htmlTemplate)

	// t1 := template.New("t1")
	// t2, err := t1.Parse(htmlTitle)

	if err != nil {
		panic(err)
	}

	browser.OpenURL(fmt.Sprintf("localhost:%v", port))

	http.HandleFunc("/", func(w http.ResponseWriter, r *http.Request) {
		// err = t.Execute(w, &gInfo)

		err = t.Execute(w, stat)
		if err != nil {
			panic(err)
		}
		go func() {
			defer os.Exit(0)
		}()
	})

	locPort := fmt.Sprintf(":%v", port)
	http.ListenAndServe(locPort, nil)
}

func checkTiTv(ref, alt string) string {
	var typeOf string

	if ref != "" && alt != "" {
		lref, lalt := strings.ToUpper(ref), strings.ToUpper(alt)
		if lref == "A" && lalt == "G" || lref == "G" && lalt == "A" || lref == "C" && lalt == "T" || lref == "T" && lalt == "C" {
			typeOf = "Ti"
		} else {
			typeOf = "Tv"
		}
	}
	return typeOf
}

func testGeneInfo(genes []geneInfo) {
	var (
		prod       string
		start, end int
		// gc, gc1, gc2, gc3 float64
	)
	for _, g := range genes {

		start = g.Start
		end = g.End
		// seq = getGeneSequence(g.Locus)
		prod = getProductByName(g.Locus)

		if start != 0 && end != 0 && g.TypeOf == "CDS" {
			// gc, gc1, gc2, gc3 = codon.GcCodonCalc(seq)
			fmt.Printf("%v:%v\t%v\t%v\n", start, end, g.Locus, prod)

		}
	}
}

func makeAltString(locus string, positions []allPositionsInGene) string {
	// var lStart, lEnd int
	var (
		seqSplit []string
		seq      string
	)
	// lStart, lEnd := getGenePosByName(locus)

	// seq = getNucFromGenome(start, end)

	seq = getGeneSequence(locus)

	for _, nuc := range seq {
		seqSplit = append(seqSplit, string(nuc))

	}

	// fmt.Println(locus, seqSplit)

	for _, val := range positions {
		// fmt.Println(val.pos, val.ref, ">", val.alt, seqSplit[val.pos-1], ">", val.alt)
		seqSplit[val.pos-1] = val.alt
		// fmt.Println(seqSplit[val.pos-1])
	}

	return strings.Join(seqSplit, "")

}

func getGenePosByName(locus string) (int, int) {
	var start, end int

	g := geneCoordinates[locus]

	if g.Start != 0 {
		start = g.Start
		end = g.End
	}

	// for _, g := range allGenesVal {

	// 	if strings.ToUpper(locus) == strings.ToUpper(g.Locus) {
	// 		start = g.Start
	// 		end = g.End
	// 		break
	// 	}
	// }

	return start, end
}

func getProductByName(locus string) string {
	var prod string

	for _, g := range allGenesVal {
		if strings.ToUpper(locus) == strings.ToUpper(g.Locus) {
			prod = g.Product
			break
		}
	}

	return prod
}

func getDirectionByName(locus string) string {
	var direction string

	for _, g := range allGenesVal {
		if strings.ToUpper(locus) == strings.ToUpper(g.Locus) {
			direction = g.Direction

			break
		}
	}

	return direction
}

func getProductByPos(start, end int) (string, string) {
	var prod, note string

	for _, g := range allGenesVal {

		if start >= g.Start && end <= g.End {
			// genesArr.WriteString(g.Locus)
			prod = g.Product
			note = g.Note

		}

	}

	return prod, note

}

func checkLocus(locus string) int {
	var locType int

	for _, g := range allGenesVal {
		if strings.ToUpper(locus) == strings.ToUpper(g.Locus) && g.TypeOf == "CDS" {
			locType = 1
			break
		} else if strings.ToUpper(locus) == strings.ToUpper(g.Locus) && g.TypeOf != "CDS" {
			locType = 2
			break

		} else if strings.ToUpper(locus) == strings.ToUpper(g.Locus) && g.TypeOf == "" {
			locType = 0
			break

		}
	}

	return locType
}

func toCircos(genes []geneInfo) {
	var (
		// seq        string
		start, end            int
		bandColor, genomeName string
		// gc         float64
		// buffer strings.Builder
	)
	// fmt.Println("---- ideogram.txt ----")
	if *statCircosBandColor != "" {
		bandColor = *statCircosBandColor
	}
	if *statCircosGenome != "" {
		genomeName = *statCircosGenome
	} else {
		genomeName = gInfo.Strain
	}

	fmt.Printf("chr -  %v caption %v %v %v\n", genomeName, gInfo.Start, gInfo.End, gInfo.Strain)

	for _, g := range genes {
		start = g.Start
		end = g.End
		// seq = getNucFromGenome(start-1, end)
		// fmt.Println(seq)
		if start != 0 && end != 0 {
			// gc, _, _, _ = codon.GcCodonCalc(seq)
			switch *statCircosTypeOf {
			case "":
				fmt.Printf("band %v %v %v %v %v %v\n", genomeName, g.Locus, g.Locus, start, end, bandColor)
			case "cds":
				if g.TypeOf == "CDS" {
					fmt.Printf("band %v %v %v %v %v %v\n", genomeName, g.Locus, g.Locus, start, end, bandColor)
				}
			case "igr":
				if g.TypeOf == "IGR" {
					fmt.Printf("band %v %v %v %v %v %v\n", genomeName, g.Locus, g.Locus, start, end, bandColor)
				}
			}

			// buffer.WriteString(fmt.Sprintf("%v %v %v %.2f\n", gInfo.Strain, start, end, gc))
		}
	}

	// fmt.Println("---- histogram.txt ----")
	// fmt.Println(buffer.String())

}

func annotatePositions(file string) {
	cols := regexp.MustCompile(`(\w+)\W+(\d+)\W+(\d+)`)
	f, err := os.Open(file)

	if err != nil {
		fmt.Println(err)

	}
	defer f.Close()
	fmt.Println("name\tstart\tend\tlocus\tlen\tgc\tseq\tproduct")
	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		for _, colsMatch := range cols.FindAllStringSubmatch(scanner.Text(), -1) {
			// fmt.Printf("c1:%v c2:%v\n", colsMatch[1], colsMatch[2])
			lStart, _ := strconv.Atoi(colsMatch[2])
			lEnd, _ := strconv.Atoi(colsMatch[3])
			gname, glen := getGeneNameByPos(lStart, lEnd)
			// prod := getProductByName(gname)

			prod, _ := getProductByPos(lStart, lEnd)
			seq := getNucFromGenome(lStart-1, lEnd)
			gc, _, _, _ := codon.GcCodonCalc(seq)
			fmt.Printf("%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", colsMatch[1], colsMatch[2], colsMatch[3], gname, glen, gc, seq, prod)

		}
	}

}

func getGeneNameByPos(start, end int) (string, int) {
	var (
		locName string
		locLen  int
	)
	locLen = (end - start) + 1

	for _, g := range allGenesVal {
		if start >= g.Start && end <= g.End {
			locName = g.Locus

		}

	}
	if locName == "" {
		locName = "IGR"
	}
	return locName, locLen
}

// func makeMatrixOld(typeof string, fileOut string) {

// 	var (
// 		AllPos     []int
// 		allLocuses []string
// 		buffer     strings.Builder
// 		headers    strings.Builder
// 		posCount   = make(map[int]int)
// 	)
// 	// var snps []snpInfo
// 	// t0 := time.Now()

// 	// var ResSeq []seqInfo

// 	// files := getListofVCF()
// 	files := &listOfFiles

// 	// if len(snpCacheMap) == 0 {
// 	// 	for _, file := range *files {
// 	// 		qSNP := &vcfQuery{File: file, OutChan: make(chan vcfInfoQuery)}
// 	// 		go qSNP.request()
// 	// 		// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
// 	// 		snpRes := <-qSNP.OutChan
// 	// 		snpCacheMap[snpRes.File] = snpRes.SnpInfo

// 	// 	}
// 	// }

// 	// fmt.Println(files)
// 	// for i, file := range *files {

// 	// 	fmt.Printf("Preparing files: Working on %v files from %v \r", i+1, len(*files))

// 	// 	snpsChan := make(chan []snpInfo)

// 	// 	go func() {
// 	// 		snpsChan <- makeSnps(file)
// 	// 	}()
// 	// 	snps = <-snpsChan

// 	// 	for _, val := range snps {
// 	// 		AllPosUnsort = append(AllPosUnsort, val.APos)
// 	// 		posCount[val.APos] = posCount[val.APos] + 1
// 	// 		if val.TypeOf == "CDS" {
// 	// 			allLocusUnsort = append(allLocusUnsort, val.Locus)
// 	// 		}

// 	// 	}
// 	// }

// 	// AllPos = unique(AllPosUnsort)
// 	// allLocuses = removeStringDuplicates(allLocusUnsort)

// 	// sort.Ints(AllPos)

// 	// fmt.Println(snpCacheMap)

// 	switch typeof {

// 	case "binary":

// 		matrixBinary(fmt.Sprintf("%v_1", fileOut))
// 		// pos := make(map[int]string)

// 		// var posFN = make(map[int][]string)
// 		// var posFreq = map[int][]string{}

// 		// headers.WriteString("Pos\t")

// 		// for i, file := range *files {
// 		// 	headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

// 		// 	fmt.Printf("Counting SNP positions: Working on %v files from %v \r", i+1, len(*files))

// 		// 	snpsChan := make(chan []snpInfo)

// 		// 	go func() {
// 		// 		snpsChan <- makeSnps(file)
// 		// 	}()
// 		// 	snps = <-snpsChan

// 		// 	for _, val := range snps {

// 		// 		pos[val.APos] = val.Alt
// 		// 		// posTest = append(posTest, posByFile{pos: pos, file: file, apos: val.APos})
// 		// 		if strings.Contains(strings.Join(posFN[val.APos], " "), file) == false {
// 		// 			posFN[val.APos] = append(posFN[val.APos], file)
// 		// 		}

// 		// 		AllPosUnsort = append(AllPosUnsort, val.APos)
// 		// 		posCount[val.APos] = posCount[val.APos] + 1
// 		// 		if val.TypeOf == "CDS" {
// 		// 			allLocusUnsort = append(allLocusUnsort, val.Locus)
// 		// 		}

// 		// 	}

// 		// }
// 		// AllPos = unique(AllPosUnsort)
// 		// allLocuses = removeStringDuplicates(allLocusUnsort)
// 		// sort.Ints(AllPos)

// 		// for _, file := range *files {
// 		// 	for _, allpos := range AllPos {

// 		// 		if strings.Contains(strings.Join(posFN[allpos], " "), file) {

// 		// 			posFreq[allpos] = append(posFreq[allpos], "1")
// 		// 		} else {
// 		// 			// fmt.Println(allpos, 0, file)
// 		// 			posFreq[allpos] = append(posFreq[allpos], "0")
// 		// 		}

// 		// 	}
// 		// }

// 		// for _, allpos := range AllPos {
// 		// 	if posCount[allpos] < len(*files) {

// 		// 		buffer.WriteString(fmt.Sprintln(allpos, "\t", strings.Join(posFreq[allpos], "\t")))

// 		// 	}
// 		// }
// 		// headers.WriteString("\n")

// 	// case "binary":
// 	// 	var posFreq = map[int][]string{}
// 	// 	pos := make(map[int]string)
// 	// 	// posarr := make(map[string][]map[int]string)
// 	// 	headers.WriteString("Pos\t")

// 	// 	// for i, file := range files {
// 	// 	// for fname, snps := range snpCacheMap {
// 	// 	// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))
// 	// 	for i, file := range *files {
// 	// 		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

// 	// 		fmt.Printf("Preparing files: Working on %v files from %v \r", i+1, len(*files))

// 	// 		snpsChan := make(chan []snpInfo)

// 	// 		go func() {
// 	// 			snpsChan <- makeSnps(file)
// 	// 		}()
// 	// 		snps = <-snpsChan

// 	// 		for _, val := range snps {

// 	// 			pos[val.APos] = val.Alt
// 	// 			// posarr[file] = append(posarr[file], pos)
// 	// 			AllPosUnsort = append(AllPosUnsort, val.APos)
// 	// 			posCount[val.APos] = posCount[val.APos] + 1
// 	// 			if val.TypeOf == "CDS" {
// 	// 				allLocusUnsort = append(allLocusUnsort, val.Locus)
// 	// 			}

// 	// 		}

// 	// 	}
// 	// 	// fmt.Println(posarr)
// 	// 	AllPos = unique(AllPosUnsort)
// 	// 	allLocuses = removeStringDuplicates(allLocusUnsort)
// 	// 	sort.Ints(AllPos)

// 	// 	for _, file := range *files {
// 	// 		for _, allpos := range AllPos {
// 	// 			if posCount[allpos] < len(*files) {

// 	// 				if pos[allpos] != "" {

// 	// 					posFreq[allpos] = append(posFreq[allpos], "1")

// 	// 				} else {

// 	// 					posFreq[allpos] = append(posFreq[allpos], "0")

// 	// 				}
// 	// 			}

// 	// 		}
// 	// 		fmt.Println(file)
// 	// 	}
// 	// 	for _, allpos := range AllPos {
// 	// 		if posCount[allpos] < len(*files) {

// 	// 			buffer.WriteString(fmt.Sprintln(allpos, "\t", strings.Join(posFreq[allpos], "\t")))

// 	// 		}
// 	// 	}
// 	// 	headers.WriteString("\n")

// 	case "table":
// 		var posFreq = map[int][]string{}
// 		pos := make(map[int]string)

// 		// headers.WriteString("Pos\t")

// 		for fname, snps := range snpCacheMap {
// 			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

// 			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

// 			// snps := parserVCF(file, false, allGenesVal)

// 			for _, val := range snps {

// 				pos[val.APos] = val.Alt

// 			}

// 			for _, allpos := range AllPos {
// 				if posCount[allpos] < len(*files) {
// 					if pos[allpos] != "" {
// 						posFreq[allpos] = append(posFreq[allpos], strconv.Itoa(allpos))

// 					} else {
// 						posFreq[allpos] = append(posFreq[allpos], "0")

// 					}
// 				}

// 			}
// 		}
// 		for _, allpos := range AllPos {
// 			if posCount[allpos] < len(*files) {

// 				buffer.WriteString(fmt.Sprintln(strings.Join(posFreq[allpos], "\t")))

// 			}
// 		}
// 		headers.WriteString("\n")

// 	case "nc":
// 		var posFreq = map[int][]string{}
// 		pos := make(map[int]string)

// 		headers.WriteString("Pos\tRef\t")

// 		// for i, file := range files {
// 		for fname, snps := range snpCacheMap {
// 			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

// 			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

// 			// snps := parserVCF(file, false, allGenesVal)

// 			for _, val := range snps {

// 				pos[val.APos] = val.Alt

// 			}

// 			for _, allpos := range AllPos {

// 				if pos[allpos] != "" {
// 					posFreq[allpos] = append(posFreq[allpos], fmt.Sprintf("%v", pos[allpos]))

// 				} else {
// 					posFreq[allpos] = append(posFreq[allpos], getNucFromGenomePos(allpos))

// 				}

// 			}
// 		}
// 		for _, allpos := range AllPos {
// 			// if buffer.Len() == 0 {
// 			buffer.WriteString(fmt.Sprintln(allpos, "\t", getNucFromGenomePos(allpos), "\t", strings.Join(posFreq[allpos], "\t")))

// 		}
// 		headers.WriteString("\n")

// 	case "locus":
// 		var (
// 			i          int
// 			locFreq    = map[string][]string{}
// 			locusCount = make(map[string]int)
// 		)
// 		headers.WriteString("Locus\t")

// 		// for i, file := range files {
// 		for fname, snps := range snpCacheMap {
// 			i++
// 			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

// 			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

// 			// snps := parserVCF(file, false, allGenesVal)

// 			for _, val := range snps {

// 				locusCount[fmt.Sprintf("%v_%v", val.Locus, i)] = locusCount[fmt.Sprintf("%v_%v", val.Locus, i)] + 1

// 			}
// 			for _, allloc := range allLocuses {

// 				locFreq[allloc] = append(locFreq[allloc], strconv.Itoa(locusCount[fmt.Sprintf("%v_%v", allloc, i)]))

// 			}

// 		}

// 		for _, allloc := range allLocuses {

// 			buffer.WriteString(fmt.Sprintln(allloc, "\t", strings.Join(locFreq[allloc], "\t")))

// 		}

// 		headers.WriteString("\n")

// 	case "freq":
// 		headers.WriteString("Pos\tFreq\n")

// 		for _, allpos := range AllPos {

// 			// headers.WriteString(fmt.Sprintf("P%v\n", allpos))
// 			buffer.WriteString(fmt.Sprintf("P%v\t%v\n", allpos, posCount[allpos]))

// 		}

// 	// case "dnds":
// 	// 	// var dnds [][]DnDsRes
// 	// 	var locDNDS = map[string][]string{}
// 	// 	// var locDNDSTest = map[string][]string{}
// 	// 	// var dndsArr []string
// 	// 	var i int

// 	// 	prompt := bufio.NewReader(os.Stdin)
// 	// 	fmt.Printf("It will take some time(~%v min). Continue?: ", len(*files))
// 	// 	yesNo, _ := prompt.ReadString('\n')

// 	// 	if yesNo == "y\n" || yesNo == "Y\n" {
// 	// 		headers.WriteString("Locus\t")

// 	// 		// for i, file := range files {
// 	// 		for fname := range snpCacheMap {
// 	// 			i++
// 	// 			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

// 	// 			fmt.Printf("Calculating Dn/DS: Working on %v from %v (%v) \r", i, len(*files), fname)

// 	// 			dnds := calcDnDsVal(fname, *gbVerbose)

// 	// 			// fmt.Println(fname, "\t", dnds)
// 	// 			for _, dndsVal := range dnds {

// 	// 				locDNDS[strings.Split(dndsVal.Locus, ":")[0]] = append(locDNDS[strings.Split(dndsVal.Locus, ":")[0]], fmt.Sprintf("%.2f", dndsVal.DNDS))

// 	// 			}

// 	// 			// for _, allloc := range allLocuses {
// 	// 			// 	// dndsArr = locDNDSTest[allloc]
// 	// 			// 	// found := 0
// 	// 			// 	for _, dndsVal := range dnds {

// 	// 			// 		if strings.Split(dndsVal.Locus, ":")[0] == allloc {
// 	// 			// 			dndsArr = append(dndsArr, fmt.Sprintf("%.2f", dndsVal.DNDS))
// 	// 			// 			// found = 1
// 	// 			// 			// } else {
// 	// 			// 			// 	dndsArr = append(dndsArr, "1")
// 	// 			// 		}

// 	// 			// 	}
// 	// 			// 	locDNDSTest[allloc] = dndsArr
// 	// 			// 	// dndsArr = nil

// 	// 			// }

// 	// 			// j++

// 	// 		}

// 	// 		// fmt.Println(locDNDSTest)
// 	// 		for _, allloc := range allLocuses {

// 	// 			// fmt.Println(allloc, locDNDS[allloc], len(locDNDS[allloc]))
// 	// 			if len(locDNDS[allloc]) != 0 {
// 	// 				rpt := strings.Repeat("1\t", len(*files)-len(locDNDS[allloc]))
// 	// 				locDNDS[allloc] = append(locDNDS[allloc], rpt)
// 	// 				buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))

// 	// 			}

// 	// 		}

// 	// 		headers.WriteString("\n")

// 	// 	}
// 	case "dnds":
// 		matrixDnDs(fileOut)
// 	// case "dnds":

// 	// 	var altPositions = make(map[string][]allPositionsInGene)
// 	// 	var locDNDS = map[string][]string{}
// 	// 	var countNbrOne = make(map[string]int)

// 	// 	i := 0
// 	// 	headers.WriteString("Locus\t")

// 	// 	for key := range geneCoordinates {
// 	// 		allLocuses = append(allLocuses, key)
// 	// 	}

// 	// 	for _, fname := range *files {
// 	// 		i++

// 	// 		snpsChan := make(chan []snpInfo)

// 	// 		go func() {
// 	// 			snpsChan <- makeSnps(fname)
// 	// 		}()
// 	// 		snps = <-snpsChan

// 	// 		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

// 	// 		t1 := time.Now()

// 	// 		fmt.Printf("Calculating Dn/DS: Working on %v from %v (%v)\r", i, len(*files), fname)
// 	// 		for _, val := range snps {
// 	// 			if val.TypeOf == "CDS" && containsPos(altPositions[val.Locus], val.PosInGene, val.Alt) == false {

// 	// 				altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, alt: val.Alt, ref: val.NucInPos, locus: val.Locus})
// 	// 			}

// 	// 		}

// 	// 		for _, allloc := range allLocuses {

// 	// 			// prod := getProductByName(allloc)
// 	// 			if len(altPositions[allloc]) > 2 {
// 	// 				dndsChan := make(chan []string)
// 	// 				go func() {

// 	// 					dndsChan <- getDnDsByLocus(allloc, altPositions[allloc])
// 	// 				}()
// 	// 				dndsRes, ok := <-dndsChan
// 	// 				if ok {
// 	// 					// fmt.Println(dndsRes)
// 	// 					locDNDS[dndsRes[0]] = append(locDNDS[dndsRes[0]], dndsRes[1])
// 	// 					if dndsRes[1] == "1" {
// 	// 						countNbrOne[allloc]++
// 	// 					}
// 	// 					close(dndsChan)
// 	// 				}

// 	// 			} else {

// 	// 				locDNDS[allloc] = append(locDNDS[allloc], "1")
// 	// 				countNbrOne[allloc]++
// 	// 			}

// 	// 		}
// 	// 		fmt.Printf("Calculating Dn/DS: Working on %v from %v (%v) \t\t E.t:%v\n", i, len(*files), fname, t1.Sub(t0))

// 	// 	}

// 	// 	for _, allloc := range allLocuses {

// 	// 		if countNbrOne[allloc] != len(*files) && *statAll == false {
// 	// 			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
// 	// 		} else if *statAll == true {
// 	// 			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
// 	// 		}

// 	// 	}
// 	// 	// }

// 	// 	headers.WriteString("\n")

// 	case "mst":
// 		var (
// 			exGenes = make(map[int]int)
// 			exSNP   = make(map[int]int)
// 		)

// 		exGenes = loadExcludeGenes(*gbExcludeGenes)
// 		exSNP = loadExcludeSNP(*gbExcludeSnp)

// 		matrixBinaryMST(fileOut, false, exGenes, exSNP)
// 	case "jw":
// 		// var dnds [][]DnDsRes
// 		var (
// 			locJW = map[string][]string{}
// 			i     int
// 		)
// 		headers.WriteString("Locus\t")

// 		// for i, file := range files {
// 		for fname := range snpCacheMap {
// 			i++
// 			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

// 			fmt.Printf("Calculating Jaro Winkler distance: Working on %v from %v (%v) \r", i, len(*files), fname)

// 			jw := calcJaroWinklerDist(fname, *gbVerbose)

// 			for _, jwVal := range jw {

// 				locJW[jwVal.Locus] = append(locJW[jwVal.Locus], fmt.Sprintf("%.3f", jwVal.JWDist))

// 			}

// 		}

// 		// for _, allloc := range allLocuses {

// 		// 	if len(locJW[allloc]) != 0 {
// 		// 		for i, res := range locJW[allloc] {
// 		// 			fmt.Println(i, res)
// 		// 		}
// 		// 	}
// 		// }

// 		for _, allloc := range allLocuses {

// 			if len(locJW[allloc]) != 0 {

// 				// fmt.Println(len(*files)-len(locJW[allloc]), locJW[allloc], strings.Repeat(" 1 ", len(*files)-len(locJW[allloc])))
// 				rpt := strings.Repeat("1\t", len(*files)-len(locJW[allloc]))
// 				locJW[allloc] = append(locJW[allloc], rpt)
// 				buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v\t", allloc), strings.Join(locJW[allloc], "\t")))

// 			}
// 		}

// 		headers.WriteString("\n")

// 	case "gc3":
// 		// calcGC3Val
// 		var (
// 			locGC3 = map[string][]string{}
// 			refGC3 = map[string]string{}
// 			i      int
// 		)
// 		headers.WriteString("Locus\tRefCG3\t")

// 		// for i, file := range files {
// 		for fname := range snpCacheMap {
// 			// fmt.Println(fname)
// 			i++
// 			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

// 			fmt.Printf("Calculating GC3 values: Working on %v from %v\t%v \r", i, len(*files), fname)
// 			// logger.Printf("%v\n", fname)
// 			// fmt.Println(snpCacheMap[fname])
// 			gc3 := calcGC3Val(snpCacheMap[fname])
// 			// logger.Printf("%v\t%v\n", fname, gc3)
// 			// fmt.Println(gc3)
// 			// fmt.Println(gc3)
// 			for _, gc3val := range gc3 {

// 				locGC3[gc3val.Locus] = append(locGC3[gc3val.Locus], fmt.Sprintf("%.2f", gc3val.GC3Alt))
// 				refGC3[gc3val.Locus] = fmt.Sprintf("%.2f", gc3val.GC3Ref)
// 				// if *flgDebug == true {
// 				// 	fmt.Printf("L:%v Alt:%v Ref:%v\n", gc3val.Locus, gc3val.GC3Alt, gc3val.GC3Ref)
// 				// }

// 			}

// 		}

// 		for _, allloc := range allLocuses {
// 			// if len(locGC3[allloc]) != 0 {
// 			if len(locGC3[allloc]) != 0 {
// 				rpt := strings.Repeat(refGC3[allloc]+"\t", len(*files)-len(locGC3[allloc]))
// 				locGC3[allloc] = append(locGC3[allloc], rpt)
// 				buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v\t%v\t", allloc, refGC3[allloc]), strings.Join(locGC3[allloc], "\t")))

// 			}
// 		}

// 		headers.WriteString("\n")

// 	}

// 	// case "gc":
// 	// 	// calcGC3Val
// 	// 	var locGC = map[string][]string{}
// 	// 	var refGC = map[string]string{}

// 	// 	headers.WriteString("Locus\tRefCG\t")

// 	// 	for i, file := range files {

// 	// 		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

// 	// 		fmt.Printf("Calculating GC values: Working on %v from %v \r", i+1, len(files))

// 	// 		gc := calcGC3Val(file)

// 	// 		for _, gcval := range gc {

// 	// 			locGC[gcval.Locus] = append(locGC[gcval.Locus], fmt.Sprintf("%.2f", gcval.GCalt))
// 	// 			refGC[gcval.Locus] = fmt.Sprintf("%.2f", gcval.GCref)

// 	// 		}

// 	// 	}

// 	// 	for _, allloc := range allLocuses {

// 	// 		if len(locGC[allloc]) != 0 {
// 	// 			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v\t%v\t", allloc, refGC[allloc]), strings.Join(locGC[allloc], "\t")))

// 	// 		}
// 	// 	}

// 	// 	headers.WriteString("\n")

// 	// }

// 	// if buffer.Len() != 0 && headers.Len() != 0 {
// 	// 	fOut, err := os.Create(fileOut)
// 	// 	if err != nil {
// 	// 		log.Fatal("Cannot create file", err)
// 	// 	}
// 	// 	defer fOut.Close()
// 	// 	fmt.Fprintf(fOut, headers.String())
// 	// 	fmt.Fprintf(fOut, buffer.String())
// 	// 	fmt.Printf("\n\nWell done!\n")
// 	// 	// t1 := time.Now()
// 	// 	// fmt.Printf("Elapsed time: %v", fmtDuration(t1.Sub(t0)))
// 	// }

// }

func makeMatrix(typeof string, fileOut string, verbose bool) {

	var (
		AllPos     []int
		allLocuses []string
		buffer     strings.Builder
		headers    strings.Builder
		posCount   = make(map[int]int)
	)

	files := &listOfFiles

	for i, file := range *files {
		// создаем запрос в виде типа vcfQuery, который передается через канал на выход <-qSNP.OutChan
		qSNP := &vcfQuery{File: file, OutChan: make(chan vcfInfoQuery)}
		go qSNP.request()
		// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
		snpRes := <-qSNP.OutChan
		snpCacheMap[snpRes.File] = snpRes.SnpInfo
		if verbose == true {
			fmt.Printf("Reading files:  %v from %v \r", i+1, len(*files))
		}

	}

	switch typeof {

	case "binary":

		matrixBinary(fmt.Sprintf("%v", fileOut))

	case "table":

		matrixTable(fileOut)

	case "nc":
		var posFreq = map[int][]string{}
		pos := make(map[int]string)

		headers.WriteString("Pos\tRef\t")
		// allPosChan := make(chan []int)
		// go func() {
		// 	allPosChan <- getAllPosFromCacheMap()
		// }()
		// AllPos = <-allPosChan
		// for i, file := range files {
		AllPos = getAllPosFromCacheMap()
		for fname, snps := range snpCacheMap {
			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

			// snps := parserVCF(file, false, allGenesVal)

			for _, val := range snps {

				pos[val.APos] = val.Alt
				// AllPosUnsort = append(AllPosUnsort, val.APos)

			}
			for _, allpos := range AllPos {

				if pos[allpos] != "" {
					posFreq[allpos] = append(posFreq[allpos], fmt.Sprintf("%v", pos[allpos]))
					// fmt.Println(posFreq[allpos])

				} else {
					posFreq[allpos] = append(posFreq[allpos], getNucFromGenomePos(allpos))

				}

			}

		}
		// AllPos = unique(AllPosUnsort)
		// // allLocuses = removeStringDuplicates(allLocusUnsort)
		// sort.Ints(AllPos)

		for _, allpos := range AllPos {
			// if buffer.Len() == 0 {
			buffer.WriteString(fmt.Sprintln(allpos, "\t", getNucFromGenomePos(allpos), "\t", strings.Join(posFreq[allpos], "\t")))
			// fmt.Println(buffer.String())

		}
		headers.WriteString("\n")
		matrixPrint(headers, buffer, fileOut)
		// if buffer.Len() != 0 && headers.Len() != 0 {
		// 	fOut, err := os.Create(fileOut)
		// 	if err != nil {
		// 		log.Fatal("Cannot create file", err)
		// 	}
		// 	defer fOut.Close()
		// 	fmt.Fprintf(fOut, headers.String())
		// 	fmt.Fprintf(fOut, buffer.String())
		// 	fmt.Printf("\n\nWell done!\n")
		// t1 := time.Now()
		// fmt.Printf("Elapsed time: %v", fmtDuration(t1.Sub(t0)))
		// }
	case "summary":
		var (
			posFreq                = map[int][]string{}
			pos                    = make(map[int]string)
			groupRegexp            = regexp.MustCompile(`^(\S+)\W+(\w+)\W+(\w+)`)
			group                  = make(map[string]string)
			label                  = make(map[string]string)
			groupHeader, groupBody string
			fileGroup              = make(map[int][]string)
			fileLabel              = make(map[int][]string)
		)

		if *statGroupFromFile != "" {
			f, err := os.Open(*statGroupFromFile) // открываем файл

			if err != nil {
				fmt.Println(err)

			}
			defer f.Close()

			scanner := bufio.NewScanner(f) //  новый сканер

			for scanner.Scan() {

				scanTxt := scanner.Text()

				for _, grpVal := range groupRegexp.FindAllStringSubmatch(scanTxt, -1) {
					group[grpVal[1]] = grpVal[2]
					label[grpVal[1]] = grpVal[3]

					// fmt.Println(grpVal)
				}
			}

		}
		// fmt.Println(group)
		if len(group) != 0 {
			groupHeader = "Group\t"
			if len(label) != 0 {
				groupHeader = fmt.Sprintf("%vLabel\t", groupHeader)
			}
		}

		headers.WriteString(fmt.Sprintf("Pos\tLocus\t%v", groupHeader))
		// allPosChan := make(chan []int)
		// go func() {
		// 	allPosChan <- getAllPosFromCacheMap()
		// }()
		// AllPos = <-allPosChan
		// for i, file := range files {
		AllPos = getAllPosFromCacheMap()
		for fname, snps := range snpCacheMap {

			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

			// snps := parserVCF(file, false, allGenesVal)

			for _, val := range snps {

				pos[val.APos] = val.Alt
				// AllPosUnsort = append(AllPosUnsort, val.APos)

			}
			for _, allpos := range AllPos {

				if pos[allpos] != "" {
					posFreq[allpos] = append(posFreq[allpos], fmt.Sprintf("%v/%v", getNucFromGenomePos(allpos), pos[allpos]))
					// fmt.Println(posFreq[allpos])
					if len(group) != 0 {
						fileGroup[allpos] = append(fileGroup[allpos], group[fname])
						if len(label) != 0 {
							fileLabel[allpos] = append(fileLabel[allpos], label[fname])
						}
						// groupArr = append(groupArr, group[fname])
						// groupBody = fmt.Sprintf("\t%v\t", group[fname])
					} else {
						groupBody = "\t"
					}

				} else {
					posFreq[allpos] = append(posFreq[allpos], ".")

				}

			}

		}
		// AllPos = unique(AllPosUnsort)
		// // allLocuses = removeStringDuplicates(allLocusUnsort)
		// sort.Ints(AllPos)

		for _, allpos := range AllPos {
			// if buffer.Len() == 0 {
			gname, _ := getGeneNameByPos(allpos, allpos)
			if len(group) != 0 {
				sort.Strings(fileGroup[allpos])
				uniq := removeStringDuplicates(fileGroup[allpos])
				groupBody = fmt.Sprintf("\t%v\t", strings.Join(uniq, ""))
				// fmt.Println(groupBody)
			}
			if len(label) != 0 {
				sort.Strings(fileLabel[allpos])
				uniq := removeStringDuplicates(fileLabel[allpos])
				groupBody = fmt.Sprintf("%v%v\t", groupBody, strings.Join(uniq, " "))
			}

			buffer.WriteString(fmt.Sprintln(allpos, "\t", gname, groupBody, strings.Join(posFreq[allpos], "\t")))
			// fmt.Println(buffer.String())

		}
		headers.WriteString("\n")
		matrixPrint(headers, buffer, fileOut)
		// if buffer.Len() != 0 && headers.Len() != 0 {
		// 	fOut, err := os.Create(fileOut)
		// 	if err != nil {
		// 		log.Fatal("Cannot create file", err)
		// 	}
		// 	defer fOut.Close()
		// 	fmt.Fprintf(fOut, headers.String())
		// 	fmt.Fprintf(fOut, buffer.String())
		// 	fmt.Printf("\n\nWell done!\n")
		// t1 := time.Now()
		// fmt.Printf("Elapsed time: %v", fmtDuration(t1.Sub(t0)))
		// }

	case "locus":
		var (
			i          int
			locFreq    = map[string][]string{}
			locusCount = make(map[string]int)
		)
		headers.WriteString("Locus\t")

		// for i, file := range files {
		for fname, snps := range snpCacheMap {
			i++
			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

			// snps := parserVCF(file, false, allGenesVal)

			for _, val := range snps {

				locusCount[fmt.Sprintf("%v_%v", val.Locus, i)] = locusCount[fmt.Sprintf("%v_%v", val.Locus, i)] + 1

			}
			for _, allloc := range allLocuses {

				locFreq[allloc] = append(locFreq[allloc], strconv.Itoa(locusCount[fmt.Sprintf("%v_%v", allloc, i)]))

			}

		}

		for _, allloc := range allLocuses {

			buffer.WriteString(fmt.Sprintln(allloc, "\t", strings.Join(locFreq[allloc], "\t")))

		}

		headers.WriteString("\n")
		matrixPrint(headers, buffer, fileOut)

	case "freq":
		headers.WriteString("Pos\tFreq\n")
		AllPos = getAllPosFromCacheMap()
		for _, snps := range snpCacheMap {
			for _, val := range snps {

				posCount[val.APos] = posCount[val.APos] + 1

			}

		}

		for _, allpos := range AllPos {

			// headers.WriteString(fmt.Sprintf("P%v\n", allpos))
			buffer.WriteString(fmt.Sprintf("%v\t%v\n", allpos, posCount[allpos]))
			// fmt.Println(posCount[allpos])
		}

		matrixPrint(headers, buffer, fileOut)

	case "dnds":
		matrixDnDs(fileOut)
	case "mst":
		//stat --db test_core -a matrix -t mst -o test.tsv --exclude-genes=exgenes.txt --exclude-snp=drugs2.txt --snp-number=15
		var (
			exGenes = make(map[int]int)
			exSNP   = make(map[int]int)
		)
		if *gbExcludeGenes != "" {
			exGenes = loadExcludeGenes(*gbExcludeGenes)
		} else if *gbExcludeSnp != "" {
			exSNP = loadExcludeSNP(*gbExcludeSnp)
		}

		matrixBinaryMST(fileOut, *gbVerbose, exGenes, exSNP, *statNonRandomize)

	case "jw":
		// var dnds [][]DnDsRes
		var (
			locJW = map[string][]string{}
			i     int
		)
		headers.WriteString("Locus\t")

		// for i, file := range files {
		for fname := range snpCacheMap {
			i++
			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

			fmt.Printf("Calculating Jaro Winkler distance: Working on %v from %v (%v) \r", i, len(*files), fname)

			jw := calcJaroWinklerDist(fname, *gbVerbose)

			for _, jwVal := range jw {

				locJW[jwVal.Locus] = append(locJW[jwVal.Locus], fmt.Sprintf("%.3f", jwVal.JWDist))

			}

		}

		for _, allloc := range allLocuses {

			if len(locJW[allloc]) != 0 {

				// fmt.Println(len(*files)-len(locJW[allloc]), locJW[allloc], strings.Repeat(" 1 ", len(*files)-len(locJW[allloc])))
				rpt := strings.Repeat("1\t", len(*files)-len(locJW[allloc]))
				locJW[allloc] = append(locJW[allloc], rpt)
				buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v\t", allloc), strings.Join(locJW[allloc], "\t")))

			}
		}

		headers.WriteString("\n")

	case "gc3":
		// calcGC3Val
		var (
			locGC3 = map[string][]string{}
			refGC3 = map[string]string{}
			i      int
		)
		headers.WriteString("Locus\tRefCG3\t")

		// for i, file := range files {
		for fname := range snpCacheMap {
			// fmt.Println(fname)
			i++
			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

			fmt.Printf("Calculating GC3 values: Working on %v from %v\t%v \r", i, len(*files), fname)
			// logger.Printf("%v\n", fname)
			// fmt.Println(snpCacheMap[fname])
			gc3 := calcGC3Val(snpCacheMap[fname])
			// logger.Printf("%v\t%v\n", fname, gc3)
			// fmt.Println(gc3)
			// fmt.Println(gc3)
			for _, gc3val := range gc3 {

				locGC3[gc3val.Locus] = append(locGC3[gc3val.Locus], fmt.Sprintf("%.2f", gc3val.GC3Alt))
				refGC3[gc3val.Locus] = fmt.Sprintf("%.2f", gc3val.GC3Ref)
				// if *flgDebug == true {
				// 	fmt.Printf("L:%v Alt:%v Ref:%v\n", gc3val.Locus, gc3val.GC3Alt, gc3val.GC3Ref)
				// }

			}

		}

		for _, allloc := range allLocuses {
			// if len(locGC3[allloc]) != 0 {
			if len(locGC3[allloc]) != 0 {
				rpt := strings.Repeat(refGC3[allloc]+"\t", len(*files)-len(locGC3[allloc]))
				locGC3[allloc] = append(locGC3[allloc], rpt)
				buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v\t%v\t", allloc, refGC3[allloc]), strings.Join(locGC3[allloc], "\t")))

			}
		}

		headers.WriteString("\n")

	}

}

func getAllPosFromCacheMap() (allPos []int) {
	var (
		AllPosUnsort []int
	)
	if len(snpCacheMap) != 0 {
		for _, snps := range snpCacheMap {

			for _, val := range snps {

				AllPosUnsort = append(AllPosUnsort, val.APos)

			}

		}
		allPos = unique(AllPosUnsort)
		// allLocuses = removeStringDuplicates(allLocusUnsort)
		sort.Ints(allPos)

	}
	return allPos
}

func matrixTable(fileOut string) {
	var (
		AllPosUnsort, AllPos []int
		allLocusUnsort       []string
		buffer               strings.Builder
		headers              strings.Builder
		posCount             = make(map[int]int)
		// snps                 []snpInfo
		posFN   = make(map[int][]string)
		posFreq = map[int][]string{}
		i       = 0
	// var ResSeq []seqInfo
	)
	// files := getListofVCF()
	files := &listOfFiles

	// fmt.Println(files)
	pos := make(map[int]string)

	headers.WriteString("Pos\t")

	for key, snps := range snpCacheMap {
		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(key, filepath.Ext(key))))

		fmt.Printf("Counting SNP positions: Working on %v files from %v \r", i+1, len(*files))

		snpsChan := make(chan []snpInfo)

		go func() {
			snpsChan <- makeSnps(key)
		}()
		snps = <-snpsChan

		for _, val := range snps {

			pos[val.APos] = val.Alt
			// posTest = append(posTest, posByFile{pos: pos, file: file, apos: val.APos})
			if strings.Contains(strings.Join(posFN[val.APos], " "), key) == false {
				posFN[val.APos] = append(posFN[val.APos], key)
			}

			AllPosUnsort = append(AllPosUnsort, val.APos)
			posCount[val.APos] = posCount[val.APos] + 1
			if val.TypeOf == "CDS" {
				allLocusUnsort = append(allLocusUnsort, val.Locus)
			}

		}

	}
	AllPos = unique(AllPosUnsort)
	// allLocuses = removeStringDuplicates(allLocusUnsort)
	sort.Ints(AllPos)

	for _, file := range *files {
		for _, allpos := range AllPos {

			if strings.Contains(strings.Join(posFN[allpos], " "), file) {

				posFreq[allpos] = append(posFreq[allpos], "1")
			} else {
				// fmt.Println(allpos, 0, file)
				posFreq[allpos] = append(posFreq[allpos], "0")
			}

		}
	}

	for _, allpos := range AllPos {
		if posCount[allpos] < len(*files) {

			buffer.WriteString(fmt.Sprintln(allpos, "\t", strings.Join(posFreq[allpos], "\t")))

		}
	}
	headers.WriteString("\n")

	// if buffer.Len() != 0 && headers.Len() != 0 {
	// 	fOut, err := os.Create(fileOut)
	// 	if err != nil {
	// 		log.Fatal("Cannot create file", err)
	// 	}
	// 	defer fOut.Close()
	// 	fmt.Fprintf(fOut, headers.String())
	// 	fmt.Fprintf(fOut, buffer.String())
	// 	fmt.Printf("\n\nWell done!\n")
	// t1 := time.Now()
	// fmt.Printf("Elapsed time: %v", fmtDuration(t1.Sub(t0)))
	// }
	matrixPrint(headers, buffer, fileOut)
}

func calcDnDs(vcfFile string) {
	// var AllPos []int
	var (
		allLocuses []string

		snps         []snpInfo
		altPositions = make(map[string][]allPositionsInGene)
	)

	for key, val := range geneCoordinates {
		if val.Type == "CDS" {
			allLocuses = append(allLocuses, key)
		}
	}

	sort.Strings(allLocuses)

	snpsChan := make(chan []snpInfo)

	go func() {
		snpsChan <- makeSnps(vcfFile)
	}()
	snps = <-snpsChan

	// fmt.Printf("Calculating Dn/DS: Working on %v from %v (%v)\r", i, len(*files), fname)
	for _, val := range snps {
		if val.TypeOf == "CDS" && containsPos(altPositions[val.Locus], val.PosInGene, val.Alt) == false {

			altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, alt: val.Alt, ref: val.NucInPos, locus: val.Locus})
		}

	}

	for _, allloc := range allLocuses {

		// prod := getProductByName(allloc)
		if len(altPositions[allloc]) > 2 {
			dndsChan := make(chan []string)
			go func() {

				dndsChan <- getDnDsByLocus(allloc, altPositions[allloc])
			}()
			dndsRes, ok := <-dndsChan
			if ok {

				fmt.Println(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), dndsRes[1])
				close(dndsChan)
			}

		} else {

		}

	}

}

func matrixPrint(headers strings.Builder, buffer strings.Builder, fileOut string) {
	if buffer.Len() != 0 && headers.Len() != 0 {
		fOut, err := os.Create(fileOut)
		if err != nil {
			log.Fatal("Cannot create file", err)
		}
		defer fOut.Close()
		fmt.Fprintf(fOut, headers.String())
		fmt.Fprintf(fOut, buffer.String())
		fmt.Printf("\n\nWell done!\n")
	}
}

func matrixDnDsOld(fileOut string) {
	// var AllPos []int
	var (
		allLocuses []string

		buffer  strings.Builder
		headers strings.Builder
		// var posCount = make(map[int]int)
		snps         []snpInfo
		altPositions = make(map[string][]allPositionsInGene)
		locDNDS      = map[string][]string{}
		countNbrOne  = make(map[string]int)
	)
	t0 := time.Now()

	files := &listOfFiles

	i := 0
	headers.WriteString("Locus\t")

	for key := range geneCoordinates {
		allLocuses = append(allLocuses, key)
	}

	for _, fname := range *files {
		i++

		snpsChan := make(chan []snpInfo)

		go func() {
			snpsChan <- makeSnps(fname)
		}()
		snps = <-snpsChan

		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

		t1 := time.Now()

		fmt.Printf("Calculating Dn/DS: Working on %v from %v (%v)\r", i, len(*files), fname)
		for _, val := range snps {
			if val.TypeOf == "CDS" && containsPos(altPositions[val.Locus], val.PosInGene, val.Alt) == false {

				altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, alt: val.Alt, ref: val.NucInPos, locus: val.Locus})
			}

		}

		for _, allloc := range allLocuses {

			// prod := getProductByName(allloc)
			if len(altPositions[allloc]) > 2 {
				dndsChan := make(chan []string)
				go func() {

					dndsChan <- getDnDsByLocus(allloc, altPositions[allloc])
				}()
				dndsRes, ok := <-dndsChan
				if ok {
					// fmt.Println(dndsRes)

					locDNDS[dndsRes[0]] = append(locDNDS[dndsRes[0]], dndsRes[1])
					if dndsRes[1] == "1" {
						countNbrOne[allloc]++
					}
					close(dndsChan)
				}

			} else {

				locDNDS[allloc] = append(locDNDS[allloc], "1")
				countNbrOne[allloc]++
			}

		}

		fmt.Printf("Calculating Dn/DS: Working on %v from %v (%v) \t\t Time:\t%v\n", i, len(*files), fname, t1.Sub(t0))

	}

	// fmt.Println(locDNDS)

	for _, allloc := range allLocuses {

		if countNbrOne[allloc] != len(*files) && *statAll == false {
			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
		} else if *statAll == true {
			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
		}

	}
	// }

	headers.WriteString("\n")

	if buffer.Len() != 0 && headers.Len() != 0 {
		fOut, err := os.Create(fileOut)
		if err != nil {
			log.Fatal("Cannot create file", err)
		}
		defer fOut.Close()
		fmt.Fprintf(fOut, headers.String())
		fmt.Fprintf(fOut, buffer.String())
		fmt.Printf("\n\nWell done!\n")
		// t1 := time.Now()
		// fmt.Printf("Elapsed time: %v", fmtDuration(t1.Sub(t0)))
	}
}

func matrixDnDs(fileOut string) {
	// var AllPos []int
	var (
		allLocuses []string

		buffer  strings.Builder
		headers strings.Builder
		// var posCount = make(map[int]int)
		// snps         []snpInfo
		altPositions = make(map[string][]allPositionsInGene)
		locDNDS      = map[string][]string{}
		countNbrOne  = make(map[string]int)
	)
	t0 := time.Now()

	files := &listOfFiles

	i := 1
	headers.WriteString("Locus\t")

	for key := range geneCoordinates {
		allLocuses = append(allLocuses, key)
	}

	for key, snps := range snpCacheMap {

		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(key, filepath.Ext(key))))

		t1 := time.Now()

		fmt.Printf("Calculating Dn/DS: Working on %v from %v (%v)\r", i, len(*files), key)
		for _, val := range snps {
			if val.TypeOf == "CDS" && containsPos(altPositions[val.Locus], val.PosInGene, val.Alt) == false {

				altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, alt: val.Alt, ref: val.NucInPos, locus: val.Locus})
			}

		}

		for _, allloc := range allLocuses {

			// prod := getProductByName(allloc)
			if len(altPositions[allloc]) > 2 {
				dndsChan := make(chan []string)
				go func() {

					dndsChan <- getDnDsByLocus(allloc, altPositions[allloc])
				}()
				dndsRes, ok := <-dndsChan
				if ok {
					// fmt.Println(dndsRes)

					locDNDS[dndsRes[0]] = append(locDNDS[dndsRes[0]], dndsRes[1])
					if dndsRes[1] == "1" {
						countNbrOne[allloc]++
					}
					close(dndsChan)
				}

			} else {

				locDNDS[allloc] = append(locDNDS[allloc], "1")
				countNbrOne[allloc]++
			}

		}

		fmt.Printf("Calculating Dn/DS: Working on %v from %v (%v) \t\t Time:\t%v\n", i, len(*files), key, t1.Sub(t0))
		i++
	}

	// fmt.Println(locDNDS)

	for _, allloc := range allLocuses {

		if countNbrOne[allloc] != len(*files) && *statAll == false {
			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
		} else if *statAll == true {
			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
		}

	}
	// }

	headers.WriteString("\n")

	// if buffer.Len() != 0 && headers.Len() != 0 {
	// 	fOut, err := os.Create(fileOut)
	// 	if err != nil {
	// 		log.Fatal("Cannot create file", err)
	// 	}
	// 	defer fOut.Close()
	// 	fmt.Fprintf(fOut, headers.String())
	// 	fmt.Fprintf(fOut, buffer.String())
	// 	fmt.Printf("\n\nWell done!\n")
	// 	// t1 := time.Now()
	// 	// fmt.Printf("Elapsed time: %v", fmtDuration(t1.Sub(t0)))
	// }
	matrixPrint(headers, buffer, fileOut)
}

func matrixBinary(fileOut string) {
	var (
		AllPosUnsort, AllPos []int
		allLocusUnsort       []string
		buffer               strings.Builder
		headers              strings.Builder
		posCount             = make(map[int]int)
		// snps                 []snpInfo
		posFN   = make(map[int][]string)
		posFreq = map[int][]string{}
	)
	// var ResSeq []seqInfo

	// files := getListofVCF()
	files := &listOfFiles

	// fmt.Println(files)
	pos := make(map[int]string)

	headers.WriteString("Pos\t")
	i := 1
	for key, snps := range snpCacheMap {
		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(key, filepath.Ext(key))))

		fmt.Printf("Counting SNP positions: Working on %v files from %v \r", i+1, len(*files))

		for _, val := range snps {

			pos[val.APos] = val.Alt
			// posTest = append(posTest, posByFile{pos: pos, file: file, apos: val.APos})
			if strings.Contains(strings.Join(posFN[val.APos], " "), key) == false {
				posFN[val.APos] = append(posFN[val.APos], key)
			}

			AllPosUnsort = append(AllPosUnsort, val.APos)
			posCount[val.APos] = posCount[val.APos] + 1
			if val.TypeOf == "CDS" {
				allLocusUnsort = append(allLocusUnsort, val.Locus)
			}

		}

	}
	AllPos = unique(AllPosUnsort)
	// allLocuses = removeStringDuplicates(allLocusUnsort)
	sort.Ints(AllPos)

	for _, file := range *files {
		for _, allpos := range AllPos {

			if strings.Contains(strings.Join(posFN[allpos], " "), file) {

				posFreq[allpos] = append(posFreq[allpos], "1")

			} else {
				// fmt.Println(allpos, 0, file)
				posFreq[allpos] = append(posFreq[allpos], "0")
			}

		}
	}

	for _, allpos := range AllPos {
		if posCount[allpos] < len(*files) {

			buffer.WriteString(fmt.Sprintln(allpos, "\t", strings.Join(posFreq[allpos], "\t")))

		}
	}
	headers.WriteString("\n")

	// if buffer.Len() != 0 && headers.Len() != 0 {
	// 	fOut, err := os.Create(fileOut)

	// 	if err != nil {
	// 		log.Fatal("Cannot create file", err)
	// 	}
	// 	defer fOut.Close()
	// 	fmt.Fprintf(fOut, headers.String())
	// 	fmt.Fprintf(fOut, buffer.String())
	// 	fmt.Printf("\n\nWell done!\n")
	// }
	matrixPrint(headers, buffer, fileOut)

}

func matrixBinaryMST(fileOut string, verbose bool, exGenes map[int]int, exSNPs map[int]int, nonrandomize bool) {

	var (
		AllPos, SelectedPos []int
		// ResSeq              []seqInfo
		// passSNP = make(map[string]int)
		uniqSNP         = make(map[int]int)
		nbrOfSNP        = 50
		headers, buffer strings.Builder
		snpCount        = make(map[int]int)
		// posCount             = make(map[int]int)
	)
	// files := getListofVCF()

	// queryChan := make(chan vcfInfoQuery)

	if *statNbrOfNSP != 0 {
		nbrOfSNP = *statNbrOfNSP
	}

	files := &listOfFiles
	for i, file := range *files {
		// создаем запрос в виде типа vcfQuery, который передается через канал на выход <-qSNP.OutChan
		qSNP := &vcfQuery{File: file, OutChan: make(chan vcfInfoQuery), Print: verbose}
		go qSNP.request()
		// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
		snpRes := <-qSNP.OutChan
		snpCacheMap[snpRes.File] = snpRes.SnpInfo
		if verbose == true {
			fmt.Printf("Reading files:  %v from %v \r", i+1, len(*files))
		}

	}

	for _, snps := range snpCacheMap {

		for _, val := range snps {
			if len(exGenes) != 0 {

				for key, value := range exGenes {
					if val.APos >= key && val.APos <= value {
						// fmt.Println(val.Locus, val.Start, val.End, val.Product)
						// passSNP["genes"] = passSNP["genes"] + 1
						uniqSNP[val.APos] = 2 //2-EXCLUDED
						continue
					} else if exSNPs[val.APos] == 1 {
						// passSNP["snp"] = passSNP["snp"] + 1
						uniqSNP[val.APos] = 2 //2-EXCLUDED
						continue

						// fmt.Println(val.APos)

					} else {
						if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
							uniqSNP[val.APos] = 1 //1-INCLUDED

						}
					}
				}
			} else {
				if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
					uniqSNP[val.APos] = 1

				}
				// AllPosUnsort = append(AllPosUnsort, val.APos)
			}
			snpCount[val.APos] = snpCount[val.APos] + 1
			// if exGenes[val.Locus] != 1 && exSNPs[string(val.APos)] != 1 {
			// 	AllPosUnsort = append(AllPosUnsort, val.APos)
			// }
		}
	}

	// fmt.Println(uniqSNP)

	for key, value := range uniqSNP {

		if value == 1 && len(*files) != snpCount[key] {

			AllPos = append(AllPos, key)
			// } else if value == 2 {
			// 	fmt.Println(key)

		}

	}
	// go process("Working...          ")
	// AllPos = unique(AllPosUnsort)

	sort.Ints(AllPos)

	rand.Seed(time.Now().UnixNano())

	if nbrOfSNP > len(AllPos) {
		nbrOfSNP = len(AllPos) - 1
	}

	if nonrandomize == false {
		for i := 1; i <= nbrOfSNP; i++ {
			rnd := rand.Intn(len(AllPos)-i) + i
			// fmt.Println(AllPos[rnd])
			SelectedPos = append(SelectedPos, AllPos[rnd])
		}
	} else {
		for i := 0; i <= nbrOfSNP; i++ {
			// rnd := rand.Intn(len(AllPos)-i) + i
			// fmt.Println(AllPos[rnd])
			SelectedPos = append(SelectedPos, AllPos[i])
		}
	}
	// fmt.Println(AllPos[0])

	// rand.Intn(max - min) + min

	sort.Ints(SelectedPos)

	// for i := 0; i < len(SelectedPos); i++ {
	// 	headers.WriteString(fmt.Sprintln(SelectedPos[i], "\t"))
	// }
	// headers.WriteString(fmt.Sprintf("%v\n", strings.Join(SelectedPos, "\t")))
	headers.WriteString("\t")
	// buffer.WriteString("\t")
	headers.WriteString(strings.Trim(strings.Join(strings.Fields(fmt.Sprint(SelectedPos)), "\t"), "[]"))
	headers.WriteString("\n")
	if verbose == true {
		fmt.Println(headers.String())
		// fmt.Println(AllPos)
	}
	for fname, snps := range snpCacheMap {
		buffer.WriteString(fname)
		// if verbose == true {
		// 	fmt.Printf("Generating sequences: Working on  %v from %v \r", i+1, len(snps.File))
		// }
		pos := make(map[int]string)
		// var buffer strings.Builder

		// buffer.WriteString(fmt.Sprintf(">%v\n", strings.ToUpper(fname)))

		for _, val := range snps {
			pos[val.APos] = val.Alt
			// if exGenes[val.Locus] == 1 {
			// 	fmt.Println(val.Locus)
			// }

		}
		for _, allpos := range SelectedPos {
			// posCount[allpos] = posCount[allpos] + 1
			// fmt.Println(snpCount[allpos])

			if pos[allpos] != "" {

				buffer.WriteString("\t1")
			} else {

				buffer.WriteString("\t0")
			}
			// }

		}
		buffer.WriteString("\n")
		// ResSeq = append(ResSeq, seqInfo{Name: fname, Seq: buffer.String(), UsedPositions: SelectedPos})
		// fmt.Println(len(buffer.String()))

		// if *gbDebug == true {
		// 	fmt.Printf("%v\t:\nThere was passed %v SNPs from exclude gene file\n And %v SNPs from exclude snp file\n", fname, passSNP["genes"], passSNP["snp"])
		// }
	}
	if verbose == true {
		fmt.Println(buffer.String())
	}
	matrixPrint(headers, buffer, fileOut)
	// return ResSeq
}

func matrixBinaryOld(fileOut string) {
	var (
		AllPosUnsort, AllPos []int
		allLocusUnsort       []string
		buffer               strings.Builder
		headers              strings.Builder
		posCount             = make(map[int]int)
		snps                 []snpInfo
		posFN                = make(map[int][]string)
		posFreq              = map[int][]string{}
	)
	// var ResSeq []seqInfo

	// files := getListofVCF()
	files := &listOfFiles

	// fmt.Println(files)
	pos := make(map[int]string)

	headers.WriteString("Pos\t")

	for i, file := range *files {
		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

		fmt.Printf("Counting SNP positions: Working on %v files from %v \r", i+1, len(*files))

		snpsChan := make(chan []snpInfo)

		go func() {
			snpsChan <- makeSnps(file)
		}()
		snps = <-snpsChan

		for _, val := range snps {

			pos[val.APos] = val.Alt
			// posTest = append(posTest, posByFile{pos: pos, file: file, apos: val.APos})
			if strings.Contains(strings.Join(posFN[val.APos], " "), file) == false {
				posFN[val.APos] = append(posFN[val.APos], file)
			}

			AllPosUnsort = append(AllPosUnsort, val.APos)
			posCount[val.APos] = posCount[val.APos] + 1
			if val.TypeOf == "CDS" {
				allLocusUnsort = append(allLocusUnsort, val.Locus)
			}

		}

	}
	AllPos = unique(AllPosUnsort)
	// allLocuses = removeStringDuplicates(allLocusUnsort)
	sort.Ints(AllPos)

	for _, file := range *files {
		for _, allpos := range AllPos {

			if strings.Contains(strings.Join(posFN[allpos], " "), file) {

				posFreq[allpos] = append(posFreq[allpos], "1")
			} else {
				// fmt.Println(allpos, 0, file)
				posFreq[allpos] = append(posFreq[allpos], "0")
			}

		}
	}

	for _, allpos := range AllPos {
		if posCount[allpos] < len(*files) {

			buffer.WriteString(fmt.Sprintln(allpos, "\t", strings.Join(posFreq[allpos], "\t")))

		}
	}
	headers.WriteString("\n")

	if buffer.Len() != 0 && headers.Len() != 0 {
		fOut, err := os.Create(fileOut)
		if err != nil {
			log.Fatal("Cannot create file", err)
		}
		defer fOut.Close()
		fmt.Fprintf(fOut, headers.String())
		fmt.Fprintf(fOut, buffer.String())
		fmt.Printf("\n\nWell done!\n")
		// t1 := time.Now()
		// fmt.Printf("Elapsed time: %v", fmtDuration(t1.Sub(t0)))
	}

}

func makeSnps(fname string) (snps []snpInfo) {
	qSNP := &vcfQuery{File: fname, OutChan: make(chan vcfInfoQuery)}
	go qSNP.request()
	// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
	snpRes := <-qSNP.OutChan
	snps = snpRes.SnpInfo

	return snps
}

func getDnDsByLocus(locus string, altPositions []allPositionsInGene) (dndsRes []string) {
	var (
		dndsLoc string
	)
	// start, end := getGenePosByName(locus)

	refS := getGeneSequence(locus)
	altS := makeAltString(locus, altPositions) // fmt.Println(val, "\n", refS)

	qDnDs := &codon.DnDsQuery{RefSeq: refS, AltSeq: altS, OutChan: make(chan codon.DnDs)}
	go qDnDs.Request()
	// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
	dnds := <-qDnDs.OutChan
	close(qDnDs.OutChan)

	if dnds.ND != 0 && dnds.NS != 0 {

		// fmt.Printf("L:%v\tdNdS:%v\t%v\n", allloc, fmt.Sprintf("%.2f", dnds.DNDS), prod)
		dndsLoc = fmt.Sprintf("%.2f", dnds.DNDS)

	} else {

		// fmt.Printf("L:%v\tdNdS:%v\t%v\n", allloc, "1", prod)
		dndsLoc = "1"
	}
	dndsRes = append(dndsRes, locus)
	dndsRes = append(dndsRes, dndsLoc)

	return dndsRes

}

func containsPos(s []allPositionsInGene, pos int, alt string) bool {
	for _, a := range s {
		if a.pos == pos && a.alt == alt {
			return true
		}
	}
	return false
}

func removeStringDuplicates(elements []string) []string {
	encountered := map[string]bool{}

	// Create a map of all unique elements.
	for v := range elements {
		encountered[elements[v]] = true
	}

	// Place all keys from the map into a slice.
	result := []string{}
	for key := range encountered {
		result = append(result, key)
	}
	sort.Strings(result)
	return result
}

func getSequenceRange(PosRange string, noseq bool) rangePosInfo {

	// var rangeParser = regexp.MustCompile(`^(\d+)\W+(\d+)|\D+\W+(\d+)\W+(\d+)`)

	var (
		rangeParser = regexp.MustCompile(`\b(\d+)\b\W+\b(\d+)\b`)

		result rangePosInfo
	)
	PosRange = strings.Replace(PosRange, ".", "", -1)

	for _, val := range rangeParser.FindAllStringSubmatch(PosRange, -1) {

		if val[1] != "" && val[2] != "" {
			startR, _ := strconv.Atoi(val[1])
			endR, _ := strconv.Atoi(val[2])
			if startR < endR {
				resStart = startR
				resEnd = endR
			} else {
				resStart = endR
				resEnd = startR
			}

		}

	}

	if resStart != 0 && resEnd != 0 {

		if noseq == true {
			gname, _ := getGeneNameByPos(resStart, resEnd)
			prod, _ := getProductByPos(resStart, resEnd)
			result = rangePosInfo{Start: resStart + 1, End: resEnd, Gname: gname, Prod: prod}

		} else if noseq == false {
			seq := getNucFromGenome(resStart-1, resEnd)
			gname, _ := getGeneNameByPos(resStart, resEnd)
			prod, _ := getProductByPos(resStart, resEnd)
			result = rangePosInfo{Start: resStart + 1, End: resEnd, Gname: gname, Len: len(seq), Seq: seq, Prod: prod}

		}
		// fmt.Println(result)
	} else {
		fmt.Println("Range positions is not valid")
	}

	return result
	// fmt.Println(positions)

}

func printSequenceRange(rangeList []rangePosInfo, web bool, port string) {

	switch web {
	case false:
		const basicAnnotation = "{{range $element := .}}" +
			"{{if .Seq}}" +
			"{{.Gname}}\t{{.Start}}\t{{.End}}\t{{.Len}}\t{{.Seq}}\t{{.Prod}}\t{{.Doubles}}\n" +
			"{{else}}" +
			"{{.Gname}}\t{{.Start}}\t{{.End}}\t{{.Len}}\t{{.Prod}}\t{{.Doubles}}\n" +
			"{{end}}" +
			"{{end}}"

		t := template.New("basic")
		t, err := t.Parse(basicAnnotation)
		err = t.Execute(os.Stdout, rangeList)
		if err != nil {
			log.Fatal("Parse: ", err)
			return
		}
	case true:

		var htmlAnnotation = `   <!DOCTYPE html>
			<html>
			<head>
			<meta charset="utf-8">			
			</head>		
			<table width="100%" cellspacing="0" cellpadding="4" border="1">
			<body>
			<tr>			
			{{range $element := .}}	
			{{if .Seq}}
			<td>{{.Gname}}</td><td>{{.Start}}</td><td>{{.End}}</td><td>{{.Len}}</td><td><p title="GC content: {{.GC}}"><textarea rows="3" style="width:400px; word-wrap:break-word;">{{.Seq}}</textarea></p></td><td>{{.Len}}</td><td><p title="{{.Note}}">{{.Prod}}</p></td><td>{{.Doubles}}</td>
			{{else}} 
			<td>{{.Gname}}</td><td>{{.Start}}</td><td>{{.End}}</td><td>{{.Len}}</td><td><p title="{{.Note}}">{{.Prod}}</p></td><td>{{.Doubles}}</td>
			{{end}}
			</tr>	
			{{end}}		
			</body>
			</table>
		
			<table width="100%" cellspacing="0" cellpadding="4" border="0">
			<tr>
			<td><a href="http://bsatool.ru" target="_blank">Created by BSATool (Bacterial Snp Annotation Tool)</a></td>
			</tr>
			</table>
			
			
`
		tH := template.New("html")

		tH, err := tH.Parse(htmlAnnotation)

		if err != nil {
			panic(err)
		}

		browser.OpenURL(fmt.Sprintf("localhost:%v", port))

		http.HandleFunc("/", func(w http.ResponseWriter, r *http.Request) {
			// err = t.Execute(w, &gInfo)

			err = tH.Execute(w, rangeList)
			if err != nil {
				panic(err)
			}
			go func() {
				defer os.Exit(0)
			}()
		})

		// if *flgPort != 8080 {
		locPort := fmt.Sprintf(":%v", port)
		http.ListenAndServe(locPort, nil)
		// } else {
		// 	http.ListenAndServe(":8080", nil)
		// }

	}
}

func getHashSNP(snp snpInfo) uint64 {
	// hasher := sha1.New()
	// hasher.Write([]byte(text))
	// return hex.EncodeToString(hasher.Sum(nil))
	hash, err := hashstructure.Hash(snp, nil)
	if err != nil {
		panic(err)
	}

	// fmt.Printf("%d", hash)
	return hash
}

// func getLocusHash(start, end int) string {
// 	lStart := strconv.Itoa(start)
// 	lEnd := strconv.Itoa(end)
// 	return getMD5Hash(fmt.Sprintf("%v%v", lStart, lEnd))
// }

func toBED() {
	var (
		genomeName string
	)

	if *statCircosGenome != "" {
		genomeName = *statCircosGenome
	} else {
		genomeName = gInfo.Version

	}

	for _, g := range allGenesVal {
		switch *statCircosTypeOf {
		case "":
			if *statShowAnnotation == true {
				fmt.Printf("%v\t%v\t%v\t%v\t%v\n", genomeName, g.Start, g.End, g.Locus, g.Product)
			} else if *statShowAnnotation == false {
				fmt.Printf("%v\t%v\t%v\t%v\n", genomeName, g.Start, g.End, g.Locus)
			}
		case "CDS":
			if g.TypeOf == "CDS" && *statShowAnnotation == true {
				fmt.Printf("%v\t%v\t%v\t%v\t%v\n", genomeName, g.Start, g.End, g.Locus, g.Product)
			} else if g.TypeOf == "CDS" && *statShowAnnotation == false {
				fmt.Printf("%v\t%v\t%v\t%v\n", genomeName, g.Start, g.End, g.Locus)
			}
		case "IGR":
			if g.TypeOf == "IGR" {
				fmt.Printf("%v\t%v\t%v\t%v\n", genomeName, g.Start, g.End, g.Locus)
			}
		}

	}

}

func getCoordRange(start, end int) {

	var coordArray []int

	coordArray = append(coordArray, start)

	for i := start; i <= end; i++ {

		for _, g := range allGenesVal {

			if g.Start == i {

				coordArray = append(coordArray, g.Start)

			} else if g.End == i {

				coordArray = append(coordArray, g.End)

			}

		}
	}
	coordArray = append(coordArray, end)
	sort.Ints(coordArray)

	var res []rangePosInfo

	for _, val := range coordArray {

		res = append(res, getSequenceRange(fmt.Sprintf("%v:%v", val, val), *gbNoSeq))

	}

	var last string
	for _, val := range res {
		if val.Gname != last {
			fmt.Println(val.Gname, val.Prod)
		}
		last = val.Gname

	}

}

func getGenomeMap() []genomeMapInfo {

	var (
		gmap    []genomeMapInfo
		gmapval genomeMapInfo
	)
	for _, g := range allGenesVal {
		gmapval = genomeMapInfo{Start: g.Start, End: g.End, Locus: g.Locus, TypeOf: g.TypeOf, Product: g.Product}
		gmap = append(gmap, gmapval)
	}

	return gmap
}

func getRangeFromFile(file string, verbose bool, noseq bool) []rangePosInfo {
	var (
		rangeParser = regexp.MustCompile(`\b(\d+)\b\W+\b(\d+)\b`)

		// var rangeArr []rangeArray
		posRange, seq string
		gc            float64
		unsorted      []rangeArray
		result        []rangePosInfo
		// var uniqRange = map[int]rangeArray{}
		// var lastHash, i, j, k int
		j, k int
	)
	f, err := os.Open(file)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	// var res []rangePosInfo
	for scanner.Scan() {

		posRange = scanner.Text()

		posRange = strings.Replace(posRange, ".", "", -1)

		for _, val := range rangeParser.FindAllStringSubmatch(posRange, -1) {

			if val[1] != "" && val[2] != "" {
				startR, _ := strconv.Atoi(val[1])
				endR, _ := strconv.Atoi(val[2])
				if startR < endR {
					resStart = startR
					resEnd = endR
				} else {
					resStart = endR
					resEnd = startR
				}

			}
			if resStart != 0 && resEnd != 0 {

				currRange := rangeArray{Start: resStart, End: resEnd}
				unsorted = append(unsorted, currRange)
				j++
				if verbose == true {
					fmt.Printf("Processed %v ranges...\r", j)
				}

			}

		}

	}
	// fmt.Println(uniqRange)
	sort.Slice(unsorted, func(i, j int) bool {
		return unsorted[i].Start < unsorted[j].Start
	})

	encountered := map[int]bool{}
	doubles := map[int]int{}
	found := []rangeArray{}

	for v := range unsorted {
		if encountered[unsorted[v].Start+unsorted[v].End] == true {
			// Do not add duplicate.
			doubles[unsorted[v].Start+unsorted[v].End] = doubles[unsorted[v].Start+unsorted[v].End] + 1
		} else {
			// Record this element as an encountered element.
			encountered[unsorted[v].Start+unsorted[v].End] = true
			// Append to result slice.
			found = append(found, unsorted[v])
		}
	}

	for _, val := range found {
		if noseq == true {
			seq = ""
		} else {
			seq = getNucFromGenome(val.Start-1, val.End)
			_, _, _, gc = codon.GcCodonCalc(seq)
		}
		gname, _ := getGeneNameByPos(val.Start, val.End)
		prod, note := getProductByPos(val.Start, val.End)
		fixedProd := strings.Replace(prod, " + ", " ", -1)
		gcRes, _ := strconv.ParseFloat(fmt.Sprintf("%.2f", gc), 64)
		// fixedProd = strings.Replace(prod, "'", " ", -1)

		res := rangePosInfo{Start: val.Start, End: val.End, Gname: gname, Prod: fixedProd, Len: val.End - val.Start + 1, Seq: seq, Doubles: doubles[val.Start+val.End] + 1, Note: note, GC: gcRes}

		result = append(result, res)
		k++
		if verbose == true {
			fmt.Printf("Annotated %v ranges from %v %v\r", k, len(found)-1, strings.Repeat(" ", 10))
		}
	}

	return result
}

func checkRuleFromFile(file string) (rules []rulesInfo) {

	var (
		res   []rulesInfo
		rinfo rulesInfo
	)
	rRule := regexp.MustCompile(`^(.*)\b=\b(\w+)$`)
	// ruleMap := map[string][]string{}
	f, err := os.Open(file) // открываем файл

	if err != nil {
		fmt.Println(err)

	}
	defer f.Close()
	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		if !strings.HasPrefix(scanner.Text(), "#") {
			for _, rrule := range rRule.FindAllStringSubmatch(scanner.Text(), -1) {
				// fmt.Println(rrule[1], rrule[2])
				// fmt.Println(strings.LastIndexAny(scanner.Text(), rrule[1]))

				convert := strings.Split(strings.ToUpper(rrule[1]), ",")
				// ruleMap[rrule[2]] = convert
				sort.Slice(convert, func(i, j int) bool {
					return convert[i] < convert[j]
				})
				rinfo = rulesInfo{Name: rrule[2], Variants: convert, Lenght: len(convert)}
				// fmt.Println(, len(convert))
			}
			res = append(res, rinfo)
		}
	}

	return res

}

// расчитывает индекс АК замен из файла. формат записи в файле: Pro tab Ala
func calcComplexIdxFromFile(file string) {
	// var1, var2 := amino.GetComplexIndex("P", "L", false)
	// 			fmt.Println(var1, var2)
	// 			fmt.Println(amino.GetShortNameAA("Pro"))
	f, err := os.Open(file) // открываем файл

	if err != nil {
		fmt.Println(err)

	}

	scanner := bufio.NewScanner(f)
	// var res []rangePosInfo
	for scanner.Scan() {
		scanTxt := scanner.Text()
		if strings.Contains(scanTxt, "\t") {
			// fmt.Println(scanTxt)
			aaArray := strings.Split(scanTxt, "\t")
			idxVal, idxRes := amino.GetComplexIndex(amino.GetShortNameAA(aaArray[0]), amino.GetShortNameAA(aaArray[1]), false)
			// fmt.Println()
			fmt.Printf("%v\t%v\t%v\t%v\n", aaArray[0], aaArray[1], idxVal, idxRes)
		}
	}
	defer f.Close()
}

func makeSeqFromSNPListFile(f string) {

	var (
		// parsedSNP    []snpCheckInfo
		snpCheckChan = make(chan snpCheckInfo)
		seq          string
	)
	file, err := os.Open(f)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {

		if !strings.HasPrefix(scanner.Text(), "#") {

			go func() {

				snpCheckChan <- checkSNPnotation(scanner.Text())
			}()
			snpCheck := <-snpCheckChan
			if snpCheck.Locus == "" {
				aposInt, _ := strconv.Atoi(snpCheck.APos)
				locName, _ := getGeneNameByPos(aposInt, aposInt)
				snpCheck.Locus = locName
			}
			// fmt.Println(snpCheck)
			if snpCheck.Locus != "" {
				// start, end := getGenePosByName(snpCheck.Locus)
				// parsedSNP = append(parsedSNP, snpCheck)

				// if start != 0 && end != 0 {

				dir := getDirectionByName(snpCheck.Locus)

				// if dir == "f" {
				// // 	seq = getNucFromGenome(start-1, end)
				// } else if dir == "r" {
				// // 	seq = getReverseComplement(getNucFromGenome(start-1, end))

				// }
				seq = getGeneSequence(snpCheck.Locus)

				prod := getProductByName(snpCheck.Locus)
				seqArr := strings.Split(seq, "")
				posInGene, _ := strconv.Atoi(snpCheck.PosInGene)
				for i := 0; i < len(seqArr); i++ {
					if *gbVerbose == true {

						seqArr[i] = fmt.Sprintf("%v->%v,", i+1, seqArr[i])
					}
					if i+1 == posInGene {
						if dir == "f" {
							seqArr[i] = fmt.Sprintf("[%v>%v]", strings.ToUpper(seqArr[i]), snpCheck.Alt)
						} else if dir == "r" {
							seqArr[i] = fmt.Sprintf("[%v>%v]", strings.ToUpper(seqArr[i]), getComplement(snpCheck.Alt))
						}
					}
				}
				// seqArr[posInGene+1] = fmt.Sprintf("(%v_%v>%v)", posInGene, strings.ToUpper(seqArr[posInGene+1]), strings.ToUpper(snp.Alt))
				altSeq := strings.Join(seqArr, "")
				// fmt.Println(strings.TrimRight(altSeq, ","))
				fmt.Printf(">%v_%v_%v>%v (%v) %v\n%v\n", snpCheck.Locus, snpCheck.PosInGene, snpCheck.Ref, snpCheck.Alt, prod, snpCheck.Name, altSeq)
				// fmt.Println("!\n", getReverseComplement(seq))
				// fmt.Println(seqArr[posInGene+1], snp)
				// fmt.Println(seq, gc, gc1, gc2, gc3, prod)
				// }

			}
		}
	}

}

func getReverseComplement(sequence string) string {
	// var (
	// 	altSeq []string
	// 	resSeq string
	// )
	// var result string
	seqArr := strings.Split(sequence, "")

	for i := 0; i < len(seqArr); i++ {
		switch seqArr[i] {
		case "a":
			seqArr[i] = "t"
		case "c":
			seqArr[i] = "g"
		case "t":
			seqArr[i] = "a"
		case "g":
			seqArr[i] = "c"
		case "A":
			seqArr[i] = "T"
		case "C":
			seqArr[i] = "G"
		case "T":
			seqArr[i] = "A"
		case "G":
			seqArr[i] = "C"

		}
	}
	complStr := strings.Join(seqArr, "")

	runes := []rune(complStr)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}

	return string(runes)

}

func getComplement(sequence string) string {
	// var (
	// 	altSeq []string
	// 	resSeq string
	// )
	// var result string
	seqArr := strings.Split(sequence, "")

	for i := 0; i < len(seqArr); i++ {
		switch seqArr[i] {
		case "a":
			seqArr[i] = "t"
		case "c":
			seqArr[i] = "g"
		case "t":
			seqArr[i] = "a"
		case "g":
			seqArr[i] = "c"
		case "A":
			seqArr[i] = "T"
		case "C":
			seqArr[i] = "G"
		case "T":
			seqArr[i] = "A"
		case "G":
			seqArr[i] = "C"

		}
	}
	complStr := strings.Join(seqArr, "")

	return complStr

}

func getGeneSequence(locus string) string {
	var seq string
	start, end := getGenePosByName(locus)
	if start != 0 && end != 0 {
		dir := getDirectionByName(locus)
		if dir == "f" {
			seq = getNucFromGenome(start-1, end)
		} else if dir == "r" {
			seq = getReverseComplement(getNucFromGenome(start-1, end))

		}
	} else {
		fmt.Println("Locus not found [getGeneSequence func]")
	}
	return seq
}

func loadExcludeGenes(file string) map[int]int {
	var (
		exGenes       = make(map[int]int)
		exLocusCoords = regexp.MustCompile(`^\w+\W+(\d+)\W+(\d+)`)
		start, end    int
	)

	f, err := os.Open(file) // открываем файл

	if err != nil {
		fmt.Println(err)

	}
	defer f.Close()

	scanner := bufio.NewScanner(f) //  новый сканер

	for scanner.Scan() {

		scanTxt := scanner.Text()

		for _, exGene := range exLocusCoords.FindAllStringSubmatch(scanTxt, -1) {
			start, _ = strconv.Atoi(exGene[1])
			end, _ = strconv.Atoi(exGene[2])

			// fmt.Println(exGene)
			// exGenes[exGene[0]] = 1
			exGenes[start] = end

		}
	}
	return exGenes
}

func loadExcludeSNP(file string) map[int]int {

	var (
		exSNPs = make(map[int]int)
		snpPos int
	)

	locSNPcheck := readSNPFromFile(file)
	for _, val := range locSNPcheck {
		snpPos, _ = strconv.Atoi(val.APos)
		exSNPs[snpPos] = 1
	}
	// fmt.Println(exSNPs)
	return exSNPs
}

func openbrowser(url string) {
	var err error

	switch runtime.GOOS {
	case "linux":
		err = exec.Command("xdg-open", url).Start()
	case "windows":
		err = exec.Command("rundll32", "url.dll,FileProtocolHandler", url).Start()
	case "darwin":
		err = exec.Command("open", url).Start()
	default:
		err = fmt.Errorf("unsupported platform")
	}
	if err != nil {
		log.Fatal(err)
	}

}
