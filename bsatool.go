package main

import (
	"bufio"
	"io/ioutil"
	"math/rand"
	"os/exec"
	"path"
	"runtime"
	"runtime/pprof"
	"time"

	"fmt"

	"log"

	"encoding/gob"

	"net/http"
	"os"

	// "reflect"
	"regexp"

	"strings"

	"compress/lzw"
	"sort"
	"strconv"

	"./amino"

	"gopkg.in/alecthomas/kingpin.v2"

	"html/template"

	"./codonPkg"

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
		`BSATool - Bacterial Snp Annotation Tool ` + "\n" +
		`      Laboratory of Social and Epidemic Infections
 Scientific Centre for Family Health and Human Reproduction Problems
     	(c) V.Sinkov, P.Khromova, O.Ogarkov, Irkutsk, Russia, 2017-2020                                   
                                                  
	`
	// list   = "list"
	ncFlag = "NC"
	aaFlag = "AA"
	tLMN   = "LMN"   // locus:Mutation:NAME
	tPMLN  = "PMLN"  // position:Mutation:locus:NAME
	tPMN   = "PMN"   // position:Mutation:NAME
	tLSAAN = "LSAAN" // locus:shortAA:codon:shortAA:name
	tLLAAN = "LLAAN" // locus:longAA:codon:longAA:name
	tLCN   = "LCN"   // locus:codon:name
	tSEN   = "SEN"   // start|end:name (any)
	// tPMNT  = "PMNL"  //position_Ref>Alt{Tab}Name;tag (position, mutation,name, tag )
	// tPOS   = "checkPOS"   // check:apos:name
	// tCODON = "checkCODON" // check:locus:codon:name
	// vcfExt = "*.vcf"
	// pbtmpl      = `{{counters . }}{{ bar . "⎜" "agtc"   "⬤" "TCAG" "⎜"}}{{percent . }}{{rtime . " ETA %s  "}}{{speed .}} `
	lmnRegExp  = `^(\w+)\W+(\d+)(\D)>(\D)\s+(.*)` // LMN Regular expression
	pmlnRegExp = `^(\d+)_(\D)>(\D)\W+(\w+)\W+(.*)`
	pmnRegExp  = `(\d+)_(\D)>(\D)\W+(.*)$`
	// pmntRegExp  = `^(\d+)_(\D)>(\D)(.*)\b(;|:)\b(.*)`
	lsaanRegExp = `^(\w+)\W+(\D{1})(\d+)(\D{1})\W+(.*)`
	llaanRegExp = `^(\w+)\W+(\w{3})\W+(\d+)\W+(\w{3})\W+(.*)`
	lcnRegExp   = `^(\w+)\W+codon(\d+)\W+(.*)`
	senRegExp   = `^(\d+)\|(\d+):(\w+)`
	// chkPosRegExp   = `^check\W+(\d+)\W+(.*)`
	// chkCodonRegExp = `^check\W+(\w+)\W+(\d+)\W+(.*)`
)

var (

	//
	// Database flags
	version string
	build   string
	// app       = kingpin.New(logo, "BSATool - Bacterial Snp Annotation Tool")
	// appAuthor = app.Author("V.Sinkov")
	// appVer    = app.Version(fmt.Sprintf("%v %v", version, build))

	gbWeb     = kingpin.Flag("web", " Open results in web browser").Short('w').Bool()
	gbXLSX    = kingpin.Flag("xlsx", " Export to XLSX format").Short('x').String()
	gbVerbose = kingpin.Flag("verbose", "Show additional information ").Short('v').Default("false").Bool()
	gbIndex   = kingpin.Flag("index", "Calculate Complex Index for amino acid changes").Default("false").Bool()
	gbPort    = kingpin.Flag("port", "Use your own localhost:port (default:8080)").Default("8080").String()
	gbNoSeq   = kingpin.Flag("noseq", "Don't show nucleotides").Default("false").Bool()
	gbDebug   = kingpin.Flag("debug", "Debug mode").Default("false").Bool()
	// gbLog            = kingpin.Flag("log", "write log file").Default("false").Bool()
	gbExcludeGenes   = kingpin.Flag("exclude-genes", "file with genes which should be excluded from mkseq").String()
	gbExcludeRegions = kingpin.Flag("exclude-regions", "file with coordinates(start end)  which should be excluded from mkseq").String()
	gbExcludeSnp     = kingpin.Flag("exclude-snp", "file with genes which should be excluded from mkseq").String()
	gbRandomize      = kingpin.Flag("randomize", "Set on randomizer").Default("false").Bool()
	mkdb             = kingpin.Command("mkdb", "Create database")
	dbName           = mkdb.Flag("out", "Name of database").Short('o').Required().String()
	dbGenbank        = mkdb.Flag("gb", "Name of genbank file").Short('i').Required().String()
	gbDP             = kingpin.Flag("dp", "dp value for filtering vcfs").Default("1").Int()
	gbIGR            = kingpin.Flag("igr", "Show igr regions in results").Default("false").Bool()
	// gbAbout        = kingpin.Flag("about", "About programm").Bool()

	// Annotation flags

	annAction = kingpin.Command("annotate", "This is command to allow annotate VCF file").Alias("annotation").Alias("ann").Alias("a")
	annDB     = annAction.Flag("db", "Name of database").Short('b').Required().String()
	annVCF    = annAction.Flag("vcf", "Input VCF file").Short('i').Required().String()
	// annWeb           = annAction.Flag("web", "").Bool()
	annMakeSeq       = annAction.Flag("mkseq", "NC, AA, Binary").Short('m').String()
	annOutFormat     = annAction.Flag("output-format", "Output format for binary data (phylip, nexus)").String()
	annMakeSeqRef    = annAction.Flag("ref", "Generate reference sequence").Short('r').Default("false").Bool()
	annWithFilenames = annAction.Flag("wfn", "Show filenames in list annotated VCF's").Short('n').Bool()
	annInDel         = annAction.Flag("indel", "indel detection").Bool()
	annBench         = annAction.Flag("annprof", "cpuprofile").String()
	annSeqLen        = annAction.Flag("len", "length of sequence").Int()
	annShowFileName  = annAction.Flag("filename", "Show filename in first column for next analysis").Default("false").Bool()
	annPosFile       = annAction.Flag("pos-file", "Create positions.txt file with information of each pos in sequence").String()
	annMinPosNbr     = annAction.Flag("min-pos-nbr", "Set minimal threshhold of presence of position in sequence").Default("1").Int()
	// annAltPosPercent     = annAction.Flag("alt-percent", "Maximal percent of presence of reference allele in sequence (99 default)").Default("1").Int()
	annAltRange = annAction.Flag("alt-range", "Min-Max range values (percent) of presence alternative nucleotide in sequence").Default("1:99").String()
	// annBench         = annAction.Flag("cpuprofile", "cpuprofile").String()

	// compute statistic options

	statAction = kingpin.Command("stat", "Calculates statistic tests").Alias("s")
	statDB     = statAction.Flag("db", "Database file").Short('b').Required().String()
	statTask   = statAction.Flag("action", "Type of action:share, snp,dnds, ginfo,bed, matrix,range, circos").Short('a').Required().String()
	// statWeb     = annAction.Flag("web", "").Short('w').Bool()
	statInFile  = statAction.Flag("in", "Input file").Short('i').String()
	statOutFile = statAction.Flag("out", "Output file").Short('o').String()
	statTypeOf  = statAction.Flag("type", "Type of matrix (binary, gc3, dnds, nc, locus. freq, jw, summary").Short('t').String()
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
	statGroupFromFile   = statAction.Flag("group", "File with filenames and their groups").String()
	statVCF             = statAction.Flag("vcf", "Input VCF file").String()
	statLocus           = statAction.Flag("locus", "locus name").String()
	statDnDsCountTh     = statAction.Flag("th", "DnDs threshold (How many files consists nonzero values (Default > 1))").Int()
	// statGStrain         = statAction.Flag("strain", "set strain name").String()
	statMkSeq    = statAction.Flag("with_seq", "make sequence").Bool()
	statCircos   = statAction.Flag("circos", "export to circos file").Bool()
	statCircosGC = statAction.Flag("gc", "export to circos file").Bool()
	// statGName           = statAction.Flag("gname", "set strain name").String()
	statFlankLeft  = statAction.Flag("flank_left", "number of left nucleotides").Int()
	statFlankRight = statAction.Flag("flank_right", "number of right nucleotides").Int()
	statAsTable    = statAction.Flag("table", "show results of check info as table").Bool()
	statReverse    = statAction.Flag("reverse", "make seq for complementnary genes in gene direction").Bool()
	// statAsBinary        = statAction.Flag("binary", "show results of check info as binary sequence").Bool()
	statBinTH    = statAction.Flag("nbr_pos", "nbr found SNP (default = nbr files-1)").Int()
	statDnDsStat = statAction.Flag("show-info", "show total statistics per genome").Bool()
	// statMkSeq   = statAction.Flag("mkseq", "").Bool()
	// statTH      = statAction.Flag("th", "").Int()

	// filterAction = kingpin.Command("filter", "filterVCFs")
	// filterDB     = filterAction.Flag("db", "Database file [Created by mkdb command]").Short('b').Required().String()

	infoAction     = kingpin.Command("info", "Get information")
	infoLocus      = infoAction.Flag("locus", "locus name").String()
	infoDB         = infoAction.Flag("db", "Database file [Created by mkdb command]").Short('b').Required().String()
	infoRanges     = infoAction.Flag("range", "Show sequence of nucleotides [start:end]").String()
	infoCodons     = infoAction.Flag("codons", "Show sequence of codons [start:end]").String()
	infoShowAs     = infoAction.Flag("showas", "Show as:direct (from left to right),gene(direction as in gene) ").String()
	infoMask       = infoAction.Flag("mask", "pos:ref:alt").String()
	infoFlankLeft  = infoAction.Flag("flank_left", "number of left nuclotides").Int()
	infoFlankRight = infoAction.Flag("flank_right", "number of right nuclotides").Int()
	infoIUPAc      = infoAction.Flag("iupac", "pos:ref:alt").Bool()

	devAction = kingpin.Command("dev", "Developer mode.").Alias("debug")
	devDB     = devAction.Flag("db", "Database file").Short('b').Required().String()
	devTask   = devAction.Flag("action", "Action...").Short('a').Required().String()
	devPwd    = devAction.Flag("pwd", "Password").Short('p').Required().String()

	betaAction = kingpin.Command("beta", "Beta test mode for testing new functions")
	betaTask   = betaAction.Flag("action", "Action...").Short('a').Required().String()
	betaDB     = betaAction.Flag("db", "Database").Short('b').Required().String()
	betaInFile = betaAction.Flag("in", "Input file").Short('i').String()
	// betaOutFile = betaAction.Flag("out", "Output file").Short('o').String()

	pileup2fasta  = kingpin.Command("parse_pileup", "Pileup parser")
	pileupInFile  = pileup2fasta.Flag("in", "Input file").Short('i').String()
	pileupTH      = pileup2fasta.Flag("th", "DnDs threshold (How many files consists nonzero values (Default > 1))").Int()
	pileupGStrain = pileup2fasta.Flag("strain", "set strain name").String()
	pileupMkSeq   = pileup2fasta.Flag("mkseq", "make sequence").Bool()
	pileupShowSeq = pileup2fasta.Flag("showseq", "make sequence").Bool()
	pileupGName   = pileup2fasta.Flag("genome", "set genome name").String()
	pileupCircos  = pileup2fasta.Flag("circos", "show as circos file format").Bool()
	pileupDB      = pileup2fasta.Flag("db", "Database file").Short('b').Required().String()

	// uniqAction = kingpin.Command("uniq", "Find out unique SNPs").Alias("unique")
	// uniqInFile = uniqAction.Flag("in", "Input file list").Short('i').Required().String()
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

	// colorReset  = "\x1b[39;49;0m"
	colorRed = "\x1b[31;1m"
	// colorGreen  = "\x1b[32;1m"
	// colorYellow = "\x1b[33;1m"
	// colorDbBlue = "\x1b[35;1m"
	// colorBlue   = "\x1b[36;1m"
	// colorWhite  = "\x1b[37;1m"
)

type (
	gCoords struct {
		Start, End int    // начало и конец гена
		Type       string // тип гена CDS
	}

	/*
		информация об гене. Имя локуса, имя гена, Направление и прочее
	*/
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
		APos, PosInGene, PosInCodonG, CodonNbrInG, GeneLen, Start, End, NucCore, TangIdxVal, DP, Indel int
		RefCodon, AltCodon, RefAA, AltAA, Locus,
		Direction, NucInPos, Product, Name,
		RefAAShort, AltAAShort, Mutation, Tang, Alt, Note, ReportType, ProteinID, GeneID, GOA, TiTv, TypeOf, ComplexIndex, FName, IndelType, IndelRef,
		IndelAlt string
		InterPro, PDB, ProSite []string
	}

	snpCheckInfo struct {
		Locus, PosInGene, CodonNbrInG, Ref, Alt, Name, TypeOf, APos, AASref, AASalt, AALref, AALalt, Raw, Tag string
		StartRange, EndRange                                                                                  int
	}

	seqInfo struct {
		Name, Seq, TypeOfSeq string //
		UsedPositions        []int
		DebugInfo            []string
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
		FileName, FoundSNP, ExactWithRule, RuleNames string
	}

	// DnDsRes is....
	// DnDsRes struct {
	// 	N, S, PN, PS, DN, DS, DNDS float64
	// 	ND, NS                     int
	// 	Locus, Product             string
	// }

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
		Start, End, Len, Doubles          int
		Seq, Prod, GeneName, Note, Genome string
		GC                                float64
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

	// ------------------------------------------------------------------------

	allPositionsInGene struct {
		pos             int
		alt, ref, locus string
	}

	// ------------------------------------------------------------------------

	altStringResult struct {
		start, end                   int
		locus, prod, altSeq, vcfFile string
	}

	// ------------------------------------------------------------------------

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
	q.OutChan <- vcfInfoQuery{File: q.File, SnpInfo: parserVCF(q.File, false, *gbDP, allGenesVal)}

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

func main() {
	// загружаем список vcf файлов
	listOfFiles = getListofVCF()
	// парсинг флагов
	defer os.Exit(0)
	// flag.Parse()
	kingpin.New(logo, "BSATool - Bacterial Snp Annotation Tool")
	kingpin.Version(fmt.Sprintf("%v %v", version, build)).Author("V.Sinkov")
	kingpin.UsageTemplate(kingpin.LongHelpTemplate)

	switch kingpin.Parse() {

	// bsatool mkdb/create/makedb -i /--genbank FILE --out /-o FILE

	case "mkdb":

		if _, err := os.Stat(*dbGenbank); os.IsNotExist(err) {
			fmt.Printf("The %v file is not exist!\n", *dbGenbank)
			os.Exit(3)
		}
		allGenesVal, genomeSeqSlice = geneBankFileParser(*dbGenbank)
		writeDB(*dbName, allGenesVal)

	case "annotate":

		allGenesVal = readDB(*annDB)
		var (
			// excludedGenes []string

			exGenes      = make(map[int]int)
			exGenesRaw   = make(map[int]int)
			exRegions    = make(map[int]int)
			exGenesBySNP = make(map[int]int)
			exSNPs       = make(map[int]int)
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

				exGenesRaw = loadExcludeGenes(*gbExcludeGenes)

			}
		}

		if *gbExcludeRegions != "" {
			if _, err := os.Stat(*gbExcludeRegions); os.IsNotExist(err) {
				fmt.Println("No file with excluded regions was found")
				os.Exit(3)

			} else {
				// exFileGenes = true

				exRegions = loadExcludeRegion(*gbExcludeRegions)
			}
		}

		if *gbExcludeSnp != "" {
			if _, err := os.Stat(*gbExcludeSnp); os.IsNotExist(err) {
				fmt.Println("No file with excluded SNPs was found")
				os.Exit(3)

			} else {
				// exSNPpos = true

				exSNPs = loadExcludeSNP(*gbExcludeSnp)
				exGenesBySNP = excludeGeneByPos(*gbExcludeSnp)

			}

		}
		if len(exGenesRaw) != 0 {
			for k, v := range exGenesRaw {
				exGenes[k] = v
			}
		}
		if len(exRegions) != 0 {
			for k, v := range exRegions {
				exGenes[k] = v
			}
		}
		if len(exGenesBySNP) != 0 {
			for k, v := range exGenesBySNP {
				exGenes[k] = v
			}
		}

		// fmt.Println(exGenes)

		if *annVCF == "list" || *annVCF == "*" || *annVCF == "all" {
			if *gbWeb == false && *annMakeSeq == "" {

				parserBulkVCF(*annWithFilenames)
				// go run bsatool.go annotate --vcf list --mkseq=NC  --db test_core

			} else if *gbWeb == true && strings.ToUpper(*annMakeSeq) == "NC" {
				createNCWebServer(*gbPort, exGenes, exSNPs)
			} else if *gbWeb == false && *annMakeSeq == "NC" && strings.ToUpper(*annOutFormat) != "NEXUS" {
				// seq := makeSeq(*annMakeSeq, *gbVerbose, *annMakeSeqRef)

				/* alt-range значит, что процент встречаемости альтернативного аллеля в последовательности будет в диапазоне --alt-range=10:80
								go run bsatool.go annotate --vcf list --mkseq=NC  --db test_core --debug --alt-range=10:80
				makeSeq
				*/

				seq := makeSeq(*annMakeSeq, *gbVerbose, *annMakeSeqRef, exGenes, exSNPs, *gbRandomize)
				// fmt.Println(len(seq))
				for i := 0; i < len(seq); i++ {
					fmt.Println(seq[i].Seq)
				}
				if *gbDebug == true {
					fmt.Println(seq[0].UsedPositions)
				}
			} else if *gbWeb == false && strings.ToUpper(*annMakeSeq) != "BINARY" && strings.ToUpper(*annOutFormat) == "NEXUS" {
				// fmt.Println(exGenes,exSNPs)
				makeSeqNex(*annMakeSeq, *gbVerbose, *annMakeSeqRef, exGenes, exSNPs, *gbRandomize)

			} else if *gbWeb == false && strings.ToUpper(*annMakeSeq) == "BINARY" {
				// go run bsatool.go annotate --vcf list --mkseq=NC  --db test_core
				// go run bsatool.go annotate --vcf list --mkseq=Binary  --db test_core
				// go run  bsatool.go annotate --db test_core --vcf list --exclude-genes=ppe.genes --exclude-snp=drugs3.
				// txt --mkseq=Binary --dp=30 --exclude-regions=exgenes.txt --min-pos-nbr=2 --output-format="nexus"

				makeSeqBinary(*annMakeSeq, *gbVerbose, *annMakeSeqRef, exGenes, exSNPs, *gbRandomize, *gbDP)
				// bsatool annotate --db sars --vcf list  --mkseq=Alignment --indel

			} else if *gbWeb == false && strings.ToUpper(*annMakeSeq) == "ALIGNMENT" {
				makeAlign(false, false)
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
					printTextResults(snpRes.SnpInfo[i], *gbIGR)
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
			err = pprof.StartCPUProfile(f)
			if err != nil {
				log.Fatal(err)
			}
			defer pprof.StopCPUProfile()
		}
		allGenesVal = readDB(*statDB)

		switch *statTask {
		case "share":
			/*
			 go run bsatool.go   stat -b test_core  -a share  -v
			*/

			getShareSNP(*gbVerbose, *gbIGR, *gbWeb, listOfFiles)

		case "snp":
			/*
				go run bsatool.go   stat -b test_core  -a snp  -v
			*/
			snpStat()
		case "dnds":

			/*
				go run bsatool.go   stat -b test_core  -a dnds  -i 161_RuU_m.vcf -v
			*/

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
			/*
				go run bsatool.go   stat -b test_core  -a ginfo
			*/
			testGeneInfo(allGenesVal)
		case "circos":
			/*
				go run bsatool.go   stat -b test_core  -a circos --typeof=[cds,igr,gene] --color=lblue --genome NC_000962
			*/
			toCircos(allGenesVal)
		case "bed":
			/*
				go run bsatool.go   stat -b test_core  -a bed
			*/

			toBED()
		case "matrix":

			/*
			 go run bsatool.go stat -a matrix --db test_core -t binary (binary, gc3, dnds, table, nc, locus. freq, jw, summary) -o test.csv
			 go run bsatool.go stat   -b test_core  -a matrix -t binary  --exclude-genes=ppe.genes --exclude-snp=drugs3.txt -o test.csv --nbr_pos 2 --debug|grep "passed"
			*/

			if *statTypeOf != "" && *statOutFile != "" {
				makeMatrix(*statTypeOf, *statOutFile, *gbVerbose)
			} else if *statTypeOf != "" && *statOutFile == "" {
				fmt.Println("Please, use key -o (--out) to save data.")
			}
		case "range":
			/*
				go run bsatool.go stat -a range --db test_core -i BWAIntersections.txt
				go run bsatool.go stat -a range --db test_core -i BWAIntersections.txt
			*/

			if *statInFile != "" {
				file := *statInFile
				res := getRangeFromFile(file, *gbVerbose, *gbNoSeq, *statCircosGenome)
				printSequenceRange(res, *gbWeb, *gbPort, *statCircos)

			}

			/*
				go run bsatool.go stat -a range --db test_core -i BWAIntersections.txt
			*/

		case "check_range":

		case "check":
			/*
				go run bsatool.go stat -a check --db test_core -i drugs2.txt -r rule.txt -w
				 go run bsatool.go stat --db test_core -a check -i drugs4.txt --table
				 go run bsatool.go stat --db test_core -a check -i drugs4.txt --mkseq
			*/

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
			/*
				go run bsatool.go stat -a check --db test_core -i drugs2.txt -r rule.txt -w
			*/

			if *statInRule != "" {
				rulesArr = checkRuleFromFile(*statInRule)

			}

		// создает fasta файл основываясь на снипах в границах указзаного диапазона
		case "coord2seq":
			/*  go run bsatool.go stat  -b test_core -a coord2seq --vcf=list -i coordsToMakeSeq.txt

			Rv0018c	21637	23181 <-structure file
			*/

			if *statInFile != "" && *statVCF == "list" {

				for _, file := range listOfFiles {
					result := coord2seq(*statInFile, file)
					fmt.Printf(">%v %v(%v:%v) %v\n%v\n", result.vcfFile, result.locus, result.start, result.end, result.prod, result.altSeq)
				}
			}

			// проверяет позицию из файла в снипах
		case "check_pos":

			if *statInFile != "" {
				// res := coord2gene(*statInFile)
				if *statVCF == "list" {

					/*
						go run bsatool.go  stat -b test_core -a check_pos -i uniq  --vcf=list --mkseq
						 go run bsatool.go  stat -b test_core -a check_pos -i uniq  --vcf=list --mkseq --reverse
					*/
					checkPosListVCF(*statInFile, *statMakeSeq)
				} else {

					checkPosList(*statInFile)
				}
			}

		case "snp2SeqByLocus":
			/*
				go run bsatool.go stat  -b test_core -a snp2SeqByLocus --vcf=list --locus=Rv1319c
			*/

			if *statLocus != "" && *statVCF == "list" {

				switch *statCircosTypeOf {
				case "":
					for _, file := range listOfFiles {
						result := snp2SeqByLocus(*statLocus, file)
						fmt.Printf(">%v %v(%v:%v) %v\n%v\n", result.vcfFile, result.locus, result.start, result.end, result.prod, result.altSeq)

					}
				case "nc":

					locus2Matrix(*statLocus, listOfFiles, "nc")
				case "nc_coded":

					locus2Matrix(*statLocus, listOfFiles, "nc_coded")
				case "binary":

					locus2Matrix(*statLocus, listOfFiles, "binary")
				}

			}
		case "countSNP":
			/*
				go run bsatool.go stat  -b test_core -a countSNP --vcf=list
			*/

			if *statVCF == "list" {

				for _, file := range listOfFiles {
					fmt.Printf("#---- %v ------#\n", file)
					countSNP(file)

				}
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
		allGenesVal = readDB(*betaDB)
		switch *betaTask {

		case "complex":

			if *betaInFile != "" {
				calcComplexIdxFromFile(*betaInFile)
			}
		case "annotateFromList":
			if len(*betaInFile) != 0 {
				annotateFromList(*betaInFile, allGenesVal)
			}
		}
	case "parse_pileup":
		allGenesVal = readDB(*pileupDB)
		if *pileupInFile != "" {
			pileup2multifasta(*pileupInFile, *pileupTH, *pileupGName, *pileupGStrain, *pileupMkSeq)
		}

	case "filter":
		if *gbDP != 0 {

			for _, file := range listOfFiles {
				filterVCF(file, *gbDP)
			}
		}

	case "info":
		/*
		 go run bsatool.go info --db test_core --locus=Rv0019c --showas=gene
		 go run bsatool.go info --db test_core --locus=Rv0019c --showas=direct
		 go run bsatool.go info --db test_core --locus=Rv0278c --codons=1:30
		 go run bsatool.go info --db test_core --locus=Rv0278c --range=1:30
		 go run bsatool.go info --db test_core --locus=Rv2629 --showas=fasta --mask=191:a:c
		 go run bsatool.go info --db test_core --locus=Rv3915 --showas=fasta --mask=1040:a:c  --iupac
		 go run bsatool.go info --db test_core --locus=Rv0294 --showas=fasta --mask=726:t:c  --iupac --flank_right=400 --flank_left=200

		*/
		allGenesVal = readDB(*infoDB)
		var (
			seq, direction              string
			asFasta                     bool
			flankSeqLeft, flankSeqRight string
			flankInfo                   strings.Builder
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
				case "fasta":
					asFasta = true
					if *infoFlankLeft != 0 {
						flankSeqLeft = getNucFromGenome((start-*infoFlankLeft)-1, end)
						flankInfo.WriteString(fmt.Sprintf("+%v_left ", *infoFlankLeft))
					}
					if *infoFlankRight != 0 {
						flankSeqRight = getNucFromGenome(end, end+*infoFlankRight)
						flankInfo.WriteString(fmt.Sprintf("+%v_right ", *infoFlankRight))
					}

					seq = getNucFromGenome(start-1, end)
					if *infoMask != "" {
						posRefAlt := strings.Split(*infoMask, ":")
						locPosInGene, _ := strconv.Atoi(posRefAlt[0])
						if locPosInGene != 0 && locPosInGene < len(seq) {
							splitSeq := strings.Split(seq, "")
							if strings.ToUpper(splitSeq[locPosInGene-1]) == strings.ToUpper(posRefAlt[1]) {
								if *infoIUPAc == false {

									splitSeq[locPosInGene-1] = fmt.Sprintf("[%v/%v]", posRefAlt[1], posRefAlt[2])
								} else {
									locUIPAC := getIUPAC(fmt.Sprintf("%v%v", posRefAlt[1], posRefAlt[2]))
									splitSeq[locPosInGene-1] = fmt.Sprintf("\n%v\n", locUIPAC)
								}

							} else {
								splitSeq[locPosInGene-1] = fmt.Sprintf("[%v/%v]", posRefAlt[1], posRefAlt[2])
								fmt.Println("Masked reference nucleotide does not match to reference!")
							}

							seq = strings.Join(splitSeq, "")

							if *infoFlankRight != 0 || *infoFlankLeft != 0 {
								seq = fmt.Sprintf("%v%v%v", flankSeqLeft, seq, flankSeqRight)
							}

						}
					}
					direction = "in forward direction"
				}

				prod := getProductByName(*infoLocus)
				note := getNoteByName(*infoLocus)
				goa := getGOAByName(*infoLocus)

				// вывод последовательности кодонов заданной ключом --codons=Start:End
				if *infoCodons != "" {

					seqByCodons := make([]string, (len(seq)/3)+1)
					codonNbr := 1
					var codonBuffer strings.Builder
					for _, c := range seq {
						codonBuffer.WriteString(string(c))
						if codonBuffer.Len() == 3 {

							seqByCodons[codonNbr] = codonBuffer.String()

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

							codonBuffer.WriteString(fmt.Sprintf("[%v]%v", i, seqByCodons[i]))
						}

					}

					fmt.Printf(">%v %v [%v:%v codons %v]\n%v\n", *infoLocus, prod, startCodon, endCodon, direction, codonBuffer.String())

				} else {

					if asFasta == false {
						gc, gc1, gc2, gc3 := codon.GcCodonCalc(seq)
						fmt.Printf(">%v %v (%v-%v %v)\n%v\n-----------------\ngc:%v gc1:%v gc2:%v gc3:%v Len:%v\n-----------------\n%v\ngoa:%v\n", *infoLocus,
							prod,
							start, end, direction, seq, gc, gc1, gc2, gc3, len(seq), note, goa)
					} else {

						fmt.Printf(">%v %v (%v-%v %v) %v \n%v\ngoa:%v\n", *infoLocus, prod, start, end, direction, flankInfo.String(), seq, goa)
					}

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
		codonVal                    string // переменная для кодона
		altCodon                    string // аналогично для альтернативного кодона
		mut                         string
		tangIdx                     string
		tangIdxVal                  int
	)
	// var trouble int
	lStart := g.Start // переменная начала гена
	lEnd := g.End
	posInGene := (apos - lStart) + 1             // позиция снипа в гене
	codonNbrInG := ((posInGene - 1) / 3) + 1     // номер кодона=номеру аминокислоты в трансляции
	posInCodonG := (codonNbrInG * 3) - posInGene // позиция в буквы в кодоне (0-первая, 1-средняя, 2-последняя)
	CPosInGene := (lEnd - apos) + 1              // комплементарный ген. позиция в гене
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
		codonVal = getNucFromGenome((posInGene+lStart)-2, ((posInGene+lStart)-1)+2)
	} else if posInCodonG == 1 {
		codonVal = getNucFromGenome((posInGene+lStart)-3, ((posInGene+lStart)-1)+1)
	} else if posInCodonG == 2 {
		codonVal = getNucFromGenome((posInGene+lStart)-4, (posInGene+lStart)-1)
	}

	/*
		ревертируем кодон для комплементарных генов

	*/
	nucG := getNucFromGenomePos((posInGene + lStart) - 1)
	typeOf = g.TypeOf

	// fmt.Println(apos, getNucFromGenomePos((posInGene+lStart)-1))

	// if apos == 2288836 {
	// 	fmt.Println(codon+"!!!", "genome:", nucG, "/", alt, codon, posInCodonG+1)
	// }
	if g.Direction == "r" {

		// alt = getComplement(alt)
		// nucG = getComplement(nucG)

		codonVal = getReverseComplement(codonVal)
		posInGene = CPosInGene
		codonNbrInG = CCodonNbrInG
		alt = getComplement(alt)
		nucG = getComplement(nucG)
		if posInCodonG == 2 {
			posInCodonG = 0
		} else if posInCodonG == 0 {
			posInCodonG = 2
		}

	}
	//

	codonPositions = strings.Split(codonVal, "")

	// if apos == 2288836 {
	// 	fmt.Println(strings.Join(codonPositions, ""), "!!!")
	// }

	codonPositions[posInCodonG] = strings.ToUpper(codonPositions[posInCodonG])

	codonVal = strings.Join(codonPositions, "")

	altCodonPositions = codonPositions

	altCodonPositions[posInCodonG] = alt

	altCodonPositions[posInCodonG] = strings.ToUpper(altCodonPositions[posInCodonG])

	altCodon = strings.Join(altCodonPositions, "")

	aaRef, aaRefShort := amino.Codon2AA(codonVal)
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
	// if apos == 2288836 {
	// 	fmt.Println(codon+"!!!", "reverse:", nucG, "/", alt, codon, posInCodonG+1, codon, altCodon, altCodonArr)
	// }

	// if strings.ToUpper(nucG) == strings.ToUpper(alt) {
	// 	fmt.Println("!!!", apos, nucG, alt)
	// }

	snp = snpInfo{APos: apos, PosInGene: posInGene, PosInCodonG: posInCodonG,
		RefCodon: codonVal, RefAA: aaRef, NucInPos: strings.ToUpper(nucG), Locus: g.Locus,
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

func geneBankFileParser(file string) (g []geneInfo, genomeSplice []string) {
	// функция для считывания файла в формате генбанка и занесения строк в массив linesFromGB

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
					scanTxt = strings.Replace(scanTxt, rfeat[1], fmt.Sprintf("*%v", rfeat[1]), -1)
				} else if rfeat[4] != "" {
					scanTxt = strings.Replace(scanTxt, rfeat[4], fmt.Sprintf("*%v", rfeat[4]), -1)

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
			changedStr = regDelSpaces.ReplaceAllString(scanTxt, "$2 ") // удаляем пробелы

			changedStr = strings.Replace(changedStr, "/note=", "!note=", -1) // меняем / на ! для дальнейшего парсинга

			changedStr = strings.Replace(changedStr, "/product=", "!product=", -1) // см выше.

			changedStr = strings.Replace(changedStr, "\"", "", -1)

			changedStr = makeAnchors.ReplaceAllString(changedStr, "!!")
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
	if gobFile != nil {
		defer gobFile.Close()
	}

	compressedGobFile := lzw.NewWriter(gobFile, lzw.LSB, 8)
	defer compressedGobFile.Close()
	gobParser := gob.NewEncoder(compressedGobFile)
	_ = gobParser.Encode(&gene)
	_ = gobParser.Encode(&genomeSeqSlice)
	_ = gobParser.Encode(&gInfo)
	_ = gobParser.Encode(&geneCoordinates)

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
	if gobFile != nil {
		defer gobFile.Close()
	}

	compressedGobFile := lzw.NewReader(gobFile, lzw.LSB, 8)
	defer compressedGobFile.Close()
	gobParser := gob.NewDecoder(compressedGobFile)
	_ = gobParser.Decode(&gene)
	_ = gobParser.Decode(&genomeSeqSlice)
	_ = gobParser.Decode(&gInfo)
	_ = gobParser.Decode(&geneCoordinates)

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

func parserVCF(f string, print bool, dpFilter int, genes []geneInfo) []snpInfo {
	var (
		vcf         = regexp.MustCompile(`^\S+\s+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)`)
		indel       = regexp.MustCompile(`^\S+\s+(\d+)\W+(\w+)\s+(\w+).*INDEL.*DP=(\d+)`)
		validateVCF = regexp.MustCompile(`(##fileformat)=VCF`)
		vcfValid    bool

		//
		snpFromVCF []snpInfo
		// locNbrSyn, locNbrNonSyn = 1, 1
		// locNbrNonsense          = 1
	)
	cpuprofile := *annBench
	if cpuprofile != "" {
		f, err := os.Create(cpuprofile)
		if err != nil {
			log.Fatal(err)

		}
		_ = pprof.StartCPUProfile(f)
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
			for _, matchindel := range indel.FindAllStringSubmatch(scanner.Text(), -1) {
				if vcfValid == true {
					apos, _ := strconv.Atoi(matchindel[1])
					ref := matchindel[2]
					alt := matchindel[3]
					dp, _ := strconv.Atoi(matchindel[4])
					// for _, g := range genes {
					for z := 0; z < len(genes); z++ {
						// g := genes[z]
						// lStart := genes[z].Start
						// lEnd := genes[z].End

						if apos >= genes[z].Start && apos <= genes[z].End && dp >= dpFilter {
							if *gbDebug == true {
								fmt.Printf("apos:%v\tref:%v\talt:%v\tdp: %v\t f:%v\n", apos, ref, alt, dp, f)
							}

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
							snp.DP = dp
							snp.Indel = 1
							snp.Alt = fmt.Sprintf("%v/%v", ref, alt)
							snp.IndelAlt = alt
							snp.IndelRef = ref
							if len(ref) > len(alt) {
								snp.IndelType = fmt.Sprintf("del%v", alt)

							} else if len(ref) < len(alt) {
								snp.IndelType = fmt.Sprintf("ins%v", alt)
							}
							// if snp.Mutation == "missense" {
							// 	locNbrNonSyn++
							// } else if snp.Mutation == "synonymous" {
							// 	locNbrSyn++
							// } else if snp.Mutation == "nonsense" {
							// 	locNbrNonsense++
							// }
							// snp.FileName = f
							// fmt.Println(g.PDB)
							// br := testing.Benchmark(snp)
							// fmt.Println(br)

							// if len(ref) == 1 && len(alt) == 1 {

							if *gbDebug == true {
								if strings.ToUpper(snp.NucInPos) == strings.ToUpper(snp.Alt) {
									snp.Product = fmt.Sprintf("%v|!ref=alt[%v]dp=%v", snp.Product, f, dp)
								}
								// fmt.Printf("apos:%v\tref:%v\talt:%v\tdp: %v\t f:%v\n", apos, ref, alt, dp, f)
							}

							if strings.ToUpper(snp.NucInPos) != strings.ToUpper(snp.Alt) {
								snpFromVCF = append(snpFromVCF, snp)
							}

							if print == true {
								// printResults(snp)
								printTextResults(snp, *gbIGR)
								// fmt.Println(snp)
							}
							// }

						}

					}
				}
			}
		}

		for _, match := range vcf.FindAllStringSubmatch(scanner.Text(), -1) {

			if vcfValid == true {
				apos, _ := strconv.Atoi(match[1])
				ref := match[2]
				alt := match[3]
				dp, _ := strconv.Atoi(match[4])

				// for _, g := range genes {
				for z := 0; z < len(genes); z++ {
					// g := genes[z]
					// lStart := genes[z].Start
					// lEnd := genes[z].End

					if apos >= genes[z].Start && apos <= genes[z].End && dp >= dpFilter {
						// if *gbDebug == true {
						// 	// fmt.Printf("apos:%v\tref:%v\talt:%v\tdp: %v\t f:%v\n", apos, ref, alt, dp, f)
						// }

						qSnpInfo := &snpInfoQuery{OutChan: make(chan snpInfo), apos: apos, g: genes[z], alt: alt, index: *gbIndex}
						go qSnpInfo.request()
						snp := <-qSnpInfo.OutChan
						// fmt.Println(snp.APos, ref, alt, snp.NucInPos, snp.Alt, snp.Direction, "!!!")
						// go qSNP.request()
						// // snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
						// snpRes := <-qSNP.OutChan
						// snpCacheMap[snpRes.File] = snpRes.SnpInfo

						// snp := getSnpInfo(apos, g, alt, *gbIndex)
						snp.InterPro = genes[z].InterPro
						snp.PDB = genes[z].PDB
						snp.ProSite = genes[z].ProSite
						snp.DP = dp
						snp.Indel = 0
						if *gbDebug == true {
							if strings.ToUpper(snp.NucInPos) == strings.ToUpper(snp.Alt) {
								snp.Product = fmt.Sprintf("%v|!ref=alt[%v]dp=%v", snp.Product, f, dp)
							}
							// fmt.Printf("apos:%v\tref:%v\talt:%v\tdp: %v\t f:%v\n", apos, ref, alt, dp, f)
						}

						// if snp.Mutation == "missense" {
						// 	locNbrNonSyn++
						// } else if snp.Mutation == "synonymous" {
						// 	locNbrSyn++
						// } else if snp.Mutation == "nonsense" {
						// 	locNbrNonsense++
						// }
						// snp.FileName = f
						// fmt.Println(g.PDB)
						// br := testing.Benchmark(snp)
						// fmt.Println(br)

						if len(ref) == 1 && len(alt) == 1 && strings.ToUpper(snp.NucInPos) != strings.ToUpper(snp.Alt) {
							snpFromVCF = append(snpFromVCF, snp)
							if print == true {
								// printResults(snp)
								printTextResults(snp, *gbIGR)
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

func filterVCF(f string, dpFilter int) {
	var (
		vcf = regexp.MustCompile(`^\S+\s+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)`)
		// var indel = regexp.MustCompile(`^^\S+\W+(\d+)\W+(\w+)\s+(\w+).*(INDEL).*DP=(\d+)`)
		validateVCF    = regexp.MustCompile(`(##fileformat)=VCF`)
		vcfSpecialInfo = regexp.MustCompile(`(##.*)`)
		// filterName     = regexp.MustCompile(`(.*).vcf`)
		vcfValid bool
		// outName  string

		//

	)

	file, err := os.Open(f)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	// for _, filtername := range filterName.FindAllStringSubmatch(f, -1) {
	// 	outName = fmt.Sprintf("%v_%v_f.vcf", filtername[1], *gbDP)
	// }

	fOut, err := os.Create("tmp")
	if err != nil {
		log.Fatal("Cannot create file", err)
	}
	defer fOut.Close()
	// var posArr []string

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {

		for _, vcfvalid := range validateVCF.FindAllStringSubmatch(scanner.Text(), -1) {

			if vcfvalid[0] == "##fileformat=VCF" {
				vcfValid = true
			}

		}

		for _, scpecInfo := range vcfSpecialInfo.FindAllStringSubmatch(scanner.Text(), -1) {

			_, _ = fmt.Fprintln(fOut, scpecInfo[0])

		}
		if vcfValid == false {
			fmt.Printf("\n%v is not VCF file!!! Check it!\n", file.Name())
			break
		}
		if *annInDel == true {

		}

		for _, match := range vcf.FindAllStringSubmatch(scanner.Text(), -1) {

			if vcfValid == true {

				dp, _ := strconv.Atoi(match[4])

				// for _, g := range genes {
				if dp >= dpFilter {
					_, _ = fmt.Fprintln(fOut, scanner.Text())
				}

			}
		}

	}

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}
	_ = os.Rename(f, fmt.Sprintf("%v.bak", f))
	_ = os.Rename("tmp", f)

}

// func about() {
// 	fmt.Println("\n", logo)
// }

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
		// if *annShowFileName == false {
		for _, val := range snps {
			val.FName = file
			printTextResults(val, *gbIGR)
		}

	}

}

func makeSeq(typeof string, verbose bool, ref bool, exGenes map[int]int, exSNPs map[int]int, randomize bool) []seqInfo {

	var (
		AllPos, SelectedPos, TempPos []int
		ResSeq                       []seqInfo
		// passSNP = make(map[string]int)
		uniqueSNP      = make(map[int]int)
		nbrOfSNP       int
		aaAltCoords    = make(map[string]map[int]string)
		aaRefCoords    = make(map[int]string)
		dpMAP          = make(map[int][]int)
		locus          = make(map[int]string)
		prod           = make(map[int]string)
		filesPOS       = make(map[int][]string)
		indel          = make(map[int]int)
		posFN          = make(map[int][]string)
		posFreq        = map[int][]string{}
		altPercent     int
		altMinMax      []string
		altMin, altMax int
		// nexusTaxa   = make(map[string][]string)
		// nChar       int
		// nTax        int

		// excludedLocus = make(map[int]string)
		// excludedProd  = make(map[int]string)
		// excludedDP    = make(map[int][]int)

		// posCount             = make(map[int]int)
	)

	// files := getListofVCF()

	// queryChan := make(chan vcfInfoQuery)

	if *annSeqLen != 0 {
		nbrOfSNP = *annSeqLen
	}

	files := &listOfFiles
	for i, file := range *files {
		aaAltCoords[file] = make(map[int]string)

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

	for fname, snps := range snpCacheMap {

		for _, val := range snps {
			// DP - глубина прочтения СНИПа
			dpMAP[val.APos] = append(dpMAP[val.APos], val.DP)
			// indel[val.APos] = val.Indel
			locus[val.APos] = val.Locus
			prod[val.APos] = val.Product
			filesPOS[val.APos] = append(filesPOS[val.APos], fname)
			indel[val.APos] = val.Indel
			if val.DP < *gbDP {
				uniqueSNP[val.APos] = 2
			}

			// exGenes  - гены, которые необходимо исключить из анализа. Загружаются в виде списка командой --exclude-genes
			if len(exGenes) != 0 {

				for key, value := range exGenes {
					if val.APos >= key && val.APos <= value {
						// значение 2- пропустить данный ген
						uniqueSNP[val.APos] = 2

						continue
					} else if exSNPs[val.APos] == 1 {
						/*
							exSNPs - снипы, которые исключаются из анализа. Загружаются командой --exclude-snp
							значение =2 пропускаем позицию
						*/
						uniqueSNP[val.APos] = 2
						continue

						// fmt.Println(val.APos)

					} else {
						if uniqueSNP[val.APos] != 2 && uniqueSNP[val.APos] != 1 && val.Indel == 0 {
							/*
								значение =1 , включаем позицию в формирование последовательности. Главное условие, что данная
								прозиция является не инделом
							*/
							uniqueSNP[val.APos] = 1
							aaAltCoords[fname][val.APos] = val.AltAAShort
							posFN[val.APos] = append(posFN[val.APos], fname)

						} else if val.Indel == 1 {
							/*
								если в позиции находится индел, то пропускаем
							*/
							uniqueSNP[val.APos] = 2
							continue
						}
					}
				}
			} else {
				if uniqueSNP[val.APos] != 2 && uniqueSNP[val.APos] != 1 && val.Indel == 0 {
					uniqueSNP[val.APos] = 1
					aaAltCoords[fname][val.APos] = val.AltAAShort
					aaRefCoords[val.APos] = val.RefAAShort

				}
				if uniqueSNP[val.APos] == 1 {
					posFN[val.APos] = append(posFN[val.APos], fname)
				}

			}

		}
	}

	// fmt.Println(uniqueSNP)

	for key, value := range uniqueSNP {

		if value == 1 {
			if len(dpMAP[key]) >= *annMinPosNbr {
				AllPos = append(AllPos, key)

			}

		}

	}
	// go process("Working...          ")
	// AllPos = unique(AllPosUnsort)

	sort.Ints(AllPos)
	// fmt.Println(AllPos)
	for _, allpos := range AllPos {
		for _, file := range *files {

			if strings.Contains(strings.Join(posFN[allpos], " "), file) {
				posFreq[allpos] = append(posFreq[allpos], "1")

			} else {
				// fmt.Println(allpos, 0, file)
				posFreq[allpos] = append(posFreq[allpos], "0")
			}

		}
	}
	// fmt.Println(posFreq)
	if nbrOfSNP == 0 {
		nbrOfSNP = len(AllPos) - 1
	}

	if nbrOfSNP > len(AllPos) {
		nbrOfSNP = len(AllPos) - 1
	}

	if randomize == true && nbrOfSNP != 0 {
		posAlreadyIs := make(map[int]int)
		rand.Seed(time.Now().UnixNano())
		for i := 1; i <= nbrOfSNP; i++ {
			rnd := rand.Intn(len(AllPos)-i) + i
			// fmt.Println(AllPos[rnd])
			if posAlreadyIs[rnd] < 1 {
				TempPos = append(TempPos, AllPos[rnd])
				posAlreadyIs[rnd] = posAlreadyIs[rnd] + 1
			}

		}
	} else if randomize == false {

		for i := 0; i <= nbrOfSNP; i++ {
			// rnd := rand.Intn(len(AllPos)-i) + i
			// fmt.Println(AllPos[rnd])

			TempPos = append(TempPos, AllPos[i])
		}
	}

	sort.Ints(TempPos)
	// fmt.Println(TempPos)

	for _, pos := range TempPos {
		var count0, count1 int

		for i := 0; i < len(posFreq[pos]); i++ {
			if posFreq[pos][i] == "0" {
				count0++
			} else if posFreq[pos][i] == "1" {
				count1++
			}
		}
		if count0 != 0 {
			// if *annMinPosNbr!=0 &&  {
			//
			// }
			altPercent = (count1 * 100) / len(posFreq[pos])
			altMinMax = strings.Split(*annAltRange, ":")

			if len(altMinMax) != 0 {
				altMin, _ = strconv.Atoi(altMinMax[0])
				altMax, _ = strconv.Atoi(altMinMax[1])
			}
			if altPercent >= altMin && altPercent <= altMax {
				SelectedPos = append(SelectedPos, pos)
				if *gbDebug {
					fmt.Printf("pos:%v isRef:%v,isAlt:%v %v legnth_array: %v AltPerc: %v \n", pos, count0, count1, posFreq[pos], len(posFreq[pos]), altPercent)
					// fmt.Println(*annTest)
				}
			}

		}

	}

	// -pos_file ФЛАГ

	sort.Ints(SelectedPos)

	// fmt.Println(SelectedPos)
	if *annPosFile != "" && typeof == ncFlag {
		fOut, err := os.Create(*annPosFile)
		fOutF, err := os.Create(fmt.Sprintf("%v_files.txt", *annPosFile))
		if err != nil {
			log.Fatal("Cannot create file", err)
		}
		defer fOut.Close()
		defer fOutF.Close()
		// var posArr []string
		_, _ = fmt.Fprintln(fOut, fmt.Sprintf("---NC method----"))
		_, _ = fmt.Fprintln(fOut, fmt.Sprintf("[pos]apos:indel:dp:loc:prod"))
		for i, value := range SelectedPos {

			_, _ = fmt.Fprintln(fOut, fmt.Sprintf("[%v]%v:%v:%v:%v:%v\t", i, value, indel[value], dpMAP[value], locus[value], prod[value]))
		}
		_, _ = fmt.Fprintln(fOutF, fmt.Sprintf("---NC method----"))
		_, _ = fmt.Fprintln(fOutF, fmt.Sprintf("[pos]apos:indel:file:loc:prod"))
		for i, value := range SelectedPos {

			_, _ = fmt.Fprintln(fOutF, fmt.Sprintf("[%v]%v:%v:%v:%v:%v\t", i, value, indel[value], filesPOS[value], locus[value], prod[value]))
		}

	}

	if ref == true {
		switch typeof {
		case ncFlag:
			var refBuffer strings.Builder
			refBuffer.WriteString(fmt.Sprintf(">%v\n", "REFERENCE"))
			for _, allpos := range SelectedPos {
				refBuffer.WriteString(getNucFromGenomePos(allpos))
			}
			ResSeq = append(ResSeq, seqInfo{Name: "reference", Seq: refBuffer.String(), UsedPositions: SelectedPos, TypeOfSeq: "NC"})

		}
	}

	// fmt.Println(AllPos)
	// fmt.Println(snpCacheMap)
	for fname, snps := range snpCacheMap {

		// if verbose == true {
		// 	fmt.Printf("Generating sequences: Working on  %v from %v \r", i+1, len(snps.File))
		// }
		pos := make(map[int]string)
		// altAA := make(map[int]string)
		// refAA := make(map[int]string)
		var buffer strings.Builder
		// var i int
		var aaSNPpos []int

		// var test strings.Builder

		buffer.WriteString(fmt.Sprintf(">%v\n", fname))
		switch typeof {
		case ncFlag:
			for _, val := range snps {
				pos[val.APos] = val.Alt
				// if exGenes[val.Locus] == 1 {
				// 	fmt.Println(val.Locus)
				// }
				// fmt.Println(pos)

			}
			// nTax = len(snpCacheMap)
			// nChar = len(SelectedPos)
			for _, allpos := range SelectedPos {
				// posCount[allpos] = posCount[allpos] + 1
				if pos[allpos] != "" {

					buffer.WriteString(pos[allpos])
					// nexusTaxa[fname] = append(nexusTaxa[fname], pos[allpos])

				} else {

					buffer.WriteString(getNucFromGenomePos(allpos))
					// nexusTaxa[fname] = append(nexusTaxa[fname], getNucFromGenomePos(allpos))
				}

				// 	var taxNbr int
				// 	// fmt.Println(nexusTaxa)
				// 	// fmt.Printf("#nexus \n\nBEGIN Taxa;\nDI
				//			}
				//
				//			// if len(nexusTaxa) != 0 {MENSIONS\nntax=%v;\nTAXLABELS\n", nTax)
				// 	for key := range nexusTaxa {
				// 		taxNbr++
				// 		// fmt.Printf("[%v]\t'%v'\n", taxNbr, key)
				// 	}
				//
				// 	// fmt.Printf(";\nEND; [Taxa]\n\nBEGIN Characters;\nDIMENSIONS nchar=%v;\nFORMAT\n\tdatatype=STANDARD\n\tmissing=?\n\tgap=-\n\tsymbols=\"01\"\n\t labels=left\n\ttranspose=no\n\tinterleave=no\n;\nMATRIX\n", nChar)
				// 	// for key, val := range nexusTaxa {
				// 	//
				// 	// 	fmt.Printf("'%v'\t%v\n", key, strings.Join(val, ""))
				// 	// }
				//
				// 	// fmt.Printf("\n;\nEnd;\n")
				// }

			}
			ResSeq = append(
				ResSeq, seqInfo{Name: fname, Seq: buffer.String(), UsedPositions: SelectedPos, TypeOfSeq: "AA"})
			// fmt.Println(ResSeq)
			// fmt.Println(buffer.Len())
			// fmt.Println(len(buffer.String()))

		case aaFlag:
			// for _, val := range snps {
			// 	pos[val.APos] = val.Alt
			// 	// if exGenes[val.Locus] == 1 {
			// 	// 	fmt.Println(val.Locus)
			// 	// }

			// }

			for _, allpos := range SelectedPos {
				if aaAltCoords[fname][allpos] == "" {
					buffer.WriteString(strings.ToLower(aaRefCoords[allpos]))
					aaSNPpos = append(aaSNPpos, allpos)
				} else if aaAltCoords[fname][allpos] != "" {
					aaSNPpos = append(aaSNPpos, allpos)
					buffer.WriteString(aaAltCoords[fname][allpos])
				}
				// posCount[allpos] = posCount[allpos] + 1
				// if pos[allpos] != "" {
				// 	fmt.Println(aaAltCoords[fname][allpos].alt)
				// 	// buffer.WriteString(pos[allpos])
				// } else {
				// 	fmt.Println(aaAltCoords[fname][allpos].ref)

				// 	// buffer.WriteString(getNucFromGenomePos(allpos))
				// }

			}
			// fmt.Println("reserved")
			// for _, val := range snps {
			// 	pos[val.APos] = val.Alt
			// 	altAA[val.APos] = val.AltAAShort
			// 	refAA[val.APos] = val.RefAAShort
			// 	i++
			// 	if i <= nbrOfSNP {
			// 		aaSNPpos = append(aaSNPpos, val.APos)
			// 	}

			// }
			// fmt.Println(allGenesVal)
			// for i := 0; i < len(SelectedPos); i++ {
			// fmt.Println(SelectedPos[i])
			// for _, g := range allGenesVal {
			// 	snp := getSnpInfo(SelectedPos[i], g, pos[SelectedPos[i]], false)
			// 	fmt.Println(snp.APos, snp.RefAAShort, snp.AltAAShort)
			// }
			// }

			// for _, allpos := range aaSNPpos {

			// 	for _, val := range geneCoordinates {
			// 		if allpos >= val.Start && allpos <= val.End {
			// 			if altAA[allpos] != "" {

			// 				buffer.WriteString(altAA[allpos])
			// 			} else {
			// 				buffer.WriteString(refAA[allpos])

			// 			}

			// 		}
			// 	}

			// 	// for z := 0; z < len(allGenesVal); z++ {
			// 	// if allpos >= allGenesVal[z].Start && allpos <= allGenesVal[z].End {
			// 	// 	fmt.Println("YES", allpos)
			// 	// 	// 		snp := getSnpInfo(allpos, allGenesVal[z], pos[allpos], false)
			// 	// 	// 		// fmt.Println(allpos, snp.RefAAShort)
			// 	// 	// 		buffer.WriteString(strings.ToLower(snp.RefAAShort))
			// 	// 	// 		// test.WriteString("," + strconv.Itoa(allpos))
			// 	// 	// 	}
			// 	// } else {
			// 	// 	fmt.Println("NO", allpos)
			// 	// }
			// 	// }
			// 	// posCount[allpos] = posCount[allpos] + 1

			// }

			// // fmt.Println(test.String())
			ResSeq = append(ResSeq, seqInfo{Name: fname, Seq: buffer.String(), UsedPositions: aaSNPpos})

		}
	}

	// }
	// if *gbDebug == true {
	// 	fmt.Printf("%v\t:\nThere was passed %v SNPs from exclude gene file\n And %v SNPs from exclude snp file\n", fname, passSNP["genes"], passSNP["snp"])
	// }
	// }
	// fmt.Println(ResSeq)
	return ResSeq
}

func makeSeqNex(typeof string, verbose bool, ref bool, exGenes map[int]int, exSNPs map[int]int, randomize bool) {

	var (
		AllPos, SelectedPos, TempPos []int
		// passSNP = make(map[string]int)
		uniqueSNP      = make(map[int]int)
		nbrOfSNP       int
		aaAltCoords    = make(map[string]map[int]string)
		aaRefCoords    = make(map[int]string)
		dpMAP          = make(map[int][]int)
		locus          = make(map[int]string)
		prod           = make(map[int]string)
		filesPOS       = make(map[int][]string)
		indel          = make(map[int]int)
		posFN          = make(map[int][]string)
		posFreq        = map[int][]string{}
		altPercent     int
		altMinMax      []string
		altMin, altMax int
		nexusTaxa      = make(map[string][]string)
		nChar          int
		nTax           int

		// excludedLocus = make(map[int]string)
		// excludedProd  = make(map[int]string)
		// excludedDP    = make(map[int][]int)

		// posCount             = make(map[int]int)
	)

	// files := getListofVCF()

	// queryChan := make(chan vcfInfoQuery)

	if *annSeqLen != 0 {
		nbrOfSNP = *annSeqLen
	}

	files := &listOfFiles
	for i, file := range *files {
		aaAltCoords[file] = make(map[int]string)

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

	for fname, snps := range snpCacheMap {

		for _, val := range snps {
			// DP - глубина прочтения СНИПа
			dpMAP[val.APos] = append(dpMAP[val.APos], val.DP)
			// indel[val.APos] = val.Indel
			locus[val.APos] = val.Locus
			prod[val.APos] = val.Product
			filesPOS[val.APos] = append(filesPOS[val.APos], fname)
			indel[val.APos] = val.Indel
			if val.DP < *gbDP {
				uniqueSNP[val.APos] = 2
			}

			// exGenes  - гены, которые необходимо исключить из анализа. Загружаются в виде списка командой --exclude-genes
			if len(exGenes) != 0 {

				for key, value := range exGenes {
					if val.APos >= key && val.APos <= value {
						// значение 2- пропустить данный ген
						uniqueSNP[val.APos] = 2

						continue
					} else if exSNPs[val.APos] == 1 {
						/*
							exSNPs - снипы, которые исключаются из анализа. Загружаются командой --exclude-snp
							значение =2 пропускаем позицию
						*/
						uniqueSNP[val.APos] = 2
						continue

						// fmt.Println(val.APos)

					} else {
						if uniqueSNP[val.APos] != 2 && uniqueSNP[val.APos] != 1 && val.Indel == 0 {
							/*
								значение =1 , включаем позицию в формирование последовательности. Главное условие, что данная
								прозиция является не инделом
							*/
							uniqueSNP[val.APos] = 1
							aaAltCoords[fname][val.APos] = val.AltAAShort
							posFN[val.APos] = append(posFN[val.APos], fname)

						} else if val.Indel == 1 {
							/*
								если в позиции находится индел, то пропускаем
							*/
							uniqueSNP[val.APos] = 2
							continue
						}
					}
				}
			} else {
				if uniqueSNP[val.APos] != 2 && uniqueSNP[val.APos] != 1 && val.Indel == 0 {
					uniqueSNP[val.APos] = 1
					aaAltCoords[fname][val.APos] = val.AltAAShort
					aaRefCoords[val.APos] = val.RefAAShort
				}
				if uniqueSNP[val.APos] == 1 {
					posFN[val.APos] = append(posFN[val.APos], fname)
				}

			}

		}
	}

	// fmt.Println(uniqueSNP)

	for key, value := range uniqueSNP {

		if value == 1 {
			if len(dpMAP[key]) >= *annMinPosNbr {
				AllPos = append(AllPos, key)

			}

		}

	}
	// go process("Working...          ")
	// AllPos = unique(AllPosUnsort)

	sort.Ints(AllPos)
	for _, allpos := range AllPos {

		for _, file := range *files {
			if strings.Contains(strings.Join(posFN[allpos], " "), file) {
				posFreq[allpos] = append(posFreq[allpos], "1")

			} else {
				// fmt.Println(allpos, 0, file)
				posFreq[allpos] = append(posFreq[allpos], "0")
			}

		}
	}

	if nbrOfSNP == 0 {
		nbrOfSNP = len(AllPos) - 1
	}

	if nbrOfSNP > len(AllPos) {
		nbrOfSNP = len(AllPos) - 1
	}

	if randomize == true && nbrOfSNP != 0 {
		posAlreadyIs := make(map[int]int)
		rand.Seed(time.Now().UnixNano())
		for i := 1; i <= nbrOfSNP; i++ {
			rnd := rand.Intn(len(AllPos)-i) + i
			// fmt.Println(AllPos[rnd])
			if posAlreadyIs[rnd] < 1 {
				TempPos = append(TempPos, AllPos[rnd])
				posAlreadyIs[rnd] = posAlreadyIs[rnd] + 1
			}

		}
	} else if randomize == false {

		for i := 0; i <= nbrOfSNP; i++ {
			// rnd := rand.Intn(len(AllPos)-i) + i
			// fmt.Println(AllPos[rnd])

			TempPos = append(TempPos, AllPos[i])
		}
	}

	sort.Ints(TempPos)
	// fmt.Println(posFreq)
	for _, pos := range TempPos {
		var count0, count1 int
		for i := 0; i < len(posFreq[pos]); i++ {
			if posFreq[pos][i] == "0" {
				count0++
			} else if posFreq[pos][i] == "1" {
				count1++
			}
		}
		if count0 != 0 {
			// if *annMinPosNbr!=0 &&  {
			//
			// }
			altPercent = (count1 * 100) / len(posFreq[pos])
			altMinMax = strings.Split(*annAltRange, ":")
			if len(altMinMax) != 0 {
				altMin, _ = strconv.Atoi(altMinMax[0])
				altMax, _ = strconv.Atoi(altMinMax[1])
			}
			// fmt.Println(altPercent,altMin,altMax)
			if altPercent >= altMin && altPercent <= altMax {
				SelectedPos = append(SelectedPos, pos)
				if *gbDebug {
					fmt.Printf("pos:%v isRef:%v,isAlt:%v %v legnth_array: %v AltPerc: %v \n", pos, count0, count1, posFreq[pos], len(posFreq[pos]), altPercent)
					// fmt.Println(*annTest)
				}
			}

		}

	}

	// -pos_file ФЛАГ

	sort.Ints(SelectedPos)
	for fname, snps := range snpCacheMap {

		pos := make(map[int]string)
		var buffer strings.Builder
		// var aaSNPpos []int

		buffer.WriteString(fmt.Sprintf(">%v\n", fname))
		switch typeof {
		case ncFlag:
			for _, val := range snps {
				pos[val.APos] = val.Alt

			}
			nTax = len(snpCacheMap)
			nChar = len(SelectedPos)
			for _, allpos := range SelectedPos {
				// posCount[allpos] = posCount[allpos] + 1
				if pos[allpos] != "" {

					// buffer.WriteString(pos[allpos])
					nexusTaxa[fname] = append(nexusTaxa[fname], pos[allpos])

				} else {

					// buffer.WriteString(getNucFromGenomePos(allpos))
					nexusTaxa[fname] = append(nexusTaxa[fname], getNucFromGenomePos(allpos))
				}

			}

		}

	}
	if len(nexusTaxa) != 0 {
		// var taxNbr int
		// fmt.Println(nexusTaxa)
		// fmt.Printf("#nexus\nBEGIN Taxa;\nDIMENSIONS\nntax=%v;\nTAXLABELS\n", nTax)
		fmt.Printf("#NEXUS\nBegin data;\nDimensions ntax=%v nchar=%v\nFormat datatype=dna missing=N gap=-;\n", nTax, nChar)
		// for key := range nexusTaxa {
		// 	taxNbr++
		// 	fmt.Printf("[%v]\t'%v'\n", taxNbr, key)
		// }
		fmt.Printf("matrix\n")
		for key, val := range nexusTaxa {

			fmt.Printf("%v  \n%v\n", key, strings.Join(val, ""))
		}
		fmt.Printf(";\nEnd;\n")
	}
	// fmt.Println(buffer.Len())
	// fmt.Println(len(buffer.String()))

	// 	case aaFlag:
	//
	// 		for _, allpos := range SelectedPos {
	// 			if aaAltCoords[fname][allpos] == "" {
	// 				buffer.WriteString(strings.ToLower(aaRefCoords[allpos]))
	// 				aaSNPpos = append(aaSNPpos, allpos)
	// 			} else if aaAltCoords[fname][allpos] != "" {
	// 				aaSNPpos = append(aaSNPpos, allpos)
	// 				buffer.WriteString(aaAltCoords[fname][allpos])
	// 			}
	//
	// 		ResSeq = append(ResSeq, seqInfo{Name: fname, Seq: buffer.String(), UsedPositions: aaSNPpos})
	//
	// 	}
	// }

	// return ResSeq
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

func printTextResults(snps snpInfo, igr bool) {

	const fullAnnotations = " {{if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T0\") (eq .Indel 0))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T1\")  (eq .Indel 0))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Tang}}({{.TangIdxVal}})\t{{.Product}}\n" +
		"{{else if (and (ne .TypeOf \"CDS\") (eq .ReportType \"T0\")  (eq .Indel 0))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.APos}}{{.NucInPos}}>{{.Alt}}\t-\t-\t-\t-\t{{.Product}}\n" +
		"{{else if (and (ne .TypeOf \"CDS\") (eq .ReportType \"T1\")  (eq .Indel 0))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.APos}}{{.NucInPos}}>{{.Alt}}\t-\t-\t-\t-\t-\t{{.Product}}\n" +
		"{{end}}"

	const fullAnnotationsWithName = " {{if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T0\"))}}" +
		"{{.FName}}\t{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T1\"))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Tang}}({{.TangIdxVal}})\t{{.Product}}\n" +
		"{{else if (and (ne .TypeOf \"CDS\") (eq .ReportType \"T0\"))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.APos}}{{.NucInPos}}>{{.Alt}}\t-\t-\t-\t-\t{{.Product}}\n" +
		"{{else if (and (ne .TypeOf \"CDS\") (eq .ReportType \"T1\"))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.APos}}{{.NucInPos}}>{{.Alt}}\t-\t-\t-\t-\t-\t{{.Product}}\n" +
		"{{end}}"

	const cdsAnnotations = "{{if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T0\")  (eq .Indel 0))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T1\")  (eq .Indel 0))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Tang}}({{.TangIdxVal}})\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T0\")  (eq .Indel 1))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.IndelType}}\t{{.Alt}}\t-\t-\t-\t-\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T1\")  (eq .Indel 1))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Tang}}({{.TangIdxVal}})\t{{.Product}}\t{{.IndelType}}\n" +
		"{{end}}"

	const cdsAnnotationsWithName = "{{if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T0\"))}}" +
		"{{.FName}}\t{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T1\"))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Tang}}({{.TangIdxVal}})\t{{.Product}}\n" +
		"{{end}}"

		// if *annHeader == true {
		// 	fmt.Printf("fname\tlocus\tapos\tpos_in_gene\tcodon_change\taa_change\tmutation\tproduct\t")
		// }

	t := template.New("report")
	switch *annShowFileName {
	case false:
		if igr == true {
			t, _ = t.Parse(fullAnnotations)
			_ = t.Execute(os.Stdout, snps)

		} else {
			t, _ = t.Parse(cdsAnnotations)
			_ = t.Execute(os.Stdout, snps)
		}
	case true:
		if igr == true {
			t, _ = t.Parse(fullAnnotationsWithName)
			_ = t.Execute(os.Stdout, snps)

		} else {

			t, _ = t.Parse(cdsAnnotationsWithName)
			_ = t.Execute(os.Stdout, snps)
		}
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

	locURL := fmt.Sprintf("localhost:%v/", port)
	_ = browser.OpenURL(locURL)

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

	// locPort := fmt.Sprintf(":%v", port)
	_ = http.ListenAndServe(":8080", nil)
}

func checkSNPnotation(inString string) snpCheckInfo {
	var (
		snpCheck snpCheckInfo
	)
	rLMN := regexp.MustCompile(lmnRegExp) // L:1 PiG:2 REF:3 ALT:4 NAME:5
	// LocusMutationName(LMN)
	rPMLN := regexp.MustCompile(pmlnRegExp) // APOS:1 REF:2 ALT:3 L:4 NAME:5
	// PositionMutationLocusName (PMLN)
	rPMN := regexp.MustCompile(pmnRegExp)     // APOS:1 REF:2 ALT:3 NAME:4
	rLSAAN := regexp.MustCompile(lsaanRegExp) // L:1 AA_REF:2 POS:3 AA_ALT:4 NAME:5
	// LocusShortAminoAcidName (LSAAN)
	rLLAAN := regexp.MustCompile(llaanRegExp) // L:1 LAA_REF:2 PiG:3 LAA_ALT:4 NAME:5
	// LocusLongAminoAcidName (LLAAN)
	rLCN := regexp.MustCompile(lcnRegExp)
	// rPMNT := regexp.MustCompile(pmntRegExp) //Pos(1)_Ref(2)>Alt(3){Tab}NAME(4):|;(5)Tag(6)
	rSEN := regexp.MustCompile(senRegExp) // start(1)|end(2):(name)

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
	for _, matchSEN := range rSEN.FindAllStringSubmatch(inString, -1) {
		startRange, err1 := strconv.Atoi(matchSEN[1])
		endRange, err2 := strconv.Atoi(matchSEN[2])
		if err1 == nil && err2 == nil {
			snpCheck = snpCheckInfo{StartRange: startRange, EndRange: endRange, Tag: matchSEN[3], TypeOf: tSEN}
			// fmt.Println(snpCheck)
		} else {
			fmt.Println(err1, err2)
		}

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

func getShareSNP(verbose bool, igr bool, web bool, files []string) {
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

			hash := getHash(fmt.Sprintf("%v%v%v%v%v%v%v", val.Locus, val.APos, val.PosInGene, val.Alt, val.Direction, val.TypeOf, val.Product))

			pos[hash] = pos[hash] + 1

			// fmt.Println(pos[val.APos+val.Start], val.Locus, upperLimit)
			// fmt.Println(val.APos, file, pos[val.APos])

			alt[hash] = val

			// fmt.Println(val.APos, hash)

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
			// if val.APos == 4284429 {
			// 	fmt.Println(hash)
			// }
		}

		for i, lPos := range pos {

			if lPos == upperLimit {

				share = append(share, i)

			}

		}

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

	//  --------------------------------------------
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

			printTextResults(res, igr)
		}
	}
	//  ---------------------------------------------

}

// func getSNPFromSnpInfo(snps []snpInfo, c chan int) {
// 	for _, val := range snps {
// 		c <- val.APos
// 	}
// }

func snpStat() {
	// var countSNPs = 1
	var (
		pos                  = make(map[int]int)
		alt                  = make(map[int]snpInfo)
		f                    = make(map[int][]string)
		positions, posUnsort []int
		// snpMask              = make(map[string]map[string][]int)
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
		// snpMask[fname]=make(map[string][]int)
		for _, val := range snps {
			if pos[val.APos] <= upperLimit {
				pos[val.APos] = pos[val.APos] + 1        // count
				alt[val.APos] = val                      // pos
				f[val.APos] = append(f[val.APos], fname) // files
				posUnsort = append(posUnsort, val.APos)  // array of positions

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
	var diffStr []string
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

	_ = browser.OpenURL(fmt.Sprintf("localhost:%v/", port))

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
	_ = http.ListenAndServe(locPort, nil)
}

func createNCWebServer(port string, exGenes map[int]int, exSNPs map[int]int) {
	/*

	 */

	seq := makeSeq(ncFlag, *gbVerbose, *annMakeSeqRef, exGenes, exSNPs, *gbRandomize)
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
	_ = browser.OpenURL(fmt.Sprintf("localhost:%v/", port))

	locPort := fmt.Sprintf(":%v", port)

	_ = http.ListenAndServe(locPort, nil)

	// fmt.Println("Для выхода из программы, нажмите Ctl+C")

}

func getListofVCF() []string {
	/*
		возвращает список VCF файлов в папке в виде массива
	*/
	// files, err := filepath.Glob(vcfExt)
	var (
		fList []string
	)
	files, err := ioutil.ReadDir(".")
	if err != nil {
		log.Fatal(err)
	}
	for _, f := range files {
		// fmt.Println(f.Name())
		if strings.ToUpper(path.Ext(f.Name())) == ".VCF" {
			fList = append(fList, f.Name())
		}

	}

	sort.Strings(fList)

	return fList
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

	isNotFound := make(map[string]int)

	snpFound := make(map[string]map[string]int)

	var (
		chkSNP        []checkSNP
		mutationsLits []string
	)

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
		snpFound[file] = make(map[string]int)

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
				// "LMN"   //locus:Mutation:NAME
				// "PMLN"  // position:Mutation:locus:NAME
				// "PMN"   // position:Mutation:NAME 			4326236_C>T	ETO_ethA_Gly413Asp
				// "LSAAN" //locus:shortAA:codon:shortAA:name
				// "LLAAN" // locus:longAA:codon:longAA:name
				// "LCN"   //locus:codon:name  -> Rv0667:codon491:RpoB_codon491
				case tLMN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && strings.ToUpper(val.Alt) == strings.ToUpper(snpFromFile.Alt) && lGpoS == snpFromFile.PosInGene {
						// fmt.Println(val.Locus, "\t", val.PosInGene, "\t", snpFromFile.PosInGene, "\t", snpFromFile.Locus)
						// chkSNP = checkSNP{FileName: file, FoundSNP: fmt.Sprintf("%v[%v:%v%v>%v]\t", val.Name, val.Locus, val.PosInGene, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt))}
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:%v%v>%v]", val.Name, val.Locus, val.PosInGene, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))
						// buffer.WriteString(fmt.Sprintf("%v_%v:%v%v>%v\t", val.Name, val.Locus, val.PosInGene, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))
						//  fmt.Println(strings.ToUpper(buffer.String()))
						mutationsLits = append(mutationsLits, fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name))
						snpFound[file][fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name)] = 1
					}
				case tPMLN:

					if lAPos == snpFromFile.APos && strings.ToUpper(val.Alt) == strings.ToUpper(snpFromFile.Alt) {
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:%v_%v>%v]", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))
						// chkSNP = checkSNP{FileName: file, FoundSNP: fmt.Sprintf("%v[%v:%v_%v>%v]\t", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt))}
						// snpArr = append(snpArr, checkSNP{FileName: file, FoundSNP: fmt.Sprintf("%v[%v:%v_%v>%v]\t", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt))})
						// buffer.WriteString(fmt.Sprintf("%v_%v:%v_%v>%v\t", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))
						mutationsLits = append(mutationsLits, fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name))
						snpFound[file][fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name)] = 1
					}
				case tLSAAN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && CodonNbrInG == snpFromFile.CodonNbrInG &&

						val.AASalt == snpFromFile.AltAAShort {
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:%v%v%v]", val.Name, val.Locus, val.AASref, CodonNbrInG, val.AASalt))
						mutationsLits = append(mutationsLits, fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name))
						snpFound[file][fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name)] = 1 // buffer.WriteString(fmt.Sprintf("%v_%v:%v%v%v\t", val.Name, val.Locus, val.AASref, CodonNbrInG, val.AASalt))

					}
				case tLLAAN:
					// if val.Locus == "Rv1908c" && snpFromFile.Locus == "Rv1908c" && CodonNbrInG == snpFromFile.CodonNbrInG {
					// 	fmt.Println(val.Locus, snpFromFile.Locus, CodonNbrInG, snpFromFile.CodonNbrInG, snpFromFile.RefAA, "->", val.AALalt, snpFromFile.AltAA)
					// }

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && CodonNbrInG == snpFromFile.CodonNbrInG &&
						val.AALalt == snpFromFile.AltAA {
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:%v%v%v]", val.Name, val.Locus, val.AALref, CodonNbrInG, val.AALalt))
						mutationsLits = append(mutationsLits, fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name))
						snpFound[file][fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name)] = 1 // buffer.WriteString(fmt.Sprintf("%v_%v:%v%v%v\t", val.Name, val.Locus, val.AALref, CodonNbrInG, val.AALalt))

					}
				case tLCN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && CodonNbrInG == snpFromFile.CodonNbrInG {
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:codon%v]", val.Name, val.Locus, CodonNbrInG))
						mutationsLits = append(mutationsLits, fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, snpFromFile.APos,
							strings.ToUpper(snpFromFile.NucInPos),
							strings.ToUpper(snpFromFile.Alt), val.Name))
						snpFound[file][fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, snpFromFile.APos,
							strings.ToUpper(snpFromFile.NucInPos),
							strings.ToUpper(snpFromFile.Alt), val.Name)] = 1 // buffer.WriteString(fmt.Sprintf("%v_%v:codon%v\t", val.Name, val.Locus, CodonNbrInG))
					}
				case tPMN:

					// "PMN"   // position:Mutation:NAME
					// 497491_C>T   N:RD105

					if lAPos == snpFromFile.APos && strings.ToUpper(val.Alt) == strings.ToUpper(snpFromFile.Alt) || lAPos == snpFromFile.APos && strings.ToUpper(getComplement(val.Alt)) == strings.ToUpper(snpFromFile.Alt) && snpFromFile.Direction == "r" {
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:%v_%v>%v]", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))
						mutationsLits = append(mutationsLits, fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name))
						snpFound[file][fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name)] = 1

					}
				case tSEN:
					if val.StartRange < snpFromFile.APos && val.EndRange > snpFromFile.APos {
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:%v_%v>%v]", snpFromFile.Locus, snpFromFile.Locus, snpFromFile.APos,
							strings.ToUpper(snpFromFile.NucInPos),
							strings.ToUpper(snpFromFile.Alt)))
						mutationsLits = append(mutationsLits, fmt.Sprintf("%v:%v_%v>%v[%v>%v](%v)", snpFromFile.Locus, snpFromFile.APos,
							strings.ToUpper(snpFromFile.NucInPos), strings.ToUpper(snpFromFile.Alt), snpFromFile.RefAAShort, snpFromFile.AltAAShort, val.Tag))
						// mutationsLits = append(mutationsLits, fmt.Sprintf("%v:%v_%v>%v_%v", fmt.Sprintf("%v[%v:%v_%v>%v]", snpFromFile.Locus, snpFromFile.Locus,
						// 	snpFromFile.APos, strings.ToUpper(snpFromFile.NucInPos),	strings.ToUpper(snpFromFile.Alt))))
						snpFound[file][fmt.Sprintf("%v:%v_%v>%v[%v>%v](%v)", snpFromFile.Locus, snpFromFile.APos, strings.ToUpper(snpFromFile.NucInPos),
							strings.ToUpper(snpFromFile.Alt), snpFromFile.RefAAShort, snpFromFile.AltAAShort,
							val.Tag)] = 1
					}
					// fmt.Println("!!!")
				}
			}

			isNotFound[file] = len(mapofSNP[file])

		}

		if isNotFound[file] == 0 {
			mapofSNP[file] = append(mapofSNP[file], "NO_SNPS_FOUND")

		}

	}

	if useRule == true {

		found := make(map[string][]string)
		foundName := make(map[string][]string)
		// fmt.Println(rulesArr)
		for fname, val := range mapofSNP {

			for _, tag := range rulesArr {
				// fmt.Println(tag.Name)

				for i := 0; i <= len(tag.Variants)-1; i++ {
					if strings.Contains(strings.ToUpper(strings.Join(val, ",")), tag.Variants[i]) == true {

						// fmt.Println(strings.Count(strings.Join(val, ","), tag.Variants[i]), found)
						found[fname] = appendIfMissing(found[fname], strings.ToUpper(tag.Variants[i]))

						if len(found[fname]) == tag.Lenght {

							foundName[fname] = appendIfMissing(foundName[fname], strings.ToUpper(tag.Name))

						}

					}
				}

			}

		}
		// for _, fname := range *files {

		// }
		for key, val := range mapofSNP {
			// sortNames := sort.Sort(foundName[key])
			sort.Slice(found[key], func(i, j int) bool {
				return found[key][i] < found[key][j]
			})

			chkSNP = append(chkSNP, checkSNP{FileName: key, FoundSNP: strings.Join(val, ","), ExactWithRule: strings.Join(found[key], ","), RuleNames: strings.Join(foundName[key], ",")})
		}

	} else {
		for key, val := range mapofSNP {
			chkSNP = append(chkSNP, checkSNP{FileName: key, FoundSNP: strings.Join(val, ",")})

		}
	}

	if *statAsTable == false {

		if web == true {

			printSNPfromFile(chkSNP, *gbPort)
		} else {

			for _, key := range chkSNP {

				fmt.Printf("%v\t%v\t%v\t%v\n", key.FileName, key.FoundSNP, key.ExactWithRule, key.RuleNames)

			}

		}
	} else {
		// удаляем дубликаты
		mutationsLits = removeStringDuplicates(mutationsLits)
		if len(mutationsLits) != 0 {
			// добавляем заголовок
			mutationsLits = append([]string{"file"}, mutationsLits...)
			results := make([]string, len(mutationsLits))
			// выводим заголовки имен мутаций
			fmt.Println(strings.Join(mutationsLits, "\t"))
			for _, file := range *files {
				for i := 0; i < len(mutationsLits); i++ {
					results[0] = file
					if snpFound[file][mutationsLits[i]] == 1 {
						results[i] = "1"
					} else if snpFound[file][mutationsLits[i]] == 0 {
						results[i] = "0"

					}
					// if i == 7 || i == 3 {
					// 	results[i] = "+"

					// } else {
					// 	results[i] = "-"
					// }
				}
				fmt.Println(strings.Join(results, "\t"))
			}
		} else {
			fmt.Println("No SNP found in range")
		}

	}
	// if *statAsBinary=true {

	// }

	// fmt.Println(snpFound)

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

	_ = browser.OpenURL(fmt.Sprintf("localhost:%v/", port))

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
	_ = http.ListenAndServe(locPort, nil)
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
		prod              string
		start, end        int
		gc, gc1, gc2, gc3 float64
	)
	for _, g := range genes {

		start = g.Start
		end = g.End

		if start != 0 && end != 0 && g.TypeOf == "CDS" {
			seq := getGeneSequence(g.Locus)
			prod = getProductByName(g.Locus)
			gc, gc1, gc2, gc3 = codon.GcCodonCalc(seq)
			fmt.Printf("%v:%v\t%v\t%v\t%.2f\t%.2f\t%.2f\t%.2f\n", start, end, g.Locus, prod, gc, gc1, gc2, gc3)

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
		if len(seqSplit) != 0 {
			seqSplit[val.pos-1] = val.alt
		}

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

func getNoteByName(locus string) string {
	var note string

	for _, g := range allGenesVal {
		if strings.ToUpper(locus) == strings.ToUpper(g.Locus) {
			note = g.Note
			break
		}
	}

	return note
}

func getGOAByName(locus string) string {
	var goa string

	for _, g := range allGenesVal {
		if strings.ToUpper(locus) == strings.ToUpper(g.Locus) {
			goa = g.GOA
			break
		}
	}

	return goa
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
	if *statCircosTypeOf == "" || *statCircosTypeOf == "cds" || *statCircosTypeOf == "igr" {
		fmt.Printf("chr -  %v caption %v %v %v\n", genomeName, gInfo.Start, gInfo.End, gInfo.Strain)
	}
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
			case "gene":
				if g.TypeOf == "CDS" {
					fmt.Printf("%v %v %v %v\n", genomeName, start, end, g.Locus)
				}
			case "all":
				fmt.Printf("%v %v %v %v\n", genomeName, start, end, g.Locus)

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

func getGeneNameCoordsByApos(pos int) (int, int) {
	var (
		gStart, gEnd int
	)

	for _, g := range allGenesVal {
		if pos >= g.Start && pos <= g.End {
			gStart = g.Start
			gEnd = g.End

		}

	}

	return gStart, gEnd
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

func nuc2IntCode(nuc string) string {

	var (
		locNuc = strings.ToUpper(nuc)
		res    string
	)

	switch locNuc {

	}

	switch locNuc {

	case "A":
		res = "0"

	case "T":
		res = "1"
	case "G":
		res = "2"
	case "C":

		res = "3"
	}
	return res
}

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

	for key, val := range geneCoordinates {
		if val.Type == "CDS" {
			allLocuses = append(allLocuses, key)
		}
	}
	sort.Strings(allLocuses)

	var (
		exGenes = make(map[int]int)
		exSNP   = make(map[int]int)
	)
	if *gbExcludeGenes != "" {
		exGenes = loadExcludeGenes(*gbExcludeGenes)
	} else if *gbExcludeSnp != "" {
		exSNP = loadExcludeSNP(*gbExcludeSnp)
	}

	switch typeof {

	case "binary":

		matrixBinary(fileOut, exGenes, exSNP)

	// case "table":

	// 	matrixTable(fileOut)

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
		AllPos = getAllPosFromCacheMap(exGenes, exSNP)
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
	case "nc_coded":
		// A, T, G,C, ---> 0, 1, 2, 3
		var posFreq = map[int][]string{}
		pos := make(map[int]string)
		posCount := make(map[int]int)
		files := &listOfFiles
		headers.WriteString("Pos\t")
		// allPosChan := make(chan []int)
		// go func() {
		// 	allPosChan <- getAllPosFromCacheMap()
		// }()
		// AllPos = <-allPosChan
		// for i, file := range files {
		AllPos = getAllPosFromCacheMap(exGenes, exSNP)
		for fname, snps := range snpCacheMap {
			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

			// snps := parserVCF(file, false, allGenesVal)

			for _, val := range snps {

				pos[val.APos] = val.Alt
				// AllPosUnsort = append(AllPosUnsort, val.APos)
				posCount[val.APos] = posCount[val.APos] + 1
			}
			for _, allpos := range AllPos {

				if pos[allpos] != "" {

					posFreq[allpos] = append(posFreq[allpos], fmt.Sprintf("%v", nuc2IntCode(pos[allpos])))
					// fmt.Println(posFreq[allpos])

				} else {
					posFreq[allpos] = append(posFreq[allpos], nuc2IntCode(getNucFromGenomePos(allpos)))

				}

			}

		}
		// AllPos = unique(AllPosUnsort)
		// // allLocuses = removeStringDuplicates(allLocusUnsort)
		// sort.Ints(AllPos)

		for _, allpos := range AllPos {
			// if buffer.Len() == 0 {
			if posCount[allpos] <= len(*files)-1 {
				buffer.WriteString(fmt.Sprintln(allpos, "\t", strings.Join(posFreq[allpos], "\t")))
				// fmt.Println(buffer.String())
			}

		}
		headers.WriteString("\n")
		matrixPrint(headers, buffer, fileOut)

	case "summary":
		var (
			posFreq = map[int][]string{}
			// pos                    = make(map[int]string)
			groupRegexp            = regexp.MustCompile(`^(\S+)\W+(\w+)\W+(\w+)`)
			group                  = make(map[string]string)
			label                  = make(map[string]string)
			groupHeader, groupBody string
			fileGroup              = make(map[int][]string)
			fileLabel              = make(map[int][]string)
			PosInGenome            = make(map[string]map[int]string)
			PosAdditionalInfo      = make(map[int]string)
			genomes                []string
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
					group[strings.ToUpper(grpVal[1])] = strings.ToUpper(grpVal[2])

					label[strings.ToUpper(grpVal[1])] = grpVal[3]

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

		headers.WriteString(fmt.Sprintf("Pos\tLocus\t%vGene\tAAchange\tDirection\tMutation\tCIndex\tCIndex_Res\t", groupHeader))
		// allPosChan := make(chan []int)
		// go func() {
		// 	allPosChan <- getAllPosFromCacheMap()
		// }()
		// AllPos = <-allPosChan
		// for i, file := range files {
		AllPos = getAllPosFromCacheMap(exGenes, exSNP)

		for fname, snps := range snpCacheMap {
			PosInGenome[fname] = make(map[int]string)
			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))
			genomes = append(genomes, fname)
			// snps := parserVCF(file, false, allGenesVal)

			for _, val := range snps {

				// pos[val.APos] = val.Alt

				PosInGenome[fname][val.APos] = fmt.Sprintf("%v/%v", val.NucInPos, val.Alt)
				idxVal, idxRes := amino.GetComplexIndex(val.RefAAShort, val.AltAAShort, false)
				// fmt.Println(idxVal, idxRes)
				if val.TypeOf == "CDS" {
					// fmt.Printf("%v\t%v%v%v\t%v\t%v\t%v\t%v\t", val.Name, val.RefAA, val.CodonNbrInG, val.AltAA, val.Direction,
					// 	val.Mutation, idxVal, idxRes)
					// fmt.Println(idxVal)
					PosAdditionalInfo[val.APos] = fmt.Sprintf("%v\t%v%v%v\t%v\t%v\t%v\t%v\t", val.Name, val.RefAA, val.CodonNbrInG, val.AltAA, val.Direction,
						val.Mutation, idxVal, idxRes)
				} else {
					PosAdditionalInfo[val.APos] = fmt.Sprintf("%v\t%v%v%v\t%v\t%v\t%v\t%v\t", val.Name, "-", "-", "-", val.Direction, "-", "-", "-")
				}
				// AllPosUnsort = append(AllPosUnsort, val.APos)

			}

			// for _, allpos := range AllPos {

			// 	if pos[allpos] != "" {
			// 		posFreq[allpos] = append(posFreq[allpos], fmt.Sprintf("%v/%v", getNucFromGenomePos(allpos), pos[allpos]))
			// 		// fmt.Println(posFreq[allpos])
			// 		if len(group) != 0 {
			// 			fileGroup[allpos] = append(fileGroup[allpos], group[fname])
			// 			if len(label) != 0 {
			// 				fileLabel[allpos] = append(fileLabel[allpos], label[fname])
			// 			}
			// 			// groupArr = append(groupArr, group[fname])
			// 			// groupBody = fmt.Sprintf("\t%v\t", group[fname])
			// 		} else {
			// 			groupBody = "\t"
			// 		}

			// 	} else {
			// 		posFreq[allpos] = append(posFreq[allpos], ".")

			// 	}

			// }

		}
		// AllPos = unique(AllPosUnsort)
		// // allLocuses = removeStringDuplicates(allLocusUnsort)
		// sort.Ints(AllPos)
		// fmt.Println(PosInGenome)

		for _, gname := range genomes {
			for _, allpos := range AllPos {
				mut := PosInGenome[gname][allpos]
				if mut != "" {
					posFreq[allpos] = append(posFreq[allpos], mut)
					if len(group) != 0 {
						fileGroup[allpos] = append(fileGroup[allpos], group[strings.ToUpper(gname)])

						if len(label) != 0 {
							fileLabel[allpos] = append(fileLabel[allpos], label[strings.ToUpper(gname)])
						}
						// fmt.Printf("%v\t%v\t%v\n", allpos, group[strings.ToUpper(gname)], strings.ToUpper(gname))

					} else {
						groupBody = "\t"
					}
				} else {
					posFreq[allpos] = append(posFreq[allpos], ".")
				}

			}

		}

		for _, allpos := range AllPos {
			// if buffer.Len() == 0 {
			gname, _ := getGeneNameByPos(allpos, allpos)
			if len(group) != 0 {
				sort.Strings(fileGroup[allpos])
				uniq := removeStringDuplicates(fileGroup[allpos])
				if len(uniq) == 1 {
					if uniq[0] == "" {
						uniq = append(uniq, "-")
					}
				}
				groupBody = fmt.Sprintf("\t%v\t", strings.Trim(strings.Join(uniq, " "), " "))

			}

			if len(label) != 0 {
				sort.Strings(fileLabel[allpos])
				uniq := removeStringDuplicates(fileLabel[allpos])

				groupBody = fmt.Sprintf("%v%v\t", groupBody, strings.Trim(strings.Join(uniq, " "), " "))

			}

			buffer.WriteString(fmt.Sprintln(allpos, "\t", gname, groupBody, PosAdditionalInfo[allpos], strings.Join(posFreq[allpos], "\t")))
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
		AllPos = getAllPosFromCacheMap(exGenes, exSNP)
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
		// stat --db test_core -a matrix -t mst -o test.tsv --exclude-genes=exgenes.txt --exclude-snp=drugs2.txt --snp-number=15
		var (
			exGenes = make(map[int]int)
			exSNP   = make(map[int]int)
		)
		if *gbExcludeGenes != "" {
			exGenes = loadExcludeGenes(*gbExcludeGenes)
		} else if *gbExcludeSnp != "" {
			exSNP = loadExcludeSNP(*gbExcludeSnp)
		}

		matrixBinaryMST(fileOut, *gbVerbose, exGenes, exSNP, *gbRandomize)

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

			// if len(locJW[allloc]) != 0 {

			// fmt.Println(len(*files)-len(locJW[allloc]), locJW[allloc], strings.Repeat(" 1 ", len(*files)-len(locJW[allloc])))
			rpt := strings.Repeat("1\t", len(*files)-len(locJW[allloc]))
			locJW[allloc] = append(locJW[allloc], rpt)
			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v\t", allloc), strings.Join(locJW[allloc], "\t")))

			// }
		}

		headers.WriteString("\n")
		matrixPrint(headers, buffer, fileOut)

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
		matrixPrint(headers, buffer, fileOut)

	}

}

func getAllPosFromCacheMap(exGenes map[int]int, exSNPs map[int]int) (allPos []int) {
	var (
		AllPosUnsort []int
		uniqSNP      = make(map[int]int)
		passedSNP    int
	)
	if len(snpCacheMap) != 0 {
		for _, snps := range snpCacheMap {

			for _, val := range snps {

				if len(exGenes) != 0 {

					for key, value := range exGenes {
						if val.APos >= key && val.APos <= value {
							// fmt.Println(val.Locus, val.Start, val.End, val.Product)
							// passSNP["genes"] = passSNP["genes"] + 1
							uniqSNP[val.APos] = 2 // 2-EXCLUDED
							// countPassSNPinGenes = countPassSNPinGenes + 1
							continue
						} else if exSNPs[val.APos] == 1 {
							// passSNP["snp"] = passSNP["snp"] + 1
							uniqSNP[val.APos] = 2 // 2-EXCLUDED
							// countPassSNPfromSNPList = countPassSNPfromSNPList + 1
							continue

							// fmt.Println(val.APos)

						} else {
							if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
								uniqSNP[val.APos] = 1 // 1-INCLUDED

							}
						}
					}
				} else {
					if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
						uniqSNP[val.APos] = 1

					}
					// AllPosUnsort = append(AllPosUnsort, val.APos)
				}

				if uniqSNP[val.APos] == 1 {
					AllPosUnsort = append(AllPosUnsort, val.APos)
				} else {
					passedSNP = passedSNP + 1
				}
			}

		}
		allPos = unique(AllPosUnsort)
		// allLocuses = removeStringDuplicates(allLocusUnsort)
		sort.Ints(allPos)

	}
	if *gbDebug {
		fmt.Println("passed : ", passedSNP, " snps")
	}
	return allPos
}

// func matrixTable(fileOut string) {
// 	var (
// 		AllPosUnsort, AllPos []int
// 		allLocusUnsort       []string
// 		buffer               strings.Builder
// 		headers              strings.Builder
// 		posCount             = make(map[int]int)
// 		// snps                 []snpInfo
// 		posFN   = make(map[int][]string)
// 		posFreq = map[int][]string{}
// 		i       = 0
// 	// var ResSeq []seqInfo
// 	)
// 	// files := getListofVCF()
// 	files := &listOfFiles

// 	// fmt.Println(files)
// 	pos := make(map[int]string)

// 	headers.WriteString("Pos\t")

// 	for key, snps := range snpCacheMap {
// 		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(key, filepath.Ext(key))))

// 		fmt.Printf("Counting SNP positions: Working on %v files from %v \r", i+1, len(*files))

// 		snpsChan := make(chan []snpInfo)

// 		go func() {
// 			snpsChan <- makeSnps(key)
// 		}()
// 		snps = <-snpsChan

// 		for _, val := range snps {

// 			pos[val.APos] = val.Alt
// 			// posTest = append(posTest, posByFile{pos: pos, file: file, apos: val.APos})
// 			if strings.Contains(strings.Join(posFN[val.APos], " "), key) == false {
// 				posFN[val.APos] = append(posFN[val.APos], key)
// 			}

// 			AllPosUnsort = append(AllPosUnsort, val.APos)
// 			posCount[val.APos] = posCount[val.APos] + 1
// 			if val.TypeOf == "CDS" {
// 				allLocusUnsort = append(allLocusUnsort, val.Locus)
// 			}

// 		}

// 	}
// 	AllPos = unique(AllPosUnsort)
// 	// allLocuses = removeStringDuplicates(allLocusUnsort)
// 	sort.Ints(AllPos)

// 	for _, file := range *files {
// 		for _, allpos := range AllPos {

// 			if strings.Contains(strings.Join(posFN[allpos], " "), file) {

// 				posFreq[allpos] = append(posFreq[allpos], "1")
// 			} else {
// 				// fmt.Println(allpos, 0, file)
// 				posFreq[allpos] = append(posFreq[allpos], "0")
// 			}

// 		}
// 	}

// 	for _, allpos := range AllPos {
// 		if posCount[allpos] < len(*files) {

// 			buffer.WriteString(fmt.Sprintln(allpos, "\t", strings.Join(posFreq[allpos], "\t")))

// 		}
// 	}
// 	headers.WriteString("\n")

// 	// if buffer.Len() != 0 && headers.Len() != 0 {
// 	// 	fOut, err := os.Create(fileOut)
// 	// 	if err != nil {
// 	// 		log.Fatal("Cannot create file", err)
// 	// 	}
// 	// 	defer fOut.Close()
// 	// 	fmt.Fprintf(fOut, headers.String())
// 	// 	fmt.Fprintf(fOut, buffer.String())
// 	// 	fmt.Printf("\n\nWell done!\n")
// 	// t1 := time.Now()
// 	// fmt.Printf("Elapsed time: %v", fmtDuration(t1.Sub(t0)))
// 	// }
// 	matrixPrint(headers, buffer, fileOut)
// }

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
		_, _ = fmt.Fprintf(fOut, headers.String())
		_, _ = fmt.Fprintf(fOut, buffer.String())
		fmt.Printf("\n\nWell done!\n")
	}
}

// func matrixDnDsOld(fileOut string) {
// 	// var AllPos []int
// 	var (
// 		allLocuses []string

// 		buffer  strings.Builder
// 		headers strings.Builder
// 		// var posCount = make(map[int]int)
// 		snps         []snpInfo
// 		altPositions = make(map[string][]allPositionsInGene)
// 		locDNDS      = map[string][]string{}
// 		countNbrOne  = make(map[string]int)
// 	)
// 	t0 := time.Now()

// 	files := &listOfFiles

// 	i := 0
// 	headers.WriteString("Locus\t")

// 	for key := range geneCoordinates {
// 		allLocuses = append(allLocuses, key)
// 	}

// 	for _, fname := range *files {
// 		i++

// 		snpsChan := make(chan []snpInfo)

// 		go func() {
// 			snpsChan <- makeSnps(fname)
// 		}()
// 		snps = <-snpsChan

// 		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

// 		t1 := time.Now()

// 		fmt.Printf("Calculating Dn/DS: Working on %v from %v (%v)\r", i, len(*files), fname)
// 		for _, val := range snps {
// 			if val.TypeOf == "CDS" && containsPos(altPositions[val.Locus], val.PosInGene, val.Alt) == false {

// 				altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, alt: val.Alt, ref: val.NucInPos, locus: val.Locus})
// 			}

// 		}

// 		for _, allloc := range allLocuses {

// 			// prod := getProductByName(allloc)
// 			if len(altPositions[allloc]) > 2 {
// 				dndsChan := make(chan []string)
// 				go func() {

// 					dndsChan <- getDnDsByLocus(allloc, altPositions[allloc])
// 				}()
// 				dndsRes, ok := <-dndsChan
// 				if ok {
// 					// fmt.Println(dndsRes)

// 					locDNDS[dndsRes[0]] = append(locDNDS[dndsRes[0]], dndsRes[1])
// 					if dndsRes[1] == "1" {
// 						countNbrOne[allloc]++
// 					}
// 					close(dndsChan)
// 				}

// 			} else {

// 				locDNDS[allloc] = append(locDNDS[allloc], "1")
// 				countNbrOne[allloc]++
// 			}

// 		}

// 		fmt.Printf("Calculating Dn/DS: Working on %v from %v (%v) \t\t Time:\t%v\n", i, len(*files), fname, t1.Sub(t0))

// 	}

// 	// fmt.Println(locDNDS)

// 	for _, allloc := range allLocuses {

// 		if countNbrOne[allloc] != len(*files) && *statAll == false {
// 			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
// 		} else if *statAll == true {
// 			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
// 		}

// 	}
// 	// }

// 	headers.WriteString("\n")

// 	if buffer.Len() != 0 && headers.Len() != 0 {
// 		fOut, err := os.Create(fileOut)
// 		if err != nil {
// 			log.Fatal("Cannot create file", err)
// 		}
// 		defer fOut.Close()
// 		fmt.Fprintf(fOut, headers.String())
// 		fmt.Fprintf(fOut, buffer.String())
// 		fmt.Printf("\n\nWell done!\n")
// 		// t1 := time.Now()
// 		// fmt.Printf("Elapsed time: %v", fmtDuration(t1.Sub(t0)))
// 	}
// }

// func matrixDnDs(fileOut string) {
// 	// var AllPos []int
// 	var (
// 		allLocuses []string

// 		buffer  strings.Builder
// 		headers strings.Builder
// 		// var posCount = make(map[int]int)
// 		// snps         []snpInfo
// 		altPositions = make(map[string][]allPositionsInGene)
// 		locDNDS      = map[string][]string{}
// 		countNbrOne  = make(map[string]int)
// 	)
// 	t0 := time.Now()

// 	files := &listOfFiles

// 	i := 1
// 	headers.WriteString("Locus\t")

// 	for key := range geneCoordinates {
// 		allLocuses = append(allLocuses, key)
// 	}

// 	for key, snps := range snpCacheMap {

// 		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(key, filepath.Ext(key))))

// 		t1 := time.Now()

// 		fmt.Printf("Calculating Dn/DS: Working on %v from %v (%v)\r", i, len(*files), key)
// 		for _, val := range snps {
// 			if val.TypeOf == "CDS" && containsPos(altPositions[val.Locus], val.PosInGene, val.Alt) == false {

// 				altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, alt: val.Alt, ref: val.NucInPos, locus: val.Locus})
// 			}

// 		}

// 		for _, allloc := range allLocuses {

// 			// prod := getProductByName(allloc)
// 			if len(altPositions[allloc]) > 2 {
// 				dndsChan := make(chan []string)
// 				go func() {

// 					dndsChan <- getDnDsByLocus(allloc, altPositions[allloc])
// 				}()
// 				dndsRes, ok := <-dndsChan
// 				if ok {
// 					// fmt.Println(dndsRes)

// 					locDNDS[dndsRes[0]] = append(locDNDS[dndsRes[0]], dndsRes[1])
// 					if dndsRes[1] == "1" {
// 						countNbrOne[allloc]++
// 					}
// 					close(dndsChan)
// 				}

// 			} else {

// 				locDNDS[allloc] = append(locDNDS[allloc], "1")
// 				countNbrOne[allloc]++
// 			}

// 		}

// 		fmt.Printf("Calculating Dn/DS: Working on %v from %v (%v) \t\t Time:\t%v\n", i, len(*files), key, t1.Sub(t0))
// 		i++
// 	}

// 	// fmt.Println(locDNDS)

// 	for _, allloc := range allLocuses {

// 		if countNbrOne[allloc] != len(*files) && *statAll == false {
// 			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
// 		} else if *statAll == true {
// 			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
// 		}

// 	}
// 	// }

// 	headers.WriteString("\n")

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
// 	matrixPrint(headers, buffer, fileOut)
// }

// func matrixBinary(fileOut string) {
// 	var (
// 		AllPosUnsort, AllPos []int
// 		allLocusUnsort       []string
// 		buffer               strings.Builder
// 		headers              strings.Builder
// 		posCount             = make(map[int]int)
// 		// snps                 []snpInfo
// 		posFN   = make(map[int][]string)
// 		posFreq = map[int][]string{}
// 	)
// 	// var ResSeq []seqInfo

// 	// files := getListofVCF()
// 	files := &listOfFiles

// 	// fmt.Println(files)
// 	pos := make(map[int]string)

// 	headers.WriteString("Pos\t")
// 	i := 1
// 	for key, snps := range snpCacheMap {
// 		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(key, filepath.Ext(key))))

// 		fmt.Printf("Counting SNP positions: Working on %v files from %v \r", i+1, len(*files))

// 		for _, val := range snps {

// 			pos[val.APos] = val.Alt
// 			// posTest = append(posTest, posByFile{pos: pos, file: file, apos: val.APos})
// 			if strings.Contains(strings.Join(posFN[val.APos], " "), key) == false {
// 				posFN[val.APos] = append(posFN[val.APos], key)
// 			}

// 			AllPosUnsort = append(AllPosUnsort, val.APos)
// 			posCount[val.APos] = posCount[val.APos] + 1
// 			if val.TypeOf == "CDS" {
// 				allLocusUnsort = append(allLocusUnsort, val.Locus)
// 			}

// 		}

// 	}
// 	AllPos = unique(AllPosUnsort)
// 	// allLocuses = removeStringDuplicates(allLocusUnsort)
// 	sort.Ints(AllPos)

// 	for _, file := range *files {
// 		for _, allpos := range AllPos {

// 			if strings.Contains(strings.Join(posFN[allpos], " "), file) {

// 				posFreq[allpos] = append(posFreq[allpos], "1")

// 			} else {
// 				// fmt.Println(allpos, 0, file)
// 				posFreq[allpos] = append(posFreq[allpos], "0")
// 			}

// 		}
// 	}

// 	for _, allpos := range AllPos {
// 		if posCount[allpos] < len(*files) {

// 			buffer.WriteString(fmt.Sprintln(allpos, "\t", strings.Join(posFreq[allpos], "\t")))

// 		}
// 	}
// 	headers.WriteString("\n")

// 	// if buffer.Len() != 0 && headers.Len() != 0 {
// 	// 	fOut, err := os.Create(fileOut)

// 	// 	if err != nil {
// 	// 		log.Fatal("Cannot create file", err)
// 	// 	}
// 	// 	defer fOut.Close()
// 	// 	fmt.Fprintf(fOut, headers.String())
// 	// 	fmt.Fprintf(fOut, buffer.String())
// 	// 	fmt.Printf("\n\nWell done!\n")
// 	// }
// 	matrixPrint(headers, buffer, fileOut)

// }

func matrixDnDs(fileOut string) {
	// var AllPos []int

	// type dndsPerGenome struct {
	// 	dnds, locus string
	// }
	var (
		allLocuses []string

		buffer  strings.Builder
		headers strings.Builder
		// var posCount = make(map[int]int)
		// snps         []snpInfo
		// altPositions        = make(map[string][]allPositionsInGene)
		altPositionsPerFile = make(map[string]map[string][]allPositionsInGene)
		// locDNDS                        = map[string][]string{}
		countNbrOne        = make(map[string]int)
		countNonZeroValues = make(map[string]int)
		countPositive      = make(map[string]int)
		countNegative      = make(map[string]int)
		countNeutral       = make(map[string]int)
		positiveGenes      = make(map[string][]string)
		geneDnDs           = make(map[string]map[string]string)
		positiveGenesList  []string
		positiveGenesCheck = make(map[string]int)

		// locInGenome                    = make(map[string]map[string]string)
		// genomes                        []string
		usedLocusesUnsort, usedLocuses []string
		dndsPerLocus                   = make(map[string][]string)
		// filePerLocus                   = make(map[string][]string)
		groupRegexp            = regexp.MustCompile(`^(\S+)\W+(\w+)\W+(\w+)`)
		group                  = make(map[string]string)
		label                  = make(map[string]string)
		groupHeader, groupBody string
		fileGroup              = make(map[string][]string)
		fileLabel              = make(map[string][]string)
		dndsPerGenomeAll       = make(map[string]map[string]string)
		t0                     = time.Now()
		th                     int
		dndsFloat              float64
	)

	if *statDnDsCountTh != 0 {
		th = *statDnDsCountTh
	} else {
		th = 1
	}
	files := &listOfFiles
	// fmt.Println(files)

	i := 1
	// headers.WriteString("Locus\t")

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
			groupHeader = fmt.Sprintf("%vLabel\tProduct\t", groupHeader)
		}

	} else {

		groupHeader = "Product\t"
	}

	headers.WriteString(fmt.Sprintf("Locus\t%v", groupHeader))

	for key, val := range geneCoordinates {
		if val.Type == "CDS" {
			allLocuses = append(allLocuses, key)
		}

	}
	// fmt.Println(geneCoordinates)
	for fname, snps := range snpCacheMap {
		// t0 = time.Now()
		// locInGenome[fname] = make(map[string]string)
		// headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(key, filepath.Ext(key))))

		// altPositions = nil

		// genomes = append(genomes, fname)

		fmt.Printf("Gathering information: Working on %v from %v (%v)\r", i, len(*files), fname)
		altPositionsPerFile[fname] = make(map[string][]allPositionsInGene)
		for _, val := range snps {

			if val.TypeOf == "CDS" && containsPos(altPositionsPerFile[fname][val.Locus], val.PosInGene, val.Alt) == false {

				// altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, alt: val.Alt, ref: val.NucInPos, locus: val.Locus})
				altPositionsPerFile[fname][val.Locus] = append(altPositionsPerFile[fname][val.Locus], allPositionsInGene{pos: val.PosInGene, alt: val.Alt, ref: val.NucInPos, locus: val.Locus})
				usedLocusesUnsort = append(usedLocusesUnsort, val.Locus)
			}

		}
	}
	// fmt.Println(altPositionsPerFile)
	// dndsPerGenomeAll = make(map[string]string)
	usedLocuses = removeStringDuplicates(usedLocusesUnsort)
	for _, fname := range *files {
		dndsPerGenomeAll[fname] = make(map[string]string)
		geneDnDs[fname] = map[string]string{}
		t1 := time.Now()
		for _, allloc := range allLocuses {

			// prod := getProductByName(allloc)
			if len(altPositionsPerFile[fname][allloc]) > 2 {
				dndsChan := make(chan []string)
				go func() {

					dndsChan <- getDnDsByLocus(allloc, altPositionsPerFile[fname][allloc])
				}()
				dndsRes, ok := <-dndsChan

				// fmt.Println(fname, dndsRes)
				if ok {

					// if allloc == "Rv3879c" {

					// 	fmt.Println(allloc, dndsRes, altPositions[allloc], fname)

					// }

					// fmt.Println(dndsRes)

					// locDNDS[dndsRes[0]] = append(locDNDS[dndsRes[0]], dndsRes[1])
					// if allloc == "Rv3854c" {
					// 	fmt.Println(allloc, dndsRes, altPositions[allloc])
					// }
					// fmt.Println(allloc, dndsRes)
					// locInGenome[fname][dndsRes[0]] = dndsRes[1]
					if dndsRes[1] == "1.00" {
						countNbrOne[allloc]++
						countNeutral[fname]++
					}
					// close(dndsChan)
					dndsPerGenomeAll[fname][dndsRes[0]] = dndsRes[1]
					dndsFloat, _ = strconv.ParseFloat(dndsRes[1], 64)
					geneDnDs[fname][allloc] = dndsRes[1]
					if dndsFloat > 1 {
						countPositive[fname]++
						positiveGenes[fname] = append(positiveGenes[fname], allloc)
						positiveGenesList = append(positiveGenesList, allloc)
						positiveGenesCheck[allloc] = 1
					} else if dndsFloat < 1 {
						countNegative[fname]++
					}
					countNonZeroValues[allloc]++
				}

			} else {

				// locDNDS[allloc] = append(locDNDS[allloc], "1")
				// locInGenome[fname][allloc] = "1.00"
				dndsPerGenomeAll[fname][allloc] = "1.00"
				geneDnDs[fname][allloc] = "1.00"
				// geneDnDs[fname][allloc] = "1.00"
				countNeutral[fname]++
				countNbrOne[allloc]++
			}

		}

		// if *gbVerbose == true {
		fmt.Printf("Calculating DN/DS: Working on %v from %v (%v) \t\t Time:\t%v\n", i, len(*files), fname, t1.Sub(t0))
		// }
		i++
		// fmt.Println(locInGenome)
		// fmt.Println(len(altPositions))

	}
	i = 1

	// sort.Strings(genomes)
	fmt.Println(usedLocuses, len(usedLocuses))
	for _, gname := range *files {
		t1 := time.Now()
		for i := 0; i < len(usedLocuses); i++ {
			// for j := 0; j < len(dndsPerGenomeAll); j++ {
			// if dndsPerGenomeAll[gname] == usedLocuses[i] {
			dndsPerLocus[usedLocuses[i]] = append(dndsPerLocus[usedLocuses[i]], dndsPerGenomeAll[gname][usedLocuses[i]])
			// filePerLocus[allloc] = append(filePerLocus[allloc], gname)

			// }
			// }
		}
		// for _, allloc := range usedLocuses {
		// 	for _, dnds := range dndsPerGenomeAll {
		// 		if dnds.genome == gname && dnds.locus == allloc {
		// 			dndsPerLocus[allloc] = append(dndsPerLocus[allloc], dnds.dnds)
		// 			// filePerLocus[allloc] = append(filePerLocus[allloc], gname)

		// 		}
		// 	}
		// }

		fmt.Printf("Making matrix: Working on %v from %v (%v) \t\t Time:\t%v\n", i, len(*files), gname, t1.Sub(t0))
		i++
	}

	i = 1

	for _, gname := range *files {
		// t1 := time.Now()
		headers.WriteString(fmt.Sprintf("%v\t", gname))
		for _, allloc := range usedLocuses {
			// dndsPerLocus[allloc] = append(dndsPerLocus[allloc], locInGenome[gname][allloc])

			if len(group) != 0 {
				sort.Strings(fileGroup[allloc])
				uniq := removeStringDuplicates(fileGroup[allloc])
				groupBody = fmt.Sprintf("\t%v\t", strings.Join(uniq, ""))
				// fmt.Println(groupBody)
			}
			if len(label) != 0 {
				sort.Strings(fileLabel[allloc])
				uniq := removeStringDuplicates(fileLabel[allloc])
				groupBody = fmt.Sprintf("%v%v\t", groupBody, strings.Join(uniq, " "))
			}

			if len(group) != 0 {
				fileGroup[allloc] = append(fileGroup[allloc], group[gname])
				if len(label) != 0 {
					fileLabel[allloc] = append(fileLabel[allloc], label[gname])
				}
			} else {
				groupBody = "\t"
			}
		}

		// fmt.Printf("Sorting data: Working on %v from %v (%v) \t\t Time:\t%v\n", i, len(*files), gname, t1.Sub(t0))
		// i++
	}

	for _, allloc := range usedLocuses {

		prod := getProductByName(allloc)

		if *statAll == true {
			buffer.WriteString(fmt.Sprintln(allloc, groupBody, prod, "\t", strings.Join(dndsPerLocus[allloc], "\t"), "\t"))

		} else if countNbrOne[allloc] < len(*files) && *statAll == false && countNonZeroValues[allloc] > th {
			// fmt.Println(countNbrOne[allloc], allloc)
			// prod := getProductByName(allloc)
			// fmt.Println(strings.Join(filePerLocus[allloc], "\t"))
			buffer.WriteString(fmt.Sprintln(allloc, groupBody, prod, "\t", strings.Join(dndsPerLocus[allloc], "\t"), "\t"))
		}
	}

	// fmt.Println(dndsPerLocus)

	// for _, gname := range genomes {

	// for _, allloc := range usedLocuses {
	// 	prod := getProductByName(allloc)
	// 	// dnds := locInGenome[gname][allloc]

	// 	if countNbrOne[allloc] != len(*files) && *statAll == false {

	// 		buffer.WriteString(fmt.Sprintln(allloc, groupBody, prod, "\t", strings.Join(dndsPerLocus[allloc], "\t"), "\t"))
	// 	} else if *statAll == true {

	// 		buffer.WriteString(fmtdndsPerLocus[allloc] = append(dndsPerLocus[allloc], locInGenome[gname][allloc]).Sprintln(allloc, groupBody, prod, "\t", strings.Join(dndsPerLocus[allloc], "\t"), "\t"))
	// 		// 			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
	// 	}

	// }

	// }

	headers.WriteString("\n")
	matrixPrint(headers, buffer, fileOut)

	// headers.WriteString("\n")
	// matrixPrint(headers, buffer, fileOut)

	// fmt.Println(locDNDS)

	// for _, allloc := range allLocuses {

	// 	if countNbrOne[allloc] != len(*files) && *statAll == false {
	// 		buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
	// 	} else if *statAll == true {
	// 		buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
	// 	}

	// }
	// // }

	// headers.WriteString("\n")

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
	// matrixPrint(headers, buffer, fileOut)
	if *statDnDsStat {
		var posGenesHeader strings.Builder
		var posGenesBody strings.Builder
		var posGenesArray = make(map[string][]string)
		uniqPosGenes := removeStringDuplicates(positiveGenesList)

		posGenesHeader.WriteString(fmt.Sprintf("genome\t%v\n", strings.Join(uniqPosGenes, "\t")))
		fOut, err := os.Create("dnds_info.txt")

		if err != nil {
			log.Fatal("Cannot create file", err)
		}
		defer fOut.Close()

		fOutPos, err := os.Create("dnds_pos_genes.txt")

		if err != nil {
			log.Fatal("Cannot create file", err)
		}
		defer fOutPos.Close()

		_, _ = fmt.Fprint(fOut, fmt.Sprintf("Genome\tNeutral\tNegative\tPositive\tGenesUnderPosSelection\n"))
		for _, fname := range *files {

			_, _ = fmt.Fprint(fOut, fmt.Sprintf("%v\t%v\t%v\t%v\t%v\n", fname, countNeutral[fname], countNegative[fname], countPositive[fname], strings.Join(positiveGenes[fname], ",")))

			for _, val := range uniqPosGenes {

				if geneDnDs[fname][val] != "" {
					posGenesArray[fname] = append(posGenesArray[fname], geneDnDs[fname][val])
					// 	posGenesBody.WriteString(fmt.Sprintf("%v\t", fname))
					// 	posGenesBody.WriteString(fmt.Sprintf("%v\t", geneDnDs[fname][val]))
					// } else {
					// 	posGenesArray[fname] = append(posGenesArray[fname], "1.0")
				} else if positiveGenesCheck[val] == 1 {
					posGenesArray[fname] = append(posGenesArray[fname], "NA")
				}
				// posGenesBody.WriteString(fmt.Sprint("\n"))
			}

			// for _, val := range uniqPosGenes {

			// 	// if geneDnDs[fname][val] != 0 {
			// 	// 	posGenesBody.WriteString(fmt.Sprintf("%v\t", fname))
			// 	// 	posGenesBody.WriteString(fmt.Sprintf("%v\t", geneDnDs[fname][val]))
			// 	// }
			// 	// posGenesBody.WriteString(fmt.Sprint("\n"))
			// }

		}

		for _, fname := range *files {
			posGenesBody.WriteString(fmt.Sprintf("%v\t%v\n", fname, strings.Join(posGenesArray[fname], "\t")))
		}

		_, _ = fmt.Fprint(fOutPos, fmt.Sprintf("%v", posGenesHeader.String()))
		_, _ = fmt.Fprint(fOutPos, fmt.Sprintf("%v", posGenesBody.String()))

		fmt.Println("DnDs statistics writen to file dnds_info.txt and dnds_pos_genes.txt")

	}
	// fmt.Println(countNegative, countPositive, countNeutral)
}

func matrixBinary(fileOut string, exGenes map[int]int, exSNP map[int]int) {

	/*
		создание бинарной матрицы

	*/

	var (
		AllPosUnsort, AllPos []int
		allLocusUnsort       []string
		buffer               strings.Builder
		headers              strings.Builder
		posCount             = make(map[int]int)
		// snps                 []snpInfo
		posFN   = make(map[int][]string)
		posFreq = map[int][]string{}
		// resPOS  = make(map[int]int)
		th int
		// countPassSNPfromSNPList, countPassSNPinGenes int
	)

	// type binMatrixStruct struct {
	// 	pos, res int
	// 	posArray []string
	// }
	// var ResSeq []seqInfo

	// files := getListofVCF()
	files := &listOfFiles

	if *statBinTH < len(*files) {
		th = *statBinTH
		// fmt.Println("!", th)
	} else if *statBinTH >= len(*files) {
		th = len(*files) - 1
		// fmt.Println("!!", th)

	}

	// fmt.Println(files)
	pos := make(map[int]string)

	headers.WriteString("\t")
	// i := 1
	for key, snps := range snpCacheMap {
		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(key, filepath.Ext(key))))
		//
		// fmt.Printf("Counting SNP positions: Working on %v files from %v \r", i+1, len(*files))

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
	// fmt.Println("Passed SNP from Genes: ", countPassSNPinGenes, " , passed from SNP list: ", countPassSNPfromSNPList)
	AllPos = getAllPosFromCacheMap(exGenes, exSNP)
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

	for _, pos := range AllPos {
		var count0, count1 int
		for i := 0; i < len(posFreq[pos]); i++ {
			if posFreq[pos][i] == "0" {
				count0++
			} else if posFreq[pos][i] == "1" {
				count1++
			}
		}
		// resPOS[pos] = n
		if *statBinTH == 0 {
			// if *gbDebug {
			// 	fmt.Println(count0, count1, th, posFreq[pos], len(posFreq[pos]))
			// }
			if count1 <= len(posFreq[pos])-1 && count0 <= len(posFreq[pos])-1 {
				buffer.WriteString(fmt.Sprintln(pos, "\t", strings.Join(posFreq[pos], "\t")))
			}
		} else {

			if count0 >= th && count1 >= th && th <= len(posFreq[pos]) && count1 < len(posFreq[pos]) && count1 < len(posFreq[pos]) {
				if *gbDebug == true {
					fmt.Printf("pos:%v isRef:%v,isAlt:%v,th:%v %v legnth_array: %v\n", pos, count0, count1, th, posFreq[pos], len(posFreq[pos]))
				}
				buffer.WriteString(fmt.Sprintln(pos, "\t", strings.Join(posFreq[pos], "\t")))
			}

		}

	}

	headers.WriteString("\n")
	matrixPrint(headers, buffer, fileOut)

}

func matrixBinaryMST(fileOut string, verbose bool, exGenes map[int]int, exSNPs map[int]int, randomize bool) {

	var (
		AllPos, SelectedPos []int
		// ResSeq              []seqInfo
		// passSNP = make(map[string]int)
		uniqSNP         = make(map[int]int)
		nbrOfSNP        = 5000
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
						uniqSNP[val.APos] = 2 // 2-EXCLUDED
						continue
					} else if exSNPs[val.APos] == 1 {
						// passSNP["snp"] = passSNP["snp"] + 1
						uniqSNP[val.APos] = 2 // 2-EXCLUDED
						continue

						// fmt.Println(val.APos)

					} else {
						if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
							uniqSNP[val.APos] = 1 // 1-INCLUDED

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

	if nbrOfSNP > len(AllPos) {
		nbrOfSNP = len(AllPos) - 1
	}

	if randomize == true {
		rand.Seed(time.Now().UnixNano())
		for i := 1; i <= nbrOfSNP; i++ {
			rnd := rand.Intn(len(AllPos)-i) + i
			// fmt.Println(AllPos[rnd])
			if len(AllPos) != 0 {
				SelectedPos = append(SelectedPos, AllPos[rnd])
			}

		}
	} else {
		for i := 0; i <= nbrOfSNP; i++ {
			// rnd := rand.Intn(len(AllPos)-i) + i
			// fmt.Println(AllPos[rnd])
			if len(AllPos) != 0 {
				SelectedPos = append(SelectedPos, AllPos[i])
			}

		}
	}
	// fmt.Println(AllPos[0])

	// rand.Intn(max - min) + min

	sort.Ints(SelectedPos)

	// for i := 0; i < len(SelectedPos); i++ {
	// 	headers.WriteString(fmt.Sprintln(SelectedPos[i], "\t"))
	// }
	// headers.WriteString(fmt.Sprintf("%v\n", strings.Join(SelectedPos, "\t")))
	headers.WriteString("FILE\t")
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

func makeSeqBinary(typeof string, verbose bool, ref bool, exGenes map[int]int, exSNPs map[int]int, randomize bool, dp int) {

	var (
		AllPos, SelectedPos []int

		// passSNP = make(map[string]int)
		uniqSNP     = make(map[int]int)
		nbrOfSNP    int
		aaAltCoords = make(map[string]map[int]string)
		aaRefCoords = make(map[int]string)
		dpMAP       = make(map[int][]int)
		locus       = make(map[int]string)
		prod        = make(map[int]string)
		indel       = make(map[int]int)
		filesPOS    = make(map[int][]string)
		nexusTaxa   = make(map[string][]string)
		nChar       int
		nTax        int
	)

	if *annSeqLen != 0 {
		nbrOfSNP = *annSeqLen
	}

	files := &listOfFiles
	for i, file := range *files {
		aaAltCoords[file] = make(map[int]string)

		// создаем запрос в виде типа vcfQuery, который передается через канал на выход <-qSNP.OutChan
		qSNP := &vcfQuery{File: file, OutChan: make(chan vcfInfoQuery), Print: verbose}
		go qSNP.request()
		snpRes := <-qSNP.OutChan
		snpCacheMap[snpRes.File] = snpRes.SnpInfo
		if verbose == true {
			fmt.Printf("Reading files:  %v from %v \r", i+1, len(*files))
		}

	}

	for fname, snps := range snpCacheMap {

		// fmt.Println(tmp)
		for _, val := range snps {
			dpMAP[val.APos] = append(dpMAP[val.APos], val.DP)
			// indel[val.APos] = val.Indel
			locus[val.APos] = val.Locus
			prod[val.APos] = val.Product
			indel[val.APos] = val.Indel
			filesPOS[val.APos] = append(filesPOS[val.APos], fname)

			if val.DP < *gbDP {
				uniqSNP[val.APos] = 2
			}
			if len(exGenes) != 0 {

				for key, value := range exGenes {
					if val.APos >= key && val.APos <= value {

						uniqSNP[val.APos] = 2

						continue
					} else if exSNPs[val.APos] == 1 {

						uniqSNP[val.APos] = 2

						continue

					} else {
						if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 && val.DP >= dp {
							uniqSNP[val.APos] = 1
							aaAltCoords[fname][val.APos] = val.AltAAShort

						}
					}
				}
			} else {
				if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
					uniqSNP[val.APos] = 1
					aaAltCoords[fname][val.APos] = val.AltAAShort
					aaRefCoords[val.APos] = val.RefAAShort

				}

			}

		}
	}

	for key, value := range uniqSNP {

		if value == 1 {

			if len(dpMAP[key]) >= *annMinPosNbr {
				AllPos = append(AllPos, key)

			}

		}

	}

	sort.Ints(AllPos)

	if nbrOfSNP == 0 {
		nbrOfSNP = len(AllPos) - 1
	}

	if nbrOfSNP > len(AllPos) {
		nbrOfSNP = len(AllPos) - 1
	}

	if randomize == true && *annSeqLen != 0 {
		rand.Seed(time.Now().UnixNano())
		for i := 1; i <= nbrOfSNP; i++ {
			rnd := rand.Intn(len(AllPos)-i) + i
			// fmt.Println(AllPos[rnd])
			if len(AllPos) != 0 {
				SelectedPos = append(SelectedPos, AllPos[rnd])
			}

		}
	} else {
		for i := 0; i <= nbrOfSNP; i++ {
			if len(AllPos) != 0 {
				SelectedPos = append(SelectedPos, AllPos[i])
			}

		}
	}

	sort.Ints(SelectedPos)

	if *annPosFile != "" {
		fOut, err := os.Create(*annPosFile)
		fOutF, err := os.Create(fmt.Sprintf("%v_files.txt", *annPosFile))
		if err != nil {
			log.Fatal("Cannot create file", err)
		}
		defer fOut.Close()
		defer fOutF.Close()
		// var posArr []string
		_, _ = fmt.Fprintln(fOut, fmt.Sprintf("---Binary method----"))
		_, _ = fmt.Fprintln(fOut, fmt.Sprintf("[pos]apos:indel:dp:loc:prod"))
		for i, value := range SelectedPos {

			_, _ = fmt.Fprintln(fOut, fmt.Sprintf("[%v]%v:%v:%v:%v:%v\t", i, value, indel[value], dpMAP[value], locus[value], prod[value]))
		}
		_, _ = fmt.Fprintln(fOutF, fmt.Sprintf("---Binary method----"))
		_, _ = fmt.Fprintln(fOutF, fmt.Sprintf("[pos]apos:indel:file:loc:prod"))
		for i, value := range SelectedPos {

			_, _ = fmt.Fprintln(fOutF, fmt.Sprintf("[%v]%v:%v:%v:%v:%v\t", i, value, indel[value], filesPOS[value], locus[value], prod[value]))
		}

	}
	nTax = len(snpCacheMap)
	nChar = len(SelectedPos)
	// fmt.Println(AllPos)
	if *annOutFormat == "phylip" {
		fmt.Printf("%v\t%v\n", nTax, nChar)
	}
	for fname, snps := range snpCacheMap {

		pos := make(map[int]string)

		var buffer strings.Builder

		switch typeof {
		case "Binary":

			if *annOutFormat == "phylip" {
				for _, val := range snps {
					pos[val.APos] = val.Alt

				}
				for _, allpos := range SelectedPos {
					// posCount[allpos] = posCount[allpos] + 1
					if pos[allpos] != "" {

						buffer.WriteString("1")
					} else {

						buffer.WriteString("0")
					}

				}

				fmt.Println(fname, "\t", buffer.String())

			} else if *annOutFormat == "nexus" {

				for _, val := range snps {
					pos[val.APos] = val.Alt

				}
				for _, allpos := range SelectedPos {
					// posCount[allpos] = posCount[allpos] + 1
					if pos[allpos] != "" {

						nexusTaxa[fname] = append(nexusTaxa[fname], "1")
					} else {

						nexusTaxa[fname] = append(nexusTaxa[fname], "0")
					}

				}

			}

		}
	}

	if len(nexusTaxa) != 0 {
		var taxNbr int
		// fmt.Println(nexusTaxa)
		fmt.Printf("#nexus \n\nBEGIN Taxa;\nDIMENSIONS\nntax=%v;\nTAXLABELS\n", nTax)
		for key := range nexusTaxa {
			taxNbr++
			fmt.Printf("[%v]\t'%v'\n", taxNbr, key)
		}
		fmt.Printf(";\nEND; [Taxa]\n\nBEGIN Characters;\nDIMENSIONS nchar=%v;\nFORMAT\n\tdatatype=STANDARD\n\tmissing=?\n\tgap=-\n\tsymbols=\"01\"\n\t labels=left\n\ttranspose=no\n\tinterleave=no\n;\nMATRIX\n", nChar)
		for key, val := range nexusTaxa {

			fmt.Printf("'%v'\t%v\n", key, strings.Join(val, ""))
		}

		fmt.Printf("\n;\nEnd;\n")
	}

}

// #nexus

// BEGIN Taxa;
// DIMENSIONS ntax=2;
// TAXLABELS
// [1] 'taxon_1'
// [2] 'taxon_2'
// ;
// END; [Taxa]

// BEGIN Characters;
// DIMENSIONS nchar=50;
// FORMAT
//         datatype=STANDARD
//         missing=?
//         gap=-
//         symbols="01"
//         labels=left
//         transpose=no
//         interleave=no
// ;
// MATRIX
// 'taxon_1'  10001010000000100010000001000001010100000010001000
// 'taxon_2'  01010001000000000000100000000010000100000000001000
// ;
// End;

// func matrixBinaryOld(fileOut string) {
// 	var (
// 		AllPosUnsort, AllPos []int
// 		allLocusUnsort       []string
// 		buffer               strings.Builder
// 		headers              strings.Builder
// 		posCount             = make(map[int]int)
// 		snps                 []snpInfo
// 		posFN                = make(map[int][]string)
// 		posFreq              = map[int][]string{}
// 	)
// 	// var ResSeq []seqInfo
//
// 	// files := getListofVCF()
// 	files := &listOfFiles
//
// 	// fmt.Println(files)
// 	pos := make(map[int]string)
//
// 	headers.WriteString("Pos\t")
//
// 	for i, file := range *files {
// 		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))
//
// 		fmt.Printf("Counting SNP positions: Working on %v files from %v \r", i+1, len(*files))
//
// 		snpsChan := make(chan []snpInfo)
//
// 		go func() {
// 			snpsChan <- makeSnps(file)
// 		}()
// 		snps = <-snpsChan
//
// 		for _, val := range snps {
//
// 			pos[val.APos] = val.Alt
// 			// posTest = append(posTest, posByFile{pos: pos, file: file, apos: val.APos})
// 			if strings.Contains(strings.Join(posFN[val.APos], " "), file) == false {
// 				posFN[val.APos] = append(posFN[val.APos], file)
// 			}
//
// 			AllPosUnsort = append(AllPosUnsort, val.APos)
// 			posCount[val.APos] = posCount[val.APos] + 1
// 			if val.TypeOf == "CDS" {
// 				allLocusUnsort = append(allLocusUnsort, val.Locus)
// 			}
//
// 		}
//
// 	}
// 	AllPos = unique(AllPosUnsort)
// 	// allLocuses = removeStringDuplicates(allLocusUnsort)
// 	sort.Ints(AllPos)
//
// 	for _, file := range *files {
// 		for _, allpos := range AllPos {
//
// 			if strings.Contains(strings.Join(posFN[allpos], " "), file) {
//
// 				posFreq[allpos] = append(posFreq[allpos], "1")
// 			} else {
// 				// fmt.Println(allpos, 0, file)
// 				posFreq[allpos] = append(posFreq[allpos], "0")
// 			}
//
// 		}
// 	}
//
// 	for _, allpos := range AllPos {
// 		if posCount[allpos] < len(*files) {
//
// 			buffer.WriteString(fmt.Sprintln(allpos, "\t", strings.Join(posFreq[allpos], "\t")))
//
// 		}
// 	}
// 	headers.WriteString("\n")
//
// 	if buffer.Len() != 0 && headers.Len() != 0 {
// 		fOut, err := os.Create(fileOut)
// 		if err != nil {
// 			log.Fatal("Cannot create file", err)
// 		}
// 		defer fOut.Close()
// 		_, _ = fmt.Fprintf(fOut, headers.String())
// 		_, _ = fmt.Fprintf(fOut, buffer.String())
// 		fmt.Printf("\n\nWell done!\n")
// 		// t1 := time.Now()
// 		// fmt.Printf("Elapsed time: %v", fmtDuration(t1.Sub(t0)))
// 	}
//
// }

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
		// altStringArray []string
	)

	// start, end := getGenePosByName(locus)

	refS := getGeneSequence(locus)
	altS := makeAltString(locus, altPositions) // fmt.Println(val, "\n", refS)

	qDnDs := &codon.DnDsQuery{RefSeq: refS, AltSeq: altS, OutChan: make(chan codon.DnDs)}
	go qDnDs.Request()
	// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
	dnds := <-qDnDs.OutChan
	// if locus == "Rv3879c" {

	// 	fmt.Println(locus, dnds, altPositions)

	// }

	close(qDnDs.OutChan)
	// fmt.Println(locus, dnds)

	if dnds.ND != 0 && dnds.NS != 0 {

		// fmt.Printf("L:%v\tdNdS:%v\t%v\n", allloc, fmt.Sprintf("%.2f", dnds.DNDS), prod)
		dndsLoc = fmt.Sprintf("%.2f", dnds.DNDS)

	} else {

		// fmt.Printf("L:%v\tdNdS:%v\t%v\n", allloc, "1", prod)
		dndsLoc = "1.00"
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
	var (
		result      []string
		encountered = map[string]bool{}
	)

	// Create a map of all unique elements.
	for v := range elements {
		encountered[elements[v]] = true
	}

	// Place all keys from the map into a slice.

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
			result = rangePosInfo{Start: resStart + 1, End: resEnd, GeneName: gname, Prod: prod}

		} else if noseq == false {
			seq := getNucFromGenome(resStart-1, resEnd)
			gname, _ := getGeneNameByPos(resStart, resEnd)
			prod, _ := getProductByPos(resStart, resEnd)
			result = rangePosInfo{Start: resStart + 1, End: resEnd, GeneName: gname, Len: len(seq), Seq: seq, Prod: prod}

		}
		// fmt.Println(result)
	} else {
		fmt.Println("Range positions is not valid")
	}

	return result
	// fmt.Println(positions)

}

func printSequenceRange(rangeList []rangePosInfo, web bool, port string, toCircos bool) {
	var basicAnnotation string
	switch web {
	case false:
		if toCircos == true && *statCircosGC == false {
			basicAnnotation = "{{range $element := .}}" +
				"{{.Genome}}\t{{.Start}}\t{{.End}}\t{{.GeneName}}\n" +
				"{{end}}"

		} else if toCircos == true && *statCircosGC == true {
			basicAnnotation = "{{range $element := .}}" +
				"{{if (ge .GC 1.0) }}" +
				"{{.Genome}}\t{{.Start}}\t{{.End}}\t{{.GC}}\n" +
				"{{end}}" +
				"{{end}}"
		} else if *statMakeSeq == true || *statMkSeq == true {

			basicAnnotation = "{{range $element := .}}" +
				">" +
				"{{.Genome}}[{{.Start}}" + ":" + "{{.End}}]{{.GeneName}}\n" +
				"{{.Seq}}\n" +
				"{{end}}"

		} else {
			basicAnnotation = "{{range $element := .}}" +
				"{{if .Seq}}" +
				"{{.GeneName}}\t{{.Start}}\t{{.End}}\t{{.Len}}\t{{.GC}}\t{{.Seq}}\t{{.Prod}}\t{{.Doubles}}\n" +
				"{{else}}" +
				"{{.GeneName}}\t{{.Start}}\t{{.End}}\t{{.Len}}\t{{.GC}}\t{{.Prod}}\t{{.Doubles}}\n" +
				"{{end}}" +
				"{{end}}"

		}

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
			<td>{{.GeneName}}</td><td>{{.Start}}</td><td>{{.End}}</td><td>{{.Len}}</td><td><p title="GC content: {{.GC}}"><textarea rows="3" style="width:400px; word-wrap:break-word;">{{.Seq}}</textarea></p></td><td>{{.Len}}</td><td><p title="{{.Note}}">{{.Prod}}</p></td><td>{{.Doubles}}</td>
			{{else}} 
			<td>{{.GeneName}}</td><td>{{.Start}}</td><td>{{.End}}</td><td>{{.Len}}</td><td><p title="{{.Note}}">{{.Prod}}</p></td><td>{{.Doubles}}</td>
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

		_ = browser.OpenURL(fmt.Sprintf("localhost:%v/", port))

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
		_ = http.ListenAndServe(locPort, nil)
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

func getHash(str string) uint64 {
	// hasher := sha1.New()
	// hasher.Write([]byte(text))
	// return hex.EncodeToString(hasher.Sum(nil))
	hash, err := hashstructure.Hash(str, nil)
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
		if val.GeneName != last {
			fmt.Println(val.GeneName, val.Prod)
		}
		last = val.GeneName

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

func getRangeFromFile(file string, verbose bool, noseq bool, genome string) []rangePosInfo {
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

	if genome == "" {
		genome = gInfo.Strain
	}
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
	var found []rangeArray

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
			if len(seq) > 3 {
				gc, _, _, _ = codon.GcCodonCalc(seq)
			}
		}
		gname, _ := getGeneNameByPos(val.Start, val.End)
		prod, note := getProductByPos(val.Start, val.End)
		fixedProd := strings.Replace(prod, " + ", " ", -1)
		gcRes, _ := strconv.ParseFloat(fmt.Sprintf("%.2f", gc), 64)
		// fixedProd = strings.Replace(prod, "'", " ", -1)

		res := rangePosInfo{Start: val.Start, End: val.End, GeneName: gname, Prod: fixedProd, Len: val.End - val.Start + 1, Seq: seq, Doubles: doubles[val.Start+val.End] + 1, Note: note, GC: gcRes, Genome: genome}

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
				// fmt.Println(rinfo)
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
			if len(aaArray[0]) == 3 && len(aaArray[1]) == 3 {
				idxVal, idxRes := amino.GetComplexIndex(amino.GetShortNameAA(aaArray[0]), amino.GetShortNameAA(aaArray[1]), false)
				fmt.Printf("%v\t%v\t%v\t%v\n", aaArray[0], aaArray[1], idxVal, idxRes)
			} else if len(aaArray[0]) == 1 && len(aaArray[1]) == 1 {
				idxVal, idxRes := amino.GetComplexIndex(aaArray[0], aaArray[1], false)
				fmt.Printf("%v\t%v\t%v\t%v\n", aaArray[0], aaArray[1], idxVal, idxRes)
			}

			// fmt.Println()

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

			if snpCheck.Locus != "" {

				dir := getDirectionByName(snpCheck.Locus)

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
		fmt.Println("Locus not found [getGeneSequence func]", locus)
	}
	return seq
}

func loadExcludeGenes(file string) map[int]int {
	var (
		exGenes    = make(map[int]int)
		exLocus    = regexp.MustCompile(`(^\w+)`)
		start, end int
	)

	f, err := os.Open(file) // открываем файл

	if err != nil {
		fmt.Println(err)

	}
	defer f.Close()

	scanner := bufio.NewScanner(f) //  новый сканер

	for scanner.Scan() {

		scanTxt := scanner.Text()

		for _, exGene := range exLocus.FindAllStringSubmatch(scanTxt, -1) {
			start, end = getGenePosByName(exGene[0])

			// fmt.Println(exGene)
			// fmt.Println(exGene)
			// exGenes[exGene[0]] = 1
			exGenes[start] = end

		}
	}
	return exGenes
}

func loadExcludeRegion(file string) map[int]int {
	var (
		exGenes     = make(map[int]int)
		coordsType1 = regexp.MustCompile(`^\w+\W+(\d+)\W+(\d+)`)
		// coordsType2 = regexp.MustCompile(`^(\d+)\W+(\d+)`)

		start, end int
	)

	f, err := os.Open(file) // открываем файл

	if err != nil {
		fmt.Println(err)

	}
	defer f.Close()

	scanner := bufio.NewScanner(f) //  новый сканер

	for scanner.Scan() {

		scanTxt := scanner.Text()

		for _, exGene := range coordsType1.FindAllStringSubmatch(scanTxt, -1) {
			start, _ = strconv.Atoi(exGene[1])
			end, _ = strconv.Atoi(exGene[2])

			// fmt.Println(exGene)
			// fmt.Println(exGene)
			// exGenes[exGene[0]] = 1
			exGenes[start] = end

		}
	}
	return exGenes
}

func excludeGeneByPos(file string) map[int]int {

	var (
		exGenes = make(map[int]int)
		snpPos  int
	)

	locSNPcheck := readSNPFromFile(file)
	for _, val := range locSNPcheck {
		snpPos, _ = strconv.Atoi(val.APos)

		start, end := getGeneNameCoordsByApos(snpPos)
		if start != 0 && end != 0 {
			// locname, _ := getGeneNameByPos(start, end)
			// prod := getProductByName(locname)
			// fmt.Println(locname, prod)
			exGenes[start] = end
		}

	}
	// fmt.Println(exSNPs)
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
		// start, end := getGeneNameCoordsByApos(snpPos)
		// if start != 0 && end != 0 {
		// 	locname, _ := getGeneNameByPos(start, end)
		// 	prod := getProductByName(locname)
		// 	fmt.Println(locname, prod)
		// }

	}
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

func coord2seq(file string, vcfFile string) altStringResult {
	var (
		coords = regexp.MustCompile(`^(\S+)\W+(\d+)\W+(\d+)`)
		// coordsTwo  = regexp.MustCompile(`^(\w+)\W+(\d+)\W+(\d+)`)
		start, end          int
		locus, prod, altseq string
		altPostitions       []allPositionsInGene
		result              altStringResult
	)

	f, err := os.Open(file) // открываем файл

	if err != nil {
		fmt.Println(err)

	}
	defer f.Close()

	scanner := bufio.NewScanner(f) //  новый сканер

	for scanner.Scan() {

		scanTxt := scanner.Text()

		for _, pos := range coords.FindAllStringSubmatch(scanTxt, -1) {
			start, _ = strconv.Atoi(pos[2])
			end, _ = strconv.Atoi(pos[3])
			// fmt.Println(pos)
			// fmt.Println(exGene)
			// exGenes[exGene[0]] = 1
			// exGenes[start] = end
			locus, _ = getGeneNameByPos(start, end)
			prod, _ = getProductByPos(start, end)
			altPostitions = getAltPositions(start, end, vcfFile)
			altseq = makeAltString(locus, altPostitions)
			result = altStringResult{start: start, end: end, locus: locus, altSeq: altseq, prod: prod, vcfFile: vcfFile}
			// fmt.Println(len(altPostitions), start, end, locus, prod, altseq)

		}

	}
	return result
}

func snp2SeqByLocus(locus string, vcfFile string) altStringResult {
	var (
		// coords = regexp.MustCompile(`^(\S+)\W+(\d+)\W+(\d+)`)
		// coordsTwo  = regexp.MustCompile(`^(\w+)\W+(\d+)\W+(\d+)`)
		start, end    int
		prod, altseq  string
		altPostitions []allPositionsInGene
		result        altStringResult
	)

	// f, err := os.Open(file) // открываем файл

	// if err != nil {
	// 	fmt.Println(err)

	// }
	// defer f.Close()

	// scanner := bufio.NewScanner(f) //  новый сканер

	// for scanner.Scan() {

	// 	scanTxt := scanner.Text()

	// 	for _, pos := range coords.FindAllStringSubmatch(scanTxt, -1) {
	// 		start, _ = strconv.Atoi(pos[2])
	// 		end, _ = strconv.Atoi(pos[3])
	// fmt.Println(pos)
	// fmt.Println(exGene)
	// exGenes[exGene[0]] = 1
	// exGenes[start] = end
	start, end = getGenePosByName(locus)
	prod, _ = getProductByPos(start, end)
	altPostitions = getAltPositions(start, end, vcfFile)
	altseq = makeAltString(locus, altPostitions)
	result = altStringResult{start: start, end: end, locus: locus, altSeq: altseq, prod: prod, vcfFile: vcfFile}
	// fmt.Println(len(altPostitions), start, end, locus, prod, altseq)

	// }

	// }
	return result
}

func getAltPositions(start int, end int, vcfFile string) []allPositionsInGene {
	var (
		// allLocuses []string
		// locus        string
		snps         []snpInfo
		altPositions []allPositionsInGene
		// count        = 1
	)

	// for key, val := range geneCoordinates {
	// 	if val.Type == "CDS" && start >= val.Start && end <= val.End {
	// 		// allLocuses = append(allLocuses, key)
	// 		// altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, alt: val.Alt, ref: val.NucInPos, locus: val.Locus})
	// 		locus = key
	// 	}
	// }
	// if locus=="" {
	// 	fmt.Println("Coor")
	// }
	// fmt.Println(allLocuses)
	// fmt.Println(start, end)
	// sort.Strings(allLocuses)

	snpsChan := make(chan []snpInfo)

	go func() {
		snpsChan <- makeSnps(vcfFile)
	}()
	snps = <-snpsChan
	// fmt.Println(locus)
	for _, val := range snps {
		// fmt.Println(val.Locus, val.PosInGene, val.Alt)
		if start >= val.Start && end <= val.End {
			altPositions = append(altPositions, allPositionsInGene{pos: val.PosInGene, alt: val.Alt, ref: val.NucInPos, locus: val.Locus})
			// locus = val.Locus
		}
	}

	// // fmt.Println(snps)
	// // fmt.Printf("Calculating Dn/DS: Working on %v from %v (%v)\r", i, len(*files), fname)
	// for _, val := range snps {
	// 	if val.TypeOf == "CDS" && containsPos(altPositions[val.Locus], val.PosInGene, val.Alt) == false {

	// 		altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, alt: val.Alt, ref: val.NucInPos, locus: val.Locus})
	// 	}

	// }
	// fmt.Println(altPositions)

	// if len(altPositions) != 0 {

	// // fmt.Println(altPositions)
	// prod := getProductByName(locus)
	// fmt.Printf(">%v (%v)\n%v\n", locus, prod, makeAltString(locus, altPositions))
	// }
	// } else {
	// 	fmt.Printf("No alterations was found in %v region (%v:%v)\n", locus, start, end)
	// }
	// fmt.Println()
	return altPositions
}

func countSNP(vcfFile string) {
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
			prod := getProductByName(allloc)
			fmt.Println(allloc, prod, len(altPositions[allloc]))
		}

	}

}

func pileup2multifasta(file string, th int, gstrain string, gname string, mkseq bool) {
	var (
		nuc                    string
		pos, first, start, end int
		positions              []int
		sequence               []string
		dp                     = 1
		//  foundEnd bool
		fileFormat1 = regexp.MustCompile(`(^\S+)\t(\d+)\W+(\d+)`)
		fileFormat2 = regexp.MustCompile(`(^\S+)\t(\d+)\W+(\w)\W+(\d+)`)

		endFlg                  bool
		resSequence             string
		karyoBuffer, bandBuffer strings.Builder
	)

	if gstrain == "" {
		gstrain = gInfo.Strain
	}

	karyoBuffer.WriteString(fmt.Sprintf("chr -  %v caption %v %v %v\n", gstrain, gInfo.Start, gInfo.End, gstrain))

	if file != "" {

		file, _ := os.Open(file)
		// if err != nil {
		// 	log.Fatal(err)
		// }
		defer file.Close()

		scanner := bufio.NewScanner(file)
		for scanner.Scan() {

			// fmt.Println(scanner.Text())

			format2Res := fileFormat2.FindAllStringSubmatch(scanner.Text(), -1)

			for _, val2 := range format2Res {
				if val2[4] != "" {
					// if gstrain == "" {
					// 	strain = val2[1]
					// } else {
					// 	strain = gstrain
					// }
					pos, _ = strconv.Atoi(val2[2])
					dp, _ = strconv.Atoi(val2[4])
					nuc = val2[3]

					if endFlg == true {
						positions = nil
						endFlg = false
					}

					if first == 0 && dp >= th {
						first = pos

					}

					if pos == first+1 && endFlg == false && dp >= th {
						first = pos
						// fmt.Println("yes", strain, end, first, pos)
						positions = append(positions, pos)
						sequence = append(sequence, nuc)
					} else {
						if len(positions) != 0 {
							// fmt.Println(positions[0], positions[len(positions)-1], strain)
							start = positions[0] - 1
							end = positions[len(positions)-1]
							switch *pileupCircos {
							case false:
								if *pileupShowSeq == true {
									resSequence = strings.Join(sequence, "")
								}
								if start < end && mkseq == false {
									fmt.Printf("%v\t%v\t%v\t%v\t%v\t\n", gstrain, start, end, fmt.Sprintf("%vL%v", gname, end-start), resSequence)

								} else if start < end && mkseq == true {
									// here is + added after start
									fmt.Printf(">%v:%v:%v\n%v\n", fmt.Sprintf("%vL%v", gname, end-start), start, end, strings.Join(sequence, ""))

								} else if start > end {
									fmt.Printf("%v\t%v\t%v\t%v\t\n", gstrain, end, start, fmt.Sprintf("%vL%v", gname, start-end))
								}
								endFlg = true
							case true:

								if start < end {

									karyoBuffer.WriteString(fmt.Sprintf("band %v %v %v %v %v lblue\n", gstrain, fmt.Sprintf("%vL%v", gname, end-start), fmt.Sprintf("%vL%v", gname, end-start), start, end))
									bandBuffer.WriteString(fmt.Sprintf("%v %v %v %v\n", gstrain, start, end, fmt.Sprintf("%vL%v", gname, end-start)))
								} else if start > end {
									karyoBuffer.WriteString(fmt.Sprintf("band %v %v %v %v %v lblue\n", gstrain, fmt.Sprintf("%vL%v", gname, start-end), fmt.Sprintf("%vL%v", gname, start-end), end, start))
									bandBuffer.WriteString(fmt.Sprintf("%v %v %v %v\n", gstrain, start, end, fmt.Sprintf("%vL%v", gname, start-end)))
								}
								endFlg = true
								// fmt.Println(buffer.String())

							}
							first = pos
							sequence = nil
						}

					}
				}

			}

			format1Res := fileFormat1.FindAllStringSubmatch(scanner.Text(), -1)

			for _, val := range format1Res {
				// if gstrain == "" {
				// 	strain = val[1]
				// } else {
				// 	strain = gstrain
				// }
				pos, _ = strconv.Atoi(val[2])
				dp, _ = strconv.Atoi(val[3])

				if endFlg == true {
					positions = nil
					endFlg = false
				}

				if first == 0 && dp >= th {
					first = pos

				}

				if pos == first+1 && endFlg == false && dp >= th {
					first = pos
					// fmt.Println("yes", strain, end, first, pos)
					positions = append(positions, pos)
				} else {
					if len(positions) != 0 {
						// fmt.Println(positions[0], positions[len(positions)-1], strain)
						start = positions[0] - 1
						end = positions[len(positions)-1]
						if start < end {
							fmt.Printf("%v\t%v\t%v\t%v\t\n", gstrain, start, end, fmt.Sprintf("%vL%v", gname, end-start))
						} else if start > end {
							fmt.Printf("%v\t%v\t%v\t%v\t\n", gstrain, end, start, fmt.Sprintf("%vL%v", gname, start-end))
						}
						endFlg = true

					}
					first = pos

				}

				// if count == pos {

				// }

				// fmt.Println(strain, first, end, fmt.Sprintf("%v(%v)L%v",gname, dp, end-first))

				// fmt.Printf("s:%v p:%v s:%v\n", first, pos, pos+1)
				//
				// if (pos == first+1) == false {
				// 	fmt.Printf("s:%v p:%v s:%v\n", first, pos, pos+1)
				// }

			}
			// first = pos

		}

		// if err := scanner.Err(); err != nil {
		// 	log.Fatal(err)
		// }

	}

	if karyoBuffer.Len() != 0 && bandBuffer.Len() != 0 {

		var (
			fKaryo = fmt.Sprintf("%v.karyo", file)
			fBand  = fmt.Sprintf("%v.band", file)
		)
		fileKaryo, err := os.Create(fKaryo)
		if err != nil {
			return
		}
		defer fileKaryo.Close()

		_, _ = fileKaryo.WriteString(karyoBuffer.String())
		fmt.Println(fKaryo, " was successful created!")

		fileBand, err := os.Create(fBand)
		if err != nil {
			return
		}
		defer fileBand.Close()

		_, _ = fileBand.WriteString(bandBuffer.String())
		fmt.Println(fBand, " was successful created!")

		// toCircos(allGenesVal)
	}

}

func locus2Matrix(locus string, listOfFiles []string, typeof string) {
	var (
		// coords = regexp.MustCompile(`^(\S+)\W+(\d+)\W+(\d+)`)
		// coordsTwo  = regexp.MustCompile(`^(\w+)\W+(\d+)\W+(\d+)`)
		start, end        int
		altseq, tabMatrix string
		altPostitions     []allPositionsInGene
		locmatrix         []string
		locPosLetter      = make(map[int][]string)
	)

	// f, err := os.Open(file) // открываем файл

	// if err != nil {
	// 	fmt.Println(err)

	// }
	// defer f.Close()

	// scanner := bufio.NewScanner(f) //  новый сканер

	// for scanner.Scan() {

	// 	scanTxt := scanner.Text()

	// 	for _, pos := range coords.FindAllStringSubmatch(scanTxt, -1) {
	// 		start, _ = strconv.Atoi(pos[2])
	// 		end, _ = strconv.Atoi(pos[3])
	// fmt.Println(pos)
	// fmt.Println(exGene)
	// exGenes[exGene[0]] = 1
	// exGenes[start] = end
	for _, file := range listOfFiles {
		start, end = getGenePosByName(locus)
		// prod, _ = getProductByPos(start, end)
		altPostitions = getAltPositions(start, end, file)
		// fmt.Println(altPostitions)
		altseq = makeAltString(locus, altPostitions)
		locmatrix = strings.Split(altseq, "")
		// fmt.Println(locmatrix)

		// for i := 0; i < len(locmatrix)-1; i++ {
		switch typeof {

		case "nc":
			for i, val := range locmatrix {
				locPosLetter[i] = append(locPosLetter[i], val)

			}
		case "nc_coded":
			for i, val := range locmatrix {
				locPosLetter[i] = append(locPosLetter[i], nuc2IntCode(val))

			}
		case "binary":
			for i, val := range locmatrix {
				if val == strings.ToUpper(val) {
					locPosLetter[i] = append(locPosLetter[i], "1")
				} else {
					locPosLetter[i] = append(locPosLetter[i], "0")
				}

			}

		}

		// tabMatrix = strings.Join(locmatrix, "\t")
	}

	// result = altStringMatrix{start: start, end: end, locus: locus, tabMatrix: tabMatrix, prod: prod, vcfFile: vcfFile}
	// fmt.Println(len(altPostitions), start, end, locus, prod, altseq)

	// }

	// }

	fmt.Println(strings.Join(listOfFiles, "\t"))

	for i := 0; i < len(locmatrix)-1; i++ {

		tabMatrix = strings.Join(locPosLetter[i], "\t")
		fmt.Println(tabMatrix)
	}

	// return tabMatrix
}

func getIUPAC(nucleotides string) string {

	var iupac string
	locNuc := strings.ToTitle(nucleotides)

	if locNuc == "AC" || locNuc == "CA" {
		iupac = "M"
	} else if locNuc == "AG" || locNuc == "GA" {
		iupac = "R"
	} else if locNuc == "CT" || locNuc == "TC" {
		iupac = "Y"
	} else if locNuc == "GC" || locNuc == "CG" {
		iupac = "S"
	} else if locNuc == "AT" || locNuc == "TA" {
		iupac = "W"
	} else if locNuc == "GT" || locNuc == "TG" {
		iupac = "K"
	}
	return iupac
}

func checkPosList(file string) {
	var (
		posInGene   int
		refNuc, seq string
		seqArr      []string
	)
	f, err := os.Open(file) // открываем файл

	if err != nil {
		fmt.Println(err)

	}
	defer f.Close()

	scanner := bufio.NewScanner(f) //  новый сканер

	for scanner.Scan() {

		scanTxt := scanner.Text()
		if pos, err := strconv.Atoi(scanTxt); err == nil {
			locus, _ := getGeneNameByPos(pos, pos)
			prod := getProductByName(locus)
			start, end := getGenePosByName(locus)

			if getDirectionByName(locus) == "f" {
				refNuc = getNucFromGenomePos(pos)
				posInGene = (pos - start) + 1
				seq = getNucFromGenome(start, end)
				seqArr = strings.Split(seq, "")
				seqArr[posInGene-1] = fmt.Sprintf("\n[%v/%v]\n", seqArr[posInGene-1], refNuc)
			} else if getDirectionByName(locus) == "r" {
				refNuc = getComplement(getNucFromGenomePos(pos))
				posInGene = (end - pos) + 1
				seq = getReverseComplement(getNucFromGenome(start, end))
				seqArr = strings.Split(seq, "")
				seqArr[posInGene-1] = fmt.Sprintf("\n%v(%v)[%v|%v]\n", pos, posInGene, seqArr[posInGene-1], refNuc)
			}
			// posInGene := posInGeneFromEnd + posInGeneFromStart
			// seqArr[posInGene+1] = fmt.Sprintf("[%v|%v]", seqArr[posInGene+1], refNuc)
			// seq = strings.Join(seqArr, "")

			seq = strings.Join(seqArr, "")
			fmt.Printf("\n>%v %v %v", locus, prod, seq)
			// fmt.Println(pos, locus, start, end, prod, strings.ToUpper(refNuc), posInGene, seq)

		}
	}
}

func checkPosListVCF(fileUniq string, makeSeq bool) {

	var (
		posCounter                                                                                                        = make(map[int]int)
		foundPos                                                                                                          = make(map[int]snpInfo)
		posInGene                                                                                                         int
		refNuc, altNuc, seq, locus, locusDescr, posfromfile, leftFlankSeq, rightFlankSeq, leftFlankDescr, rightFlankDescr string
		seqArr                                                                                                            []string
		posFromFile                                                                                                       = regexp.MustCompile(`(\d+)`)
	)

	files := &listOfFiles

	f, err := os.Open(fileUniq) // открываем файл

	if err != nil {
		fmt.Println(err)

	}
	defer f.Close()

	scanner := bufio.NewScanner(f) //  новый сканер

	for scanner.Scan() {

		scanTxt := scanner.Text()
		for _, getPos := range posFromFile.FindAllStringSubmatch(scanTxt, -1) {
			posfromfile = getPos[0]
		}

		if pos, err := strconv.Atoi(posfromfile); err == nil {
			for _, file := range *files {

				snpsChan := make(chan []snpInfo)

				go func() {
					snpsChan <- makeSnps(file)
				}()
				snps := <-snpsChan

				// if *annShowFileName == false {
				for _, val := range snps {

					if val.APos == pos {
						posCounter[pos] = posCounter[pos] + 1
						if len(*files) == posCounter[pos] {
							foundPos[pos] = val
						}
						// fmt.Println(val)
					}

				}

			}
		}

	}

	var keys []int
	for key := range foundPos {
		keys = append(keys, key)
	}
	sort.Ints(keys)

	if makeSeq == false {

		for _, key := range keys {

			printTextResults(foundPos[key], false)

		}
	} else {
		for _, key := range keys {

			if foundPos[key].TypeOf == "CDS" || foundPos[key].TypeOf == "repeat_region" {

				direction := foundPos[key].Direction
				locus = foundPos[key].Locus
				start, end := getGenePosByName(locus)
				prod := getProductByName(locus)

				if direction == "f" {
					if *statFlankLeft != 0 {
						leftPos := start - *statFlankLeft
						leftFlankSeq = getNucFromGenome(leftPos, end)
						leftFlankDescr = fmt.Sprintf("+left=%v", *statFlankLeft)
					}
					if *statFlankRight != 0 {
						rightPos := end + *statFlankLeft
						rightFlankSeq = getNucFromGenome(start, rightPos)
						rightFlankDescr = fmt.Sprintf("+right=%v", *statFlankRight)
					}
					refNuc = foundPos[key].NucInPos
					posInGene = (key - start) + 1
					altNuc = foundPos[key].Alt
					// posInGeneT := foundPos[key].PosInGene
					seq = getNucFromGenome(start, end)
					seqArr = strings.Split(seq, "")
					seqArr[posInGene-1] = fmt.Sprintf("\n[%v/%v]\n", refNuc, altNuc)
					locusDescr = fmt.Sprintf("%v [apos:%v,posInG:%v,iupac:%v, %v %v] ", locus, key, posInGene, getIUPAC(fmt.Sprintf("%v%v", refNuc, altNuc)), leftFlankDescr, rightFlankDescr)
				} else if direction == "r" && *statReverse == true {
					if *statFlankLeft != 0 {
						rightPos := end + *statFlankLeft
						rightFlankSeq = getReverseComplement(getNucFromGenome(end, rightPos))
						leftFlankDescr = fmt.Sprintf("+c_left=%v", *statFlankLeft)
					}
					if *statFlankRight != 0 {
						leftPos := start - *statFlankLeft
						leftFlankSeq = getReverseComplement(getNucFromGenome(leftPos, start))
						rightFlankDescr = fmt.Sprintf("+c_right=%v", *statFlankRight)
					}
					refNuc = foundPos[key].NucInPos
					altNuc = foundPos[key].Alt
					posInGene = (end - key) + 1
					seq = getReverseComplement(getNucFromGenome(start, end))
					seqArr = strings.Split(seq, "")

					seqArr[posInGene-1] = fmt.Sprintf("\n[%v/%v]\n", refNuc, altNuc)
					locusDescr = fmt.Sprintf("%v [apos:%v,posInG:%v,iupac:%v, %v %v ] ", locus, key, posInGene, getIUPAC(fmt.Sprintf("%v%v", refNuc, altNuc)), leftFlankDescr, rightFlankDescr)

				} else if direction == "r" && *statReverse == false {
					if *statFlankLeft != 0 {
						leftPos := start - *statFlankLeft
						leftFlankSeq = getNucFromGenome(leftPos, end)
						leftFlankDescr = fmt.Sprintf("+left=%v", *statFlankLeft)
					}
					if *statFlankRight != 0 {
						rightPos := end + *statFlankLeft
						rightFlankSeq = getNucFromGenome(start, rightPos)
						rightFlankDescr = fmt.Sprintf("+right=%v", *statFlankRight)
					}
					refNuc = getReverseComplement(foundPos[key].NucInPos)
					posInGene = (end - key) + 1
					altNuc = getReverseComplement(foundPos[key].Alt)
					// posInGeneT := foundPos[key].PosInGene
					seq = getNucFromGenome(start, end)
					seqArr = strings.Split(seq, "")
					seqArr[posInGene-1] = fmt.Sprintf("\n[%v/%v]\n", refNuc, altNuc)
					locusDescr = fmt.Sprintf("%v [apos:%v,posInG:%v,iupac:%v, %v %v] ", locus, key, posInGene, getIUPAC(fmt.Sprintf("%v%v", refNuc, altNuc)), leftFlankDescr, rightFlankDescr)
				}
				// posInGene := posInGeneFromEnd + posInGeneFromStart
				// seqArr[posInGene+1] = fmt.Sprintf("[%v|%v]", seqArr[posInGene+1], refNuc)
				// seq = strings.Join(seqArr, "")

				seq = strings.Join(seqArr, "")
				if *statFlankLeft != 0 {
					seq = fmt.Sprintf("%v%v", leftFlankSeq, seq)
					leftFlankSeq = ""
				}
				if *statFlankRight != 0 {
					seq = fmt.Sprintf("%v%v", seq, rightFlankSeq)
					rightFlankSeq = ""
				}
				fmt.Printf("\n>%v %v \n%v\n", locusDescr, prod, seq)
			}
		}
	}
}

func reverse(data []string) []string {
	for i := 0; i < len(data)/2; i++ {
		j := len(data) - i - 1
		data[i], data[j] = data[j], data[i]
	}
	return data
}

func getUniqSNP(snps []snpInfo, exGenes map[int]int, exSNPs map[int]int) map[int]int {
	var uniqSNP = make(map[int]int)
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

						// fmt.Println(fname, val.APos, val.AltAAShort, val.RefAAShort)
						// fmt.Println(aaAltCoords)

					}
				}
			}
		} else {
			if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
				uniqSNP[val.APos] = 1

				// fmt.Println(fname, val.APos, val.AltAAShort, val.RefAAShort)
				// fmt.Println(aaAltCoords)

			}
			// AllPosUnsort = append(AllPosUnsort, val.APos)
		}

		// if exGenes[val.Locus] != 1 && exSNPs[string(val.APos)] != 1 {
		// 	AllPosUnsort = append(AllPosUnsort, val.APos)
		// }

	}
	return uniqSNP
}

func annotateFromList(file string, genes []geneInfo) {

	var (
		posAlt = regexp.MustCompile(`^(\d+)\W+(\w)\W+(\w)`)
		snp    snpInfo
	)
	f, err := os.Open(file) // открываем файл

	if err != nil {
		fmt.Println(err)

	}
	defer f.Close()

	scanner := bufio.NewScanner(f) //  новый сканер

	for scanner.Scan() {
		for _, findPosAlt := range posAlt.FindAllStringSubmatch(scanner.Text(), -1) {

			if len(findPosAlt) == 4 {

				for z := 0; z < len(genes); z++ {
					// g := genes[z]
					// lStart := genes[z].Start
					// lEnd := genes[z].End
					apos, _ := strconv.Atoi(findPosAlt[1])
					alt := findPosAlt[3]
					if apos >= genes[z].Start && apos <= genes[z].End {
						qSnpInfo := &snpInfoQuery{OutChan: make(chan snpInfo), apos: apos, g: genes[z], alt: alt, index: true}
						go qSnpInfo.request()
						snp = <-qSnpInfo.OutChan
						printTextResults(snp, false)
						fmt.Println(snp.TangIdxVal)
					}
				}

			}

		}
	}

}

func makeAlign(verbose bool, ref bool) {

	type posPerFile struct {
		Fname string
		Pos   []map[int]string
	}

	var (
		// AllPos, SelectedPos, TempPos []int
		// ResSeq                       []seqInfo
		// passSNP = make(map[string]int)
		uniqueSNP = make(map[int]int)
		// nbrOfSNP       int
		aaAltCoords = make(map[string]map[int]string)
		// aaRefCoords    = make(map[int]string)
		dpMAP    = make(map[int][]int)
		locus    = make(map[int]string)
		prod     = make(map[int]string)
		filesPOS = make(map[int][]string)
		indel    = make(map[int]int)

		// seq strings.Builder
		lineSeq []string

	// posPerFile = make(map[string]map[int]string)
	// wg           sync.WaitGroup
	// lenCount =1
	)

	// files := getListofVCF()

	// queryChan := make(chan vcfInfoQuery)

	// if *annSeqLen != 0 {
	// 	nbrOfSNP = *annSeqLen
	// }

	files := &listOfFiles

	for i, file := range *files {
		aaAltCoords[file] = make(map[int]string)

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
	// var wg sync.WaitGroup
	// wg.Add(len(*files))
	// go func() {
	// 	defer wg.Done()
	for fname, snps := range snpCacheMap {
		fmt.Printf("\n>%v\n", fname)
		snpPos := make(map[int]string)
		for _, val := range snps {
			// DP - глубина прочтения СНИПа
			dpMAP[val.APos] = append(dpMAP[val.APos], val.DP)
			// indel[val.APos] = val.Indel
			locus[val.APos] = val.Locus
			prod[val.APos] = val.Product
			filesPOS[val.APos] = append(filesPOS[val.APos], fname)
			indel[val.APos] = val.Indel
			if val.DP < *gbDP {
				uniqueSNP[val.APos] = 2
			}
			if val.Indel == 0 {
				// fmt.Println(val.APos,val.Alt)
				snpPos[val.APos] = val.Alt

			} else {
				// fmt.Println(val.APos,val.IndelAlt)
				snpPos[val.APos] = val.IndelAlt
			}

		}

		for i := gInfo.Start; i < gInfo.End; i++ {

			// lenCount++
			// fmt.Println(i,getNucFromGenomePos(i),snpPos[i])
			if len(snpPos[i]) != 0 {
				// lineSeq.WriteString(snpPos[i])
				lineSeq = append(lineSeq, snpPos[i])
				// fmt.Println(snpPos[i])
			} else {
				lineSeq = append(lineSeq, getNucFromGenomePos(i))
				// fmt.Println(getNucFromGenomePos(i))
			}
			if len(lineSeq) == 70 {
				// seq.WriteString(lineSeq.String())
				// fmt.Println(lineSeq.String())
				seq := strings.Join(lineSeq, "")
				fmt.Printf("%v\n", seq)
				lineSeq = nil
				// seq.WriteString("\n")

			}
		}
	}

	// }()
	// wg.Wait()

	// fmt.Println(posPerFile)
}
