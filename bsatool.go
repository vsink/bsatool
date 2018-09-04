package main

import (
	"bufio"
	// "context"
	// "testing"
	// "compress/gzip"
	"fmt"
	// "encoding/json"
	"log"
	// "github.com/CrowdSurge/banner"
	// "io"
	// "encoding/binary"
	"encoding/gob"
	// "encoding/hex"
	"net/http"
	"os"
	// "os/exec"
	// "reflect"
	// "math/rand"
	"regexp"
	// "encoding/json"
	// "flag"
	"strings"
	// "time"
	// "gopkg.in/alecthomas/kingpin.v2"
	// "encoding/gob"
	// "bufio"
	"compress/lzw"
	// "io"
	"sort"
	"strconv"

	"./amino"
	"./gene"
	// "./tang"
	kingpin "gopkg.in/alecthomas/kingpin.v2"
	// "text/template"
	// "net/http"
	// "github.com/davecgh/go-spew/spew"
	// tm "github.com/buger/goterm"
	// "gopkg.in/cheggaaa/pb.v2"
	"html/template"
	// "text/template"

	"./codon"

	// "github.com/fatih/color"

	textdistance "github.com/masatana/go-textdistance"
	"github.com/mitchellh/hashstructure"
	"github.com/pkg/browser"
	// "io/ioutil"
	// "encoding/csv"
	// "crypto/md5"

	"path/filepath"
	// "text/template"
	// "sync"
	// "arg"
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
                                                                          
 BSATool - Bacterial Snp Annotation Tool ver.0.1.220518alpha
       Laboratory of Social and Epidemic Infections
 Scientific Centre for Family Health and Human Reproduction Problems
     (c) V.Sinkov & O.Ogarkov,Irkutsk, Russia, 2018                                   
                                                  
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
	// url         = "localhost:8080"

)

// var url = fmt.Sprintf("localhost:%v", *flgPort)

// var args struct {
// 	MkDB, Debug, GenomeMap, Tang, Web, Share, InDel, Verbose, WithRef, Dev bool
// 	GBFile, Out, DB, VCF, MakeSeq, SNP, TestMode, InFile, Type, Action     string
// 	Pos, Port                                                              int
// }

// type SNP interface {
//
// }

// var linesFromGB []string             // массив строк считанных из файла genbank
// var genomeAnchors = map[int]string{} // карта резльтатов парсинга генбанк-файла
// var totalCDSFound int // количество найденных CDS
// var allGenesVal []gene.Gene          // массив информации об генах, загруженный из файла базы данных
// var (
// 	flgMkDB   = flag.Bool("mkdb", false, "Создание файла локальной базы данных") // флаг создания базы данных
// 	flgGBfile = flag.String("gb", "", "Genbank файл")                            // генбанк файл
// 	flgabout  = flag.Bool("v", false, "О программе")                             // о программе

// 	flgOut       = flag.String("out", "", "Имя локальной базы данных") // выходной файл базы данных при парсинге файла генбанка
// 	flgDB        = flag.String("db", "", "Файл локальной базы данных") // файл основной базы данных
// 	flgVCF       = flag.String("vcf", "", "vcf файл")                  // файл vcf, с параметром=list все файлы в директории
// 	flgPos       = flag.Int("pos", 0, "SNP-позиция")
// 	flgDebug     = flag.Bool("debug", false, "Режим отладки") // режим отладки
// 	flgGenomeMap = flag.Bool("gmap", false, "Карта генома")
// 	flgmakeSeq   = flag.String("mkseq", "", "NC или AA")   // сгенерировать синтетическую последовательность
// 	flgTang      = flag.Bool("tang", false, "Индек Танга") // индекс танга
// 	// var flgDevMode = flag.Bool("dev", false, "")

// 	// var flgReportType = flag.String("rep", "", "")
// 	flgSNP         = flag.String("snp", "", "Файл с известными SNP, которые необходимо проверить")
// 	flgWeb         = flag.Bool("web", false, "Режим вывода результатов в Веб браузер")             // вывод результатов в браузер
// 	flgShare       = flag.Bool("share", false, "Поиск общих, для всех vcf файлов в папке, СНИПОВ") // поиск общих снипов для всех Vcf файлов в папке
// 	flgIndel       = flag.Bool("indel", false, "Поиск ИнДелов")
// 	flgVerbose     = flag.Bool("verb", false, "Вывод дополнительной информации")
// 	flgWithRef     = flag.Bool("ref", false, "Создавать референсную последовательность при использовании команды -mkseq")
// 	flgTestMode    = flag.String("test", "", "test mode")
// 	flgInFile      = flag.String("in", "", "input")
// 	flgType        = flag.String("type", "", "type of matrix")
// 	flgDev         = flag.Bool("dev", false, "dev mode")
// 	flgPort        = flag.Int("port", 8080, "set localhost:port")
// 	flgAction      = flag.String("action", "", "select action of work")
// 	flgTestOptions = flag.String("opt", "", "Test Options")
// 	flgTH          = flag.Int("th", 0, "set treshold to range len sequence")
// 	flgNoSeq       = flag.Bool("noseq", false, "don't show sequence in range output")
// )

// var flgParam =flag.String("param", "", "Test Options")
var (

	//
	//Database flags
	app       = kingpin.New(logo, "BSATool - Bacterial Snp Annotation Tool")
	appAuthor = app.Author("V.Sinkov && O. Ogarkov")
	appVer    = app.Version("0.1.0")

	gbWeb     = kingpin.Flag("web", " Open results in web browser").Short('w').Bool()
	gbVerbose = kingpin.Flag("verbose", "Show additional information ").Short('v').Default("false").Bool()
	gbIndex   = kingpin.Flag("index", "Calculate Complex Index for amino acid changes").Default("false").Bool()
	gbPort    = kingpin.Flag("port", "Use your own localhost:port (default:8080)").Default("8080").String()
	gbNoSeq   = kingpin.Flag("noseq", "Don't show nucleotides").Default("false").Bool()
	gbDebug   = kingpin.Flag("debug", "Debug mode").Default("false").Bool()

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

	// compute statistic options

	statAction = kingpin.Command("stat", "Calculates statistic tests").Alias("s")
	statDB     = statAction.Flag("db", "Database file").Short('b').Required().String()
	statTask   = statAction.Flag("action", "Type of action:share, snp,dnds, ginfo,bed, matrix,range").Short('a').Required().String()
	// statWeb     = annAction.Flag("web", "").Short('w').Bool()
	statInFile  = statAction.Flag("in", "Input file").Short('i').String()
	statOutFile = statAction.Flag("out", "Output file").Short('o').String()
	statTypeOf  = statAction.Flag("type", "Type of matrix (binary, gc3, dnds, table, nc, locus. freq, jw").Short('t').String()
	statInRule  = statAction.Flag("rule", "Input rule file").Short('r').String()
	// statMkSeq   = statAction.Flag("mkseq", "").Bool()
	// statTH      = statAction.Flag("th", "").Int()

	devAction = kingpin.Command("dev", "Developer mode.").Alias("debug")
	devDB     = devAction.Flag("db", "Database file").Short('b').Required().String()
	devTask   = devAction.Flag("action", "Action...").Short('a').Required().String()
	devPwd    = devAction.Flag("pwd", "Password").Short('p').Required().String()

	betaAction  = kingpin.Command("beta", "Beta test mode for testing new functions")
	betaTask    = betaAction.Flag("action", "Action...").Short('a').Required().String()
	betaDB      = betaAction.Flag("db", "Database").Short('b').Required().String()
	betaInFile  = betaAction.Flag("in", "Input file").Short('i').String()
	betaOutFile = betaAction.Flag("out", "Output file").Short('o').String()
)
var gInfo genomeInfo

// var featAllGenesVal []featInfo

// var allPos []allPositionsInGene

type seqInfo struct {
	Name, Seq string //
	// Len       int
}
type genomeInfo struct {
	Organism, Strain, Version string
	Start, End, NucCore       int
}

type statInfo struct {
	Pos, Count, Perc        int
	FilesWith, FilesWithout string
}
type checkSNP struct {
	FileName, FoundSNP, ExactWithRule string
}

type allPositionsInGene struct {
	pos int
	nuc string
}

// DnDsRes is....
type DnDsRes struct {
	N, S, PN, PS, DN, DS, DNDS float64
	ND, NS                     int
	Locus, Product             string
}

// JaroWinklerInfo is....
type JaroWinklerInfo struct {
	Locus, Product string
	JWDist         float64
}

// GC3Type is....
type GC3Type struct {
	GC3Alt, GC3Ref, GCalt, GCref float64
	Locus                        string
}

type rangePosInfo struct {
	Start, End, Len, Doubles int
	Seq, Prod, Gname, Note   string
	GC                       float64
}

// g.Start, g.End, g.Locus, g.Product, g.Direction, g.TypeOf

type genomeMapInfo struct {
	Start, End             int
	Locus, Product, TypeOf string
}

type rangeArray struct {
	Start, End int
}

type rulesInfo struct {
	Name     string
	Variants []string
	Lenght   int
}

// type debugInfo struct {

// }

// type genomeCoordInfo struct {
// 	Start, End  int
// 	Locus, Prod string
// }

// type featInfo struct {
// 	Start, End, Direction, Name, Product, Note, TypeOf, Locus string
// }

// type genomeCoords struct {
// 	Start, End                              int
// 	Locus, Product, Direction, TypeOf, Note string
// }

// type allPositionsInGeneArray struct {
// 	items []allPositionsInGene
// }

// type locusInfo struct {
// 	start, end int
// 	alt        []allPositionsInGene
// }

var genomeSeqSlice []string // информация об генах, загруженная из базы
// var genomeCoordinates []genomeCoordInfo
var allGenesVal []gene.Gene
var rulesArr []rulesInfo

func main() {

	// парсинг флагов
	defer os.Exit(0)
	// flag.Parse()
	kingpin.New("BSATool", "BSATool - Bacterial Snp Annotation Tool")
	kingpin.Version("0.1.0").Author("V.Sinkov")

	switch kingpin.Parse() {

	// bsatool mkdb/create/makedb -i /--genbank FILE --out /-o FILE
	case "mkdb":

		if _, err := os.Stat(*dbGenbank); os.IsNotExist(err) {
			fmt.Printf("The %v file is not exist!\n", *dbGenbank)
			os.Exit(3)
		}
		allGenesVal, genomeSeqSlice = readGBFile(*dbGenbank, *gbVerbose)
		writeDB(*dbName, allGenesVal)

	case "annotate":

		allGenesVal = readDB(*annDB)
		if *annVCF == "list" || *annVCF == "*" || *annVCF == "all" {
			if *gbWeb == false && *annMakeSeq == "" {
				listOfVCFFiles(*annWithFilenames)
			} else if *gbWeb == true && *annMakeSeq == "NC" {
				createNCWebServer(*gbPort)
			} else if *gbWeb == false && *annMakeSeq == "NC" {
				seq := makeSeq(*annMakeSeq, *gbVerbose, *annMakeSeqRef)
				for _, val := range seq {
					fmt.Println(val.Seq)
				}

			}
		} else {
			if _, err := os.Stat(*annVCF); os.IsNotExist(err) {
				fmt.Printf("The %v file is not exist!\n", *annVCF)
				os.Exit(3)
			}
			if *gbWeb == false && *annMakeSeq == "" {
				parserVCF(*annVCF, true, allGenesVal)
			} else if *gbWeb == true && *annMakeSeq == "" {
				snps := parserVCF(*annVCF, false, allGenesVal)
				printWebResults(snps, *gbPort)
			}
		}
	case "stat":
		allGenesVal = readDB(*statDB)
		switch *statTask {
		case "share":

			getShareSNP(*gbVerbose, *gbWeb, getListofVCF())

		case "snp":

			snpStat()
		case "dnds":
			if _, err := os.Stat(*statInFile); os.IsNotExist(err) {
				fmt.Println("No input file found")
				os.Exit(3)
			}
			calcDnDsVal(*statInFile, *gbVerbose)
		case "ginfo":
			testGeneInfo(allGenesVal)
		case "circos":
			toCircos(allGenesVal)
		case "bed":
			toBED()
		case "matrix":
			if *statTypeOf != "" && *statOutFile != "" {
				makeMatrix(*statTypeOf, *statOutFile)
			} else if *statTypeOf != "" && *statOutFile == "" {
				fmt.Println("Please, use key -o (--out) to save data.")
			}
		case "range":
			// go run bsatool.go stat -a range --db test_core -i BWAIntersections.txt

			// for _, val := range *statOptions {
			// var opt []string
			// opt =strings.Split(*statOptions, ",")

			if *statInFile != "" {
				file := *statInFile
				res := getRangeFromFile(file, *gbVerbose, *gbNoSeq)
				printSequenceRange(res, *gbWeb, *gbPort)
			}
		case "check":
			// go run bsatool.go stat -a check --db test_core -i drugs2.txt -r rule.txt -w

			var useRule bool
			if *statInRule != "" {
				useRule = true
				rulesArr = checkRuleFromFile(*statInRule)
			}
			if *statInFile != "" {

				checkSNPfromFile(*statInFile, *gbVerbose, *gbWeb, useRule)
			}
		case "rule":
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
			}
		}
		fmt.Println("DevMode")
	case "beta":
		allGenesVal = readDB(*betaDB)
		fmt.Println("Beta mode")
		switch *betaTask {
		// case "rule":
		// 	if *betaInFile != "" {

		// 	}
		}

	}

	// о программе
	// if *flgabout == true {
	// 	about()
	// }

	//  файл базы данных
	// if *flgMkDB == true {
	// 	if *flgGBfile != "" && *flgOut != "" {
	// 		if _, err := os.Stat(*flgGBfile); os.IsNotExist(err) {
	// 			fmt.Printf("The %v file is not exist!\n", *flgGBfile)
	// 			os.Exit(3)
	// 		}
	// 		//
	// 		// if *flgTestMode == true {

	// 		// 	allGenesVal, genomeSeqSlice = readGBFile(*flgGBfile)
	// 		// 	// readGBFile(*flgGBfile)
	// 		// 	writeDB(*flgOut, allGenesVal)
	// 		// } else {
	// 		// allGenesVal, genomeSeqSlice = readGBFile(*flgGBfile)
	// 		allGenesVal, genomeSeqSlice = readGBFile(*flgGBfile)
	// 		writeDB(*flgOut, allGenesVal)
	// 		// }
	// 		// fmt.Printf("%v\n", genomeSeqSlice)
	// 	} else {
	// 		fmt.Println("Не указаны необходимые аргументы (-mkdb -gb=GENBANK_FILE -out=OUTPUT_DB_FILE)")
	// 	}
	// }

	// if *flgDB != "" {
	// 	// g := Gene
	// 	if _, err := os.Stat(*flgDB); os.IsNotExist(err) {
	// 		fmt.Printf("The %v file is not exist!\n", *flgDB)
	// 		os.Exit(3)
	// 	}

	// 	allGenesVal = readDB(*flgDB)

	// 	if *flgTestMode == "stat" {
	// 		snpStat()
	// 	} else if *flgTestMode == "dnds" {
	// 		calcDnDsVal("13_MoU_m.vcf", true)
	// 	} else if *flgTestMode == "ginfo" {
	// 		testGeneInfo(allGenesVal)
	// 		// fmt.Println(gInfo)
	// 	} else if *flgTestMode == "circos" {
	// 		toCircos(allGenesVal)
	// 	} else if *flgTestMode == "bed" || *flgTestMode == "2bed" || *flgTestMode == "ToBed" {
	// 		toBED()
	// 		// } else if *flgTestMode == "range" && *flgTestOptions != "" && *flgInFile == "" {
	// 		// 	var res []rangePosInfo

	// 		// 	res = append(res, getSequenceRange(*flgTestOptions))

	// 		// 	printSequenceRange(res)

	// 	} else if *flgTestMode == "range" && *flgInFile != "" {

	// 		// file, err := os.Open(*flgInFile)
	// 		// if err != nil {
	// 		// 	log.Fatal(err)
	// 		// }
	// 		// defer file.Close()

	// 		// scanner := bufio.NewScanner(file)
	// 		// var res []rangePosInfo
	// 		// for scanner.Scan() {

	// 		// 	res = append(res, getSequenceRange(scanner.Text()))
	// 		// }
	// 		// printSequenceRange(res)
	// 		file := *flgInFile
	// 		res := getRangeFromFile(file, false)
	// 		printSequenceRangeTest(res, *flgWeb)

	// 		// getSequenceRange(*flgTestOptions)
	// 	} else if *flgTestMode == "annotate" && *flgInFile != "" {
	// 		annotatePositions()
	// 	} else if *flgTestMode == "matrix" && *flgOut != "" && *flgType != "" {
	// 		makeMatrix(*flgType, *flgOut)

	// 	} else if *flgTestMode == "dump" {
	// 		fmt.Println(allGenesVal)

	// 	}

	// 	switch {
	// 	case *flgVCF != list && *flgVCF != "" && *flgWeb == false && *flgmakeSeq == "":
	// 		if _, err := os.Stat(*flgVCF); os.IsNotExist(err) {
	// 			fmt.Printf("The %v file is not exist!\n", *flgVCF)
	// 			os.Exit(3)
	// 		}
	// 		parserVCF(*flgVCF, true, allGenesVal)
	// 	case *flgVCF != list && *flgVCF != "" && *flgWeb == true && *flgmakeSeq == "":
	// 		if _, err := os.Stat(*flgVCF); os.IsNotExist(err) {
	// 			fmt.Printf("The %v file is not exist!\n", *flgVCF)
	// 			os.Exit(3)
	// 		}
	// 		snps := parserVCF(*flgVCF, false, allGenesVal)
	// 		printWebResults(snps)
	// 	case *flgVCF == list && *flgWeb == false && *flgmakeSeq == "" && *flgShare == false && *flgSNP == "":
	// 		listOfVCFFiles(false)
	// 	case *flgWeb == false && *flgmakeSeq == ncFlag && *flgShare == false:
	// 		seq := makeSeq(ncFlag)
	// 		for _, val := range seq {
	// 			fmt.Println(val.Seq)
	// 		}
	// 	case *flgWeb == true && *flgmakeSeq == ncFlag && *flgShare == false:
	// 		createNCWebServer()

	// 	case *flgWeb == false && *flgmakeSeq == aaFlag && *flgShare == false:
	// 		makeSeq(aaFlag)
	// 	case *flgShare == true:
	// 		// getShareSNP()
	// 		getShareSNP(*flgDebug, *flgWeb)
	// 	case *flgSNP != "":
	// 		if _, err := os.Stat(*flgSNP); os.IsNotExist(err) {
	// 			fmt.Printf("The %v file is not exist!\n", *flgSNP)
	// 			os.Exit(3)
	// 		}

	// 		checkSNPfromFile(*flgSNP)

	// 		// case *flgDebug == true:
	// 		// 	fmt.Printf("%v\n", allGenesVal)
	// 	}

	// 	if *flgGenomeMap == true && *flgTestOptions != "coord" {
	// 		gm := getGenomeMap()
	// 		for _, val := range gm {
	// 			if *flgNoSeq == true {
	// 				fmt.Println(val.Locus, val.Start, val.End, val.Product, val.TypeOf)
	// 			} else {
	// 				seq := getNucFromGenome(val.Start-1, val.End)
	// 				fmt.Printf(">%v %v:%v:%v\n%v\n", val.Locus, val.Start, val.End, val.Product, seq)
	// 			}
	// 		}
	// 		// 	var igensS []int
	// 		// 	var igensE []int
	// 		// 	var leftGene []string

	// 		// 	for i, g := range allGenesVal {

	// 		// 		igensS = append(igensS, g.Start)
	// 		// 		igensE = append(igensE, g.End)
	// 		// 		leftGene = append(leftGene, g.Locus)

	// 		// 		if i >= 1 && i < len(allGenesVal)-1 {
	// 		// 			// 	// igenEnd, _ := strconv.Atoi(g.Start)
	// 		// 			checkEndOfCDS := (igensE[i-1] + 1) - g.Start
	// 		// 			if checkEndOfCDS < 0 {
	// 		// 				fmt.Printf("IGR: s%v e:%v  %v<--->%v\n", igensE[i-1]+1, igensS[i]-1, leftGene[i-1], g.Locus)
	// 		// 			}
	// 		// 		}
	// 		// 		fmt.Printf("s:%v e:%v n:%v p:%v d:%v t:%v \n", g.Start, g.End, g.Locus, g.Product, g.Direction, g.TypeOf)

	// 		// }

	// 	} else if *flgGenomeMap == true && *flgTestOptions == "coord" {

	// 		getCoordRange(1000, 20000)
	// 		// fmt.Println("!!!")
	// 	}

	// }

}

func getSNPInfo(apos int, g gene.Gene, alt string, flgTang bool) (gene.SNPinfo, int) {
	var snp gene.SNPinfo           // структура SNP
	var codonPositions []string    // срез для разбивки кодона побуквенно
	var altCodonPositions []string // срез для разбивки кодона побуквенно альтернативным нуклеотидом
	var locReportType, typeOf, titv string
	var geneLen int

	var trouble int
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

	var codon string    // переменная для кодона
	var altCodon string // аналогично для альтернативного кодона
	var mut string
	var tangIdx string
	var tangIdxVal int
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
		alt = codonReverse(alt)
		nucG = codonReverse(nucG)
		codon = codonReverse(codon)
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
		tangIdx = "----"
	} else if aaRefShort != aaAltShort && aaAltShort != "X" {
		mut = "missense"

		// tangIdx = strconv.FormatFloat(amino.GetTangInx(aaRefShort, aaAltShort), 'f', 2, 64)
		tangIdx, tangIdxVal = amino.GetComplexIndex(aaRefShort, aaAltShort, *gbVerbose)
	} else if aaRefShort != aaAltShort && aaAltShort == "X" {
		mut = "nonsense"
		tangIdx = "----"
	}

	if strings.ToUpper(alt) == strings.ToUpper(nucG) {
		trouble = 1
	}

	if flgTang == true {
		locReportType = "T1"
	} else if flgTang == false {
		locReportType = "T0"
	}
	// amino.GetTangInx(aaRefShort, aaAltShort)
	titv = checkTiTv(nucG, alt)
	// fmt.Println(lStart, lEnd, g.Direction, geneLen, g.Locus, g.Product)

	snp = gene.SNPinfo{APos: apos, PosInGene: posInGene, PosInCodonG: posInCodonG,
		RefCodon: codon, RefAA: aaRef, NucInPos: strings.ToUpper(nucG), Locus: g.Locus,
		Direction: g.Direction, Name: g.Locus, Product: g.Product,
		Start: g.Start, End: g.End, CodonNbrInG: codonNbrInG, AltCodon: altCodon,
		AltAA: aaAlt, RefAAShort: aaRefShort, AltAAShort: aaAltShort,
		Mutation: mut, Tang: tangIdx, TangIdxVal: tangIdxVal, Alt: alt,
		Note: g.Note, ReportType: locReportType, ProteinID: g.ProteinID,
		GeneID: g.GeneID, GOA: g.GOA, GeneLen: geneLen, TiTv: titv, TypeOf: typeOf}

	return snp, trouble
}

func getNucFromGenome(start int, end int) string {
	/*
		возвращает последовательность нуклеотидов из генома:

	*/
	var result string
	var slice []string
	slice = genomeSeqSlice[start:end]
	result = strings.Join(slice, "")
	return result

}

func getNucFromGenomePos(pos int) string {
	/*
		возвращает последовательность нуклеотидов из генома:

	*/
	var result string
	var slice []string
	slice = genomeSeqSlice[pos-1 : pos]
	result = strings.Join(slice, "")
	return result

}

// func readGBFile(file string) (g []gene.Gene, genomeSplice []string) {
// 	//функция для считывания файла в формате генбанка и занесения строк в массив linesFromGB

// 	var qenomeSeq = regexp.MustCompile(`\s+(\d+)\s+(\w{6}.*)`)
// 	var regDelSpaces = regexp.MustCompile(`(\s+)(.*)`)
// 	var checkOrigin = regexp.MustCompile(`(\d+\s\w{10}\s\w{10}\s\w{10})`)
// 	// var checkCDS = regexp.MustCompile(`(CDS\s+complement\W\d+\W+\d+\W)|(CDS\s+\d+\W+\d+)`)
// 	// var endOfCDS = regexp.MustCompile(`(gene\s+\d+\W+\d+)|(gene\s+complement\W\d+\W+\d+)`)
// 	var startOrigin = regexp.MustCompile(`(ORIGIN\s+)|(translation=)`)
// 	var feature = regexp.MustCompile(`^\s{5}(\w+)\s+complement\W(\d)+\W+(\d)+\W|^\s{5}(\w+)\s+(\d)+\W+(\d+)`)
// 	var makeAnchors = regexp.MustCompile(`^(/)`)
// 	var genomeSource = regexp.MustCompile(`organism=\W(.*)\W`)
// 	var genomeStrain = regexp.MustCompile(`strain=\W(.*)\W`)
// 	var genomeStartEnd = regexp.MustCompile(`source\W+(\d+)..(\d+)`)

// 	var gStart, gEnd int

// 	//  source          1..4411532
// 	//                  /organism="Mycobacterium tuberculosis H37Rv"
// 	//                  /mol_type="genomic DNA"
// 	//                  /strain="H37Rv"
// 	//                  /db_xref="taxon:83332"
// 	// var replaceCDS = regexp.MustCompile(`(CDS\s+.*)`)

// 	var resString []string

// 	var extractedData, extractedDataFeat []string
// 	var splitedGenome []string
// 	var originBlock, CDSblock, firstCDS, nonCDSFeatStart, nonCDSFeatEnd int

// 	var changedStr, organismName, organismStrain, changedStrFeat string
// 	var noteBuffer, noteBufferFeat strings.Builder

// 	// if _, err := os.Stat(file); os.IsNotExist(err) {
// 	// 	fmt.Printf("The %v file is not exist!\n", file)
// 	// 	os.Exit(3)
// 	// }
// 	fmt.Println(logo)
// 	f, err := os.Open(file) // открываем файл

// 	if err != nil {
// 		fmt.Println(err)

// 	}
// 	defer f.Close()
// 	extractedData = append(extractedData, "START_BLOCK") // начало блока записи
// 	extractedDataFeat = append(extractedData, "START_BLOCK")
// 	scanner := bufio.NewScanner(f) //  новый сканер

// 	for scanner.Scan() {

// 		scanTxt := scanner.Text()
// 		// размер генома
// 		for _, gStartEnd := range genomeStartEnd.FindAllStringSubmatch(scanTxt, -1) {
// 			gStart, _ = strconv.Atoi(gStartEnd[1])
// 			gEnd, _ = strconv.Atoi(gStartEnd[2])
// 		}
// 		//  организм
// 		for _, gName := range genomeSource.FindAllStringSubmatch(scanTxt, -1) {
// 			organismName = gName[1]
// 			fmt.Printf("Organism:%v\n", organismName)
// 		}
// 		// штамм
// 		for _, gStrain := range genomeStrain.FindAllStringSubmatch(scanTxt, -1) {
// 			organismStrain = gStrain[1]
// 			fmt.Printf("Strain:%v\n%v%v...\r", organismStrain, "Working on ", file)
// 		}

// 		// поиск начала CDS блока. если найден, флаг = 1
// 		// for _, cdsFound := range checkCDS.FindAllStringSubmatch(scanTxt, -1) {
// 		// 	if len(cdsFound[1]) != 0 || len(cdsFound[2]) != 0 {
// 		// 		CDSblock = 1
// 		// 		// nonCDSFeatStart = 0

// 		// 		// fmt.Println(scanTxt)
// 		// 	}
// 		// }
// 		// // выход из CDS блока, найден  флаг = 0 (определяется как обнаружение начала блока gene)
// 		// for _, cdsEnd := range endOfCDS.FindAllStringSubmatch(scanTxt, -1) {
// 		// 	if len(cdsEnd[1]) != 0 || len(cdsEnd[2]) != 0 {
// 		// 		CDSblock = 0
// 		// 		nonCDSFeatEnd = 1
// 		// 		nonCDSFeatStart = 0
// 		// 		// fmt.Println(scanTxt)
// 		// 	}
// 		// }
// 		// начало блока Origin (нуклеотидной последовательности)
// 		for _, sOrigin := range startOrigin.FindAllStringSubmatch(scanTxt, -1) {
// 			if len(sOrigin[1]) != 0 || len(sOrigin[2]) != 0 {
// 				CDSblock = 0
// 			}
// 		}

// 		for _, rfeat := range feature.FindAllStringSubmatch(scanTxt, -1) {
// 			if rfeat[1] != "gene" && rfeat[4] != "gene" && rfeat[1] != "CDS" && rfeat[4] != "CDS" && strings.Contains(rfeat[0], "source") == false {
// 				CDSblock = 0
// 				nonCDSFeatEnd = 0
// 				nonCDSFeatStart = 1
// 				// if strings.Contains(rfeat[0], "source") {
// 				// 	fmt.Println(rfeat[0])
// 				// }
// 				if rfeat[1] != "" {
// 					scanTxt = strings.Replace(string(scanTxt), rfeat[1], fmt.Sprintf("*%v", rfeat[1]), -1)
// 				} else if rfeat[4] != "" {
// 					scanTxt = strings.Replace(string(scanTxt), rfeat[4], fmt.Sprintf("*%v", rfeat[4]), -1)

// 				}

// 			} else if rfeat[1] == "CDS" || rfeat[4] == "CDS" && strings.Contains(rfeat[0], "source") == false {
// 				// fmt.Println(rfeat[1], rfeat[4], "!!!")
// 				CDSblock = 1
// 			} else if rfeat[1] == "gene" || rfeat[4] == "gene" && strings.Contains(rfeat[0], "source") == false {
// 				// fmt.Println(rfeat[1], rfeat[4], "???")
// 				CDSblock = 0
// 				nonCDSFeatEnd = 1
// 				nonCDSFeatStart = 0
// 			}
// 		}

// 		for _, origninFound := range checkOrigin.FindAllStringSubmatch(scanTxt, -1) {
// 			if len(origninFound[1]) != 0 {
// 				originBlock = 1

// 			} else {
// 				originBlock = 0
// 			}
// 		}

// 		if nonCDSFeatStart == 1 && nonCDSFeatEnd == 0 {
// 			// changedStrFeat = strings.TrimSuffix(changedStrFeat, "\n")
// 			// changedStrFeat = replaceCodonStart.ReplaceAllString(string(scanTxt), "")               // удаляем пробелы
// 			changedStrFeat = regDelSpaces.ReplaceAllString(string(scanTxt), "$2 ")                 // удаляем пробелы
// 			changedStrFeat = strings.Replace(string(changedStrFeat), "/note=", "!note=", -1)       // меняем / на ! для дальнейшего парсинга
// 			changedStrFeat = strings.Replace(string(changedStrFeat), "/product=", "!product=", -1) // см выше.
// 			changedStrFeat = strings.Replace(string(changedStrFeat), "\"", "", -1)

// 			changedStrFeat = makeAnchors.ReplaceAllString(string(changedStrFeat), "!!")
// 			// if strings.Index(changedStrFeat, "!!") == -1 && strings.Index(changedStrFeat, "!") == -1 && strings.Index(changedStrFeat, "*") == -1 {
// 			// 	noteBufferFeat.WriteString(changedStrFeat)
// 			// } else {
// 			// 	if noteBufferFeat.Len() != 0 {
// 			// 		fmt.Println(noteBufferFeat.String())
// 			// 		noteBufferFeat.Reset()
// 			// 	}
// 			// }

// 			if strings.Index(changedStrFeat, "!!") != 0 && strings.Index(changedStrFeat, "*") != 0 {

// 				noteBufferFeat.WriteString(changedStrFeat)
// 				// fmt.Println(changedStr)

// 			} else {
// 				if noteBufferFeat.Len() != 0 {

// 					// fmt.Printf("%v\n", strings.TrimPrefix(noteBufferFeat.String(), "!"))
// 					// extractedData = append(extractedData, strings.TrimPrefix(noteBuffer.String(), "!"))

// 					// fmt.Println()
// 					extractedDataFeat = append(extractedDataFeat, strings.TrimPrefix(strings.Replace(noteBufferFeat.String(), "!", "\n!", -1), "!"))
// 					noteBufferFeat.Reset()
// 					// fmt.Println(strings.Repeat("*", 30))
// 				}

// 			}

// 			if strings.Index(changedStrFeat, "!!") == 0 {
// 				extractedDataFeat = append(extractedDataFeat, strings.TrimPrefix(changedStrFeat, "!!"))

// 				// fmt.Println(changedStrFeat)
// 			} else if strings.Index(changedStrFeat, "*") == 0 {
// 				extractedDataFeat = append(extractedDataFeat, strings.TrimPrefix(changedStrFeat, "*"))

// 				// fmt.Println(changedStrFeat)
// 			}

// 			// fmt.Println(changedStrFeat)
// 		}
// 		// удаляем перенос строки
// 		changedStr = strings.TrimSuffix(changedStr, "\n")
// 		// changedStr = replaceCodonStart.ReplaceAllString(string(scanTxt), "") // удаляем пробелы
// 		changedStr = regDelSpaces.ReplaceAllString(string(scanTxt), "$2 ") // удаляем пробелы

// 		switch originBlock {
// 		case 0:
// 			if CDSblock == 1 {
// 				if firstCDS == 0 {

// 					firstCDS = 1

// 				}
// 				changedStr = strings.TrimSuffix(changedStr, "\n")
// 				changedStr = strings.Replace(string(changedStr), "/note=", "!note=", -1)       // меняем / на ! для дальнейшего парсинга
// 				changedStr = strings.Replace(string(changedStr), "/product=", "!product=", -1) // см выше.
// 				changedStr = strings.Replace(string(changedStr), "\"", "", -1)

// 				changedStr = makeAnchors.ReplaceAllString(string(changedStr), "!!")
// 				// fmt.Println(changedStr)

// 			}

// 		case 1:

// 			for _, genomeMatch := range qenomeSeq.FindAllStringSubmatch(scanTxt, -1) {

// 				resString = append(resString, strings.Replace(genomeMatch[2], " ", "", -1))
// 			}
// 		}

// 		if CDSblock == 1 {
// 			// fmt.Println(changedStr)
// 			if strings.Index(changedStr, "!!") != 0 && strings.Index(changedStr, "CDS  ") != 0 {
// 				noteBuffer.WriteString(changedStr)
// 				// fmt.Println(changedStr)

// 			} else {
// 				if noteBuffer.Len() != 0 {

// 					// fmt.Printf("%v\n", tmpStr)
// 					extractedData = append(extractedData, strings.TrimPrefix(noteBuffer.String(), "!"))
// 					// fmt.Println()
// 					noteBuffer.Reset()
// 				}

// 			}

// 			if strings.Index(changedStr, "!!") == 0 {
// 				extractedData = append(extractedData, strings.TrimPrefix(changedStr, "!!"))
// 				// fmt.Println(changedStr)

// 			} else if strings.Index(changedStr, "CDS  ") == 0 {

// 				// if firstCDS
// 				if firstCDS == 1 {
// 					extractedData = append(extractedData, "START_OF")
// 					firstCDS = 2
// 				} else if firstCDS == 2 {
// 					extractedData = append(extractedData, "END_OF")
// 					extractedData = append(extractedData, "START_OF")
// 				}

// 				extractedData = append(extractedData, strings.TrimPrefix(changedStr, "!!"))
// 			}

// 		}

// 	}
// 	extractedData = append(extractedData, "END_OF")
// 	extractedData = append(extractedData, "END_BLOCK")
// 	splitedGenome = strings.SplitAfter(strings.Join(resString, ""), "")
// 	// for i, val := range extractedData {
// 	// 	fmt.Printf("%v %v", i, val)
// 	// }
// 	var lStart, lEnd, lDir, lLoc, lName, lProd, lNote, gID, pID, lGOA string
// 	var nucCore string
// 	var cdsStEnd = regexp.MustCompile(`^CDS\s+(\d+)\W+(\d+)|^CDS\s+complement\W(\d+)\W+(\d+)`)
// 	var cdsLocus = regexp.MustCompile(`locus_tag=(.*)`)
// 	var cdsName = regexp.MustCompile(`gene=(.*)`)
// 	var cdsProd = regexp.MustCompile(`product=(.*)`)
// 	var cdsNote = regexp.MustCompile(`note=(.*)`)
// 	var cdsgID = regexp.MustCompile(`db_xref=GeneID:(\d+)`)
// 	var cdsprotID = regexp.MustCompile(`protein_id=(.*)`)
// 	var cdsGOA = regexp.MustCompile(`db_xref=GOA:(.*)`)

// 	var startOfBlock, endOfBlock, cdsStart, cdsEnd int
// 	// var cdsOpen, cdsClosed int
// 	cdsCount := 1
// 	// fmt.Println(extractedData)

// 	for _, val := range extractedData {

// 		if strings.Index(val, "START_BLOCK") != -1 {
// 			startOfBlock = 1
// 			endOfBlock = 0
// 		} else if strings.Index(val, "END_BLOCK") != -1 {
// 			endOfBlock = 1
// 			startOfBlock = 0
// 		}
// 		if strings.Index(val, "START_OF") != -1 {
// 			cdsStart = 1
// 			cdsEnd = 0
// 			cdsCount++
// 		} else if strings.Index(val, "END_OF") != -1 {
// 			cdsStart = 0
// 			cdsEnd = 1
// 		}

// 		if startOfBlock == 1 && endOfBlock == 0 {

// 			if cdsStart == 1 && cdsEnd == 0 {

// 				for _, cdsStEndMatch := range cdsStEnd.FindAllStringSubmatch(val, -1) {

// 					if len(cdsStEndMatch[1]) != 0 {
// 						lStart, lEnd, lDir = strings.TrimSpace(cdsStEndMatch[1]), strings.TrimSpace(cdsStEndMatch[2]), "f"
// 						// cdsOpen = 1

// 					} else if len(cdsStEndMatch[3]) != 0 {

// 						lStart, lEnd, lDir = strings.TrimSpace(cdsStEndMatch[3]), strings.TrimSpace(cdsStEndMatch[4]), "r"
// 						// cdsOpen = 1

// 					}

// 				}

// 				for _, cdsLocusMatch := range cdsLocus.FindAllStringSubmatch(val, -1) {
// 					lLoc = strings.Replace(cdsLocusMatch[1], " ", "", -1)
// 					// fmt.Println(lLoc)
// 					// cdsOpen = 1
// 					// cdsClosed = 0
// 				}

// 				// cdsLocusMatch := cdsLocus.FindAllStringSubmatch(val, -1)
// 				// if len(cdsLocusMatch) != 0 {
// 				// 	lLoc = cdsLocusMatch[0][1]
// 				// }

// 				for _, cdsNameMatch := range cdsName.FindAllStringSubmatch(val, -1) {
// 					lName = strings.Replace(cdsNameMatch[1], " ", "", -1)
// 					// cdsOpen = 1
// 					// cdsClosed = 0
// 					// fmt.Println(lName)
// 				}
// 				for _, cdsProdMatch := range cdsProd.FindAllStringSubmatch(val, -1) {
// 					lProd = cdsProdMatch[1]
// 					// cdsOpen = 1
// 					// cdsClosed = 0
// 					// fmt.Println(lProd)
// 				}
// 				for _, cdsNoteMatch := range cdsNote.FindAllStringSubmatch(val, -1) {
// 					lNote = cdsNoteMatch[1]

// 					// fmt.Println(lNote)
// 				}
// 				for _, cdsgIDMatch := range cdsgID.FindAllStringSubmatch(val, -1) {
// 					gID = strings.Replace(cdsgIDMatch[1], " ", "", -1)
// 					// cdsOpen = 1
// 					// cdsClosed = 0
// 					// fmt.Println(gID)

// 				}
// 				for _, cdsgprotIDMatch := range cdsprotID.FindAllStringSubmatch(val, -1) {
// 					pID = strings.Replace(cdsgprotIDMatch[1], " ", "", -1)
// 					// cdsOpen = 1
// 					// cdsClosed = 0
// 					// fmt.Println(pID)

// 				}
// 				for _, cdsGOAMatch := range cdsGOA.FindAllStringSubmatch(val, -1) {
// 					lGOA = strings.Replace(cdsGOAMatch[1], " ", "", -1)

// 					// fmt.Println(lGOA)
// 					// cdsOpen = 1
// 					// cdsClosed = 0
// 				}

// 				if *flgDev == true {
// 					fmt.Println(val)
// 				}

// 				if len(gID) == 0 {
// 					nucCore = "yes"

// 				} else {
// 					nucCore = "no"
// 				}
// 				// fmt.Println(nucCore)

// 				if lName == "" {
// 					lName = lLoc
// 				}

// 			} else if cdsStart == 0 && cdsEnd == 1 {
// 				g := gene.Gene{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir,
// 					Product: lProd, Name: lName, GeneID: gID, ProteinID: pID, Note: lNote, GOA: lGOA}

// 				allGenesVal = append(allGenesVal, g)
// 				gcoordStart, _ := strconv.Atoi(lStart)
// 				gcoordEnd, _ := strconv.Atoi(lEnd)

// 				gcoords := genomeCoords{Start: gcoordStart, End: gcoordEnd, Name: lLoc, Product: lProd, Direction: lDir, TypeOf: "CDS"}

// 				genomeCoordinates = append(genomeCoordinates, gcoords)

// 				// fmt.Printf("g:%v b:%v %v\n", cdsOpen, cdsClosed, val)
// 				if *flgDebug == true {

// 					fmt.Printf("l:%v s:%v e:%v d:%v p:%v gId:%v pId:%v n:%v GOA:%v\n", lLoc, lStart, lEnd, lDir, lProd, gID, pID, lNote, lGOA)
// 					// }
// 				} else {
// 					fmt.Printf("Found %v CDS, nucCore:%v                   \r", cdsCount, nucCore)
// 				}
// 				lName, lStart, lEnd, lLoc, lDir, lProd, gID, pID, lGOA, lNote = "", "", "", "", "", "", "", "", "", ""

// 				// fmt.Printf("s:%v e:%v d:%v o:%v i:%v val:%v\n", lStart, lEnd, lDir, cdsOpen, i, val)

// 			}
// 		}

// 		// if cdsOpen != 0 {
// 		// 	g := gene.Gene{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir,
// 		// 		Product: lProd, Name: lName, GeneID: gID, ProteinID: pID, Note: lNote, GOA: lGOA}

// 		// 	allGenesVal = append(allGenesVal, g)
// 		// 	// fmt.Printf("g:%v b:%v %v\n", cdsOpen, cdsClosed, val)
// 		// 	if *flgDebug == true {

// 		// 		fmt.Printf("l:%v s:%v e:%v d:%v p:%v gId:%v pId:%v n:%v GOA:%v\n", lLoc, lStart, lEnd, lDir, lProd, gID, pID, lNote, lGOA)
// 		// 	}
// 		// 	// } else {
// 		// 	// 	fmt.Printf("\rFound %v CDS, nucCore:%v", cdsCount-1, nucCore)
// 		// 	// }
// 		// 	lName, lStart, lEnd, lLoc, lDir, lProd, gID, pID, lGOA, lNote = "", "", "", "", "", "", "", "", "", ""

// 		// }

// 		// fmt.Println(val, cdsOpen, cdsClosed)

// 		// if cdsOpen == 1 && cdsClosed == 1 && lStart != "" && lEnd != "" {
// 		// 	// fmt.Println(lStart, lEnd, lDir, lLoc, lProd, cdsOpen, cdsClosed)
// 		// 	// cdsCount++
// 		// 	g := gene.Gene{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir,
// 		// 		Product: lProd, Name: lName, GeneID: gID, ProteinID: pID, Note: lNote, GOA: lGOA}

// 		// 	allGenesVal = append(allGenesVal, g)
// 		// 	// fmt.Printf("g:%v b:%v %v\n", cdsOpen, cdsClosed, val)
// 		// 	if *flgDebug == true {

// 		// 		fmt.Printf("l:%v s:%v e:%v d:%v p:%v gId:%v pId:%v n:%v GOA:%v\n", lLoc, lStart, lEnd, lDir, lProd, gID, pID, lNote, lGOA)
// 		// 	}
// 		// 	// } else {
// 		// 	// 	fmt.Printf("\rFound %v CDS, nucCore:%v", cdsCount-1, nucCore)
// 		// 	// }
// 		// 	lName, lStart, lEnd, lLoc, lDir, lProd, gID, pID, lGOA, lNote = "", "", "", "", "", "", "", "", "", ""

// 		// }

// 	}

// 	var featStEnd = regexp.MustCompile(`^(\w+)\s+(\d+)\W+(\d+)|^(\w+)\s+complement\W(\d+)\W+(\d+)`)
// 	var featLocus = regexp.MustCompile(`locus_tag=(.*)`)
// 	var featName = regexp.MustCompile(`gene=(.*)`)
// 	var lType string
// 	var featProd = regexp.MustCompile(`product=(.*)`)
// 	var featNote = regexp.MustCompile(`note=(.*)`)
// 	var featCount = map[string]int{}
// 	for _, feat := range extractedDataFeat {
// 		for _, stend := range featStEnd.FindAllStringSubmatch(feat, -1) {
// 			lType = stend[1]

// 			if stend[2] != "" {
// 				lStart, lEnd, lDir = strings.TrimSpace(stend[2]), strings.TrimSpace(stend[3]), "f"
// 				// fmt.Println(lStart, lEnd, lDir)
// 			} else if stend[5] != "" {
// 				lStart, lEnd, lDir = strings.TrimSpace(stend[5]), strings.TrimSpace(stend[6]), "r"
// 				// fmt.Println(lStart, lEnd, lDir)
// 			}
// 		}

// 		for _, gene := range featName.FindAllStringSubmatch(feat, -1) {
// 			if gene[1] != "" {
// 				lName = gene[1]
// 				// if lName == "" {
// 				// 	lName = lType
// 				// }
// 				// fmt.Println(lName)
// 			}
// 		}

// 		for _, loctag := range featLocus.FindAllStringSubmatch(feat, -1) {
// 			if loctag[1] != "" {
// 				lLoc = loctag[1]
// 				// if lLoc == "" {
// 				// 	lLoc = lName
// 				// }
// 				// fmt.Println(lLoc)
// 			}
// 		}

// 		for _, prod := range featProd.FindAllStringSubmatch(feat, -1) {
// 			if prod[1] != "" {
// 				lProd = prod[1]
// 			}
// 		}
// 		for _, note := range featNote.FindAllStringSubmatch(feat, -1) {
// 			if note[1] != "" {
// 				lNote = note[1]
// 			}
// 		}
// 		// if lLoc == "" && lName == "" {
// 		// 	lLoc = "-"
// 		// 	lType = "-"
// 		// }

// 		if lType != "" {
// 			featCount[lType] = featCount[lType] + 1
// 			// featG := featInfo{TypeOf: lType, Start: lStart, End: lEnd, Name: lName, Product: lProd, Note: lNote, Locus: lLoc}
// 			// featAllGenesVal = append(featAllGenesVal, featG)
// 			resName := lName
// 			if resName == "" {
// 				resName = lLoc
// 				if lLoc == "" {
// 					resName = lType
// 				}
// 			}
// 			gcoordStart, _ := strconv.Atoi(lStart)
// 			gcoordEnd, _ := strconv.Atoi(lEnd)

// 			gcoords := genomeCoords{Start: gcoordStart, End: gcoordEnd, Name: resName, Product: lProd, Direction: lDir, TypeOf: lType}
// 			fmt.Println(gcoords)
// 			genomeCoordinates = append(genomeCoordinates, gcoords)
// 			sort.Slice(genomeCoordinates, func(i, j int) bool {
// 				return genomeCoordinates[i].Start < genomeCoordinates[j].Start
// 			})
// 			// fmt.Printf("%v\t%v\t%v\t%v\t%v\t%v\t%v\n", lType, lStart, lEnd, lName, lLoc, lProd, lNote)
// 			lType, lStart, lEnd, lName, lLoc, lProd, lNote = "", "", "", "", "", "", ""
// 		}

// 	}
// 	for key, val := range featCount {
// 		fmt.Printf("Found %v %v\n", val, key)
// 	}
// 	gInfo = genomeInfo{TypeOf: nucCore, Organism: organismName, Start: gStart, End: gEnd, Strain: organismStrain}
// 	// go process("Анализ файла закончен!")

// 	return allGenesVal, splitedGenome

// }

func readGBFile(file string, verbose bool) (g []gene.Gene, genomeSplice []string) {
	//функция для считывания файла в формате генбанка и занесения строк в массив linesFromGB

	var qenomeSeq = regexp.MustCompile(`\s+(\d+)\s+(\w{6}.*)`)
	var regDelSpaces = regexp.MustCompile(`(\s+)(.*)`)
	var checkOrigin = regexp.MustCompile(`(\d+\s\w{10}\s\w{10}\s\w{10})`)
	// var checkCDS = regexp.MustCompile(`(CDS\s+complement\W\d+\W+\d+\W)|(CDS\s+\d+\W+\d+)`)
	// var endOfCDS = regexp.MustCompile(`(gene\s+\d+\W+\d+)|(gene\s+complement\W\d+\W+\d+)`)
	var startOrigin = regexp.MustCompile(`(ORIGIN\s+)|(translation=)`)
	var feature = regexp.MustCompile(`^\s{5}(\w+)\s+complement\W(\d)+\W+(\d)+\W|^\s{5}(\w+)\s+(\d)+\W+(\d+)`)
	var makeAnchors = regexp.MustCompile(`^(/)`)
	var genomeSource = regexp.MustCompile(`organism=\W(.*)\W`)
	var genomeStrain = regexp.MustCompile(`strain=\W(.*)\W`)
	var genomeVersion = regexp.MustCompile(`^VERSION\s+(.*)`)
	var genomeStartEnd = regexp.MustCompile(`source\W+(\d+)..(\d+)`)

	var gStart, gEnd int
	var gVersion string

	//  source          1..4411532
	//                  /organism="Mycobacterium tuberculosis H37Rv"
	//                  /mol_type="genomic DNA"
	//                  /strain="H37Rv"
	//                  /db_xref="taxon:83332"
	// var replaceCDS = regexp.MustCompile(`(CDS\s+.*)`)

	var resString []string

	var extractedData, extractedDataFeat []string
	var splitedGenome []string
	var originBlock, CDSblock, firstCDS int

	var changedStr, organismName, organismStrain string
	var noteBuffer strings.Builder

	// if _, err := os.Stat(file); os.IsNotExist(err) {
	// 	fmt.Printf("The %v file is not exist!\n", file)
	// 	os.Exit(3)
	// }

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
			fmt.Printf("Vsersion:%v\n", gVersion)
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

	var lDir, lLoc, lName, lProd, lNote, gID, pID, lGOA, lType string
	var lPDBArr, lInterProArr, lProSiteArr []string
	var nucCore int
	var cdsStEnd = regexp.MustCompile(`^(\w+)\s+(\d+)\W+(\d+)|^(\w+)\s+complement\W(\d+)\W+(\d+)`)
	var cdsLocus = regexp.MustCompile(`locus_tag=(.*)`)
	var cdsName = regexp.MustCompile(`gene=(.*)`)
	var cdsProd = regexp.MustCompile(`product=(.*)`)
	var cdsNote = regexp.MustCompile(`note=(.*)`)
	var cdsgID = regexp.MustCompile(`db_xref=GeneID:(\d+)`)
	var cdsprotID = regexp.MustCompile(`protein_id=(.*)`)
	var cdsGOA = regexp.MustCompile(`db_xref=GOA:(.*)`)
	var cdsPDB = regexp.MustCompile(`db_xref=PDB:(.*)`)
	var cdsInterPro = regexp.MustCompile(`db_xref=InterPro:(.*)`)
	var cdsProSite = regexp.MustCompile(`inference=protein motif:PROSITE:(.*)`)

	var startOfBlock, endOfBlock, cdsStart, cdsEnd, lStart, lEnd int
	// var cdsOpen, cdsClosed int
	var cdsCount = map[string]int{}
	// fmt.Println(extractedData)

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
						lDir = "f"

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

				if lName == "" {
					lName = lLoc
				}

			} else if cdsStart == 0 && cdsEnd == 1 {

				if lLoc == "" && lName != "" {
					lLoc = lName
				} else if lLoc == "" && lName == "" {
					lLoc = fmt.Sprintf("%v_%v_%v", lType, lStart, lEnd)

				}

				if lProd == "" && lNote != "" {
					lProd = lNote
				} else if lProd == "" && lNote == "" {
					lProd = lType
				}

				g := gene.Gene{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir,
					Product: lProd, Name: lName, GeneID: gID, ProteinID: pID, Note: lNote, GOA: lGOA, TypeOf: lType, PDB: lPDBArr, InterPro: lInterProArr, ProSite: lProSiteArr}

				allGenesVal = append(allGenesVal, g)

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

	var igensS []int
	var igensE []int
	var leftGene []string
	var igrCount int

	for i, g := range allGenesVal {

		igensS = append(igensS, g.Start)
		igensE = append(igensE, g.End)
		leftGene = append(leftGene, g.Locus)

		if i >= 1 && i < len(allGenesVal)-1 {

			igrCount++
			checkEndOfCDS := (igensE[i-1] + 1) - g.Start

			if checkEndOfCDS < 0 {
				igr := gene.Gene{Start: igensE[i-1] + 1, End: igensS[i] - 1, Locus: fmt.Sprintf("IGR_%v_%v", igensE[i-1]+1, igensS[i]-1), Direction: "f", Product: "Intergenic region (IGR)", Name: "IGR", TypeOf: "IGR"}
				allGenesVal = append(allGenesVal, igr)
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

func writeDB(file string, gene []gene.Gene) {

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
	// gobParser.Encode(&featAllGenesVal)
	// gobParser.Encode(&genomeCoordinates)

	fmt.Println(file, " was created successfully.")

}

func readDB(file string) []gene.Gene {
	var gene []gene.Gene
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
	// gobParser.Decode(&featAllGenesVal)
	// gobParser.Decode(&genomeCoordinates)
	return gene

}

func codonReverse(codon string) string {
	var lCodonSplit []string
	var result string
	lCodonSplit = strings.Split(codon, "")
	// fmt.Println(lCodonSplit)
	for _, n := range lCodonSplit {

		switch n {
		case "a":
			result = "t" + result
		case "t":
			result = "a" + result
		case "c":
			result = "g" + result
		case "g":
			result = "c" + result
		case "A":
			result = "T" + result
		case "T":
			result = "A" + result
		case "C":
			result = "G" + result
		case "G":
			result = "C" + result
		}

	}
	return result

}

//

func parserVCF(f string, print bool, genes []gene.Gene) []gene.SNPinfo {
	var vcf = regexp.MustCompile(`^\S+\s+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)`)
	// var indel = regexp.MustCompile(`^^\S+\W+(\d+)\W+(\w+)\s+(\w+).*(INDEL).*DP=(\d+)`)
	var validateVCF = regexp.MustCompile(`(##fileformat)=VCF`)
	var vcfValid bool

	//
	var snpFromVCF []gene.SNPinfo

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
			// for _, indelMatch := range indel.FindAllStringSubmatch(scanner.Text(), -1) {
			// 	if vcfValid == true && indelMatch[4] == "INDEL" {
			// 		fmt.Println(indelMatch[1], indelMatch[2], indelMatch[3])

			// 	}

			// }
			// fmt.Println("Under construction")
		}

		for _, match := range vcf.FindAllStringSubmatch(scanner.Text(), -1) {

			if vcfValid == true {
				apos, _ := strconv.Atoi(match[1])
				ref := match[2]
				alt := match[3]

				for _, g := range genes {

					lStart := g.Start
					lEnd := g.End

					if apos >= lStart && apos <= lEnd {
						snp, err := getSNPInfo(apos, g, alt, *gbIndex)
						snp.InterPro = g.InterPro
						snp.PDB = g.PDB
						snp.ProSite = g.ProSite
						// fmt.Println(g.PDB)
						// br := testing.Benchmark(snp)
						// fmt.Println(br)

						if len(ref) == 1 && len(alt) == 1 && err != 1 {
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

func process(str string) {
	fmt.Printf("%v\r", str)
}

func listOfVCFFiles(withFilenames bool) {

	files := getListofVCF()
	for _, file := range files {
		if strings.Contains(file, ".vcf") {
			if withFilenames == true {
				fmt.Printf("\n\n%v:\n\n", file)
			}

			parserVCF(file, true, allGenesVal)

		}
	}
}

func makeSeq(typeof string, verbose bool, ref bool) []seqInfo {

	var AllPosUnsort, AllPos []int
	var ResSeq []seqInfo

	files := getListofVCF()

	for i, file := range files {

		if verbose == true {
			fmt.Printf("Reading files: Read %v from %v \r", i+1, len(files))
		}
		snps := parserVCF(file, false, allGenesVal)
		switch typeof {
		case ncFlag:
			for _, val := range snps {
				AllPosUnsort = append(AllPosUnsort, val.APos)

			}
		case aaFlag:

		}

	}
	// go process("Working...          ")
	AllPos = unique(AllPosUnsort)
	sort.Ints(AllPos)
	if ref == true {
		switch typeof {
		case ncFlag:
			var refBuffer strings.Builder
			refBuffer.WriteString(fmt.Sprintf(">%v\n", "REFERENCE"))
			for _, allpos := range AllPos {
				refBuffer.WriteString(getNucFromGenomePos(allpos))
			}
			ResSeq = append(ResSeq, seqInfo{Name: "reference", Seq: refBuffer.String()})

		}
	}
	for i, file := range files {
		fmt.Printf("Generating sequences: Working on  %v from %v \r", i+1, len(files))
		pos := make(map[int]string)
		var buffer strings.Builder

		buffer.WriteString(fmt.Sprintf(">%v\n", strings.ToUpper(file)))

		snps := parserVCF(file, false, allGenesVal)
		switch typeof {
		case ncFlag:
			for _, val := range snps {
				pos[val.APos] = val.Alt

			}
			for _, allpos := range AllPos {

				if pos[allpos] != "" {

					buffer.WriteString(pos[allpos])
				} else {

					buffer.WriteString(getNucFromGenomePos(allpos))
				}
			}
			ResSeq = append(ResSeq, seqInfo{Name: file, Seq: buffer.String()})

		case aaFlag:

			fmt.Println("Reserved")

		}

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
	var igensS []int
	var igensE []int
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

func printTextResults(snps gene.SNPinfo, verbose bool) {

	const fullAnnotations = " {{if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T0\"))}}" +
		"{{.Locus}}\t{{.APos}}\tc.{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\tp.{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T1\"))}}" +
		"{{.Locus}}\t{{.APos}}\tc.{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\tp.{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Tang}}({{.TangIdxVal}})\t{{.Product}}\n" +
		"{{else if (and (ne .TypeOf \"CDS\") (eq .ReportType \"T0\"))}}" +
		"{{.Locus}}\t{{.APos}}\tg.{{.APos}}{{.NucInPos}}>{{.Alt}}\t-\t-\t-\t-\t{{.Product}}\n" +
		"{{else if (and (ne .TypeOf \"CDS\") (eq .ReportType \"T1\"))}}" +
		"{{.Locus}}\t{{.APos}}\tg.{{.APos}}{{.NucInPos}}>{{.Alt}}\t-\t-\t-\t-\t-\t{{.Product}}\n" +
		"{{end}}"

	const cdsAnnotations = "{{if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T0\"))}}" +
		"{{.Locus}}\t{{.APos}}\tc.{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\tp.{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T1\"))}}" +
		"{{.Locus}}\t{{.APos}}\tc.{{.PosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\tp.{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Tang}}({{.TangIdxVal}})\t{{.Product}}\n" +
		"{{end}}"

	t := template.New("report")
	if verbose == true {
		t, _ := t.Parse(fullAnnotations)
		t.Execute(os.Stdout, snps)

	} else {
		t, _ := t.Parse(cdsAnnotations)
		t.Execute(os.Stdout, snps)
	}

}

func printWebResults(snps []gene.SNPinfo, port string) {

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
				<td bgcolor="#CECEF6"><p title="Complex Index: {{.Tang}}({{.TangIdxVal}})">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}" target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ebi.ac.uk/QuickGO/GProtein?ac={{.GOA}}" target="_blank">{{.GOA}}</a>
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

func getSNPNotation(f string) []gene.SNPcheck {
	var snpCheck gene.SNPcheck
	var parsedSNP []gene.SNPcheck
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
	file, err := os.Open(f)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {

		for _, matchLMN := range rLMN.FindAllStringSubmatch(scanner.Text(), -1) {
			snpCheck = gene.SNPcheck{Locus: strings.ToUpper(matchLMN[1]), PosInGene: matchLMN[2], Ref: matchLMN[3], Alt: matchLMN[4], Name: matchLMN[5], TypeOf: tLMN, Raw: matchLMN[0]}
			parsedSNP = append(parsedSNP, snpCheck)
			// fmt.Printf("%v\n", snpCheck)
			// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
		}
		for _, matchPMLN := range rPMLN.FindAllStringSubmatch(scanner.Text(), -1) {
			snpCheck = gene.SNPcheck{APos: matchPMLN[1], Ref: matchPMLN[2], Alt: matchPMLN[3], Locus: matchPMLN[4], Name: matchPMLN[5], TypeOf: tPMLN, Raw: matchPMLN[0]}
			parsedSNP = append(parsedSNP, snpCheck)
			// fmt.Printf("%v\n", snpCheck)
			// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
		}
		for _, matchPMN := range rPMN.FindAllStringSubmatch(scanner.Text(), -1) {

			snpCheck = gene.SNPcheck{APos: matchPMN[1], Ref: matchPMN[2], Alt: matchPMN[3], Name: matchPMN[4], TypeOf: tPMN, Raw: matchPMN[0]}

			parsedSNP = append(parsedSNP, snpCheck)

			// fmt.Printf("%v\n", snpCheck)
			// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
		}

		// for _, matchPMNT := range rPMNT.FindAllStringSubmatch(scanner.Text(), -1) {

		// 	snpCheck = gene.SNPcheck{APos: matchPMNT[1], Ref: matchPMNT[2], Alt: matchPMNT[3], Name: matchPMNT[4], TypeOf: tPMNT, Raw: matchPMNT[0], Tag: matchPMNT[6]}

		// 	parsedSNP = append(parsedSNP, snpCheck)

		// 	// fmt.Println("!!!")
		// 	// fmt.Printf("%v\n", snpCheck.Tag)
		// 	// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
		// }

		for _, matchLSAAN := range rLSAAN.FindAllStringSubmatch(scanner.Text(), -1) {
			snpCheck = gene.SNPcheck{Locus: strings.ToUpper(matchLSAAN[1]), AASref: matchLSAAN[2], CodonNbrInG: matchLSAAN[3], AASalt: matchLSAAN[4], Name: matchLSAAN[5], TypeOf: tLSAAN, Raw: matchLSAAN[0]}
			parsedSNP = append(parsedSNP, snpCheck)
			// fmt.Printf("%v\n", snpCheck)
			// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
		}
		for _, matchLLAAN := range rLLAAN.FindAllStringSubmatch(scanner.Text(), -1) {
			snpCheck = gene.SNPcheck{Locus: strings.ToUpper(matchLLAAN[1]), AALref: matchLLAAN[2], CodonNbrInG: matchLLAAN[3], AALalt: matchLLAAN[4], Name: matchLLAAN[5], TypeOf: tLLAAN, Raw: matchLLAAN[0]}
			parsedSNP = append(parsedSNP, snpCheck)
			// fmt.Printf("%v\n", snpCheck)
			// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
		}
		for _, matchLCN := range rLCN.FindAllStringSubmatch(scanner.Text(), -1) {
			snpCheck = gene.SNPcheck{Locus: strings.ToUpper(matchLCN[1]), CodonNbrInG: matchLCN[2], Name: matchLCN[3], TypeOf: tLCN, Raw: matchLCN[0]}
			parsedSNP = append(parsedSNP, snpCheck)
			// fmt.Printf("%v\n", snpCheck)
			// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
		}

	}

	return parsedSNP
}

func getShareSNP(verbose bool, web bool, files []string) {
	// var countSNPs = 1
	var pos = map[uint64]int{}
	var alt = map[uint64]gene.SNPinfo{}
	var share []uint64
	// var g gene.Gene
	var snpToWeb []gene.SNPinfo
	var snpToConsole []gene.SNPinfo

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
		// bar.Increment()
		snps := parserVCF(file, false, allGenesVal)
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
			// 	g = gene.Gene{Locus: val.Locus, Start: val.Start, End: val.End, Name: val.Name, Product: val.Product, Direction: val.Direction, GeneI } else if lPos > upperLimit {
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

// func getShareSNPTest(verbose bool, web bool, files []string) {
// 	// c := make(chan int)
// 	// var countSNPs = 1
// 	// var pos = map[int]int{}
// 	// var alt = map[int]gene.SNPinfo{}
// 	// var share []int
// 	// var g gene.Gene
// 	// var snpToWeb []gene.SNPinfo
// 	// var snpToConsole []gene.SNPinfo

// 	// upperLimit := len(files)
// 	// bar := pb.StartNew(len(files))
// 	// bar := pb.ProgressBarTemplate(pbtmpl).Start(len(files))
// 	// bar.SetWidth(90)
// 	// tmp := 0
// 	// var allSNPs = map[string][]gene.SNPinfo{}

// 	for i, file := range files {
// 		// tmp++
// 		// fmt.Printf("processed %v files\r", cnt+1)
// 		if verbose == true {
// 			fmt.Printf("Pass: %v from %v \r", i+1, len(files))
// 		}
// 		// bar.Increment()
// 		snps := parserVCF(file, false, allGenesVal)

// 		// getSNPFromSNPInfo(snps, c)
// 		// x := <-c
// 		for _, val := range snps {
// 			if val.APos == 2163790 {
// 				fmt.Println(val, file)
// 			}
// 			// fmt.Println(val.APos)
// 		}

// 		// allSNPs[file] = snps
// 	}
// 	// for _, val := range allSNPs {
// 	// 	for _, snp := range val {
// 	// 		fmt.Println(snp.APos)
// 	// 	}
// 	// }
// 	// fmt.Println(allSNPs)
// 	// for key := range allSNPs {
// 	// 	for i := 0; i <= len(allSNPs[key])-1; i++ {
// 	// 		fmt.Println(allSNPs[key][i].APos)
// 	// 	}

// 	// }
// 	// fmt.Println(tmp, "2163790")
// 	// 	for _, val := range snps {
// 	// 		pos[val.APos] = pos[val.APos] + 1
// 	// 		fmt.Println(val.APos, file, pos[val.APos])

// 	// 		alt[val.APos] = val
// 	// 		// if pos[val.APos] > upperLimit {
// 	// 		// 	fmt.Println(pos[val.APos], val.APos)
// 	// 		// }
// 	// 		// 	g = gene.Gene{Locus: val.Locus, Start: val.Start, End: val.End, Name: val.Name, Product: val.Product, Direction: val.Direction, GeneID: val.GeneID, ProteinID: val.ProteinID}
// 	// 		// }
// 	// 	}

// 	// 	for i, lPos := range pos {
// 	// 		if lPos == upperLimit {
// 	// 			share = append(share, i)
// 	// 		}
// 	// 	}
// 	// 	sort.Ints(share)
// 	// 	for _, sharePos := range share {
// 	// 		countSNPs++
// 	// 		// fmt.Println(alt[sharePos])

// 	// 		if web == true {
// 	// 			snpToWeb = append(snpToWeb, alt[sharePos])
// 	// 			// printResults(alt[sharePos])

// 	// 		} else {
// 	// 			// printResults(alt[sharePos])
// 	// 			snpToConsole = append(snpToConsole, alt[sharePos])
// 	// 		}
// 	// 	}
// 	// }
// 	// // fmt.Println()
// 	// // if verbose == true {
// 	// // 	fmt.Printf("f:%v snp: %v\n%v\n", len(files), countSNPs, files)
// 	// // }
// 	// if web == true && len(snpToWeb) != 0 {
// 	// 	printWebResults(snpToWeb, *gbPort)
// 	// } else if web == false && len(snpToConsole) != 0 {
// 	// 	for _, res := range snpToConsole {
// 	// 		// printResults(res)
// 	// 		printTextResults(res, verbose)
// 	// 	}
// 	// }

// }

func getSNPFromSNPInfo(snps []gene.SNPinfo, c chan int) {
	for _, val := range snps {
		c <- val.APos
	}
}

func snpStat() {
	// var countSNPs = 1
	var pos = make(map[int]int)
	var alt = make(map[int]gene.SNPinfo)
	var f = make(map[int][]string)
	var positions, posUnsort []int

	files := getListofVCF()
	upperLimit := len(files)

	for i, file := range files {

		fmt.Printf("Reading: %v (%v from %v)%v \r", file, i+1, len(files), strings.Repeat(" ", 60))
		// }
		// bar.Increment()
		snps := parserVCF(file, false, allGenesVal)

		for _, val := range snps {
			if pos[val.APos] <= upperLimit {
				pos[val.APos] = pos[val.APos] + 1       //count
				alt[val.APos] = val                     // pos
				f[val.APos] = append(f[val.APos], file) //files
				posUnsort = append(posUnsort, val.APos) //array of positions
			}

		}
		positions = unique(posUnsort)
		sort.Ints(positions)

	}
	var stat []statInfo

	for _, p := range positions {
		perc := (pos[p] * 100) / upperLimit
		filesNotInList := compareSlices(files, f[p])
		stat = append(stat, statInfo{Pos: p, Count: pos[p], Perc: perc, FilesWith: strings.Join(f[p], ",\n"), FilesWithout: strings.Join(filesNotInList, ",\n")})

		// }
	}
	printWebStat(stat, *gbPort)
	// fmt.Println(stat)
}

func calcDnDsVal(file string, print bool) []DnDsRes {

	var altPositions = make(map[string][]allPositionsInGene)
	var validData []string
	var dndsArray []DnDsRes
	var i int
	snps := parserVCF(file, false, allGenesVal)

	for _, val := range snps {
		// fmt.Println(val.Locus, val.PosInGene, val.Alt)
		if val.TypeOf == "CDS" {
			altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, nuc: val.Alt})
		}
		// start, end := getGenePosByName(val.Locus)

		// codonBias[strings.ToUpper(val.AltCodon)] = codonBias[strings.ToUpper(val.AltCodon)] + 1

		// if locusName == val.Locus {

		// 	// switch val.Mutation {
		// 	// case "missense":
		// 	// 	mCount++
		// 	// case "synonymous":
		// 	// 	sCount++

		// 	// }
		// 	allPos = append(allPos, allPositionsInGene{pos: val.PosInGene, nuc: val.Alt})
		// 	// mutRes["m"] = mCount
		// 	// mutRes["s"] = sCount
		// fmt.Println()

		// if amino.CDACodon(val.RefCodon) < 0 && amino.CDACodon(val.AltCodon) > 0 || amino.CDACodon(val.RefCodon) > 0 && amino.CDACodon(val.AltCodon) < 0 {

		// 	fmt.Printf("cda: %v(%v) -> %v(%v)%v %v %v %v %v \n", amino.CDACodon(val.RefCodon), val.RefCodon, amino.CDACodon(val.AltCodon), val.AltCodon, val.Tang, val.Mutation, val.Locus, val.APos, val.CodonNbrInG)

		// }
		//

		// } else {

		// 	locusName = val.Locus
		// 	allPos = nil

		// 	// mCount, sCount = 0, 0

		// }

	}

	// fmt.Println(codonBias)
	// fmt.Println(altPositions)

	for key, val := range altPositions {
		// fmt.Println(key, len(val))
		if len(val) > 1 {
			validData = append(validData, key)
		}
	}
	sort.Strings(validData)
	for _, val := range validData {

		altS := makeAltString(val, altPositions[val])
		start, end := getGenePosByName(val)
		refS := getNucFromGenome(start-1, end)
		dnds := codon.CalcDnDs(refS, altS)
		// fmt.Println(val, "\n", refS)
		prod := getProductByName(val)
		// 	N, S, PN, PS, DN, DS, DNDS float64
		// ND, NS                     int

		if dnds.ND != 0 && dnds.NS != 0 {
			i++
			if print == true {
				// 	fmt.Printf("\nl:%v\tdn/ds:%.2f (%v:%v)\nREF:%v\nALT:%v\n", val, dnds.DNDS, start, end, refS, altS)
				// 	// fmt.Printf("Calculated %v dn/ds from %v\r", i, len(validData))
				// } else if print == true && *gbNoSeq == true {
				fmt.Printf("\nl:%v\tdn/ds:%.2f (%v:%v)\t%v", val, dnds.DNDS, start, end, prod)
			}

			dndsArray = append(dndsArray, DnDsRes{Locus: val, Product: prod, N: dnds.N, S: dnds.S, ND: dnds.ND, NS: dnds.NS, PN: dnds.PN, PS: dnds.PS, DN: dnds.DN, DS: dnds.DS, DNDS: dnds.DNDS})
		}

	}

	return dndsArray

}

func calcJaroWinklerDist(file string, print bool) []JaroWinklerInfo {

	var altPositions = make(map[string][]allPositionsInGene)
	var validData []string
	var altSequences []string
	var jwRes []JaroWinklerInfo
	// var dndsArray []DnDsRes

	var jarwinkl float64
	snps := parserVCF(file, false, allGenesVal)

	for _, val := range snps {
		// fmt.Println(val.Locus, val.PosInGene, val.Alt)
		if val.TypeOf == "CDS" {
			altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, nuc: val.Alt})
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

		altS := makeAltString(val, altPositions[val])
		start, end := getGenePosByName(val)
		altSequences = append(altSequences, altS)
		refS := getNucFromGenome(start-1, end)
		jarwinkl, _ = textdistance.JaroDistance(refS, altS)
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

func calcGC3Val(file string) []GC3Type {
	// var wg sync.WaitGroup
	var altPositions = make(map[string][]allPositionsInGene)
	// var validData []string
	var gcArray []GC3Type
	// wg.Add(1)

	var allLocusUnsort, allLocuses []string

	snps := parserVCF(file, false, allGenesVal)

	for _, val := range snps {

		altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, nuc: val.Alt})

		allLocusUnsort = append(allLocusUnsort, val.Locus)
	}
	allLocuses = removeStringDuplicates(allLocusUnsort)

	for _, val := range allLocuses {

		// defer wg.Done()
		altS := makeAltString(val, altPositions[val])

		start, end := getGenePosByName(val)

		refS := getNucFromGenome(start-1, end)
		// if len(altS) == 0 {
		// 	altS = refS
		// 	fmt.Println("!!!!!")
		// }

		gcAlt, _, _, gc3Alt := codon.GcCodonCalc(altS)
		gcRef, _, _, gc3Ref := codon.GcCodonCalc(refS)

		gcArray = append(gcArray, GC3Type{Locus: val, GC3Alt: gc3Alt, GC3Ref: gc3Ref, GCalt: gcAlt, GCref: gcRef})

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

func createNCWebServer(port string) {
	/*

	 */
	seq := makeSeq(ncFlag, *gbVerbose, *annMakeSeqRef)
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
	browser.OpenURL(fmt.Sprintf("localhost:%v", port))
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
	mapOfVCF := make(map[string][]gene.SNPinfo)
	mapofSNP := make(map[string][]string)
	// var snpArr []checkSNP
	files := getListofVCF()
	locSNPcheck := getSNPNotation(f)
	var chkSNP []checkSNP

	// var buffer strings.Builder
	// bar := pb.ProgressBarTemplate(pbtmpl).Start(len(files))
	// bar.SetWidth(90)
	for i, file := range files {
		if verbose == true {
			fmt.Printf("Reading %v (%v:%v)%v\r", file, i+1, len(files), strings.Repeat(" ", 60))
		}
		snps := parserVCF(file, false, allGenesVal)
		mapOfVCF[file] = snps

		// bar.Increment()
		// i++
		// fmt.Printf("processed %v files\r", i)
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

func testGeneInfo(genes []gene.Gene) {
	var seq string
	var start, end int
	var gc, gc1, gc2, gc3 float64

	for _, g := range genes {

		start = g.Start
		end = g.End
		seq = getNucFromGenome(start-1, end)

		if start != 0 && end != 0 {
			gc, gc1, gc2, gc3 = codon.GcCodonCalc(seq)
			fmt.Printf("%v\t%.2f\t%.2f\t%.2f\t%.2f\n", g.Locus, gc, gc1, gc2, gc3)

		}
	}
}

func makeAltString(locus string, positions []allPositionsInGene) string {
	// var lStart, lEnd int
	var seqSplit []string
	var seq string

	lStart, lEnd := getGenePosByName(locus)

	seq = getNucFromGenome(lStart-1, lEnd)

	for _, nuc := range seq {
		seqSplit = append(seqSplit, string(nuc))
	}

	for _, val := range positions {

		seqSplit[val.pos-1] = val.nuc

	}

	return strings.Join(seqSplit, "")

}

func getGenePosByName(locus string) (int, int) {
	var start, end int

	for _, g := range allGenesVal {
		if strings.ToUpper(locus) == strings.ToUpper(g.Locus) {
			start = g.Start
			end = g.End
			break
		}
	}

	return start - 1, end
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

func toCircos(genes []gene.Gene) {
	var seq string
	var start, end int
	var gc float64
	var buffer strings.Builder
	fmt.Println("---- ideogram.txt ----")

	fmt.Printf("chr -  %v caption %v %v %v\n", gInfo.Strain, gInfo.Start, gInfo.End, gInfo.Strain)

	for _, g := range genes {
		start = g.Start
		end = g.End
		seq = getNucFromGenome(start-1, end)
		// fmt.Println(seq)
		if start != 0 && end != 0 {
			gc, _, _, _ = codon.GcCodonCalc(seq)

			fmt.Printf("band %v %v %v %v %v\n", gInfo.Strain, g.Locus, g.Locus, start, end)

			buffer.WriteString(fmt.Sprintf("%v %v %v %.2f\n", gInfo.Strain, start, end, gc))
		}
	}

	fmt.Println("---- histogram.txt ----")
	fmt.Println(buffer.String())

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
	var locName string
	var locLen int

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

func makeMatrix(typeof string, fileOut string) {

	var AllPosUnsort, AllPos []int
	var allLocusUnsort, allLocuses []string
	var buffer strings.Builder
	var headers strings.Builder
	var posCount = make(map[int]int)

	// var ResSeq []seqInfo

	files := getListofVCF()
	// fmt.Println(files)
	for i, file := range files {

		// if *flgVerbose == true {
		fmt.Printf("Preparing files: Working on %v files from %v \r", i+1, len(files))
		// }
		snps := parserVCF(file, false, allGenesVal)
		for _, val := range snps {
			AllPosUnsort = append(AllPosUnsort, val.APos)
			posCount[val.APos] = posCount[val.APos] + 1
			allLocusUnsort = append(allLocusUnsort, val.Locus)

		}

	}

	AllPos = unique(AllPosUnsort)
	allLocuses = removeStringDuplicates(allLocusUnsort)

	sort.Ints(AllPos)

	switch typeof {
	case "binary":
		var posFreq = map[int][]string{}
		pos := make(map[int]string)

		headers.WriteString("Pos\t")

		for i, file := range files {
			fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

			snps := parserVCF(file, false, allGenesVal)

			for _, val := range snps {

				pos[val.APos] = val.Alt

			}

			for _, allpos := range AllPos {
				if posCount[allpos] < len(files) {
					if pos[allpos] != "" {
						posFreq[allpos] = append(posFreq[allpos], "1")

					} else {
						posFreq[allpos] = append(posFreq[allpos], "0")

					}
				}

			}
		}
		for _, allpos := range AllPos {
			if posCount[allpos] < len(files) {

				buffer.WriteString(fmt.Sprintln(allpos, "\t", strings.Join(posFreq[allpos], "\t")))

			}
		}
		headers.WriteString("\n")

	case "table":
		var posFreq = map[int][]string{}
		pos := make(map[int]string)

		// headers.WriteString("Pos\t")

		for i, file := range files {
			fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

			snps := parserVCF(file, false, allGenesVal)

			for _, val := range snps {

				pos[val.APos] = val.Alt

			}

			for _, allpos := range AllPos {
				if posCount[allpos] < len(files) {
					if pos[allpos] != "" {
						posFreq[allpos] = append(posFreq[allpos], strconv.Itoa(allpos))

					} else {
						posFreq[allpos] = append(posFreq[allpos], "0")

					}
				}

			}
		}
		for _, allpos := range AllPos {
			if posCount[allpos] < len(files) {

				buffer.WriteString(fmt.Sprintln(strings.Join(posFreq[allpos], "\t")))

			}
		}
		headers.WriteString("\n")

	case "nc":
		var posFreq = map[int][]string{}
		pos := make(map[int]string)

		headers.WriteString("Pos\tRef\t")

		for i, file := range files {
			fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

			snps := parserVCF(file, false, allGenesVal)

			for _, val := range snps {

				pos[val.APos] = val.Alt

			}

			for _, allpos := range AllPos {

				if pos[allpos] != "" {
					posFreq[allpos] = append(posFreq[allpos], fmt.Sprintf("%v", pos[allpos]))

				} else {
					posFreq[allpos] = append(posFreq[allpos], getNucFromGenomePos(allpos))

				}

			}
		}
		for _, allpos := range AllPos {
			// if buffer.Len() == 0 {
			buffer.WriteString(fmt.Sprintln(allpos, "\t", getNucFromGenomePos(allpos), "\t", strings.Join(posFreq[allpos], "\t")))

		}
		headers.WriteString("\n")

	case "locus":

		var locFreq = map[string][]string{}

		var locusCount = make(map[string]int)

		headers.WriteString("Locus\t")

		for i, file := range files {
			fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

			snps := parserVCF(file, false, allGenesVal)

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

	case "freq":
		headers.WriteString("Pos\tFreq\n")

		for _, allpos := range AllPos {

			// headers.WriteString(fmt.Sprintf("P%v\n", allpos))
			buffer.WriteString(fmt.Sprintf("P%v\t%v\n", allpos, posCount[allpos]))

		}

	case "dnds":
		// var dnds [][]DnDsRes
		var locDNDS = map[string][]string{}

		prompt := bufio.NewReader(os.Stdin)
		fmt.Printf("It will take some time(~%v min). Continue?: ", len(files))
		yesNo, _ := prompt.ReadString('\n')

		if yesNo == "y\n" || yesNo == "Y\n" {
			headers.WriteString("Locus\t")

			for i, file := range files {
				headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

				fmt.Printf("Calculating Dn/DS: Working on %v from %v (%v) \r", i+1, len(files), file)

				dnds := calcDnDsVal(file, *gbVerbose)

				for _, dndsVal := range dnds {

					locDNDS[dndsVal.Locus] = append(locDNDS[dndsVal.Locus], fmt.Sprintf("%.2f", dndsVal.DNDS))

				}

			}

			for _, allloc := range allLocuses {
				if len(locDNDS[allloc]) != 0 {

					buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))

				}
			}

			headers.WriteString("\n")

		}
	case "jw":
		// var dnds [][]DnDsRes
		var locJW = map[string][]string{}

		headers.WriteString("Locus\t")

		for i, file := range files {
			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

			fmt.Printf("Calculating Jaro Winkler distance: Working on %v from %v (%v) \r", i+1, len(files), file)

			jw := calcJaroWinklerDist(file, *gbVerbose)

			for _, jwVal := range jw {

				locJW[jwVal.Locus] = append(locJW[jwVal.Locus], fmt.Sprintf("%.3f", jwVal.JWDist))

			}

		}

		for _, allloc := range allLocuses {
			if len(locJW[allloc]) != 0 {

				buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v\t", allloc), strings.Join(locJW[allloc], "\t")))

			}
		}

		headers.WriteString("\n")

	case "gc3":
		// calcGC3Val
		var locGC3 = map[string][]string{}
		var refGC3 = map[string]string{}

		headers.WriteString("Locus\tRefCG3\t")

		for i, file := range files {

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

			fmt.Printf("Calculating GC3 values: Working on %v from %v \r", i+1, len(files))

			gc3 := calcGC3Val(file)

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

			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v\t%v\t", allloc, refGC3[allloc]), strings.Join(locGC3[allloc], "\t")))

			// }
		}

		headers.WriteString("\n")

	}

	// case "gc":
	// 	// calcGC3Val
	// 	var locGC = map[string][]string{}
	// 	var refGC = map[string]string{}

	// 	headers.WriteString("Locus\tRefCG\t")

	// 	for i, file := range files {

	// 		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

	// 		fmt.Printf("Calculating GC values: Working on %v from %v \r", i+1, len(files))

	// 		gc := calcGC3Val(file)

	// 		for _, gcval := range gc {

	// 			locGC[gcval.Locus] = append(locGC[gcval.Locus], fmt.Sprintf("%.2f", gcval.GCalt))
	// 			refGC[gcval.Locus] = fmt.Sprintf("%.2f", gcval.GCref)

	// 		}

	// 	}

	// 	for _, allloc := range allLocuses {

	// 		if len(locGC[allloc]) != 0 {
	// 			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v\t%v\t", allloc, refGC[allloc]), strings.Join(locGC[allloc], "\t")))

	// 		}
	// 	}

	// 	headers.WriteString("\n")

	// }

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

var resStart, resEnd int

func getSequenceRange(PosRange string, noseq bool) rangePosInfo {

	// var rangeParser = regexp.MustCompile(`^(\d+)\W+(\d+)|\D+\W+(\d+)\W+(\d+)`)
	var rangeParser = regexp.MustCompile(`\b(\d+)\b\W+\b(\d+)\b`)

	var result rangePosInfo

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

// func printSequenceRange(rangeList []rangePosInfo, noseq bool, typeOf int, th string) {
// 	var posCount = map[int]int{}
// 	var posShown = map[int]int{}
// 	// var printOptions int
// 	locTH, _ := strconv.Atoi(th)
// 	for _, value := range rangeList {
// 		posCount[value.Start+value.End] = posCount[value.Start+value.End] + 1
// 	}
// 	for _, value := range rangeList {
// 		if noseq == true {
// 			value.Seq = ""
// 		}
// 		switch {
// 		// case typeOf != "uniq" && *flgTestOptions != "makeseq" && *flgTestOptions != "mkseq":
// 		case typeOf == 1:
// 			if locTH != 0 {

// 				if value.End-value.Start <= locTH {
// 					fmt.Printf("%v\t%v\t%v\t%v\t%v\t%v\t%v\n", value.Gname, value.Start-1, value.End, value.Len, value.Seq, value.Prod, posCount[value.Start+value.End])
// 					// locPosInfo := PosInfo{start: resStart, end: resEnd}
// 					// positions = append(positions, locPosInfo)
// 				}

// 			} else {
// 				fmt.Printf("%v\t%v\t%v\t%v\t%v\t%v\t%v\n", value.Gname, value.Start-1, value.End, value.Len, value.Seq, value.Prod, posCount[value.Start+value.End])
// 			}
// 		case typeOf == 2:
// 			if locTH != 0 {

// 				if value.End-value.Start <= locTH && posShown[value.Start+value.End] == 0 {
// 					fmt.Printf("%v\t%v\t%v\t%v\t%v\t%v\t%v\n", value.Gname, value.Start-1, value.End, value.Len, value.Seq, value.Prod, posCount[value.Start+value.End])
// 					// locPosInfo := PosInfo{start: resStart, end: resEnd}
// 					// positions = append(positions, locPosInfo)
// 				}

// 			} else {
// 				if posShown[value.Start+value.End] == 0 {
// 					fmt.Printf("%v\t%v\t%v\t%v\t%v\t%v\t%v\n", value.Gname, value.Start-1, value.End, value.Len, value.Seq, value.Prod, posCount[value.Start+value.End])
// 				}
// 			}
// 		// case *flgTestOptions == "makeseq" || *flgTestOptions == "mkeseq":
// 		case typeOf == 3:

// 			if locTH != 0 {

// 				if value.End-value.Start <= locTH && value.Gname != "-" && posShown[value.Start+value.End] == 0 {
// 					fmt.Printf(">%v:%v:%v:%v:%v\n%v\n", value.Gname, value.Start, value.End, value.Len, posCount[value.Start+value.End], value.Seq)
// 				}

// 			} else if value.Gname != "-" && posShown[value.Start+value.End] == 0 {
// 				fmt.Printf(">%v:%v:%v:%v:%v\n%v\n", value.Gname, value.Start, value.End, value.Len, posCount[value.Start+value.End], value.Seq)
// 			}
// 		}
// 		if posCount[value.Start+value.End] > 0 {
// 			posShown[value.Start+value.End] = 1
// 		}

// 	}
// }

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

func getHashSNP(snp gene.SNPinfo) uint64 {
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
	for _, g := range allGenesVal {
		fmt.Printf("%v\t%v\t%v\t%v\n", gInfo.Version, g.Start, g.End, g.Locus)
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

	var gmap []genomeMapInfo
	var gmapval genomeMapInfo

	for _, g := range allGenesVal {
		gmapval = genomeMapInfo{Start: g.Start, End: g.End, Locus: g.Locus, TypeOf: g.TypeOf, Product: g.Product}
		gmap = append(gmap, gmapval)
	}

	return gmap
}

func getRangeFromFile(file string, verbose bool, noseq bool) []rangePosInfo {
	var rangeParser = regexp.MustCompile(`\b(\d+)\b\W+\b(\d+)\b`)

	// var rangeArr []rangeArray
	var posRange, seq string
	var gc float64
	var unsorted []rangeArray
	var result []rangePosInfo
	// var uniqRange = map[int]rangeArray{}
	// var lastHash, i, j, k int
	var j, k int
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
	var res []rulesInfo
	var rinfo rulesInfo
	rRule := regexp.MustCompile(`^(.*)\b=\b(\w+)$`)
	// ruleMap := map[string][]string{}
	f, err := os.Open(file) // открываем файл

	if err != nil {
		fmt.Println(err)

	}
	defer f.Close()
	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
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

	// for _, val := range res {
	// 	for _, v := range val.Variants {
	// 		if v == "RIFAMPICIN" {
	// 			fmt.Println("YES", val.Name)
	// 		}
	// 	}

	// }
	// if len(res) == 0 {
	// 	fmt.Println("No rules was found! Please check and correct it.")
	// 	os.Exit(3)
	// }
	// fmt.Println(res)
	return res
	// fmt.Println(ruleMap)

}

// func getHash(str string) uint64 {
// 	h := xxhash.New64()
// 	r := strings.NewReader(str)
// 	io.Copy(h, r)
// 	res := h.Sum64()
// 	return res
// }
