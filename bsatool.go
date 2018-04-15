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
	"net/http"
	"os"
	// "os/exec"
	// "reflect"
	// "math/rand"
	"regexp"
	// "encoding/json"
	"flag"
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
	"./tang"
	// "text/template"
	// "net/http"
	// "github.com/davecgh/go-spew/spew"
	// tm "github.com/buger/goterm"
	// "gopkg.in/cheggaaa/pb.v2"
	"html/template"

	"./codon"
	"github.com/pkg/browser"
	// "io/ioutil"
	// "encoding/csv"
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
                                                                          
 BSATool - Bacterial Snp Annotation Tool ver.0.5alpha
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
	vcfExt = "*.vcf"
	// pbtmpl      = `{{counters . }}{{ bar . "⎜" "agtc"   "⬤" "TCAG" "⎜"}}{{percent . }}{{rtime . " ETA %s  "}}{{speed .}} `
	lmnRegExp   = `^(\w+)\W+(\d+)(\D)>(\D)\s+(.*)` // LMN Regular expression
	pmlnRegExp  = `^(\d+)_(\D)>(\D)\W+(\w+)\W+(.*)`
	pmnRegExp   = `^(\d+)_(\D)>(\D)\W+N:(.*)`
	lsaanRegExp = `^(\w+)\W+(\D{1})(\d+)(\D{1})\W+(.*)`
	llaanRegExp = `^(\w+)\W+(\D{3})(\d+)(\D{3})\W+(.*)`
	lcnRegExp   = `^(\w+)\W+codon(\d+)\W+(.*)`
	url         = "localhost:8080"
)

// type SNP interface {
//
// }

// var linesFromGB []string             // массив строк считанных из файла genbank
// var genomeAnchors = map[int]string{} // карта резльтатов парсинга генбанк-файла
// var totalCDSFound int // количество найденных CDS
// var gene.AllGenesVal []gene.Gene          // массив информации об генах, загруженный из файла базы данных

var flgMkDB = flag.Bool("mkdb", false, "Создание файла локальной базы данных") // флаг создания базы данных
var flgGBfile = flag.String("gb", "", "Genbank файл")                          // генбанк файл
var flgabout = flag.Bool("v", false, "О программе")                            // о программе

var flgOut = flag.String("out", "", "Имя локальной базы данных") // выходной файл базы данных при парсинге файла генбанка
var flgDB = flag.String("db", "", "Файл локальной базы данных")  // файл основной базы данных
var flgVCF = flag.String("vcf", "", "vcf файл")                  // файл vcf, с параметром=list все файлы в директории
var flgPos = flag.Int("pos", 0, "SNP-позиция")
var flgDebug = flag.Bool("debug", false, "Режим отладки") // режим отладки
var flgGenomeMap = flag.Bool("gmap", false, "Карта генома")
var flgmakeSeq = flag.String("mkseq", "", "NC или AA") // сгенерировать синтетическую последовательность
var flgTang = flag.Bool("tang", false, "Индек Танга")  // индекс танга
// var flgDevMode = flag.Bool("dev", false, "")

// var flgReportType = flag.String("rep", "", "")
var flgSNP = flag.String("snp", "", "Файл с известными SNP, которые необходимо проверить")
var flgWeb = flag.Bool("web", false, "Режим вывода результатов в Веб браузер")               // вывод результатов в браузер
var flgShare = flag.Bool("share", false, "Поиск общих, для всех vcf файлов в папке, СНИПОВ") // поиск общих снипов для всех Vcf файлов в папке
var flgIndel = flag.Bool("indel", false, "Поиск ИнДелов")
var flgVerbose = flag.Bool("verb", false, "Вывод дополнительной информации")
var flgWithRef = flag.Bool("ref", false, "Создавать референсную последовательность при использовании команды -mkseq")
var flgTestMode = flag.String("test", "", "test mode")
var flgInFile = flag.String("in", "", "input")
var flgType = flag.String("type", "", "type of matrix")
var flgDev = flag.Bool("dev", false, "dev mode")
var flgPort = flag.Int("port", 8080, "set localhost:port")
var flgAction = flag.String("action", "", "select action of work")

var gInfo genomeInfo

// var allPos []allPositionsInGene

type seqInfo struct {
	Name, Seq string //
	// Len       int
}
type genomeInfo struct {
	TypeOf, Organism, Strain string
	Start, End               int
}

type statInfo struct {
	Pos, Count, Perc        int
	FilesWith, FilesWithout string
}
type checkSNP struct {
	FileName, FoundSNP string
}

type allPositionsInGene struct {
	pos int
	nuc string
}

type DnDsRes struct {
	N, S, PN, PS, DN, DS, DNDS float64
	ND, NS                     int
	Locus, Product             string
}

// type allPositionsInGeneArray struct {
// 	items []allPositionsInGene
// }

// type locusInfo struct {
// 	start, end int
// 	alt        []allPositionsInGene
// }

var genomeSeqSlice []string // информация об генах, загруженная из базы

func main() {

	// парсинг флагов
	defer os.Exit(0)
	flag.Parse()
	// о программе
	if *flgabout == true {
		about()
	}

	//  файл базы данных
	if *flgMkDB == true {
		if *flgGBfile != "" && *flgOut != "" {
			if _, err := os.Stat(*flgGBfile); os.IsNotExist(err) {
				fmt.Printf("The %v file is not exist!\n", *flgGBfile)
				os.Exit(3)
			}
			//
			// if *flgTestMode == true {

			// 	gene.AllGenesVal, genomeSeqSlice = readGBFileTest(*flgGBfile)
			// 	// readGBFileTest(*flgGBfile)
			// 	writeDB(*flgOut, gene.AllGenesVal)
			// } else {
			gene.AllGenesVal, genomeSeqSlice = readGBFile(*flgGBfile)
			writeDB(*flgOut, gene.AllGenesVal)
			// }
			// fmt.Printf("%v\n", genomeSeqSlice)
		} else {
			fmt.Println("Не указаны необходимые аргументы (-mkdb -gb=GENBANK_FILE -out=OUTPUT_DB_FILE)")
		}
	}

	if *flgDB != "" {
		// g := Gene
		if _, err := os.Stat(*flgDB); os.IsNotExist(err) {
			fmt.Printf("The %v file is not exist!\n", *flgDB)
			os.Exit(3)
		}

		gene.AllGenesVal = readDB(*flgDB)

		if *flgTestMode == "stat" {
			testSNPstat()
		} else if *flgTestMode == "dnds" {
			calcDnDsVal("ERR067585.vcf", true)
		} else if *flgTestMode == "ginfo" {
			testGeneInfo(gene.AllGenesVal)
			// fmt.Println(gInfo)
		} else if *flgTestMode == "circos" {
			toCircos(gene.AllGenesVal)

		} else if *flgTestMode == "annotate" && *flgInFile != "" {
			annotatePositions()
		} else if *flgTestMode == "matrix" && *flgOut != "" && *flgType != "" {
			makeMatrix()
		}

		switch {
		case *flgVCF != list && *flgVCF != "" && *flgWeb == false && *flgmakeSeq == "":
			if _, err := os.Stat(*flgVCF); os.IsNotExist(err) {
				fmt.Printf("The %v file is not exist!\n", *flgVCF)
				os.Exit(3)
			}
			parserVCF(*flgVCF, true, gene.AllGenesVal)
		case *flgVCF != list && *flgVCF != "" && *flgWeb == true && *flgmakeSeq == "":
			if _, err := os.Stat(*flgVCF); os.IsNotExist(err) {
				fmt.Printf("The %v file is not exist!\n", *flgVCF)
				os.Exit(3)
			}
			snps := parserVCF(*flgVCF, false, gene.AllGenesVal)
			printWebResults(snps)
		case *flgVCF == list && *flgWeb == false && *flgmakeSeq == "" && *flgShare == false && *flgSNP == "":
			listOfVCFFiles()
		case *flgWeb == false && *flgmakeSeq == ncFlag && *flgShare == false:
			seq := makeSeq(ncFlag)
			for _, val := range seq {
				fmt.Println(val.Seq)
			}
		case *flgWeb == true && *flgmakeSeq == ncFlag && *flgShare == false:
			createNCWebServer()

		case *flgWeb == false && *flgmakeSeq == aaFlag && *flgShare == false:
			makeSeq(aaFlag)
		case *flgShare == true:
			// getShareSNP()
			getShareSNP()
		case *flgSNP != "":
			if _, err := os.Stat(*flgSNP); os.IsNotExist(err) {
				fmt.Printf("The %v file is not exist!\n", *flgSNP)
				os.Exit(3)
			}

			checkSNPfromFile(*flgSNP)

		case *flgDebug == true:
			fmt.Printf("%v\n", gene.AllGenesVal)
		}

		if *flgGenomeMap == true {
			var igensS []string
			var igensE []string
			for i, g := range gene.AllGenesVal {

				igensS = append(igensS, g.Start)
				igensE = append(igensE, g.End)
				// igenStart, _ := strconv.Atoi(g.End)
				if i >= 1 && i < len(gene.AllGenesVal)-1 {
					// 	// igenEnd, _ := strconv.Atoi(g.Start)
					fmt.Printf("igen:s%v e:%v\n", igensS[i], igensE[i-1])
				}
				fmt.Printf("l:%v s:%v e:%v d:%v %v\n", g.Locus, g.Start, g.End, g.Direction, i)

			}
		}

		// for _, g := range gene.AllGenesVal {
		// 	if *flgPos != 0 {

		// 		lStart, _ := strconv.Atoi(g.Start)
		// 		lEnd, _ := strconv.Atoi(g.End)

		// 		if *flgPos >= lStart && *flgPos <= lEnd {

		// 			snp := getSNPInfo(*flgPos, g, "A")
		// 			fmt.Printf("apos:%v, gpos:%v, pic:%v, AA:%v, ref:%v, alt:%v, dir:%v, aa:%v\nL:%v,▶%v…%v◀, P:%v\n",
		// 				snp.APos, snp.PosInGene, snp.PosInCodonG, snp.CodonNbrInG, snp.NucInPos, snp.Alt, g.Direction,
		// 				snp.RefAA, snp.Locus, snp.Start, snp.End, snp.Product)

		// 		}

		// 	}
		// }

	}

}

func getSNPInfo(apos int, g gene.Gene, alt string) (gene.SNPinfo, int) {
	var snp gene.SNPinfo           // структура SNP
	var codonPositions []string    // срез для разбивки кодона побуквенно
	var altCodonPositions []string // срез для разбивки кодона побуквенно альтернативным нуклеотидом
	var locReportType string
	var geneLen int
	var titv string
	var trouble int
	lStart, _ := strconv.Atoi(g.Start) // переменная начала гена
	lEnd, _ := strconv.Atoi(g.End)
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
		tangIdx = "-"
	} else if aaRefShort != aaAltShort && aaAltShort != "X" {
		mut = "missense"
		tangIdx = tang.GetTangInx(aaRefShort, aaAltShort)
	} else if aaRefShort != aaAltShort && aaAltShort == "X" {
		mut = "nonsense"
	}

	if strings.ToUpper(alt) == strings.ToUpper(nucG) {
		trouble = 1
	}

	if *flgTang == true {
		locReportType = "T1"
	} else if *flgTang == false {
		locReportType = "T0"
	}
	tang.GetTangInx(aaRefShort, aaAltShort)
	titv = checkTiTv(nucG, alt)
	// fmt.Println(lStart, lEnd, g.Direction, geneLen, g.Locus, g.Product)

	snp = gene.SNPinfo{APos: apos, PosInGene: posInGene, PosInCodonG: posInCodonG,
		RefCodon: codon, RefAA: aaRef, NucInPos: strings.ToUpper(nucG), Locus: g.Locus,
		Direction: g.Direction, Name: g.Name, Product: g.Product,
		Start: g.Start, End: g.End, CodonNbrInG: codonNbrInG, AltCodon: altCodon,
		AltAA: aaAlt, RefAAShort: aaRefShort, AltAAShort: aaAltShort,
		Mutation: mut, Tang: tangIdx, Alt: alt,
		Note: g.Note, ReportType: locReportType, ProteinID: g.ProteinID,
		GeneID: g.GeneID, GOA: g.GOA, GeneLen: geneLen, TiTv: titv}

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

func readGBFile(file string) (g []gene.Gene, genomeSplice []string) {
	//функция для считывания файла в формате генбанка и занесения строк в массив linesFromGB

	var qenomeSeq = regexp.MustCompile(`\s+(\d+)\s+(\w{6}.*)`)
	var regDelSpaces = regexp.MustCompile(`(\s+)(.*)`)
	var checkOrigin = regexp.MustCompile(`(\d+\s\w{10}\s\w{10}\s\w{10})`)
	var checkCDS = regexp.MustCompile(`(CDS\s+complement\W\d+\W+\d+\W)|(CDS\s+\d+\W+\d+)`)
	var endOfCDS = regexp.MustCompile(`(gene\s+\d+\W+\d+)|(gene\s+complement\W\d+\W+\d+)`)
	var startOrigin = regexp.MustCompile(`(ORIGIN\s+)|(translation=)`)
	var repeatRegion = regexp.MustCompile(`(repeat_region\s+complement\W\d+\W+\d+\W)|(repeat_region\s+\d+\W+\d+)`)
	var makeAnchors = regexp.MustCompile(`^(/)`)
	var genomeSource = regexp.MustCompile(`organism=\W(.*)\W`)
	var genomeStrain = regexp.MustCompile(`strain=\W(.*)\W`)
	var genomeStartEnd = regexp.MustCompile(`source\W+(\d+)..(\d+)`)

	var gStart, gEnd int

	//  source          1..4411532
	//                  /organism="Mycobacterium tuberculosis H37Rv"
	//                  /mol_type="genomic DNA"
	//                  /strain="H37Rv"
	//                  /db_xref="taxon:83332"
	// var replaceCDS = regexp.MustCompile(`(CDS\s+.*)`)

	var resString []string

	var querySplice []string
	var splitedGenome []string
	var originBlock, CDSblock, firstCDS int

	var changedStr, organismName, organismStrain string
	var noteBuffer strings.Builder

	// if _, err := os.Stat(file); os.IsNotExist(err) {
	// 	fmt.Printf("The %v file is not exist!\n", file)
	// 	os.Exit(3)
	// }

	f, err := os.Open(file)

	if err != nil {
		fmt.Println(err)

	}
	defer f.Close()
	querySplice = append(querySplice, "START_BLOCK")
	scanner := bufio.NewScanner(f)
	go process("Working...")
	for scanner.Scan() {

		scanTxt := scanner.Text() // Println will add back the final '\n'
		for _, gStartEnd := range genomeStartEnd.FindAllStringSubmatch(scanTxt, -1) {
			gStart, _ = strconv.Atoi(gStartEnd[1])
			gEnd, _ = strconv.Atoi(gStartEnd[2])
		}
		for _, gName := range genomeSource.FindAllStringSubmatch(scanTxt, -1) {
			organismName = gName[1]
			fmt.Printf("Organism:%v\n", organismName)
		}
		for _, gStrain := range genomeStrain.FindAllStringSubmatch(scanTxt, -1) {
			organismStrain = gStrain[1]
			fmt.Printf("Strain:%v\n", organismStrain)
		}

		for _, cdsFound := range checkCDS.FindAllStringSubmatch(scanTxt, -1) {
			if len(cdsFound[1]) != 0 || len(cdsFound[2]) != 0 {
				CDSblock = 1

				// fmt.Println(scanTxt)
			}
		}
		for _, cdsEnd := range endOfCDS.FindAllStringSubmatch(scanTxt, -1) {
			if len(cdsEnd[1]) != 0 || len(cdsEnd[2]) != 0 {
				CDSblock = 0
				// fmt.Println(scanTxt)
			}
		}
		for _, sOrigin := range startOrigin.FindAllStringSubmatch(scanTxt, -1) {
			if len(sOrigin[1]) != 0 || len(sOrigin[2]) != 0 {
				CDSblock = 0
			}
		}
		for _, rRegion := range repeatRegion.FindAllStringSubmatch(scanTxt, -1) {
			if len(rRegion[1]) != 0 || len(rRegion[2]) != 0 {
				CDSblock = 0
			}
		}
		for _, origninFound := range checkOrigin.FindAllStringSubmatch(scanTxt, -1) {
			if len(origninFound[1]) != 0 {
				originBlock = 1

			} else {
				originBlock = 0
			}
		}

		changedStr = strings.TrimSuffix(changedStr, "\n")
		// changedStr = replaceCodonStart.ReplaceAllString(string(scanTxt), "") // удаляем пробелы
		changedStr = regDelSpaces.ReplaceAllString(string(scanTxt), "$2 ") // удаляем пробелы

		switch originBlock {
		case 0:
			if CDSblock == 1 {
				if firstCDS == 0 {

					firstCDS = 1

				}
				changedStr = strings.TrimSuffix(changedStr, "\n")
				changedStr = strings.Replace(string(changedStr), "/note=", "!note=", -1)       // меняем / на ! для дальнейшего парсинга
				changedStr = strings.Replace(string(changedStr), "/product=", "!product=", -1) // см выше.
				changedStr = strings.Replace(string(changedStr), "\"", "", -1)

				changedStr = makeAnchors.ReplaceAllString(string(changedStr), "!!")
				// fmt.Println(changedStr)

			}

		case 1:

			for _, genomeMatch := range qenomeSeq.FindAllStringSubmatch(scanTxt, -1) {

				resString = append(resString, strings.Replace(genomeMatch[2], " ", "", -1))
			}
		}

		if CDSblock == 1 {
			// fmt.Println(changedStr)
			if strings.Index(changedStr, "!!") != 0 && strings.Index(changedStr, "CDS  ") != 0 {
				noteBuffer.WriteString(changedStr)
				// fmt.Println(changedStr)

			} else {
				if noteBuffer.Len() != 0 {

					// fmt.Printf("%v\n", tmpStr)
					querySplice = append(querySplice, strings.TrimPrefix(noteBuffer.String(), "!"))
					// fmt.Println()
					noteBuffer.Reset()
				}

			}

			if strings.Index(changedStr, "!!") == 0 {
				querySplice = append(querySplice, strings.TrimPrefix(changedStr, "!!"))
				// fmt.Println(changedStr)

			} else if strings.Index(changedStr, "CDS  ") == 0 {

				// if firstCDS
				if firstCDS == 1 {
					querySplice = append(querySplice, "START_OF")
					firstCDS = 2
				} else if firstCDS == 2 {
					querySplice = append(querySplice, "END_OF")
					querySplice = append(querySplice, "START_OF")
				}

				querySplice = append(querySplice, strings.TrimPrefix(changedStr, "!!"))
			}

		}

	}
	querySplice = append(querySplice, "END_OF")
	querySplice = append(querySplice, "END_BLOCK")
	splitedGenome = strings.SplitAfter(strings.Join(resString, ""), "")
	// for i, val := range querySplice {
	// 	fmt.Printf("%v %v", i, val)
	// }
	var lStart, lEnd, lDir, lLoc, lName, lProd, lNote, gID, pID, lGOA string
	var nucCore string
	var cdsStEnd = regexp.MustCompile(`^CDS\s+(\d+)\W+(\d+)|^CDS\s+complement\W(\d+)\W+(\d+)`)
	var cdsLocus = regexp.MustCompile(`locus_tag=(.*)`)
	var cdsName = regexp.MustCompile(`gene=(.*)`)
	var cdsProd = regexp.MustCompile(`product=(.*)`)
	var cdsNote = regexp.MustCompile(`note=(.*)`)
	var cdsgID = regexp.MustCompile(`db_xref=GeneID:(\d+)`)
	var cdsprotID = regexp.MustCompile(`protein_id=(.*)`)
	var cdsGOA = regexp.MustCompile(`db_xref=GOA:(.*)`)

	var startOfBlock, endOfBlock, cdsStart, cdsEnd int
	// var cdsOpen, cdsClosed int
	cdsCount := 1
	// fmt.Println(querySplice)

	for _, val := range querySplice {

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
			cdsCount++
		} else if strings.Index(val, "END_OF") != -1 {
			cdsStart = 0
			cdsEnd = 1
		}

		if startOfBlock == 1 && endOfBlock == 0 {

			if cdsStart == 1 && cdsEnd == 0 {

				for _, cdsStEndMatch := range cdsStEnd.FindAllStringSubmatch(val, -1) {

					if len(cdsStEndMatch[1]) != 0 {
						lStart, lEnd, lDir = strings.TrimSpace(cdsStEndMatch[1]), strings.TrimSpace(cdsStEndMatch[2]), "f"
						// cdsOpen = 1

					} else if len(cdsStEndMatch[3]) != 0 {

						lStart, lEnd, lDir = strings.TrimSpace(cdsStEndMatch[3]), strings.TrimSpace(cdsStEndMatch[4]), "r"
						// cdsOpen = 1

					}

				}

				for _, cdsLocusMatch := range cdsLocus.FindAllStringSubmatch(val, -1) {
					lLoc = strings.Replace(cdsLocusMatch[1], " ", "", -1)
					// fmt.Println(lLoc)
					// cdsOpen = 1
					// cdsClosed = 0
				}

				// cdsLocusMatch := cdsLocus.FindAllStringSubmatch(val, -1)
				// if len(cdsLocusMatch) != 0 {
				// 	lLoc = cdsLocusMatch[0][1]
				// }

				for _, cdsNameMatch := range cdsName.FindAllStringSubmatch(val, -1) {
					lName = strings.Replace(cdsNameMatch[1], " ", "", -1)
					// cdsOpen = 1
					// cdsClosed = 0
					// fmt.Println(lName)
				}
				for _, cdsProdMatch := range cdsProd.FindAllStringSubmatch(val, -1) {
					lProd = cdsProdMatch[1]
					// cdsOpen = 1
					// cdsClosed = 0
					// fmt.Println(lProd)
				}
				for _, cdsNoteMatch := range cdsNote.FindAllStringSubmatch(val, -1) {
					lNote = cdsNoteMatch[1]

					// fmt.Println(lNote)
				}
				for _, cdsgIDMatch := range cdsgID.FindAllStringSubmatch(val, -1) {
					gID = strings.Replace(cdsgIDMatch[1], " ", "", -1)
					// cdsOpen = 1
					// cdsClosed = 0
					// fmt.Println(gID)

				}
				for _, cdsgprotIDMatch := range cdsprotID.FindAllStringSubmatch(val, -1) {
					pID = strings.Replace(cdsgprotIDMatch[1], " ", "", -1)
					// cdsOpen = 1
					// cdsClosed = 0
					// fmt.Println(pID)

				}
				for _, cdsGOAMatch := range cdsGOA.FindAllStringSubmatch(val, -1) {
					lGOA = strings.Replace(cdsGOAMatch[1], " ", "", -1)

					// fmt.Println(lGOA)
					// cdsOpen = 1
					// cdsClosed = 0
				}

				if *flgDev == true {
					fmt.Println(val)
				}

				if len(gID) == 0 {
					nucCore = "yes"

				} else {
					nucCore = "no"
				}
				// fmt.Println(nucCore)

				if lName == "" {
					lName = lLoc
				}

			} else if cdsStart == 0 && cdsEnd == 1 {
				g := gene.Gene{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir,
					Product: lProd, Name: lName, GeneID: gID, ProteinID: pID, Note: lNote, GOA: lGOA}

				gene.AllGenesVal = append(gene.AllGenesVal, g)
				// fmt.Printf("g:%v b:%v %v\n", cdsOpen, cdsClosed, val)
				if *flgDebug == true {

					fmt.Printf("l:%v s:%v e:%v d:%v p:%v gId:%v pId:%v n:%v GOA:%v\n", lLoc, lStart, lEnd, lDir, lProd, gID, pID, lNote, lGOA)
					// }
				} else {
					fmt.Printf("\rFound %v CDS, nucCore:%v", cdsCount, nucCore)
				}
				lName, lStart, lEnd, lLoc, lDir, lProd, gID, pID, lGOA, lNote = "", "", "", "", "", "", "", "", "", ""

				// fmt.Printf("s:%v e:%v d:%v o:%v i:%v val:%v\n", lStart, lEnd, lDir, cdsOpen, i, val)

			}
		}

		// if cdsOpen != 0 {
		// 	g := gene.Gene{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir,
		// 		Product: lProd, Name: lName, GeneID: gID, ProteinID: pID, Note: lNote, GOA: lGOA}

		// 	gene.AllGenesVal = append(gene.AllGenesVal, g)
		// 	// fmt.Printf("g:%v b:%v %v\n", cdsOpen, cdsClosed, val)
		// 	if *flgDebug == true {

		// 		fmt.Printf("l:%v s:%v e:%v d:%v p:%v gId:%v pId:%v n:%v GOA:%v\n", lLoc, lStart, lEnd, lDir, lProd, gID, pID, lNote, lGOA)
		// 	}
		// 	// } else {
		// 	// 	fmt.Printf("\rFound %v CDS, nucCore:%v", cdsCount-1, nucCore)
		// 	// }
		// 	lName, lStart, lEnd, lLoc, lDir, lProd, gID, pID, lGOA, lNote = "", "", "", "", "", "", "", "", "", ""

		// }

		// fmt.Println(val, cdsOpen, cdsClosed)

		// if cdsOpen == 1 && cdsClosed == 1 && lStart != "" && lEnd != "" {
		// 	// fmt.Println(lStart, lEnd, lDir, lLoc, lProd, cdsOpen, cdsClosed)
		// 	// cdsCount++
		// 	g := gene.Gene{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir,
		// 		Product: lProd, Name: lName, GeneID: gID, ProteinID: pID, Note: lNote, GOA: lGOA}

		// 	gene.AllGenesVal = append(gene.AllGenesVal, g)
		// 	// fmt.Printf("g:%v b:%v %v\n", cdsOpen, cdsClosed, val)
		// 	if *flgDebug == true {

		// 		fmt.Printf("l:%v s:%v e:%v d:%v p:%v gId:%v pId:%v n:%v GOA:%v\n", lLoc, lStart, lEnd, lDir, lProd, gID, pID, lNote, lGOA)
		// 	}
		// 	// } else {
		// 	// 	fmt.Printf("\rFound %v CDS, nucCore:%v", cdsCount-1, nucCore)
		// 	// }
		// 	lName, lStart, lEnd, lLoc, lDir, lProd, gID, pID, lGOA, lNote = "", "", "", "", "", "", "", "", "", ""

		// }

	}

	gInfo = genomeInfo{TypeOf: nucCore, Organism: organismName, Start: gStart, End: gEnd, Strain: organismStrain}
	// go process("Анализ файла закончен!")
	return gene.AllGenesVal, splitedGenome

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

	fmt.Println("\nWell done! ", file, " was created.")

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
	return gene

}

// func ProcessDNA(dnaString string) map[rune]int {
// 	counts := make(map[rune]int)
// 	allowedNucleotides := "ATGC"
// 	for _, letter := range dnaString {
// 		if strings.ContainsRune(allowedNucleotides, letter) {
// 			counts[letter]++
// 		}
// 	}
// 	return counts

// }

// func processReverseCodon(codon string, codonPosition int) string {
// 	var lCodons []string
// 	var result string
// 	lCodons = strings.Split(strings.ToLower(codon), "")
// 	for i, val := range lCodons {
// 		if i == codonPosition {
// 			val = strings.ToUpper(val)
// 		}
// 		result = result + val
// 		// alt_res = alt_res + alt_val
//
// 	}
// 	return result
// }

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
	var indel = regexp.MustCompile(`^^\S+\W+(\d+)\W+(\w+)\s+(\w+).*(INDEL).*DP=(\d+)`)
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
		if *flgIndel == true {
			for _, indelMatch := range indel.FindAllStringSubmatch(scanner.Text(), -1) {
				if vcfValid == true && indelMatch[4] == "INDEL" {
					fmt.Println(indelMatch[1], indelMatch[2], indelMatch[3])

				}

			}
		}

		for _, match := range vcf.FindAllStringSubmatch(scanner.Text(), -1) {

			if vcfValid == true {
				apos, _ := strconv.Atoi(match[1])
				ref := match[2]
				alt := match[3]

				for _, g := range genes {

					lStart, _ := strconv.Atoi(g.Start)
					lEnd, _ := strconv.Atoi(g.End)

					if apos >= lStart && apos <= lEnd {
						snp, err := getSNPInfo(apos, g, alt)
						// br := testing.Benchmark(snp)
						// fmt.Println(br)

						if len(ref) == 1 && len(alt) == 1 && err != 1 {
							snpFromVCF = append(snpFromVCF, snp)
							if print == true {
								printResults(snp)
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

func listOfVCFFiles() {
	// files, err := ioutil.ReadDir(".")
	// if err != nil {
	// 	log.Fatal(err)
	// }

	files := getListofVCF()
	for _, file := range files {
		if strings.Contains(file, ".vcf") {
			if *flgDebug == true {
				fmt.Printf("\n\n%v:\n\n", file)
			}
			parserVCF(file, true, gene.AllGenesVal)

			// fmt.Println(file.Name())
		}
	}
}

func makeSeq(typeof string) []seqInfo {

	var AllPosUnsort, AllPos []int
	var ResSeq []seqInfo

	files := getListofVCF()

	for i, file := range files {

		if *flgVerbose == true {
			fmt.Printf("Reading files: Read %v from %v \r", i+1, len(files))
		}
		snps := parserVCF(file, false, gene.AllGenesVal)
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
	if *flgWithRef == true {
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
		fmt.Printf("Generating sequences: Made %v from %v \r", i+1, len(files))
		pos := make(map[int]string)
		var buffer strings.Builder

		buffer.WriteString(fmt.Sprintf(">%v\n", strings.ToUpper(file)))

		snps := parserVCF(file, false, gene.AllGenesVal)
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

// func removeDuplicates(elements []int) []int {
// 	// Use map to record duplicates as we find them.
// 	encountered := map[int]bool{}
// 	result := []int{}

// 	for v := range elements {
// 		if encountered[elements[v]] == true {
// 			// Do not add duplicate.
// 		} else {
// 			// Record this element as an encountered element.
// 			encountered[elements[v]] = true
// 			// Append to result slice.
// 			result = append(result, elements[v])
// 		}
// 	}
// 	// Return the new slice.
// 	return result
// }

func unique(list []int) []int {
	unique_set := make(map[int]bool, len(list))
	for _, x := range list {
		unique_set[x] = true
	}
	result := make([]int, 0, len(unique_set))
	for x := range unique_set {
		result = append(result, x)
	}
	return result
}

func getInterGen(pos int) {
	var igensS []string
	var igensE []string
	for i, g := range gene.AllGenesVal {

		igensS = append(igensS, g.Start)
		igensE = append(igensE, g.End)
		// igenStart, _ := strconv.Atoi(g.End)
		if i >= 1 && i < len(gene.AllGenesVal)-1 {
			// 	// igenEnd, _ := strconv.Atoi(g.Start)
			fmt.Printf("igen:s%v e:%v\n", igensS[i], igensE[i-1])
		}
		fmt.Printf("l:%v s:%v e:%v d:%v %v\n", g.Locus, g.Start, g.End, g.Direction, i)

	}
}

func printResults(snp gene.SNPinfo) {
	// resFlg
	// if snp.Mutation == "missense" {
	// spew.Dump(snp)
	switch snp.ReportType {

	case "T1":
		fmt.Printf("%v\t%v\t%v%v>%v\t%v/%v\t%v%v%v\t%v\t%v\t%v\n", snp.Locus, snp.APos, snp.PosInGene, snp.NucInPos, snp.Alt, snp.RefCodon, snp.AltCodon,
			snp.RefAAShort, snp.CodonNbrInG, snp.AltAAShort, snp.Mutation,
			snp.Tang, snp.Product)

	case "T0":
		fmt.Printf("%v\t%v\t%v%v>%v\t%v/%v\t%v%v%v\t%v\t%v\n", snp.Locus, snp.APos, snp.PosInGene, snp.NucInPos, snp.Alt, snp.RefCodon, snp.AltCodon,
			snp.RefAAShort, snp.CodonNbrInG, snp.AltAAShort, snp.Mutation,
			snp.Product)

	}

}

func printWebResults(snps []gene.SNPinfo) {

	var htmlTitle = `   <!DOCTYPE html>
				<html>

			<table width="100%" cellspacing="0" cellpadding="4" border="1">
			<tbody>
			<tr>
			
			{{if eq .TypeOf "yes"}}
			
			<td>Locus</td><td>Gene</td><td>Pos.</td><td>Mutation</td><td>Codons</td>
			<td>AA</td><td>Type</td><td>Product</td><td>GOA</td><td>ProteinID</td>
			{{else}}
			<td>Locus</td><td>Gene</td><td>Pos.</td><td>Mutation</td><td>Codons</td>
			<td>AA</td><td>Type</td><td>Product</td><td>GeneID</td><td>ProteinID</td>
			
			{{end}}
			
			</tr>
			
`
	var htmlTemplate = `
			
			{{range $element := .}}
			
			{{if .GeneID}}
			<tr>
			<td><p title="{{.Note}}">{{.Locus}}</p></td><td>{{.Name}}</td><td>{{.APos}}</td><td>{{.PosInGene}}{{.NucInPos}}>{{.Alt}}</td>
			<td>{{.RefCodon}}/{{.AltCodon}}</td><td><p title="{{.RefAA}}{{.CodonNbrInG}}{{.AltAA}}">{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}</p></td>
			<td><p title="U-Index: {{.Tang}}">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}"target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ncbi.nlm.nih.gov/gene/{{.GeneID}}={{.GeneID}}" target="_blank">{{.GeneID}}</a>
			</td><td><a href="http://www.uniprot.org/uniprot/?query={{.ProteinID}}&sort=score">{{.ProteinID}}</td>
			</tr>
			{{else}}
			<tr>
			<td><p title="{{.Note}}">{{.Locus}}</p></td><td>{{.Name}}</td><td>{{.APos}}</td><td>{{.PosInGene}}{{.NucInPos}}>{{.Alt}}</td>
			<td>{{.RefCodon}}/{{.AltCodon}}</td><td><p title="{{.RefAA}}{{.CodonNbrInG}}{{.AltAA}}">{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}</p></td>
			{{if eq .Mutation "missense"}}
			<td bgcolor="#CECEF6"><p title="Tang Index: {{.Tang}}">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}" target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ebi.ac.uk/QuickGO/GProtein?ac={{.GOA}}" target="_blank">{{.GOA}}</a>
			{{else if eq .Mutation "nonsense"}}
			<td bgcolor="#F78181"><p title="Tang Index: {{.Tang}}">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}" target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ebi.ac.uk/QuickGO/GProtein?ac={{.GOA}}" target="_blank">{{.GOA}}</a>
			{{else}}
			<td bgcolor="#ddffcc"><p title="Tang Index: {{.Tang}}">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}" target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ebi.ac.uk/QuickGO/GProtein?ac={{.GOA}}" target="_blank">{{.GOA}}</a>
			{{end}}
			</td><td><a href="http://www.uniprot.org/uniprot/?query={{.ProteinID}}&sort=score" target="_blank">{{.ProteinID}}</td>
			</tr>
			{{end}}
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

	browser.OpenURL(url)

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
	if *flgPort != 0 {
		locPort := fmt.Sprintf(":%v", *flgPort)
		http.ListenAndServe(locPort, nil)
	} else {
		http.ListenAndServe(":8080", nil)
	}
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
	// fmt.Printf("%v\n", snpCheck)
	// for _, val := range parsedSNP {
	// 	fmt.Printf("%v\n", val)
	// }
	return parsedSNP
}

func getShareSNP() {
	var countSNPs = 1
	var pos = make(map[int]int)
	var alt = make(map[int]gene.SNPinfo)
	var share []int
	// var g gene.Gene
	var snpToWeb []gene.SNPinfo
	var snpToConsole []gene.SNPinfo
	files := getListofVCF()
	upperLimit := len(files)
	// bar := pb.StartNew(len(files))
	// bar := pb.ProgressBarTemplate(pbtmpl).Start(len(files))
	// bar.SetWidth(90)
	for i, file := range files {
		// fmt.Printf("processed %v files\r", cnt+1)
		if *flgVerbose == true {
			fmt.Printf("Pass: %v from %v \r", i+1, len(files))
		}
		// bar.Increment()
		snps := parserVCF(file, false, gene.AllGenesVal)

		for _, val := range snps {
			pos[val.APos] = pos[val.APos] + 1
			alt[val.APos] = val
			// 	g = gene.Gene{Locus: val.Locus, Start: val.Start, End: val.End, Name: val.Name, Product: val.Product, Direction: val.Direction, GeneID: val.GeneID, ProteinID: val.ProteinID}
			// }
		}

		for i, lPos := range pos {
			if lPos == upperLimit {
				share = append(share, i)
			}
		}
		sort.Ints(share)
		for _, sharePos := range share {
			countSNPs++
			// fmt.Println(alt[sharePos])

			if *flgWeb == true {
				snpToWeb = append(snpToWeb, alt[sharePos])
				// printResults(alt[sharePos])

			} else {
				// printResults(alt[sharePos])
				snpToConsole = append(snpToConsole, alt[sharePos])
			}
		}
	}
	// fmt.Println()
	if *flgWeb == true && len(snpToWeb) != 0 {
		printWebResults(snpToWeb)
	} else if *flgWeb == false && len(snpToConsole) != 0 {
		for _, res := range snpToConsole {
			printResults(res)
		}
	}
	if *flgDebug == true {
		fmt.Printf("f:%v snp: %v\n%v\n", len(files), countSNPs, files)
	}
}

func testSNPstat() {
	// var countSNPs = 1
	var pos = make(map[int]int)
	var alt = make(map[int]gene.SNPinfo)
	var f = make(map[int][]string)
	var positions, posUnsort []int
	// var g gene.Gene
	// var snpToWeb []gene.SNPinfo
	// var snpToConsole []gene.SNPinfo
	files := getListofVCF()
	upperLimit := len(files)
	// bar := pb.StartNew(len(files))
	// bar := pb.ProgressBarTemplate(pbtmpl).Start(len(files))
	// bar.SetWidth(90)
	for i, file := range files {
		// fmt.Printf("Re %v files\r", cnt+1)
		// if *flgVerbose == true {
		fmt.Printf("Reading: %v (%v from %v)%v \r", file, i+1, len(files), strings.Repeat(" ", 60))
		// }
		// bar.Increment()
		snps := parserVCF(file, false, gene.AllGenesVal)

		for _, val := range snps {
			if pos[val.APos] <= upperLimit {
				pos[val.APos] = pos[val.APos] + 1       //count
				alt[val.APos] = val                     // pos
				f[val.APos] = append(f[val.APos], file) //files
				posUnsort = append(posUnsort, val.APos) //array of positions
			}
			// g = gene.Gene{Locus: val.Locus, Start: val.Start, End: val.End, Name: val.Name, Product: val.Product, Direction: val.Direction, GeneID: val.GeneID, ProteinID: val.ProteinID}
			// }
		}
		positions = unique(posUnsort)
		sort.Ints(positions)

	}
	var stat []statInfo

	for _, p := range positions {
		perc := (pos[p] * 100) / upperLimit
		filesNotInList := compareSlices(files, f[p])
		// fmt.Println(filesNotInList)
		// if perc > 5 && perc < 100 {
		// fmt.Printf("%v\t%v\t%v%v>%v\t%v/%v\t%v%v%v\t%v\t%v\n", alt[p].Locus, alt[p].APos, alt[p].PosInGene, alt[p].NucInPos, alt[p].Alt, alt[p].RefCodon, alt[p].AltCodon,
		// 	alt[p].RefAAShort, alt[p].CodonNbrInG, alt[p].AltAAShort, alt[p].Mutation,
		// 	alt[p].Product)

		// fmt.Printf("----\npos:%v\ncount:%v perc:%v%% \ncontain:%v\n", p, pos[p], perc, strings.Join(f[p], ","))
		stat = append(stat, statInfo{Pos: p, Count: pos[p], Perc: perc, FilesWith: strings.Join(f[p], ",\n"), FilesWithout: strings.Join(filesNotInList, ",\n")})

		// }
	}
	printWebStat(stat)
	// fmt.Println(stat)
}

func calcDnDsVal(file string, print bool) []DnDsRes {

	var altPositions = make(map[string][]allPositionsInGene)
	var validData []string
	var dndsArray []DnDsRes

	snps := parserVCF(file, false, gene.AllGenesVal)

	for _, val := range snps {
		altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, nuc: val.Alt})

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
		refS := getNucFromGenome(start, end)
		dnds := codon.CalcDnDs(refS, altS)

		prod := getProductByName(val)
		// 	N, S, PN, PS, DN, DS, DNDS float64
		// ND, NS                     int

		if dnds.ND != 0 && dnds.NS != 0 {
			if print == true {
				fmt.Printf("l:%v (%v)  dn/ds:%.2f\n", val, prod, dnds.DNDS)
			}
			// fmt.Printf("l:%v p:%v N:%.4f S:%.4f Nd:%v Ns:%v pN:%.4f pS:%.4f dN:%.4f dS:%.4f dN/dS:%.4f \n", val, prod, dnds.N, dnds.S, dnds.ND, dnds.NS, dnds.PN, dnds.PS, dnds.DN, dnds.DS, dnds.DNDS)
			dndsArray = append(dndsArray, DnDsRes{Locus: val, Product: prod, N: dnds.N, S: dnds.S, ND: dnds.ND, NS: dnds.NS, PN: dnds.PN, PS: dnds.PS, DN: dnds.DN, DS: dnds.DS, DNDS: dnds.DNDS})
		}

	}

	return dndsArray

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

func printWebStat(stat []statInfo) {
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

	browser.OpenURL(url)

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
	if *flgPort != 0 {
		locPort := fmt.Sprintf(":%v", *flgPort)
		http.ListenAndServe(locPort, nil)
	} else {
		http.ListenAndServe(":8080", nil)
	}
}

func createNCWebServer() {
	/*

	 */
	seq := makeSeq(ncFlag)
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
	browser.OpenURL(url)
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

	if *flgPort != 0 {
		locPort := fmt.Sprintf(":%v", *flgPort)
		http.ListenAndServe(locPort, nil)
	} else {
		http.ListenAndServe(":8080", nil)
	}

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

func checkSNPfromFile(f string) {
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
		if *flgVerbose == true {
			fmt.Printf("Reading %v (%v:%v)%v\r", file, i+1, len(files), strings.Repeat(" ", 60))
		}
		snps := parserVCF(file, false, gene.AllGenesVal)
		mapOfVCF[file] = snps
		// bar.Increment()
		// i++
		// fmt.Printf("processed %v files\r", i)
	}
	// bar.Finish()
	// fmt.Println()
	for file, snp := range mapOfVCF {

		// buffer.WriteString(fmt.Sprintf("%v ", file))
		for _, val := range locSNPcheck {

			lGpoS, _ := strconv.Atoi(val.PosInGene)
			CodonNbrInG, _ := strconv.Atoi(val.CodonNbrInG)
			lAPos, _ := strconv.Atoi(val.APos)

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

				}
			}
			// fmt.Println(mapofSNP)

		}

		// fmt.Println((buffer.String()))
		// fmt.Println(chkSNP)
	}
	for key, val := range mapofSNP {
		chkSNP = append(chkSNP, checkSNP{FileName: key, FoundSNP: strings.Join(val, ",")})
	}
	if *flgWeb == true {

		printSNPfromFile(chkSNP)
	} else {
		for _, key := range chkSNP {
			fmt.Printf("%v %v\n", key.FileName, key.FoundSNP)
		}
	}

}

func printSNPfromFile(stat []checkSNP) {
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
			<td>File</td><td>Mutations</td>
			</tr>
			
			{{range $element := .}}
			<tr>		
			<td>{{.FileName}}</td><td>{{.FoundSNP}}</td>
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

	browser.OpenURL(url)

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
	if *flgPort != 0 {
		locPort := fmt.Sprintf(":%v", *flgPort)
		http.ListenAndServe(locPort, nil)
	} else {
		http.ListenAndServe(":8080", nil)
	}
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

// func gcCodonCalc(seq string) (float32, float32, float32, float32) {
// 	var codonsMap = make(map[int]string)
// 	var i, gcCount1, gcCount2, gcCount3, gcCount int
// 	var gc1, gc2, gc3, gc float32
// 	// var codons []string
// 	var buffer strings.Builder
// 	for _, val := range seq {

// 		// if i < 2 {
// 		// codons = append(codons, string(val))

// 		buffer.WriteString(strings.ToUpper(string(val)))
// 		if buffer.Len() == 3 {
// 			i++
// 			// codons = append(codons, buffer.String())
// 			codonsMap[i] = buffer.String()
// 			// fmt.Println(buffer.String())
// 			buffer.Reset()
// 		}

// 	}

// 	for _, val := range codonsMap {

// 		for pos, nuc := range val {

// 			if pos == 0 && nuc == 'C' || pos == 0 && nuc == 'G' {
// 				gcCount1 = gcCount1 + 1

// 			} else if pos == 1 && nuc == 'C' || pos == 1 && nuc == 'G' {
// 				gcCount2 = gcCount2 + 1

// 			} else if pos == 2 && nuc == 'C' || pos == 2 && nuc == 'G' {
// 				gcCount3 = gcCount3 + 1

// 			}
// 			if nuc == 'C' || nuc == 'G' {
// 				gcCount = gcCount + 1

// 			}
// 			// fmt.Println(pos, string(nuc))
// 		}
// 	}
// 	gc = (float32(gcCount) / float32(len(codonsMap)*3)) * 100
// 	gc1 = (float32(gcCount1) / float32(len(codonsMap)*3)) * 100
// 	gc2 = (float32(gcCount2) / float32(len(codonsMap)*3)) * 100
// 	gc3 = (float32(gcCount3) / float32(len(codonsMap)*3)) * 100
// 	// fmt.Printf("gc:%.2f%% gc1:%.2f%% gc2:%.2f%% gc3:%.2f%%\n", gc, gc1, gc2, gc3)
// 	// fmt.Println(codonsMap)
// 	return gc, gc1, gc2, gc3

// }

func testGeneInfo(genes []gene.Gene) {
	var seq string
	var start, end int
	var gc, gc1, gc2, gc3 float32
	// fmt.Println(genes)
	for _, g := range genes {

		start, _ = strconv.Atoi(g.Start)
		end, _ = strconv.Atoi(g.End)
		seq = getNucFromGenome(start, end)
		// fmt.Println(seq)
		if start != 0 && end != 0 {
			gc, gc1, gc2, gc3 = codon.GcCodonCalc(seq)
			fmt.Printf("%v\t%.2f\t%.2f\t%.2f\t%.2f\n", g.Locus, gc, gc1, gc2, gc3)
			// if g.Locus == "Rv0796" {
			// 	fmt.Println("!")
			// }
		}
	}
}

func makeAltString(locus string, positions []allPositionsInGene) string {
	// var lStart, lEnd int
	var seqSplit []string
	var seq string

	// lStart, _ = strconv.Atoi(start)
	// lEnd, _ = strconv.Atoi(end)

	lStart, lEnd := getGenePosByName(locus)
	seq = getNucFromGenome(lStart, lEnd)

	for _, nuc := range seq {
		seqSplit = append(seqSplit, string(nuc))
	}
	for _, val := range positions {
		seqSplit[val.pos-1] = strings.ToUpper(val.nuc)
	}
	// seqSplit[9] = "*"
	// fmt.Println(strings.Join(seqSplit, ""))
	return strings.Join(seqSplit, "")

}

func getGenePosByName(locus string) (int, int) {
	var start, end int
	for _, g := range gene.AllGenesVal {
		if strings.ToUpper(locus) == strings.ToUpper(g.Locus) {
			start, _ = strconv.Atoi(g.Start)
			end, _ = strconv.Atoi(g.End)
			break
		}
	}
	return start, end
}

func getProductByName(locus string) string {
	var prod string
	for _, g := range gene.AllGenesVal {
		if strings.ToUpper(locus) == strings.ToUpper(g.Locus) {
			prod = g.Product
			break
		}
	}
	return prod
}

// func (allpos *allPositionsInGeneArray) AddItem(item allPositionsInGene) []allPositionsInGene {
// 	allpos.items = append(allpos.items, item)
// 	return allpos.items
// }

func toCircos(genes []gene.Gene) {
	var seq string
	var start, end int
	var gc float32
	var buffer strings.Builder
	fmt.Println("---- ideogram.txt ----")

	fmt.Printf("chr -  %v caption %v %v %v\n", gInfo.Strain, gInfo.Start, gInfo.End, gInfo.Strain)

	for _, g := range genes {
		start, _ = strconv.Atoi(g.Start)
		end, _ = strconv.Atoi(g.End)
		seq = getNucFromGenome(start, end)
		// fmt.Println(seq)
		if start != 0 && end != 0 {
			gc, _, _, _ = codon.GcCodonCalc(seq)

			fmt.Printf("band %v %v %v %v %v\n", gInfo.Strain, g.Locus, g.Locus, start, end)
			// if g.Locus == "Rv0796" {
			// 	fmt.Println("!")
			// }
			// fmt.Printf("chr -  %v %v %v %v %v\n", gInfo.Strain, gInfo.Strain, start, end, "lgreen")

			buffer.WriteString(fmt.Sprintf("%v %v %v %.2f\n", gInfo.Strain, start, end, gc))
		}
	}
	// fmt.Printf("band %v 0 0 %v %v\n", gInfo.Strain, gInfo.End-1000, gInfo.End)

	fmt.Println("---- histogram.txt ----")
	fmt.Println(buffer.String())
	// fmt.Println("---- genome.txt ----")
	// for pos, nuc := range genomeSeqSlice {
	// 	fmt.Printf("H37Rv %v %v %v\n", pos+1, pos+1, strings.ToUpper(nuc))
	// }
	// fmt.Println("---- genome.txt ----")
}

func annotatePositions() {
	cols := regexp.MustCompile(`(\w+)\W+(\d+)\W+(\d+)`)
	f, err := os.Open(*flgInFile)

	if err != nil {
		fmt.Println(err)

	}
	defer f.Close()
	fmt.Println("name\tstart\tend\tlocus\tlen\tgc\tseq\tproduct")
	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		for _, colsMatch := range cols.FindAllStringSubmatch(scanner.Text(), -1) {
			// fmt.Printf("c1:%v c2:%v\n", colsMatch[1], colsMatch[2])
			gname, glen := getGeneNameByPos(colsMatch[2], colsMatch[3])
			prod := getProductByName(gname)
			lStart, _ := strconv.Atoi(colsMatch[2])
			lEnd, _ := strconv.Atoi(colsMatch[3])
			seq := getNucFromGenome(lStart-1, lEnd)
			gc, _, _, _ := codon.GcCodonCalc(seq)
			fmt.Printf("%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", colsMatch[1], colsMatch[2], colsMatch[3], gname, glen, gc, seq, prod)

		}
	}

}

func getGeneNameByPos(start, end string) (string, int) {
	var locName string
	var locLen int
	// var genesArr strings.Builder
	locStart, _ := strconv.Atoi(start)
	locEnd, _ := strconv.Atoi(end)
	locLen = (locEnd - locStart) + 1
	for _, g := range gene.AllGenesVal {
		gStart, _ := strconv.Atoi(g.Start)
		gEnd, _ := strconv.Atoi(g.End)
		if locStart >= gStart && locEnd <= gEnd {
			// genesArr.WriteString(g.Locus)
			locName = g.Locus
			// break
			// fmt.Println(locName)
		}
	}
	if locName == "" {
		locName = "-"
		// } else {
		// 	locName = genesArr.String()
	}

	return locName, locLen
}

func makeMatrix() {

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
		fmt.Printf("Working on %v files from %v \r", i+1, len(files))
		// }
		snps := parserVCF(file, false, gene.AllGenesVal)
		for _, val := range snps {
			AllPosUnsort = append(AllPosUnsort, val.APos)
			posCount[val.APos] = posCount[val.APos] + 1
			allLocusUnsort = append(allLocusUnsort, val.Locus)

		}

	}

	AllPos = unique(AllPosUnsort)
	allLocuses = removeStringDuplicates(allLocusUnsort)

	sort.Ints(AllPos)

	switch *flgType {
	case "binary":
		var posFreq = map[int][]string{}
		pos := make(map[int]string)

		headers.WriteString("Pos\t")

		for i, file := range files {
			fmt.Printf("Generating matrix: Made %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

			snps := parserVCF(file, false, gene.AllGenesVal)

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
		// fOut, err := os.Create(*flgOut)
		// if err != nil {
		// 	log.Fatal("Cannot create file", err)
		// }
		// defer fOut.Close()
		// fmt.Fprintf(fOut, headers.String())
		// fmt.Fprintf(fOut, buffer.String())
		// fmt.Println("\n\nWell done!\n")

	case "nuc":
		var posFreq = map[int][]string{}
		pos := make(map[int]string)

		headers.WriteString("Pos\t")

		for i, file := range files {
			fmt.Printf("Generating matrix: Made %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

			snps := parserVCF(file, false, gene.AllGenesVal)

			for _, val := range snps {

				pos[val.APos] = val.Alt

			}

			for _, allpos := range AllPos {

				if pos[allpos] != "" {
					posFreq[allpos] = append(posFreq[allpos], pos[allpos])

				} else {
					posFreq[allpos] = append(posFreq[allpos], strings.ToUpper(getNucFromGenomePos(allpos)))

				}

			}
		}
		for _, allpos := range AllPos {
			// if buffer.Len() == 0 {
			buffer.WriteString(fmt.Sprintln(allpos, "\t", strings.Join(posFreq[allpos], "\t")))

		}
		headers.WriteString("\n")
		// fOut, err := os.Create(*flgOut)
		// if err != nil {
		// 	log.Fatal("Cannot create file", err)
		// }
		// defer fOut.Close()
		// fmt.Fprintf(fOut, headers.String())
		// fmt.Fprintf(fOut, buffer.String())
		// fmt.Println("\n\nWell done!\n")

	case "locus":

		var locFreq = map[string][]string{}

		var locusCount = make(map[string]int)

		headers.WriteString("Locus\t")

		for i, file := range files {
			fmt.Printf("Generating matrix: Made %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(file, filepath.Ext(file))))

			snps := parserVCF(file, false, gene.AllGenesVal)

			for _, val := range snps {

				locusCount[fmt.Sprintf("%v_%v", val.Locus, i)] = locusCount[fmt.Sprintf("%v_%v", val.Locus, i)] + 1

			}
			for _, allloc := range allLocuses {

				locFreq[allloc] = append(locFreq[allloc], strconv.Itoa(locusCount[fmt.Sprintf("%v_%v", allloc, i)]))

			}

		}

		// for i := range files {
		for _, allloc := range allLocuses {

			buffer.WriteString(fmt.Sprintln(allloc, "\t", strings.Join(locFreq[allloc], "\t")))

		}
		// }

		headers.WriteString("\n")
		// fmt.Println(locFreq)
		// fOut, err := os.Create(*flgOut)
		// if err != nil {
		// 	log.Fatal("Cannot create file", err)
		// }
		// defer fOut.Close()
		// fmt.Fprintf(fOut, headers.String())
		// fmt.Fprintf(fOut, buffer.String())
		// fmt.Println("\n\nWell done!\n")

	case "freq":
		headers.WriteString("Pos\tFreq\n")

		for _, allpos := range AllPos {

			// headers.WriteString(fmt.Sprintf("P%v\n", allpos))
			buffer.WriteString(fmt.Sprintf("P%v\t%v\n", allpos, posCount[allpos]))

		}

		// case "dnds":
		// 	// dnds := DnDsRes
		// 	for i, file := range files {
		// 		fmt.Printf("Calculating Dn/DS: Working on %v from %v \r", i+1, len(files))
		// 		calcDnDsVal(file, true)
		// 	}

	}

	if buffer.Len() != 0 && headers.Len() != 0 {
		fOut, err := os.Create(*flgOut)
		if err != nil {
			log.Fatal("Cannot create file", err)
		}
		defer fOut.Close()
		fmt.Fprintf(fOut, headers.String())
		fmt.Fprintf(fOut, buffer.String())
		fmt.Println("\n\nWell done!\n")
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
