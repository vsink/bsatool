package main

import (
	"bufio"
	"bytes"
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
	"strings"
	// "time"
	"flag"
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
	"github.com/pkg/browser"
	"html/template"
	// "io/ioutil"
	"path/filepath"
)

const (
	version = "snpMiner3 0.30.11a"
	list    = "list"
	ncFlag  = "NC"
	aaFlag  = "AA"
	tLMN    = "LMN"   //locus:Mutation:NAME
	tPMLN   = "PMLN"  // position:Mutation:locus:NAME
	tPMN    = "PMN"   // position:Mutation:NAME
	tLSAAN  = "LSAAN" //locus:shortAA:codon:shortAA:name
	tLLAAN  = "LLAAN" // locus:longAA:codon:longAA:name
	tLCN    = "LCN"   //locus:codon:name
	vcfExt  = "*.vcf"
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

var linesFromGB []string             // массив строк считанных из файла genbank
var genomeAnchors = map[int]string{} // карта резльтатов парсинга генбанк-файла
var totalCDSFound int                // количество найденных CDS
var allGenesVal []gene.Gene          // массив информации об генах, загруженный из файла базы данных

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
var flgTestMode = flag.Bool("test", false, "test mode")
var flgDev = flag.Bool("dev", false, "dev mode")

var gInfo genomeInfo

type seqInfo struct {
	Name, Seq string //
	// Len       int
}
type genomeInfo struct {
	TypeOf string
}

var genomeSeqSlice []string // информация об генах, загруженная из базы

func main() {

	// парсинг флагов
	flag.Parse()
	// о программе
	if *flgabout == true {
		about()
	}
	//  файл базы данных
	if *flgMkDB == true {
		if *flgGBfile != "" && *flgOut != "" {
			go process("Working...")
			//
			// if *flgTestMode == true {

			// 	allGenesVal, genomeSeqSlice = readGBFileTest(*flgGBfile)
			// 	// readGBFileTest(*flgGBfile)
			// 	writeDB(*flgOut, allGenesVal)
			// } else {
			allGenesVal, genomeSeqSlice = readGBFile(*flgGBfile)
			writeDB(*flgOut, allGenesVal)
			// }
			// fmt.Printf("%v\n", genomeSeqSlice)
		} else {
			fmt.Println("Не указаны необходимые аргументы (-mkdb -gb=GENBANK_FILE -out=OUTPUT_DB_FILE)")
		}
	}

	if *flgDB != "" {
		// g := Gene

		allGenesVal = readDB(*flgDB)

		switch {
		case *flgVCF != list && *flgVCF != "" && *flgWeb == false && *flgmakeSeq == "":
			parserVCF(*flgVCF, true, allGenesVal)
		case *flgVCF != list && *flgVCF != "" && *flgWeb == true && *flgmakeSeq == "":
			snps := parserVCF(*flgVCF, false, allGenesVal)
			printWebResults(snps)
		case *flgVCF == list && *flgWeb == false && *flgmakeSeq == "" && *flgShare == false && *flgSNP == "":
			listOfVCFFiles()
		case *flgVCF == list && *flgWeb == false && *flgmakeSeq == ncFlag && *flgShare == false:
			seq := makeSeq(ncFlag)
			for _, val := range seq {
				fmt.Println(val.Seq)
			}
		case *flgVCF == list && *flgWeb == true && *flgmakeSeq == ncFlag && *flgShare == false:
			createNCWebServer()

		case *flgVCF == list && *flgWeb == false && *flgmakeSeq == aaFlag && *flgShare == false:
			makeSeq(aaFlag)
		case *flgVCF == list && *flgShare == true:
			// getShareSNP()
			getShareSNP()
		case *flgVCF == list && *flgSNP != "":

			checkSNPfromFile(*flgSNP)

		case *flgDebug == true:
			fmt.Printf("%v\n", allGenesVal)
		}

		if *flgGenomeMap != false {
			var igensS []string
			var igensE []string
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

		for _, g := range allGenesVal {
			if *flgPos != 0 {

				lStart, _ := strconv.Atoi(g.Start)
				lEnd, _ := strconv.Atoi(g.End)

				if *flgPos >= lStart && *flgPos <= lEnd {

					snp := getSNPInfo(*flgPos, g, "A")
					fmt.Printf("apos:%v, gpos:%v, pic:%v, AA:%v, ref:%v, alt:%v, dir:%v, aa:%v\nL:%v,▶%v…%v◀, P:%v\n",
						snp.APos, snp.PosInGene, snp.PosInCodonG, snp.CodonNbrInG, snp.NucInPos, snp.Alt, g.Direction,
						snp.RefAA, snp.Locus, snp.Start, snp.End, snp.Product)

				}

			}
		}

	}

}

func getSNPInfo(apos int, g gene.Gene, alt string) gene.SNPinfo {
	var snp gene.SNPinfo           // структура SNP
	var codonPositions []string    // срез для разбивки кодона побуквенно
	var altCodonPositions []string // срез для разбивки кодона побуквенно альтернативным нуклеотидом
	var locReportType string
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
	} else if aaRefShort != aaAltShort && aaAltShort != "X" {
		mut = "missense"
		tangIdx = tang.GetTangInx(aaRefShort, aaAltShort)
	} else if aaRefShort != aaAltShort && aaAltShort == "X" {
		mut = "nonsense"
	}

	if *flgTang == true {
		locReportType = "T1"
	} else if *flgTang == false {
		locReportType = "T0"
	}
	tang.GetTangInx(aaRefShort, aaAltShort)
	snp = gene.SNPinfo{APos: apos, PosInGene: posInGene, PosInCodonG: posInCodonG,
		RefCodon: codon, RefAA: aaRef, NucInPos: strings.ToUpper(nucG), Locus: g.Locus,
		Direction: g.Direction, Name: g.Name, Product: g.Product,
		Start: g.Start, End: g.End, CodonNbrInG: codonNbrInG, AltCodon: altCodon,
		AltAA: aaAlt, RefAAShort: aaRefShort, AltAAShort: aaAltShort,
		Mutation: mut, Tang: tangIdx, Alt: alt,
		Note: g.Note, ReportType: locReportType, ProteinID: g.ProteinID, GeneID: g.GeneID, GOA: g.GOA}

	return snp
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
// 	// (?m:).*?\s+(.*)\s+(.*)  (?m:)(.*)\s+
// 	// var regexpAnchors = regexp.MustCompile(`(?ms)\s+(.*)`)
// 	var CDS = regexp.MustCompile(`CDS\s+(\d+)\W+(\d+)`)
// 	var CDScompl = regexp.MustCompile(`CDS\s+complement\W(\d+)\W+(\d+)`)
// 	var locus = regexp.MustCompile(`;locus_tag=\W(.*?)\W;`)
// 	var product = regexp.MustCompile(`;product=\W(.*?)\W;`)
// 	var geneName = regexp.MustCompile(`;gene=\W(.*?)\W;`)
// 	// var gseq = regexp.MustCompile(`^\W+(\d+)\s+(\w{6}.*)`)
// 	var qenomeSeq = regexp.MustCompile(`\s+(\d+)\s+(\w{6}.*)`)
// 	var geneID = regexp.MustCompile(`;db_xref=\WGeneID:(\d+)`)
// 	var geneGOA = regexp.MustCompile(`;db_xref=\WGOA:(\w+)"`)
// 	var proteinID = regexp.MustCompile(`;protein_id=\W(.*?)\W;`)
// 	var Note = regexp.MustCompile(`;note=\W(.*?)\W;`)
// 	var resString []string
// 	var querySplice []string
// 	var splitedGenome []string
// 	var originBlock int
// 	var lStart, lEnd, lLoc, lProd, lDir, lName, gID, pID, lNote, lGOA string
// 	// var CDSfound bool
// 	var CDSblock int
// 	var buffer bytes.Buffer
// 	var nucCore = "no"

// 	f, err := os.Open(file)
// 	if err != nil {
// 		fmt.Println(err)

// 	}
// 	// go process("Чтение файла...")
// 	gb := bufio.NewReader(f)
// 	for {
// 		line, err := gb.ReadBytes('\n')
// 		if err == io.EOF {
// 			break
// 		}
// 		if err != nil {
// 			panic(err)
// 		}

// 		if strings.Contains(string(line), "  CDS ") {
// 			CDSblock = 1

// 		} else if strings.Contains(string(line), "     gene") {
// 			CDSblock = 0
// 		}

// 		if strings.Contains(string(line), "ORIGIN      ") {
// 			originBlock = 1
// 		} else if strings.Contains(string(line), "//") {
// 			originBlock = 0
// 		}

// 		if originBlock == 1 {

// 			for _, genomeMatch := range qenomeSeq.FindAllStringSubmatch(string(line), -1) {

// 				resString = append(resString, strings.Replace(genomeMatch[2], " ", "", -1))
// 			}

// 		}

// 		if CDSblock == 1 {
// 			// fmt.Printf("cds%v org%v\n", CDSblock, originBlock)
// 			// fmt.Printf("%v\n", string(line))
// 			changedStr := strings.TrimSuffix(strings.Replace(string(line), "  ", "", -1), "\n")

// 			buffer.WriteString(strings.Replace(changedStr, " /", ";", -1))

// 			// fmt.Println(string(line))
// 		} else if CDSblock == 0 && buffer.Len() != 0 {

// 			querySplice = append(querySplice, buffer.String())

// 			buffer.Reset()
// 		}

// 	}

// 	splitedGenome = strings.SplitAfter(strings.Join(resString, ""), "")
// 	// fmt.Println(querySplice)
// 	defer f.Close()
// 	// fmt.Printf("%v\n", querySplice)
// 	for i, val := range querySplice {

// 		for _, cdsMatch := range CDS.FindAllStringSubmatch(val, -1) {

// 			lStart, lEnd, lDir = cdsMatch[1], cdsMatch[2], "f"

// 		}
// 		for _, cdscMatch := range CDScompl.FindAllStringSubmatch(val, -1) {

// 			lStart, lEnd, lDir = cdscMatch[1], cdscMatch[2], "r"

// 		}
// 		for _, locusMatch := range locus.FindAllStringSubmatch(val, -1) {

// 			lLoc = locusMatch[1]

// 		}
// 		for _, gnameMatch := range geneName.FindAllStringSubmatch(val, -1) {

// 			lName = gnameMatch[1]

// 		}
// 		for _, productMatch := range product.FindAllStringSubmatch(val, -1) {

// 			lProd = productMatch[1]
// 		}

// 		for _, geneIDMatch := range geneID.FindAllStringSubmatch(val, -1) {

// 			gID = geneIDMatch[1]
// 		}

// 		for _, protIDMatch := range proteinID.FindAllStringSubmatch(val, -1) {

// 			pID = protIDMatch[1]
// 		}
// 		for _, geneGOAmatch := range geneGOA.FindAllStringSubmatch(val, -1) {

// 			lGOA = geneGOAmatch[1]
// 			if len(lGOA) != 0 {
// 				nucCore = "yes"
// 			}
// 		}
// 		for _, noteMatch := range Note.FindAllStringSubmatch(val, -1) {

// 			lNote = noteMatch[1]
// 			// fmt.Println(lNote)
// 		}
// 		if lName == "" {
// 			lName = lLoc
// 		}

// 		g := gene.Gene{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir,
// 			Product: lProd, Name: lName, GeneID: gID, ProteinID: pID, Note: lNote, GOA: lGOA}

// 		allGenesVal = append(allGenesVal, g)

// 		if *flgDebug == true {

// 			fmt.Printf("l:%v s:%v e:%v p:%v gId:%v pId:%v n:%v GOA:%v\n", lLoc, lStart, lEnd, lProd, gID, pID, lNote, lGOA)
// 		} else {
// 			fmt.Printf("\rобнаружено %v кодирующих участка, nucCore:%v", i-1, nucCore)
// 		}
// 		// обнуление
// 		lName, lStart, lEnd, lLoc, lDir, lProd, gID, pID, lGOA, lNote = "", "", "", "", "", "", "", "", "", ""

// 		// )
// 	}
// 	gInfo = genomeInfo{TypeOf: nucCore}
// 	// go process("Анализ файла закончен!")
// 	return allGenesVal, splitedGenome

// }

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
	// var getNote = regexp.MustCompile(`(?s)!!note=(.*)`)
	var resString []string

	var querySplice []string
	var splitedGenome []string
	var originBlock int
	// var lStart, lEnd, lLoc, lProd, lDir, lName, gID, pID, lNote, lGOA string
	// var CDSfound bool
	var CDSblock int
	// var buffer bytes.Buffer
	// var nucCore = "no"
	// var catchNote = regexp.MustCompile(`(?s)/note=\W(.*)`)
	var changedStr string
	var noteBuffer bytes.Buffer

	f, err := os.Open(file)

	if err != nil {
		fmt.Println(err)

	}
	defer f.Close()

	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		scanTxt := scanner.Text() // Println will add back the final '\n'
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
		changedStr = regDelSpaces.ReplaceAllString(string(scanTxt), "$2 ") // удаляем пробелы

		switch originBlock {
		case 0:
			if CDSblock == 1 {
				changedStr = strings.Replace(string(changedStr), "/note=", "!note=", -1)       // меняем / на ! для дальнейшего парсинга
				changedStr = strings.Replace(string(changedStr), "/product=", "!product=", -1) // см выше.
				changedStr = strings.Replace(string(changedStr), "\"", "", -1)
				// changedStr = strings.Replace(string(changedStr), "/", "!!", -1)
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

			if strings.Index(changedStr, "!!") == 0 || strings.Index(changedStr, "CDS  ") == 0 {
				querySplice = append(querySplice, strings.TrimPrefix(changedStr, "!!"))
				// fmt.Println(changedStr)

			}

		}

	}
	splitedGenome = strings.SplitAfter(strings.Join(resString, ""), "")
	// for _, val := range querySplice {
	// 	fmt.Println(val)
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
	var cdsStartFrame, bad int
	var cdsCount = 1

	for _, val := range querySplice {

		// bad = 1

		for _, cdsStEndMatch := range cdsStEnd.FindAllStringSubmatch(val, -1) {
			// cdsBlock = 1
			// fmt.Printf("%v %v\n", cdsBlock, val)

			if len(cdsStEndMatch[1]) != 0 {
				lStart, lEnd, lDir = strings.Replace(cdsStEndMatch[1], " ", "", -1), strings.Replace(cdsStEndMatch[2], " ", "", -1), "f"

			} else {
				lStart, lEnd, lDir = strings.Replace(cdsStEndMatch[3], " ", "", -1), strings.Replace(cdsStEndMatch[4], " ", "", -1), "r"

			}

			cdsStartFrame = 1
			bad = 1

			// fmt.Println(lStart, lEnd, lDir)

		}

		for _, cdsLocusMatch := range cdsLocus.FindAllStringSubmatch(val, -1) {
			lLoc = strings.Replace(cdsLocusMatch[1], " ", "", -1)
			// fmt.Println(lLoc)
			cdsStartFrame = 1
			bad = 0
		}

		// cdsLocusMatch := cdsLocus.FindAllStringSubmatch(val, -1)
		// if len(cdsLocusMatch) != 0 {
		// 	lLoc = cdsLocusMatch[0][1]
		// }

		for _, cdsNameMatch := range cdsName.FindAllStringSubmatch(val, -1) {
			lName = strings.Replace(cdsNameMatch[1], " ", "", -1)
			cdsStartFrame = 1
			bad = 0
			// fmt.Println(lName)
		}
		for _, cdsProdMatch := range cdsProd.FindAllStringSubmatch(val, -1) {
			lProd = cdsProdMatch[1]
			cdsStartFrame = 1
			bad = 0
			// fmt.Println(lProd)
		}
		for _, cdsNoteMatch := range cdsNote.FindAllStringSubmatch(val, -1) {
			lNote = cdsNoteMatch[1]

			// fmt.Println(lNote)
		}
		for _, cdsgIDMatch := range cdsgID.FindAllStringSubmatch(val, -1) {
			gID = strings.Replace(cdsgIDMatch[1], " ", "", -1)
			cdsStartFrame = 1
			bad = 0
			// fmt.Println(gID)

		}
		for _, cdsgprotIDMatch := range cdsprotID.FindAllStringSubmatch(val, -1) {
			pID = strings.Replace(cdsgprotIDMatch[1], " ", "", -1)
			cdsStartFrame = 1
			bad = 0
			// fmt.Println(pID)

		}
		for _, cdsGOAMatch := range cdsGOA.FindAllStringSubmatch(val, -1) {
			lGOA = strings.Replace(cdsGOAMatch[1], " ", "", -1)

			// fmt.Println(lGOA)
			cdsStartFrame = 1
			bad = 0
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
		if cdsStartFrame == 1 && bad == 1 {
			cdsCount++
			g := gene.Gene{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir,
				Product: lProd, Name: lName, GeneID: gID, ProteinID: pID, Note: lNote, GOA: lGOA}

			allGenesVal = append(allGenesVal, g)
			// fmt.Printf("g:%v b:%v %v\n", cdsStartFrame, bad, val)
			if *flgDebug == true {

				fmt.Printf("l:%v s:%v e:%v d:%v p:%v gId:%v pId:%v n:%v GOA:%v\n", lLoc, lStart, lEnd, lDir, lProd, gID, pID, lNote, lGOA)
			} else {
				fmt.Printf("\rFound %v CDS, nucCore:%v", cdsCount-1, nucCore)
			}
			lName, lStart, lEnd, lLoc, lDir, lProd, gID, pID, lGOA, lNote = "", "", "", "", "", "", "", "", "", ""

		}

	}

	gInfo = genomeInfo{TypeOf: nucCore}
	// go process("Анализ файла закончен!")
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
						snp := getSNPInfo(apos, g, alt)
						// br := testing.Benchmark(snp)
						// fmt.Println(br)

						if len(ref) == 1 && len(alt) == 1 {
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
	fmt.Println(strings.Repeat(".", 50), "\n", version, "\n", strings.Repeat(".", 50), "\n\u00A9 V.Sinkov & O.Ogarkov,Irkutsk, Russia, 2017")
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
			fmt.Printf("\n\n%v:\n\n", file)
			parserVCF(file, true, allGenesVal)

			// fmt.Println(file.Name())
		}
	}
}

func makeSeq(typeof string) []seqInfo {

	var AllPos []int
	var ResSeq []seqInfo

	files := getListofVCF()

	for i, file := range files {

		if *flgVerbose == true {
			fmt.Printf("Pass: %v from %v \r", i+1, len(files))
		}
		snps := parserVCF(file, false, allGenesVal)
		switch typeof {
		case ncFlag:
			for _, val := range snps {
				AllPos = append(AllPos, val.APos)

			}
		case aaFlag:

		}

	}

	sort.Ints(removeDuplicates(AllPos))
	if *flgWithRef == true {
		switch typeof {
		case ncFlag:
			var refBuffer bytes.Buffer
			refBuffer.WriteString(fmt.Sprintf(">%v\n", "REFERENCE"))
			for _, allpos := range AllPos {
				refBuffer.WriteString(getNucFromGenomePos(allpos))
			}
			ResSeq = append(ResSeq, seqInfo{Name: "reference", Seq: refBuffer.String()})

		}
	}
	for _, file := range files {
		pos := make(map[int]string)
		var buffer bytes.Buffer

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

func removeDuplicates(elements []int) []int {
	// Use map to record duplicates as we find them.
	encountered := map[int]bool{}
	result := []int{}

	for v := range elements {
		if encountered[elements[v]] == true {
			// Do not add duplicate.
		} else {
			// Record this element as an encountered element.
			encountered[elements[v]] = true
			// Append to result slice.
			result = append(result, elements[v])
		}
	}
	// Return the new slice.
	return result
}

func getInterGen(pos int) {
	var igensS []string
	var igensE []string
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
			<td>{{.RefCodon}}/{{.AltCodon}}</td><td>{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}</td>
			<td><p title="Tang Index: {{.Tang}}">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ncbi.nlm.nih.gov/gene/{{.GeneID}}={{.GeneID}}">{{.GeneID}}</a>
			</td><td><a href="http://www.uniprot.org/uniprot/?query={{.ProteinID}}&sort=score">{{.ProteinID}}</td>
			</tr>
			{{else}}
			<tr>
			<td><p title="{{.Note}}">{{.Locus}}</p></td><td>{{.Name}}</td><td>{{.APos}}</td><td>{{.PosInGene}}{{.NucInPos}}>{{.Alt}}</td>
			<td>{{.RefCodon}}/{{.AltCodon}}</td><td>{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}</td>
			<td><p title="Tang Index: {{.Tang}}">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ebi.ac.uk/QuickGO/GProtein?ac={{.GOA}}">{{.GOA}}</a>
			</td><td><a href="http://www.uniprot.org/uniprot/?query={{.ProteinID}}&sort=score">{{.ProteinID}}</td>
			</tr>
			{{end}}
			{{end}}
			</tbody>
			</table>
			<table width="100%" cellspacing="0" cellpadding="4" border="0">
			<tr>
			<td><a href="http://snpminer.ru">Created by snpMiner3</a></td>
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

		defer os.Exit(0)
	})
	browser.OpenURL(url)
	// fmt.Println("Для выхода из программы, нажмите Ctl+C")
	http.ListenAndServe(":8080", nil)
	// log.Fatal(http.ListenAndServe(":8080", nil))
	// os.Exit(0)
	// os.Exit(1)
	// for val := range tmp {
	// 	fmt.Printf("%v\n", tmp[val].Locus)
	// }
	// os.Exit(0)

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
		snps := parserVCF(file, false, allGenesVal)

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

		// os.Exit(0)
	})
	browser.OpenURL(url)
	http.ListenAndServe(":8080", nil)

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
	files := getListofVCF()
	locSNPcheck := getSNPNotation(f)

	// bar := pb.ProgressBarTemplate(pbtmpl).Start(len(files))
	// bar.SetWidth(90)
	for i, file := range files {
		if *flgVerbose == true {
			fmt.Printf("Pass: %v from %v \r", i+1, len(files))
		}
		snps := parserVCF(file, false, allGenesVal)
		mapOfVCF[file] = snps
		// bar.Increment()
		// i++
		// fmt.Printf("processed %v files\r", i)
	}
	// bar.Finish()
	fmt.Println()
	for file, snp := range mapOfVCF {
		var buffer bytes.Buffer

		buffer.WriteString(fmt.Sprintf("%v ", file))
		for _, val := range locSNPcheck {

			lGpoS, _ := strconv.Atoi(val.PosInGene)
			CodonNbrInG, _ := strconv.Atoi(val.CodonNbrInG)
			lAPos, _ := strconv.Atoi(val.APos)

			for _, snpFromFile := range snp {

				switch val.TypeOf {
				case tLMN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && strings.ToUpper(val.Alt) == strings.ToUpper(snpFromFile.Alt) && lGpoS == snpFromFile.PosInGene {
						// fmt.Println(val.Locus, "\t", val.PosInGene, "\t", snpFromFile.PosInGene, "\t", snpFromFile.Locus)
						buffer.WriteString(fmt.Sprintf("%v[%v:%v%v>%v]\t", val.Name, val.Locus, val.PosInGene, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))
						// fmt.Println(strings.ToUpper(buffer.String()))
					}
				case tPMLN:

					if lAPos == snpFromFile.APos && strings.ToUpper(val.Alt) == strings.ToUpper(snpFromFile.Alt) {
						buffer.WriteString(fmt.Sprintf("%v[%v:%v_%v>%v]\t", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))

					}
				case tLSAAN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && CodonNbrInG == snpFromFile.CodonNbrInG &&

						val.AASalt == snpFromFile.AltAAShort {

						buffer.WriteString(fmt.Sprintf("%v[%v:%v%v%v]\t", val.Name, val.Locus, val.AASref, CodonNbrInG, val.AASalt))
					}
				case tLLAAN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && CodonNbrInG == snpFromFile.CodonNbrInG &&
						val.AALalt == snpFromFile.AltAA {
						buffer.WriteString(fmt.Sprintf("%v[%v:%v%v%v]\t", val.Name, val.Locus, val.AALref, CodonNbrInG, val.AALalt))

					}
				case tLCN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && CodonNbrInG == snpFromFile.CodonNbrInG {
						buffer.WriteString(fmt.Sprintf("%v[%v:codon%v]\t", val.Name, val.Locus, CodonNbrInG))
					}

				}
			}

		}
		fmt.Println((buffer.String()))
	}
}
