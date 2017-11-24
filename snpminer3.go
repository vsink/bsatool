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
	"io"
	"sort"
	"strconv"

	"./amino"
	"./gene"
	"./tang"
	// "text/template"
	// "net/http"
	// "github.com/davecgh/go-spew/spew"
	// tm "github.com/buger/goterm"
	"html/template"
	"path/filepath"
)

const (
	version = "snpMiner3 0.17.11a"
	list    = "list"
	ncFlag  = "NC"
	aaFlag  = "AA"
	tLMN    = "LMN"
	tPMLN   = "PMLN"
	tLSAAN  = "LSAAN"
	tLLAAN  = "LLAAN"
	vcfExt  = "*.vcf"
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
var flgWeb = flag.Bool("web", false, "Режим вывода результатов в Веб браузер")             // вывод результатов в браузер
var flgUniq = flag.Bool("uniq", false, "Поиск общих, для всех vcf файлов в папке, СНИПОВ") // поиск общих снипов для всех Vcf файлов в папке
var flgIndel = flag.Bool("indel", false, "Поиск ИнДелов")
var flgWithRef = flag.Bool("ref", false, "Создавать референсную последовательность при использовании команды -mkseq")

type seqInfo struct {
	Name, Seq string
	// Len       int
}

type uniqInfo struct {
	Apos, Count int
	Alt         string
}
type uniqInfoTest struct {
}

// var GenomeAccessNbr string

// var GenomeSeq string
var genomeSeqSlice []string

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
			go process("Идет обработка...")
			//
			allGenesVal, genomeSeqSlice = readGBFile(*flgGBfile)
			writeDB(*flgOut, allGenesVal)
			// fmt.Printf("%v\n", genomeSeqSlice)
		} else {
			fmt.Println("Не указаны необходимые аргументы (-mkdb -gb=GENBANK_FILE -out=OUTPUT_DB_FILE)")
		}
	}

	if *flgDB != "" {
		// g := Gene

		allGenesVal = readDB(*flgDB)

		// if *flgmakeSeq != false && *flgOut != "" {

		// 	fmt.Println("Reserved!")
		// }
		switch {
		case *flgVCF != list && *flgVCF != "" && *flgWeb == false && *flgmakeSeq == "":
			parserVCF(*flgVCF, true, allGenesVal)
		case *flgVCF != list && *flgVCF != "" && *flgWeb == true && *flgmakeSeq == "":
			snps := parserVCF(*flgVCF, false, allGenesVal)
			printWebResults(snps)
		case *flgVCF == list && *flgWeb == false && *flgmakeSeq == "" && *flgUniq == false && *flgSNP == "":
			listOfVCFFiles()
		case *flgVCF == list && *flgWeb == false && *flgmakeSeq == ncFlag && *flgUniq == false:
			seq := makeSeq(ncFlag)
			for _, val := range seq {
				fmt.Println(val.Seq)
			}
		case *flgVCF == list && *flgWeb == true && *flgmakeSeq == ncFlag && *flgUniq == false:
			createNCWebServer()
			// for _, val := range seq {
			// 	fmt.Println(val.Seq)
			// }
		case *flgVCF == list && *flgWeb == false && *flgmakeSeq == aaFlag && *flgUniq == false:
			makeSeq(aaFlag)
		case *flgVCF == list && *flgUniq == true:
			// getUniqSNP()
			getUniqSNP()
		case *flgVCF == list && *flgSNP != "":
			// locSNPcheck := getSNPNotation(*flgSNP)
			// for _, val := range locSNPcheck {
			// 	fmt.Printf("TypeOf:%v %v \n", val.TypeOf, val.Raw)
			// }
			checkSNPfromFile(*flgSNP)
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
		// var tmp_slited []string
		for _, g := range allGenesVal {
			// fmt.Printf("%v %v \n", i.Locus, getNucFromGenome(genomeSeqSlice, i.Start, i.End))
			if *flgPos != 0 {

				lStart, _ := strconv.Atoi(g.Start)
				lEnd, _ := strconv.Atoi(g.End)

				if *flgPos >= lStart && *flgPos <= lEnd {
					// fmt.Printf("s:%v e:%v L:%v D:%v P:%v\n", g.Start, g.End, g.Locus, g.Direction, g.Product)
					if *flgDebug == true {
						// fmt.Printf(getSNPInfo(*flgPos, g))
						snp := getSNPInfo(*flgPos, g, "A")
						fmt.Printf("apos:%v, gpos:%v, pic:%v, AA:%v, ref:%v, nuc:%v, dir:%v, aa:%v\nL:%v, S:%v, E:%v, P:%v\n",
							snp.Apos, snp.Gpos, snp.PosInCodonG, snp.CodonNbrInG, snp.RefCodon, snp.NucInPos, g.Direction,
							snp.RefAA, snp.Locus, snp.Start, snp.End, snp.Product)

					}

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
		// fmt.Printf("alt:%v codon:%v alt_r:%v rev_codon:%v\n", alt, codon, codonReverse(alt), codonReverse(codon))
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
	// tmp := strconv.Itoa(apos)
	// if tmp == "386432" {
	// if g.Product == "" {
	// 	fmt.Printf("!!!%v\t%v\t%v\t%v\t%v\t%v\t%v\n", g.Locus, apos, nucG, alt, g.Direction, codon, altCodon)
	// }
	// }
	if *flgTang == true && g.Direction == "r" {
		locReportType = "T1r"
	} else if *flgTang == true && g.Direction == "f" {
		locReportType = "T1f"
	} else if *flgTang == false && g.Direction == "r" {
		locReportType = "T0r"
	} else if *flgTang == false && g.Direction == "f" {
		locReportType = "T0f"
	}
	tang.GetTangInx(aaRefShort, aaAltShort)
	snp = gene.SNPinfo{Apos: apos, Gpos: posInGene, PosInCodonG: posInCodonG,
		RefCodon: codon, RefAA: aaRef, NucInPos: strings.ToUpper(nucG), Locus: g.Locus,
		Direction: g.Direction, Name: g.Name, Product: g.Product,
		Start: g.Start, End: g.End, CodonNbrInG: codonNbrInG, AltCodon: altCodon,
		AltAA: aaAlt, RefAAShort: aaRefShort, AltAAShort: aaAltShort,
		Mutation: mut, CCodonNbrInG: CCodonNbrInG, CPosInGene: CPosInGene, Tang: tangIdx, Alt: alt,
		Note: g.Note, ReportType: locReportType, ProteinID: g.ProteinID, GeneID: g.GeneID}
	// fmt.Printf("%v\t%v\n", aaRefShort, aaAltShort)
	// fmt.Printf("pig:%v, pic:%v, codon:%v %v%v\n", posInGene, posInCodonG, codon, aaRef, codonNbrInG)
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

func readGBFile(file string) (g []gene.Gene, genomeSplice []string) {
	//функция для считывания файла в формате генбанка и занесения строк в массив linesFromGB
	// (?m:).*?\s+(.*)\s+(.*)  (?m:)(.*)\s+
	// var regexpAnchors = regexp.MustCompile(`(?ms)\s+(.*)`)
	var CDS = regexp.MustCompile(`CDS\s+(\d+)\W+(\d+)`)
	var CDScompl = regexp.MustCompile(`CDS\s+complement\W(\d+)\W+(\d+)`)
	var locus = regexp.MustCompile(`;locus_tag=\W(.*?)\W;`)
	var product = regexp.MustCompile(`;product=\W(.*?)\W;`)
	var geneName = regexp.MustCompile(`;gene=\W(.*?)\W;`)
	// var gseq = regexp.MustCompile(`^\W+(\d+)\s+(\w{6}.*)`)
	var qenomeSeq = regexp.MustCompile(`\s+(\d+)\s+(\w{6}.*)`)
	var geneID = regexp.MustCompile(`;db_xref=\WGeneID:(\d+)`)
	var proteinID = regexp.MustCompile(`;protein_id=\W(.*?)\W;`)
	var Note = regexp.MustCompile(`;note=\W(.*?)\W;`)
	var resString []string
	var querySplice []string
	var splitedGenome []string
	var originBlock int
	var lStart, lEnd, lLoc, lProd, lDir, lName, gID, pID, lNote string
	// var CDSfound bool
	var CDSblock int
	var buffer bytes.Buffer
	// var queryString string
	// var buffer bytes.Buffer
	f, err := os.Open(file)
	if err != nil {
		fmt.Println(err)

	}
	// go process("Чтение файла...")
	gb := bufio.NewReader(f)
	for {
		line, err := gb.ReadBytes('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}

		if strings.Contains(string(line), "  CDS ") {
			CDSblock = 1
		} else if strings.Contains(string(line), "translation=") {
			CDSblock = 0
		}

		if strings.Contains(string(line), "ORIGIN      ") {
			originBlock = 1
		} else if strings.Contains(string(line), "//") {
			originBlock = 0
		}

		if originBlock == 1 {
			for _, genomeMatch := range qenomeSeq.FindAllStringSubmatch(string(line), -1) {
				// fmt.Printf("seq:%v\n", genomeMatch[1])
				resString = append(resString, strings.Replace(genomeMatch[2], " ", "", -1))
			}
		}

		if CDSblock == 1 {
			changedStr := strings.TrimSuffix(strings.Replace(string(line), "  ", "", -1), "\n")
			buffer.WriteString(strings.Replace(changedStr, " /", ";", -1))
			// fmt.Printf("%v\n", buffer.String())
		} else if CDSblock == 0 && buffer.Len() != 0 {
			// fmt.Printf("b:%v\n", strings.TrimSuffix(, "\n\n"))
			querySplice = append(querySplice, buffer.String())
			buffer.Reset()
		}

	}
	splitedGenome = strings.SplitAfter(strings.Join(resString, ""), "")

	defer f.Close()

	for i, val := range querySplice {
		// fmt.Printf("\t\t\t%v\n", val)
		for _, cdsMatch := range CDS.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("s:%v\te:%v d:%v\n", cdsMatch[1], cdsMatch[2], "f")
			lStart, lEnd, lDir = cdsMatch[1], cdsMatch[2], "f"

		}
		for _, cdscMatch := range CDScompl.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("s:%v\te:%v d:%v\n", cdscMatch[1], cdscMatch[2], "r")
			lStart, lEnd, lDir = cdscMatch[1], cdscMatch[2], "r"

		}
		for _, locusMatch := range locus.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("l:%v\n", locusMatch[1])
			lLoc = locusMatch[1]
			// if lLoc == "Rv0068" {
			// 	fmt.Println("!!!!!!\n")
			// }
		}
		for _, gnameMatch := range geneName.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("g:%v\n", gnameMatch[1])
			lName = gnameMatch[1]

		}
		for _, productMatch := range product.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("p:%v\n", productMatch[1])
			lProd = productMatch[1]
		}

		for _, geneIDMatch := range geneID.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("p:%v\n", productMatch[1])
			gID = geneIDMatch[1]
		}

		for _, protIDMatch := range proteinID.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("p:%v\n", productMatch[1])
			pID = protIDMatch[1]
		}
		for _, noteMatch := range Note.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("p:%v\n", productMatch[1])
			lNote = noteMatch[1]
		}
		if lName == "" {
			lName = lLoc
		}

		g := gene.Gene{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir,
			Product: lProd, Name: lName, GeneID: gID, ProteinID: pID, Note: lNote}
		// if *flgDebug == true {
		// fmt.Printf("%v\n", g)
		// }
		allGenesVal = append(allGenesVal, g)
		if *flgDebug == true {

			fmt.Printf("l:%v s:%v e:%v p:%v gId:%v pId:%v n:%v\n", lLoc, lStart, lEnd, lProd, gID, pID, lNote)
		} else {
			fmt.Printf("\rобнаружено %v кодирующих участка", i-1)
		}
		// обнуление
		lName, lStart, lEnd, lLoc, lDir, lProd, gID, pID = "", "", "", "", "", "", "", ""

		// )
	}
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
	fmt.Println("\nРезультаты сохранены в", file, "файл.")
	// taskFile.Write(jsonParser)

	// fmt.Println(strings.Repeat("\nБаза успешно создана!")
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

func processReverseCodon(codon string, codonPosition int) string {
	var lCodons []string
	var result string
	lCodons = strings.Split(strings.ToLower(codon), "")
	for i, val := range lCodons {
		if i == codonPosition {
			val = strings.ToUpper(val)
		}
		result = result + val
		// alt_res = alt_res + alt_val

	}
	return result
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

		// result = string(n) + result

	}
	return result
	// fmt.Println(result, "\n")
	// for i, letter := range lCodonSplit {
	// 	fmt.Printf(">> %v\n", letter[i])
	// }
}

//

func parserVCF(f string, print bool, genes []gene.Gene) []gene.SNPinfo {
	var vcf = regexp.MustCompile(`^\S+\s+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)`)
	var indel = regexp.MustCompile(`^^\S+\W+(\d+)\W+(\w+)\s+(\w+).*(INDEL).*DP=(\d+)`)
	var validateVCF = regexp.MustCompile(`(##fileformat)=VCF`)
	var vcfValid bool
	// var vcf = regexp.MustCompile(`(^\S+)\W(\d+)\W+(\w+)\W+(\w+)`)

	//
	var snpFromVCF []gene.SNPinfo
	// x := make(map[string][]gene.SNPinfo)

	file, err := os.Open(f)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		for _, vcfvalid := range validateVCF.FindAllStringSubmatch(scanner.Text(), -1) {
			// fmt.Println(vcfvalid[0])
			if vcfvalid[0] == "##fileformat=VCF" {
				vcfValid = true
			}

		}
		if vcfValid == false {
			fmt.Printf("%v is wrong VCF file!!! Check it!\n", file.Name())
			break
		}
		if *flgIndel == true {
			for _, indelMatch := range indel.FindAllStringSubmatch(scanner.Text(), -1) {
				if vcfValid == true && indelMatch[4] == "INDEL" {
					fmt.Println(indelMatch[1])

				}

			}
		}

		for _, match := range vcf.FindAllStringSubmatch(scanner.Text(), -1) {

			if vcfValid == true {
				apos, _ := strconv.Atoi(match[1])
				ref := match[2]
				alt := match[3]
				// fmt.Printf("%v\n", genes)
				// orgInfo = GenomeInfo{Organism: match[1]}
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
	fmt.Println(str)
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
	var i int
	// files, err := filepath.Glob(vcfExt)
	files := getListofVCF()

	// if err != nil {
	// 	log.Fatal(err)
	// }
	for _, file := range files {
		i++
		fmt.Printf("processed %v files\r", i)
		snps := parserVCF(file, false, allGenesVal)
		switch typeof {
		case ncFlag:
			for _, val := range snps {
				// buffer.WriteString(val.Alt)
				AllPos = append(AllPos, val.Apos)
				// fmt.Printf("%v\n", val.Alt)
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
			// fmt.Println(ResSeq)
		}
	}
	for _, file := range files {
		pos := make(map[int]string)
		var buffer bytes.Buffer

		buffer.WriteString(fmt.Sprintf(">%v\n", strings.ToUpper(file)))
		// var buffer bytes.Buffer
		// fmt.Printf(">%v\n", strings.ToUpper(file))
		// buffer.WriteString("\n>" + strings.ToUpper(file) + "\n")
		snps := parserVCF(file, false, allGenesVal)
		switch typeof {
		case ncFlag:
			for _, val := range snps {
				pos[val.Apos] = val.Alt

			}
			for _, allpos := range AllPos {

				// fmt.Println(i)
				if pos[allpos] != "" {
					// fmt.Println(pos[allpos])
					buffer.WriteString(pos[allpos])
				} else {
					// fmt.Println(getNucFromGenomePos(allpos))
					buffer.WriteString(getNucFromGenomePos(allpos))
				}
			}
			ResSeq = append(ResSeq, seqInfo{Name: file, Seq: buffer.String()})
			// ResSeq = append(ResSeq, seqInfo{Name: "REFERENCE", Seq: refBuffer.String()})

		case aaFlag:
			// for _, val := range snps {
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

	case "T1r":
		fmt.Printf("%v\t%v\t%v%v>%v\t%v/%v\t%v%v%v\t%v\t%v\t%v\n", snp.Locus, snp.Apos, snp.CPosInGene, snp.NucInPos, snp.Alt, snp.RefCodon, snp.AltCodon,
			snp.RefAAShort, snp.CCodonNbrInG, snp.AltAAShort, snp.Mutation,
			snp.Tang, snp.Product)

	case "T1f":
		fmt.Printf("%v\t%v\t%v%v>%v\t%v/%v\t%v%v%v\t%v\t%v\t%v\n", snp.Locus, snp.Apos, snp.Gpos, snp.NucInPos, snp.Alt, snp.RefCodon, snp.AltCodon,
			snp.RefAAShort, snp.CodonNbrInG, snp.AltAAShort, snp.Mutation,
			snp.Tang, snp.Product)
	case "T0r":
		fmt.Printf("%v\t%v\t%v%v>%v\t%v/%v\t%v%v%v\t%v\t%v\n", snp.Locus, snp.Apos, snp.CPosInGene, snp.NucInPos, snp.Alt, snp.RefCodon, snp.AltCodon,
			snp.RefAAShort, snp.CCodonNbrInG, snp.AltAAShort, snp.Mutation,
			snp.Product)
	case "T0f":
		fmt.Printf("%v\t%v\t%v%v>%v\t%v/%v\t%v%v%v\t%v\t%v\n", snp.Locus, snp.Apos, snp.Gpos, snp.NucInPos, snp.Alt, snp.RefCodon, snp.AltCodon,
			snp.RefAAShort, snp.CodonNbrInG, snp.AltAAShort, snp.Mutation,
			snp.Product)

	}

}

func printWebResults(snps []gene.SNPinfo) {

	var htmlTemplate = `
				<!DOCTYPE html>
				<html>


				<table width="100%" cellspacing="0" cellpadding="4" border="1">
				<tbody>
			<tr>
			<td>Locus</td><td>Gene</td><td>Pos.</td><td>Mutation</td><td>Codons</td>
			<td>AA</td><td>Type</td><td>Product</td><td>GeneID</td><td>ProteinID</td>
			</tr>
			{{range $element := .}}
			<tr>
			<td>{{.Locus}}</td><td>{{.Name}}</td><td>{{.Apos}}</td><td>{{.CPosInGene}}{{.NucInPos}}>{{.Alt}}</td>
			<td>{{.RefCodon}}/{{.AltCodon}}</td><td>{{.RefAAShort}}{{.CCodonNbrInG}}{{.AltAAShort}}</td>
			<td>{{.Mutation}}</td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}">{{.Product}}</a></td><td><a href="https://www.ncbi.nlm.nih.gov/gene/{{.GeneID}}">{{.GeneID}}</a>
			</td><td><a href="http://www.uniprot.org/uniprot/?query={{.ProteinID}}&sort=score">{{.ProteinID}}</td>
			</tr>
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
	t, err := t.Parse(htmlTemplate)
	if err != nil {
		panic(err)
	}

	http.HandleFunc("/", func(w http.ResponseWriter, r *http.Request) {

		err = t.Execute(w, snps)
		if err != nil {
			panic(err)
		}

	})
	fmt.Println("Откройте в интернет браузере адрес localhost:8080 или 127.0.0.1:8080 для просмотра результатов.\nДля выхода из программы, нажмите Ctl+C")
	http.ListenAndServe(":8080", nil)
	// for val := range tmp {
	// 	fmt.Printf("%v\n", tmp[val].Locus)
	// }

}

func getSNPNotation(f string) []gene.SNPcheck {
	var snpCheck gene.SNPcheck
	var parsedSNP []gene.SNPcheck
	rLMN := regexp.MustCompile(`^(\w+)\W+(\d+)(\D)>(\D)\s+(.*)`) // L:1 PiG:2 REF:3 ALT:4 NAME:5
	// LocusMutationName(LMN)
	rPMLN := regexp.MustCompile(`^(\d+)_(\D)>(\D)\W+(\w+)\W+(.*)`) // APOS:1 REF:2 ALT:3 L:4 NAME:5
	//PositionMutationLocusName (PMLN)
	rLSAAN := regexp.MustCompile(`^(\w+)\W+(\D{1})(\d+)(\D{1})\W(.*)`) // L:1 AA_REF:2 POS:3 AA_ALT:4 NAME:5
	// LocusShortAminoAcidName (LSAAN)
	rLLAAN := regexp.MustCompile(`^(\w+)\W+(\w{3})(\d+)(\w{3})\W+(.*)`) // L:1 LAA_REF:2 PiG:3 LAA_ALT:4 NAME:5
	// LocusLongAminoAcidName (LLAAN)
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
			snpCheck = gene.SNPcheck{Apos: matchPMLN[1], Ref: matchPMLN[2], Alt: matchPMLN[3], Locus: matchPMLN[4], Name: matchPMLN[5], TypeOf: tPMLN, Raw: matchPMLN[0]}
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

	}
	// fmt.Printf("%v\n", snpCheck)
	// for _, val := range parsedSNP {
	// 	fmt.Printf("%v\n", val)
	// }
	return parsedSNP
}

func getUniqSNP() {
	var countSNPs = 1
	var pos = make(map[int]int)
	var alt = make(map[int]gene.SNPinfo)
	var uniq []int
	// var g gene.Gene
	var snpToWeb []gene.SNPinfo
	var snpToConsole []gene.SNPinfo
	files := getListofVCF()
	for cnt, file := range files {
		fmt.Printf("processed %v files\r", cnt+1)
		snps := parserVCF(file, false, allGenesVal)

		for _, val := range snps {
			pos[val.Apos] = pos[val.Apos] + 1
			alt[val.Apos] = val
			// 	g = gene.Gene{Locus: val.Locus, Start: val.Start, End: val.End, Name: val.Name, Product: val.Product, Direction: val.Direction, GeneID: val.GeneID, ProteinID: val.ProteinID}
			// }
		}

		for i, lPos := range pos {
			if lPos == len(files) {
				uniq = append(uniq, i)
			}
		}
		sort.Ints(uniq)
		for _, uniqPos := range uniq {
			countSNPs++
			// fmt.Println(alt[uniqPos])

			if *flgWeb == true {
				snpToWeb = append(snpToWeb, alt[uniqPos])
				// printResults(alt[uniqPos])

			} else {
				// printResults(alt[uniqPos])
				snpToConsole = append(snpToConsole, alt[uniqPos])
			}
		}
	}
	fmt.Println()
	if *flgWeb == true && len(snpToWeb) != 0 {
		printWebResults(snpToWeb)
	} else if *flgWeb == false && len(snpToConsole) != 0 {
		for _, res := range snpToConsole {
			printResults(res)
		}
	}
	if *flgDebug == true {
		fmt.Printf("files:%v snps: %v\n%v\n", len(files), countSNPs, files)
	}
}

// func getUniqSNP() {
// 	var countSNPs int
// 	var snpToWeb []gene.SNPinfo
// 	var uniq []uniqInfo
// 	pos := make(map[int]int)
// 	// files, err := filepath.Glob(vcfExt)
// 	files := getListofVCF()
// 	// if err != nil {
// 	// 	log.Fatal(err)
// 	// }
//
// 	for _, file := range files {
//
// 		snps := parserVCF(file, false, allGenesVal)
//
// 		for _, val := range snps {
//
// 			pos[val.Apos] = pos[val.Apos] + 1
// 			uniq = append(uniq, uniqInfo{Apos: val.Apos, Count: pos[val.Apos], Alt: val.Alt})
//
// 			if pos[val.Apos] == len(files) && *flgWeb == true {
// 				snpToWeb = append(snpToWeb, val)
//
// 			}
//
// 		}
// 	}
//
// 	sort.Slice(uniq, func(i int, j int) bool {
// 		return uniq[i].Apos < uniq[j].Apos
// 	})
// 	// spew.Dump(uniq)
// 	sort.Slice(snpToWeb, func(i int, j int) bool {
// 		return snpToWeb[i].Apos < snpToWeb[j].Apos
// 	})
//
// 	for _, g := range allGenesVal {
//
// 		lStart, _ := strconv.Atoi(g.Start)
// 		lEnd, _ := strconv.Atoi(g.End)
// 		for _, val := range uniq {
// 			if val.Count == len(files) {
//
// 				if val.Apos >= lStart && val.Apos <= lEnd {
// 					countSNPs++
// 					snp := getSNPInfo(val.Apos, g, val.Alt)
//
// 					if *flgWeb == true {
//
// 						printWebResults(snpToWeb)
// 					} else {
//
// 						printResults(snp)
// 					}
// 				}
//
// 			}
// 		}
//
// 	}
//
// 	if *flgDebug == true {
// 		fmt.Printf("files:%v snps: %v\n%v\n", len(files), countSNPs, files)
// 	}
// }

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

	})
	fmt.Println("Откройте в интернет браузере адрес localhost:8080 или 127.0.0.1:8080 для просмотра результатов.\nДля выхода из программы, нажмите Ctl+C")
	http.ListenAndServe(":8080", nil)

}

func getListofVCF() []string {
	/*
		возвращает список VCF файлов в папке в виде массива
	*/
	files, err := filepath.Glob(vcfExt)

	if err != nil {
		log.Fatal(err)
	}
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
	var i int

	// flgFirstElement bool
	for _, file := range files {
		snps := parserVCF(file, false, allGenesVal)
		mapOfVCF[file] = snps
		i++
		fmt.Printf("processed %v files\r", i)
	}
	for file, snp := range mapOfVCF {
		var buffer bytes.Buffer

		buffer.WriteString(fmt.Sprintf("%v", file))
		for _, val := range locSNPcheck {
			lGpoS, _ := strconv.Atoi(val.PosInGene)

			for _, snpFromFile := range snp {
				switch val.TypeOf {
				case tLMN:
					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && strings.ToUpper(val.Alt) == strings.ToUpper(snpFromFile.Alt) && lGpoS == snpFromFile.Gpos {
						// fmt.Println(val.Locus, "\t", val.PosInGene, "\t", snpFromFile.Gpos, "\t", snpFromFile.Locus)
						buffer.WriteString(fmt.Sprintf(",%v[%v_%v%v>%v]", val.Name, val.Locus, val.PosInGene, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))
						fmt.Println(strings.ToUpper(buffer.String()))
					}
				case tPMLN:
					lApos, _ := strconv.Atoi(val.Apos)
					if lApos == snpFromFile.Apos && strings.ToUpper(val.Alt) == strings.ToUpper(snpFromFile.Alt) {
						buffer.WriteString(fmt.Sprintf(",%v[%v_%v>%v]", snpFromFile.Locus, lApos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))
						fmt.Println(strings.ToUpper(buffer.String()))
					}
				case tLSAAN:
					// lApos, _ := strconv.Atoi(val.Apos)
					var lCodonNbrInG int

					CodonNbrInG, _ := strconv.Atoi(val.CodonNbrInG)
					if snpFromFile.Direction == "r" {
						lCodonNbrInG = snpFromFile.CCodonNbrInG
					} else {
						lCodonNbrInG = snpFromFile.CodonNbrInG
					}

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && CodonNbrInG == lCodonNbrInG {
						fmt.Printf("%v%v%v %v %v%v%v->%v\t", val.AASref, CodonNbrInG, val.AASalt, val.Locus, snpFromFile.RefAAShort, lCodonNbrInG, snpFromFile.AltAAShort, snpFromFile.Alt)
						// spew.Dump(snpFromFile)
					}
					// if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && lGpoS == snpFromFile.Gpos && val.AASalt == snpFromFile.AltAAShort {
					// 	buffer.WriteString(fmt.Sprintf(",%v[%v%v%v]", snpFromFile.Locus, snpFromFile.RefAAShort, lGpoS, snpFromFile.AltAAShort))
					// fmt.Println(strings.ToUpper(buffer.String()))
					// }
				}
			}

			// fmt.Println(strings.ToUpper(buffer.String()))
		}

	}
}
