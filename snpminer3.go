package main

import (
	"bufio"
	"bytes"
	// "compress/gzip"
	"fmt"
	"io/ioutil"
	// "encoding/json"
	"log"
	// "github.com/CrowdSurge/banner"
	// "io"
	// "encoding/binary"
	"encoding/gob"
	// "net/http"
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
	"./amino"
	"./gene"
	"./tang"
	"compress/lzw"
	"io"
	// "sort"
	"strconv"
	// "text/template"
)

var linesFromGB []string
var genome_anchors = map[int]string{}
var TotalCDSFound int
var AllGenesVal []gene.Gene

var flgMkDB = flag.Bool("mkdb", false, "Create database from genbank file")
var flgGBfile = flag.String("gb", "", "Genbank file")
var flgAbout = flag.Bool("v", false, "About")
var flgOut = flag.String("out", "", "database filename")
var flgDB = flag.String("db", "", "database file")
var flgVCF = flag.String("vcf", "", "vcf file")
var flgPos = flag.Int("pos", 0, "Position of SNP")
var flgDebug = flag.Bool("debug", false, "")
var flgGenomeMap = flag.Bool("gmap", false, "")
var flgMakeNCSeq = flag.Bool("mkseq", false, "")
var flgTang = flag.Bool("tang", false, "")
var flgDevMode = flag.Bool("dev", false, "")
var flgReportType = flag.String("rep", "", "")
var flgSNP = flag.String("snp", "", "flag snp")

// var GenomeAccessNbr string

// var GenomeSeq string
var GenomeSeqSlice []string

const PrintResult = `
Locus {{.Locus}} 
Abs. position {{.Apos}}\t{{.CPosInGene}}{{.NucInPos}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CCodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}{{.Product}}\n`

// snp.Locus, snp.Apos, snp.CPosInGene, snp.NucInPos, snp.Alt, snp.RefCodon, snp.AltCodon,
// 									snp.RefAAShort, snp.CCodonNbrInG, snp.AltAAShort, snp.Mutation,
// 									snp.Tang, snp.Product

// %v\t%v\t%v%v>%v\t%v/%v\t%v%v%v\t%v\t%v\t%v\n

// type SNPinfo struct {
// 	Apos, Gpos, PosInCodonG, CodonNbrInG, CPosInGene, CCodonNbrInG int
// 	RefCodon, AltCodon, RefAA, AltAA, Locus,
// 	Direction, NucInPos, Product, Name, Start,
// 	RefAAShort, AltAAShort, End, Mutation, Tang, Alt, Note string
// }

// type ResValues struct {
//  Locus, apos, snp.CPosInGene, snp.NucInPos, snp.Alt, snp.RefCodon, snp.AltCodon,
// 									snp.RefAAShort, snp.CCodonNbrInG, snp.AltAAShort, snp.Mutation,
// 									snp.Tang, snp.Product
// }

func main() {

	// flgListOfFiles := flag.Bool("list", false, "")

	// fmt.Println("gb", *flgGBfile)
	// if *flgMkDB == true {
	// if *flgGBfile != "" {
	flag.Parse()

	// if *flgMkDB == true {
	// fmt.Println(*flgPos)
	// }
	// fmt.Printf("%v %v %v \n", *flgMkDB, *flgGBfile, *flgOut)
	if *flgAbout == true {
		About()
	}

	if *flgMkDB == true {
		if *flgGBfile != "" && *flgOut != "" {
			go process("Идет обработка...")
			AllGenesVal, GenomeSeqSlice = ReadGBFile(*flgGBfile)
			WriteDB(*flgOut, AllGenesVal)
			// fmt.Printf("%v\n", GenomeSeqSlice)
		} else {
			fmt.Println("Не указаны необходимые аргументы (-mkdb -gb=GENBANK_FILE -out=OUTPUT_DB_FILE)")
		}
	}

	if *flgDB != "" {
		// g := Gene

		AllGenesVal = ReadDB(*flgDB)

		if *flgMakeNCSeq != false && *flgOut != "" {

			fmt.Println("Reserved!")
		}

		if *flgVCF != "" && *flgVCF != "list" {
			ParserVCF(*flgVCF, AllGenesVal)
		} else if *flgVCF != "" && *flgVCF == "list" {
			ListOfVCFFiles()
		}

		if *flgSNP != "" {
			locSNPcheck := ParseSNP(*flgSNP)
			for _, val := range locSNPcheck {
				fmt.Printf("%v %v\n", val.Locus, val.TypeOf)
			}

		}

		if *flgGenomeMap != false {
			var igensS []string
			var igensE []string
			for i, g := range AllGenesVal {

				igensS = append(igensS, g.Start)
				igensE = append(igensE, g.End)
				// igenStart, _ := strconv.Atoi(g.End)
				if i >= 1 && i < len(AllGenesVal)-1 {
					// 	// igenEnd, _ := strconv.Atoi(g.Start)
					fmt.Printf("igen:s%v e:%v\n", igensS[i], igensE[i-1])
				}
				fmt.Printf("l:%v s:%v e:%v d:%v %v\n", g.Locus, g.Start, g.End, g.Direction, i)

			}
		}
		// var tmp_slited []string
		for _, g := range AllGenesVal {
			// fmt.Printf("%v %v \n", i.Locus, GetNucFromGenome(GenomeSeqSlice, i.Start, i.End))
			if *flgPos != 0 {

				lStart, _ := strconv.Atoi(g.Start)
				lEnd, _ := strconv.Atoi(g.End)

				if *flgPos >= lStart && *flgPos <= lEnd {
					// fmt.Printf("s:%v e:%v L:%v D:%v P:%v\n", g.Start, g.End, g.Locus, g.Direction, g.Product)
					if *flgDebug == true {
						// fmt.Printf(GetSNPInfo(*flgPos, g))
						snp := GetSNPInfo(*flgPos, g, "A")
						fmt.Printf("apos:%v, gpos:%v, pic:%v, AA:%v, ref:%v, nuc:%v, dir:%v, aa:%v\nL:%v, S:%v, E:%v, P:%v\n",
							snp.Apos, snp.Gpos, snp.PosInCodonG, snp.CodonNbrInG, snp.RefCodon, snp.NucInPos, g.Direction,
							snp.RefAA, snp.Locus, snp.Start, snp.End, snp.Product)

					}

				}

			}
		}

	}

}

func GetSNPInfo(apos int, g gene.Gene, alt string) gene.SNPinfo {
	var snp gene.SNPinfo             // структура SNP
	var codon_positions []string     // срез для разбивки кодона побуквенно
	var alt_codon_positions []string // срез для разбивки кодона побуквенно альтернативным нуклеотидом
	var locReportType string
	lStart, _ := strconv.Atoi(g.Start) // переменная начала гена
	lEnd, _ := strconv.Atoi(g.End)
	posInGene := ((apos - lStart) + 1)           // позиция снипа в гене
	codonNbrInG := ((posInGene - 1) / 3) + 1     // номер кодона=номеру аминокислоты в трансляции
	posInCodonG := (codonNbrInG * 3) - posInGene // позиция в буквы в кодоне (0-первая, 1-средняя, 2-последняя)
	CPosInGene := ((lEnd - apos) + 1)
	CCodonNbrInG := ((CPosInGene - 1) / 3) + 1
	// финт, который делал в snpMiner2, сдвиг на 1 букву. взял оттуда

	if posInCodonG == 2 {
		posInCodonG = 0
	} else if posInCodonG == 0 {
		posInCodonG = 2
	}

	var codon string    // переменная для кодона
	var altCodon string // аналогично для альтернативного кодона
	var mut string
	var tang_idx string
	/*
		определяем границы кодона, в зависимости от положения нуклеотида
		0- три буквы справа
		1-одна справа, одна слева
		2-три буквы слева
	*/

	if posInCodonG == 0 {
		codon = GetNucFromGenome((posInGene+lStart)-2, ((posInGene+lStart)-1)+2)
	} else if posInCodonG == 1 {
		codon = GetNucFromGenome((posInGene+lStart)-3, ((posInGene+lStart)-1)+1)
	} else if posInCodonG == 2 {
		codon = GetNucFromGenome((posInGene+lStart)-4, ((posInGene + lStart) - 1))
	}

	/*
		ревертируем кодон для комплементарных генов

	*/
	nucG := GetNucFromGenomePos((posInGene + lStart) - 1)
	if g.Direction == "r" {
		alt = CodonReverse(alt)
		nucG = CodonReverse(nucG)
		codon = CodonReverse(codon)
		// fmt.Printf("alt:%v codon:%v alt_r:%v rev_codon:%v\n", alt, codon, CodonReverse(alt), CodonReverse(codon))
		if posInCodonG == 2 {
			posInCodonG = 0
		} else if posInCodonG == 0 {
			posInCodonG = 2
		}

	}
	//

	codon_positions = strings.Split(codon, "")
	codon_positions[posInCodonG] = strings.ToUpper(codon_positions[posInCodonG])
	codon = strings.Join(codon_positions, "")
	alt_codon_positions = codon_positions
	alt_codon_positions[posInCodonG] = alt
	alt_codon_positions[posInCodonG] = strings.ToUpper(alt_codon_positions[posInCodonG])
	altCodon = strings.Join(alt_codon_positions, "")
	aaRef, aaRefShort := amino.Codon2AA(codon)
	aaAlt, aaAltShort := amino.Codon2AA(altCodon)
	if aaRefShort == aaAltShort {
		mut = "synonymous"
	} else if aaRefShort != aaAltShort && aaAltShort != "X" {
		mut = "missense"
		tang_idx = tang.GetTangInx(aaRefShort, aaAltShort)
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
		Mutation: mut, CCodonNbrInG: CCodonNbrInG, CPosInGene: CPosInGene, Tang: tang_idx, Alt: alt,
		Note: g.Note, ReportType: locReportType}
	// fmt.Printf("%v\t%v\n", aaRefShort, aaAltShort)
	// fmt.Printf("pig:%v, pic:%v, codon:%v %v%v\n", posInGene, posInCodonG, codon, aaRef, codonNbrInG)
	return snp
}

func GetNucFromGenome(start int, end int) string {
	/*
		возвращает последовательность нуклеотидов из генома:

	*/
	var result string
	var slice []string
	slice = GenomeSeqSlice[start:end]
	result = strings.Join(slice, "")
	return result

}

func GetNucFromGenomePos(pos int) string {
	/*
		возвращает последовательность нуклеотидов из генома:

	*/
	var result string
	var slice []string
	slice = GenomeSeqSlice[pos-1 : pos]
	result = strings.Join(slice, "")
	return result

}

func ReadGBFile(file string) (g []gene.Gene, genomeSplice []string) {
	//функция для считывания файла в формате генбанка и занесения строк в массив linesFromGB
	// (?m:).*?\s+(.*)\s+(.*)  (?m:)(.*)\s+
	// var regexpAnchors = regexp.MustCompile(`(?ms)\s+(.*)`)
	var CDS = regexp.MustCompile(`CDS\s+(\d+)\W+(\d+)`)
	var CDScompl = regexp.MustCompile(`CDS\s+complement\W(\d+)\W+(\d+)`)
	var locus = regexp.MustCompile(`;locus_tag=\W(.*?)\W;`)
	var product = regexp.MustCompile(`;product=\W(.*?)\W;`)
	var gene_name = regexp.MustCompile(`;gene=\W(.*?)\W;`)
	// var gseq = regexp.MustCompile(`^\W+(\d+)\s+(\w{6}.*)`)
	var qenomeSeq = regexp.MustCompile(`\s+(\d+)\s+(\w{6}.*)`)
	var geneID = regexp.MustCompile(`;db_xref=\WGeneID:(\d+)`)
	var proteinID = regexp.MustCompile(`;protein_id=\W(.*?)\W;`)
	var Note = regexp.MustCompile(`;note=\W(.*?)\W;`)
	var res_string []string
	var query_splice []string
	var splited_genome []string
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
			for _, genome_match := range qenomeSeq.FindAllStringSubmatch(string(line), -1) {
				// fmt.Printf("seq:%v\n", genome_match[1])
				res_string = append(res_string, strings.Replace(genome_match[2], " ", "", -1))
			}
		}

		if CDSblock == 1 {
			changedStr := strings.TrimSuffix(strings.Replace(string(line), "  ", "", -1), "\n")
			buffer.WriteString(strings.Replace(changedStr, " /", ";", -1))
			// fmt.Printf("%v\n", buffer.String())
		} else if CDSblock == 0 && buffer.Len() != 0 {
			// fmt.Printf("b:%v\n", strings.TrimSuffix(, "\n\n"))
			query_splice = append(query_splice, buffer.String())
			buffer.Reset()
		}

	}
	splited_genome = strings.SplitAfter(strings.Join(res_string, ""), "")

	defer f.Close()

	for i, val := range query_splice {
		// fmt.Printf("\t\t\t%v\n", val)
		for _, cds_match := range CDS.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("s:%v\te:%v d:%v\n", cds_match[1], cds_match[2], "f")
			lStart, lEnd, lDir = cds_match[1], cds_match[2], "f"

		}
		for _, cdsc_match := range CDScompl.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("s:%v\te:%v d:%v\n", cdsc_match[1], cdsc_match[2], "r")
			lStart, lEnd, lDir = cdsc_match[1], cdsc_match[2], "r"

		}
		for _, locus_match := range locus.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("l:%v\n", locus_match[1])
			lLoc = locus_match[1]
			// if lLoc == "Rv0068" {
			// 	fmt.Println("!!!!!!\n")
			// }
		}
		for _, gname_match := range gene_name.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("g:%v\n", gname_match[1])
			lName = gname_match[1]
		}
		for _, product_match := range product.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("p:%v\n", product_match[1])
			lProd = product_match[1]
		}

		for _, geneID_match := range geneID.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("p:%v\n", product_match[1])
			gID = geneID_match[1]
		}

		for _, protID_match := range proteinID.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("p:%v\n", product_match[1])
			pID = protID_match[1]
		}
		for _, note_match := range Note.FindAllStringSubmatch(val, -1) {
			// fmt.Printf("p:%v\n", product_match[1])
			lNote = note_match[1]
		}

		g := gene.Gene{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir,
			Product: lProd, Name: lName, GeneID: gID, ProteinID: pID, Note: lNote}
		// if *flgDebug == true {
		// fmt.Printf("%v\n", g)
		// }
		AllGenesVal = append(AllGenesVal, g)
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
	return AllGenesVal, splited_genome

}

func WriteDB(file string, gene []gene.Gene) {

	gobFile, err := os.Create(file)
	if err != nil {
		log.Println(err.Error())
	}
	defer gobFile.Close()
	compressedGobFile := lzw.NewWriter(gobFile, lzw.LSB, 8)
	defer compressedGobFile.Close()
	gobParser := gob.NewEncoder(compressedGobFile)
	gobParser.Encode(&gene)
	gobParser.Encode(&GenomeSeqSlice)
	fmt.Println("\nРезультаты сохранены в", file, "файл.")
	// taskFile.Write(jsonParser)

	// fmt.Println(strings.Repeat("\nБаза успешно создана!")
}

func ReadDB(file string) []gene.Gene {
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
	gobParser.Decode(&GenomeSeqSlice)
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

func CodonReverse(codon string) string {
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

func ParserVCF(f string, genes []gene.Gene) []gene.SNPinfo {
	var vcf = regexp.MustCompile(`(^\S+)\W(\d+)\W+(\w+)\W+(\w+)`)
	var snpFromVCF []gene.SNPinfo
	// x := make(map[string][]gene.SNPinfo)

	// var result string
	//^\S+\W(\d+)\W+(\w+)\W+(\w+)
	// var check_vcf = regexp.MustCompile(`##fileformat=(VCF)`)

	file, err := os.Open(f)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		// for _, match_check := range check_vcf.FindAllStringSubmatch(scanner.Text(), -1) {
		// 	fmt.Printf("%v\n", match_check[1])
		// }
		// // fmt.Println(scanner.Text())
		for _, match := range vcf.FindAllStringSubmatch(scanner.Text(), -1) {
			apos, _ := strconv.Atoi(match[2])
			ref := match[3]
			alt := match[4]

			// orgInfo = GenomeInfo{Organism: match[1]}
			for _, g := range genes {
				lStart, _ := strconv.Atoi(g.Start)
				lEnd, _ := strconv.Atoi(g.End)

				if apos >= lStart && apos <= lEnd {
					snp := GetSNPInfo(apos, g, alt)
					snpFromVCF = append(snpFromVCF, snp)

					if len(ref) == 1 && len(alt) == 1 {
						PrintResults(snp)
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

func About() {
	fmt.Println("snpMiner3 (alpha) 0.10.17\n(c)V.Sinkov & O.Ogarkov,Irkutsk, Russia, 2017")
}

func process(str string) {
	fmt.Println(str)
}

func ListOfVCFFiles() {
	files, err := ioutil.ReadDir(".")
	if err != nil {
		log.Fatal(err)
	}

	for _, file := range files {
		if strings.Contains(file.Name(), ".vcf") {
			fmt.Printf("\n\n%v:\n\n", file.Name())
			ParserVCF(file.Name(), AllGenesVal)

			// fmt.Println(file.Name())
		}
	}
}

func GetInterGen(pos int) {
	var igensS []string
	var igensE []string
	for i, g := range AllGenesVal {

		igensS = append(igensS, g.Start)
		igensE = append(igensE, g.End)
		// igenStart, _ := strconv.Atoi(g.End)
		if i >= 1 && i < len(AllGenesVal)-1 {
			// 	// igenEnd, _ := strconv.Atoi(g.Start)
			fmt.Printf("igen:s%v e:%v\n", igensS[i], igensE[i-1])
		}
		fmt.Printf("l:%v s:%v e:%v d:%v %v\n", g.Locus, g.Start, g.End, g.Direction, i)

	}
}

func PrintResults(snp gene.SNPinfo) {
	// resFlg
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

	// loc := string []gene.SNPinfo

	// res[locus] = append(res[locus], snp)
}

func ParseSNP(f string) []gene.SNPcheck {
	var snpCheck gene.SNPcheck
	var parsedSNP []gene.SNPcheck
	rLMN := regexp.MustCompile(`^(\w+)\t(\d+)(\w)>(\w)\t(\w+)`)
	// LocusMutationName(LMN)
	rPMLN := regexp.MustCompile(`^(\d+)_(\w)>(\w)\W+(\w+)\W+(\w+)`)
	//PositionMutationLocusName (PLMN)
	rLSAAN := regexp.MustCompile(`^(\w+)\W+(\w)(\d+)(\w)\W+(\w+)`)
	// LocusShortAminoAcidName (LSAAN)
	rLLAAN := regexp.MustCompile(`^(\w+)\W+(\w{3})(\d+)(\w{3})\W+(\w+)`)
	// LocusLongAminoAcidName (LLAAN)
	file, err := os.Open(f)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {

		for _, matchLMN := range rLMN.FindAllStringSubmatch(scanner.Text(), -1) {
			snpCheck = gene.SNPcheck{Locus: matchLMN[1], PosInGene: matchLMN[2], Ref: matchLMN[3], Alt: matchLMN[4], Name: matchLMN[5], TypeOf: "LMN"}
			parsedSNP = append(parsedSNP, snpCheck)
			// fmt.Printf("%v\n", snpCheck)
			// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
		}
		for _, matchPLMN := range rPMLN.FindAllStringSubmatch(scanner.Text(), -1) {
			snpCheck = gene.SNPcheck{Apos: matchPLMN[1], Ref: matchPLMN[2], Alt: matchPLMN[3], Locus: matchPLMN[4], Name: matchPLMN[5], TypeOf: "PLMN"}
			parsedSNP = append(parsedSNP, snpCheck)
			// fmt.Printf("%v\n", snpCheck)
			// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
		}
		for _, matchLSAAN := range rLSAAN.FindAllStringSubmatch(scanner.Text(), -1) {
			snpCheck = gene.SNPcheck{Locus: matchLSAAN[1], AASref: matchLSAAN[2], PosInGene: matchLSAAN[3], AASalt: matchLSAAN[4], Name: matchLSAAN[5], TypeOf: "LSAAN"}
			parsedSNP = append(parsedSNP, snpCheck)
			// fmt.Printf("%v\n", snpCheck)
			// fmt.Printf("%v %v %v %v %v\n", strings.ToUpper(match[1]), match[2], strings.ToUpper(match[3]), strings.ToUpper(match[4]), match[5])
		}
		for _, matchLLAAN := range rLLAAN.FindAllStringSubmatch(scanner.Text(), -1) {
			snpCheck = gene.SNPcheck{Locus: matchLLAAN[1], AALref: matchLLAAN[2], PosInGene: matchLLAAN[3], AALalt: matchLLAAN[4], Name: matchLLAAN[5], TypeOf: "LSAAN"}
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
