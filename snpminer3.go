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
	"./tang"
	"compress/lzw"
	"io"
	// "sort"
	"strconv"
)

var linesFromGB []string
var genome_anchors = map[int]string{}
var TotalCDSFound int
var AllGenesVal []Gene

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

// var GenomeAccessNbr string

// var GenomeSeq string
var GenomeSeqSlice []string

// type Gene struct {
// 	start, end, locus, name, direction string
// }

type Gene struct {
	Locus string `json:"locus"`
	// TypeOf    string `json:"typeof"`
	Start     string `json:"start"`
	End       string `json:"end"`
	Name      string `json:"name"`
	Product   string `json:"product"`
	Direction string `json:"direction"`
	GeneID    string `json:"gene_id"`
	ProteinID string `json:"protein_id"`
	Note      string `json:"note"`
}

type GenomeInfo struct {
	Lenght   string `json:"lenght"`
	Strain   string `json:"strain"`
	Organism string `json:"organism"`
	MolType  string `json:"moltype"`
}

type Options struct {
	MakeDB bool
}

type SNPinfo struct {
	Apos, Gpos, PosInCodonG, CodonNbrInG, cPosInGene, cCodonNbrInG int
	RefCodon, AltCodon, RefAA, AltAA, Locus,
	Direction, NucInPos, Product, Name, Start,
	RefAAShort, AltAAShort, End, Mutation, Tang, Alt string
}

type AA struct {
	LName, SName string
}

// type GenomeAnchors struct {
// 	anchore_pos  int
// 	anchore_type string
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

func GetSNPInfo(apos int, g Gene, alt string) SNPinfo {
	var snp SNPinfo                    // структура SNP
	var codon_positions []string       // срез для разбивки кодона побуквенно
	var alt_codon_positions []string   // срез для разбивки кодона побуквенно альтернативным нуклеотидом
	lStart, _ := strconv.Atoi(g.Start) // переменная начала гена
	lEnd, _ := strconv.Atoi(g.End)
	posInGene := ((apos - lStart) + 1)           // позиция снипа в гене
	codonNbrInG := ((posInGene - 1) / 3) + 1     // номер кодона=номеру аминокислоты в трансляции
	posInCodonG := (codonNbrInG * 3) - posInGene // позиция в буквы в кодоне (0-первая, 1-средняя, 2-последняя)
	cPosInGene := ((lEnd - apos) + 1)
	cCodonNbrInG := ((cPosInGene - 1) / 3) + 1
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
	tang.GetTangInx(aaRefShort, aaAltShort)
	snp = SNPinfo{Apos: apos, Gpos: posInGene, PosInCodonG: posInCodonG,
		RefCodon: codon, RefAA: aaRef, NucInPos: strings.ToUpper(nucG), Locus: g.Locus,
		Direction: g.Direction, Name: g.Name, Product: g.Product,
		Start: g.Start, End: g.End, CodonNbrInG: codonNbrInG, AltCodon: altCodon,
		AltAA: aaAlt, RefAAShort: aaRefShort, AltAAShort: aaAltShort,
		Mutation: mut, cCodonNbrInG: cCodonNbrInG, cPosInGene: cPosInGene, Tang: tang_idx, Alt: alt}
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

func ReadGBFile(file string) (g []Gene, genomeSplice []string) {
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

		g := Gene{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir,
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

// func MakeGenomeSeq() []string {
// 	// Функция создает полную последовательность нуклеотидов генома
// 	var re = regexp.MustCompile(`^\W+(\d+)\s+(\w{6}.*)`)
// 	// флаг определяющий структуру : 123455   aaaaagggggggg gggggggccccccc
// 	var IsStrNC bool
// 	// массив со строкми, после фильтрации
// 	var res_string []string
// 	var splited_genome []string
// 	// поиск подстрок,  удаление нумерации, пробелов и сохранение в массиве
// 	for _, gb_line := range linesFromGB {
// 		IsStrNC = re.MatchString(gb_line)
// 		if IsStrNC == true {
// 			// NCSeq = append(NCSeq, sline)
// 			for _, match := range re.FindAllStringSubmatch(gb_line, -1) {

// 				res_string = append(res_string, strings.Replace(match[2], " ", "", -1))

// 			}

// 		}
// 	}
// 	splited_genome = strings.SplitAfter(strings.Join(res_string, ""), "")
// 	return splited_genome
// 	// fmt.Println(getNucFromGenome(strings.Join(res_string, ""), 1, 1524))

// }

// func GetGenomeAnchors() {
// 	/*
// 		Ищет в каждой строке ключевые фразы и заносит значения в буфер
// 		(?m:).*?\s+(.*)\s+(.*)

// 	*/
// 	var ii int
// 	for i, val := range linesFromGB {

// 		if strings.Contains(val, "    CDS   ") ||
// 			// strings.Contains(val, "    gene   ") ||
// 			strings.Contains(val, "    tRNA   ") ||
// 			// strings.Contains(val, "    misc_feature   ") ||
// 			// strings.Contains(val, "    mobile_element   ") ||
// 			strings.Contains(val, "DEFINITION") ||
// 			strings.Contains(val, "ACCESSION") ||
// 			strings.Contains(val, "ACCESSION") ||
// 			strings.Contains(val, "SOURCE") {
// 			// fmt.Println(i+1, "\tCDS")
// 			count := i + 1
// 			genome_anchors[count] = val

// 		}
// 		if strings.Contains(val, "       /") {
// 			// fmt.Println(i+1, "\tCDS")
// 			count := i + 1
// 			genome_anchors[count] = strings.Replace(val, " /", "", -1)
// 		}

// 		if strings.Contains(val, " CDS ") {
// 			ii++
// 			TotalCDSFound = ii
// 		}

// 	}
// }

// func Gb2DB() {

// 	/*
// 		Удаляем лишние значения

// 	*/

// 	keys := []int{}
// 	for key := range genome_anchors {
// 		if strings.Contains(genome_anchors[key], "db_xref=") {
// 			delete(genome_anchors, key)
// 		} else if strings.Contains(genome_anchors[key], "experiment=") {
// 			delete(genome_anchors, key)
// 		} else if strings.Contains(genome_anchors[key], "translation=") {
// 			delete(genome_anchors, key)
// 		} else if strings.Contains(genome_anchors[key], "inference=") {
// 			delete(genome_anchors, key)
// 		} else if strings.Contains(genome_anchors[key], "transl_table=") {
// 			delete(genome_anchors, key)
// 		} else if strings.Contains(genome_anchors[key], "codon_start=") {
// 			delete(genome_anchors, key)
// 		} else if strings.Contains(genome_anchors[key], "protein_id=") {
// 			delete(genome_anchors, key)
// 		} else if strings.Contains(genome_anchors[key], "ACCESSION") {
// 			// fmt.Println(genome_anchors[key])
// 		}

// 	}

// 	for key := range genome_anchors {

// 		keys = append(keys, key)
// 		// fmt.Printf("%v\n", keys[key])
// 	}
// 	sort.Ints(keys)

// 	// for _, line := range linesFromGB {
// 	var rStartEnd = regexp.MustCompile(`^\s+(\w+)\s+(\d+)..(\d+)`)
// 	var rStartEndc = regexp.MustCompile(`^\s+(\w+)\s+complement\W(\d+)\W+(\d+)`)
// 	var rLocTag = regexp.MustCompile(`locus_tag=\W(\w+)`)
// 	var rGeneName = regexp.MustCompile(`gene=\W(\w+)`)
// 	var rProduct = regexp.MustCompile(`product=\W(.*\W+.*)`)
// 	//product=\W(.*\W+.*) product=\W(.*)\W
// 	//

// 	//
// 	var rOrgName = regexp.MustCompile(`organism=\W(.*)\W`)
// 	var rOrgStrain = regexp.MustCompile(`strain=\W(.*)\W`)
// 	var rOrgMolType = regexp.MustCompile(`mol_type=\W(.*)\W`)
// 	var rOrgLenght = regexp.MustCompile(`source\W+\d+\W+(\d+)`)

// 	var lStart, lEnd, lLoc, lProd, lDir, lName string
// 	// var lGOrg, lGStrain, lGMolType,

// 	// var rNote = regexp.MustCompile(`\/note=(".*)`)
// 	// var st, en string
// 	// var g []string

// 	for _, val := range keys {

// 		// fmt.Printf("%v\n", genome_anchors[val])

// 		for _, match := range rOrgName.FindAllStringSubmatch(genome_anchors[val], -1) {

// 			// orgInfo = GenomeInfo{Organism: match[1]}
// 			fmt.Printf("Organism: %v\nTotal CDS:%v \n", match[1], TotalCDSFound)

// 		}
// 		for _, match := range rOrgStrain.FindAllStringSubmatch(genome_anchors[val], -1) {

// 			// orgInfo = GenomeInfo{Strain: match[1]}
// 			fmt.Printf("Strain: %v\n", match[1])
// 		}
// 		for _, match := range rOrgMolType.FindAllStringSubmatch(genome_anchors[val], -1) {

// 			fmt.Printf("Mol.type: %v\n", match[1])
// 			// orgInfo = GenomeInfo{MolType: match[1]}

// 		}
// 		for _, match := range rOrgLenght.FindAllStringSubmatch(genome_anchors[val], -1) {

// 			// orgInfo = GenomeInfo{Lenght: match[1]}
// 			fmt.Printf("Lenght: %v nc\n", match[1])

// 		}

// 		for _, match := range rStartEnd.FindAllStringSubmatch(genome_anchors[val], -1) {

// 			lStart, lEnd, lDir = match[2], match[3], "f"
// 			// fmt.Printf("%v %v\n", lStart, lEnd)

// 		}
// 		for _, match_c := range rStartEndc.FindAllStringSubmatch(genome_anchors[val], -1) {

// 			lStart, lEnd, lDir = match_c[2], match_c[3], "r"
// 			// fmt.Printf("%v %v\n", lStart, lEnd)

// 		}
// 		for _, matchLocTag := range rLocTag.FindAllStringSubmatch(genome_anchors[val], -1) {

// 			lLoc = matchLocTag[1]
// 			// fmt.Printf("%v\n", lLoc)
// 		}
// 		for _, matchGeneName := range rGeneName.FindAllStringSubmatch(genome_anchors[val], -1) {

// 			lName = matchGeneName[1]
// 		}
// 		for _, matchProduct := range rProduct.FindAllStringSubmatch(genome_anchors[val], -1) {

// 			lProd = matchProduct[1]
// 			if *flgDebug == true {
// 				fmt.Printf("%v\n", genome_anchors[val])
// 			}
// 		}

// 		if lStart != "" && lEnd != "" && lLoc != "" && lProd != "" {

// 			// lSeq := GetNucFromGenome(GenomeSeqSlice, lStart, lEnd)
// 			g := Gene{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir, Product: lProd, Name: lName}
// 			if *flgDebug == true {
// 				fmt.Printf("%v\n", g)
// 			}
// 			AllGenesVal = append(AllGenesVal, g)
// 			// fmt.Printf("l:%v s:%v e:%v p:%v %v \n", lLoc, lStart, lEnd, lProd, val)
// 			// обнуление
// 			lName, lStart, lEnd, lLoc, lDir, lProd = "", "", "", "", "", ""

// 		}
// 	}
// 	// fmt.Printf("%v\n", AllGenesVal)

// }

func WriteDB(file string, gene []Gene) {

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

func ReadDB(file string) []Gene {
	var gene []Gene
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

// aa["TCA"] = "Ser"
// aa["TCC"] = "Ser"
// aa["TCG"] = "Ser"
// aa["TCT"] = "Ser"

// aa["TTC"] = "Phe"
// aa["TTT"] = "Phe"

// aa["TTA"] = "Leu"
// aa["TTG"] = "Leu"

// aa["TAC"] = "Tyr"
// aa["TAT"] = "Tyr"

// aa["TAA"] = "STOP"
// aa["TAG"] = "STOP"

// aa["TGC"] = "Cys"
// aa["TGT"] = "Cys"

// aa["TGA"] = "STOP"

// aa["TGG"] = "Trp"

// aa["CTA"] = "Leu"
// aa["CTC"] = "Leu"
// aa["CTG"] = "Leu"
// aa["CTT"] = "Leu"

// aa["CCA"] = "Pro"
// aa["CCC"] = "Pro"
// aa["CCG"] = "Pro"
// aa["CCT"] = "Pro"

// aa["CAC"] = "His"
// aa["CAT"] = "His"

// aa["CAA"] = "Gln"
// aa["CAG"] = "Gln"

// aa["CGA"] = "Arg"
// aa["CGC"] = "Arg"
// aa["CGG"] = "Arg"
// aa["CGT"] = "Arg"

// aa["ATA"] = "Ile"
// aa["ATC"] = "Ile"
// aa["ATT"] = "Ile"

// aa["ATG"] = "Met"

// aa["ACA"] = "Thr"
// aa["ACC"] = "Thr"
// aa["ACG"] = "Thr"
// aa["ACT"] = "Thr"

// aa["AAC"] = "Asn"
// aa["AAT"] = "Asn"

// aa["AAA"] = "Lys"
// aa["AAG"] = "Lys"

// aa["AGC"] = "Ser"
// aa["AGT"] = "Ser"

// aa["AGA"] = "Arg"
// aa["AGG"] = "Arg"

// aa["GTA"] = "Val"
// aa["GTC"] = "Val"
// aa["GTG"] = "Val"
// aa["GTT"] = "Val"

// aa["GCA"] = "Ala"
// aa["GCC"] = "Ala"
// aa["GCG"] = "Ala"
// aa["GCT"] = "Ala"

// aa["GAC"] = "Asp"
// aa["GAT"] = "Asp"

// aa["GAA"] = "Glu"
// aa["GAG"] = "Glu"

// aa["GGA"] = "Gly"
// aa["GGC"] = "Gly"
// aa["GGG"] = "Gly"
// aa["GGT"] = "Gly"
// return aa[strings.ToUpper(codon)]

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

func ParserVCF(f string, genes []Gene) {
	var vcf = regexp.MustCompile(`(^\S+)\W(\d+)\W+(\w+)\W+(\w+)`)
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
					// comp := strings.ToUpper(snp.NucInPos) == strings.ToUpper(ref)
					// fmt.Printf("apos:%v, gpos:%v, pic:%v, AA:%v, ref:%v, nuc:%v, dir:%v, aa:%v L:%v, S:%v, E:%v, P:%v\n",
					// snp.Apos, snp.Gpos, snp.PosInCodonG, snp.CodonNbrInG, snp.RefCodon, snp.NucInPos, g.Direction,
					// snp.RefAA, snp.Locus, snp.Start, snp.End, snp.Product)
					switch *flgTang {
					case true:
						if len(ref) == 1 && len(alt) == 1 {
							if g.Direction == "r" {
								fmt.Printf("%v\t%v\t%v%v>%v\t%v/%v\t%v%v%v\t%v\t%v\t%v\n", g.Locus, apos, snp.cPosInGene, snp.NucInPos, snp.Alt, snp.RefCodon, snp.AltCodon,
									snp.RefAAShort, snp.cCodonNbrInG, snp.AltAAShort, snp.Mutation,
									snp.Tang, snp.Product)
							} else {
								fmt.Printf("%v\t%v\t%v%v>%v\t%v/%v\t%v%v%v\t%v\t%v\t%v\n", g.Locus, apos, snp.Gpos, snp.NucInPos, snp.Alt, snp.RefCodon, snp.AltCodon,
									snp.RefAAShort, snp.CodonNbrInG, snp.AltAAShort, snp.Mutation,
									snp.Tang, snp.Product)
							}

						}
					case false:
						if len(ref) == 1 && len(alt) == 1 {
							if g.Direction == "r" {
								fmt.Printf("%v\t%v\t%v%v>%v\t%v/%v\t%v%v%v\t%v\t%v\n", g.Locus, apos, snp.cPosInGene, snp.NucInPos, snp.Alt, snp.RefCodon, snp.AltCodon,
									snp.RefAAShort, snp.cCodonNbrInG, snp.AltAAShort, snp.Mutation,
									snp.Product)
							} else {
								fmt.Printf("%v\t%v\t%v%v>%v\t%v/%v\t%v%v%v\t%v\t%v\n", g.Locus, apos, snp.Gpos, snp.NucInPos, snp.Alt, snp.RefCodon, snp.AltCodon,
									snp.RefAAShort, snp.CodonNbrInG, snp.AltAAShort, snp.Mutation,
									snp.Product)
							}

						}
					}
				}

				// fmt.Printf("acc:%v pos:%v ref:%v alt:%v \n", match[1], match[2], match[3], match[4])
			}
		}

	}

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}
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
