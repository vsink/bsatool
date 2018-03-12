package amino

import (
	// "fmt"
	"strings"
)

func Codon2AA(codon string) (string, string) {
	/*


	 */
	var Lname, Sname, codonRes string
	codonRes = strings.ToUpper(codon)
	AA := map[string]map[string]string{
		"TCA": {
			"LName": "Ser",
			"SName": "S",
			"CDA":   "0.5",
		},
		"TCC": {
			"LName": "Ser",
			"SName": "S",
			"CDA":   "-0.5",
		},
		"TCG": {
			"LName": "Ser",
			"SName": "S",
			"CDA":   "0.5",
		},
		"TCT": {
			"LName": "Ser",
			"SName": "S",
			"CDA":   "0",
		},
		"TTC": {
			"LName": "Phe",
			"SName": "F",
			"CDA":   "0.5",
		},
		"TTT": {
			"LName": "Phe",
			"SName": "F",
			"CDA":   "0",
		},
		"TTA": {
			"LName": "Leu",
			"SName": "L",
			"CDA":   "1",
		},
		"TTG": {
			"LName": "Leu",
			"SName": "L",
			"CDA":   "1",
		},
		"TAC": {
			"LName": "Tyr",
			"SName": "Y",
			"CDA":   "0.5",
		},
		"TAT": {
			"LName": "Tyr",
			"SName": "Y",
		},
		"TAA": {
			"LName": "X",
			"SName": "X",
		},
		"TAG": {
			"LName": "X",
			"SName": "X",
		},
		"TGC": {
			"LName": "Cys",
			"SName": "C",
		},
		"TGT": {
			"LName": "Cys",
			"SName": "C",
		},
		"TGA": {
			"LName": "X",
			"SName": "X",
		},
		"TGG": {
			"LName": "Trp",
			"SName": "W",
		},
		"CTA": {
			"LName": "Leu",
			"SName": "L",
		},
		"CTC": {
			"LName": "Leu",
			"SName": "L",
		},
		"CTG": {
			"LName": "Leu",
			"SName": "L",
		},
		"CTT": {
			"LName": "Leu",
			"SName": "L",
		},
		"CCA": {
			"LName": "Pro",
			"SName": "P",
		},
		"CCC": {
			"LName": "Pro",
			"SName": "P",
		},
		"CCG": {
			"LName": "Pro",
			"SName": "P",
		},
		"CCT": {
			"LName": "Pro",
			"SName": "P",
		},

		"CAC": {
			"LName": "His",
			"SName": "H",
		},
		"CAT": {
			"LName": "His",
			"SName": "H",
		},
		"CAA": {
			"LName": "Gln",
			"SName": "Q",
		},
		"CAG": {
			"LName": "Gln",
			"SName": "Q",
		},
		"CGA": {
			"LName": "Arg",
			"SName": "R",
		},
		"CGC": {
			"LName": "Arg",
			"SName": "R",
		},
		"CGG": {
			"LName": "Arg",
			"SName": "R",
		},
		"CGT": {
			"LName": "Arg",
			"SName": "R",
		},
		"ATA": {
			"LName": "Ile",
			"SName": "I",
		},
		"ATC": {
			"LName": "Ile",
			"SName": "I",
		},
		"ATT": {
			"LName": "Ile",
			"SName": "I",
		},
		"ATG": {
			"LName": "Met",
			"SName": "M",
		},
		"ACA": {
			"LName": "Thr",
			"SName": "T",
		},
		"ACC": {
			"LName": "Thr",
			"SName": "T",
		},
		"ACG": {
			"LName": "Thr",
			"SName": "T",
		},
		"ACT": {
			"LName": "Thr",
			"SName": "T",
		},
		"AAC": {
			"LName": "Asn",
			"SName": "N",
		},
		"AAT": {
			"LName": "Asn",
			"SName": "N",
		},
		"AAA": {
			"LName": "Lys",
			"SName": "K",
		},
		"AAG": {
			"LName": "Lys",
			"SName": "K",
		},
		"AGC": {
			"LName": "Ser",
			"SName": "S",
		},
		"AGT": {
			"LName": "Ser",
			"SName": "S",
		},
		"AGA": {
			"LName": "Arg",
			"SName": "R",
		},
		"AGG": {
			"LName": "Arg",
			"SName": "R",
		},
		"GTA": {
			"LName": "Val",
			"SName": "V",
		},
		"GTC": {
			"LName": "Val",
			"SName": "V",
		},
		"GTG": {
			"LName": "Val",
			"SName": "V",
		},
		"GTT": {
			"LName": "Val",
			"SName": "V",
		},
		"GCA": {
			"LName": "Ala",
			"SName": "A",
		},
		"GCC": {
			"LName": "Ala",
			"SName": "A",
		},
		"GCG": {
			"LName": "Ala",
			"SName": "A",
		},
		"GCT": {
			"LName": "Ala",
			"SName": "A",
		},
		"GAC": {
			"LName": "Asp",
			"SName": "D",
		},
		"GAT": {
			"LName": "Asp",
			"SName": "D",
		},
		"GAA": {
			"LName": "Glu",
			"SName": "E",
		},
		"GAG": {
			"LName": "Glu",
			"SName": "E",
		},
		"GGA": {
			"LName": "Gly",
			"SName": "G",
		},
		"GGC": {
			"LName": "Gly",
			"SName": "G",
		},
		"GGG": {
			"LName": "Gly",
			"SName": "G",
		},
		"GGT": {
			"LName": "Gly",
			"SName": "G",
		},
	}

	if aa, ok := AA[codonRes]; ok {
		// AA[codon]["Lname"], AA[codon]["SName"]
		Lname = aa["LName"]
		Sname = aa["SName"]
		// fmt.Println(Lname, Sname)

	}
	return Lname, Sname
}

func CDACodon(codon string) float32 {
	/*

	 */
	var codonRes string
	var cda float32
	codonRes = strings.ToUpper(codon)
	AA := map[string]map[string]float32{
		"TCA": {"CDA": 0.5},
		"TCC": {"CDA": -0.5},
		"TCG": {"CDA": 0.5},
		"TCT": {"CDA": 0},
		"TTC": {"CDA": 0.5},
		"TTT": {"CDA": 0},
		"TTA": {"CDA": 1},
		"TTG": {"CDA": 1},
		"TAC": {"CDA": 0.5},
		"TAT": {"CDA": 0},
		"TAA": {"CDA": -1},
		"TAG": {"CDA": -0.5},
		"TGC": {"CDA": 0.5},
		"TGT": {"CDA": 0},
		"TGA": {"CDA": -0.5},
		"TGG": {"CDA": -1},
		"CTA": {"CDA": 0.5},
		"CTC": {"CDA": 0},
		"CTG": {"CDA": 0.5},
		"CTT": {"CDA": -1},
		"CCA": {"CDA": 1},
		"CCC": {"CDA": 0},
		"CCG": {"CDA": 1},
		"CCT": {"CDA": 0.5},

		"CAC": {"CDA": 0},
		"CAT": {"CDA": 0.5},
		"CAA": {"CDA": -1},
		"CAG": {"CDA": 0.5},
		"CGA": {"CDA": -0.5},
		"CGC": {"CDA": 0},
		"CGG": {"CDA": -1},
		"CGT": {"CDA": -0.5},
		"ATA": {"CDA": 0},
		"ATC": {"CDA": -0.5},
		"ATT": {"CDA": -1},
		"ATG": {"CDA": -0.5},
		"ACA": {"CDA": 0},
		"ACC": {"CDA": -1},
		"ACG": {"CDA": 0.5},
		"ACT": {"CDA": -0.5},
		"AAC": {"CDA": -1},
		"AAT": {"CDA": 1},
		"AAA": {"CDA": 0},
		"AAG": {"CDA": 0.5},
		"AGC": {"CDA": 0.5},
		"AGT": {"CDA": 0.5},
		"AGA": {"CDA": 0},
		"AGG": {"CDA": -0.5},
		"GTA": {"CDA": 0.5},
		"GTC": {"CDA": -0.5},
		"GTG": {"CDA": 0},
		"GTT": {"CDA": -1},
		"GCA": {"CDA": -0.5},
		"GCC": {"CDA": -1},
		"GCG": {"CDA": 0},
		"GCT": {"CDA": -0.5},
		"GAC": {"CDA": 0.5},
		"GAT": {"CDA": 0.5},
		"GAA": {"CDA": -0.5},
		"GAG": {"CDA": 0},
		"GGA": {"CDA": 0.5},
		"GGC": {"CDA": 1},
		"GGG": {"CDA": 0},
		"GGT": {"CDA": 1},
	}

	if aa, ok := AA[codonRes]; ok {
		// AA[codon]["Lname"], AA[codon]["SName"]
		cda = aa["CDA"]
		// fmt.Println(Lname, Sname)

	}
	return cda
}

func CodonVolatility(codon string) float32 {
	var codonRes string
	var volty float32
	codonRes = strings.ToUpper(codon)
	AA := map[string]map[string]float32{
		"TTT": {"V": 0.889},
		"TTC": {"V": 0.889},
		"TTA": {"V": 0.714},
		"TTG": {"V": 0.750},
		"CTT": {"V": 0.667},
		"CTC": {"V": 0.667},
		"CTA": {"V": 0.556},
		"CTG": {"V": 0.556},
		"ATT": {"V": 0.778},
		"ATC": {"V": 0.778},
		"ATA": {"V": 0.778},
		"ATG": {"V": 1.000},
		"GTT": {"V": 0.667},
		"GTC": {"V": 0.667},
		"GTA": {"V": 0.667},
		"GTG": {"V": 0.667},
		"TCT": {"V": 0.667},
		"TCC": {"V": 0.667},
		"TCA": {"V": 0.571},
		"TCG": {"V": 0.625},
		"CCT": {"V": 0.667},
		"CCC": {"V": 0.667},
		"CCA": {"V": 0.667},
		"CCG": {"V": 0.667},
		"ACT": {"V": 0.667},
		"ACC": {"V": 0.667},
		"ACA": {"V": 0.667},
		"ACG": {"V": 0.667},
		"GCT": {"V": 0.667},
		"GCC": {"V": 0.667},
		"GCA": {"V": 0.667},
		"GCG": {"V": 0.667},
		"TAT": {"V": 0.857},
		"TAC": {"V": 0.857},
		"CAT": {"V": 0.889},
		"CAC": {"V": 0.889},
		"CAA": {"V": 0.875},
		"CAG": {"V": 0.875},
		"AAT": {"V": 0.889},
		"AAC": {"V": 0.889},
		"AAA": {"V": 0.875},
		"AAG": {"V": 0.875},
		"GAT": {"V": 0.889},
		"GAC": {"V": 0.889},
		"GAA": {"V": 0.875},
		"GAG": {"V": 0.875},
		"TGT": {"V": 0.875},
		"TGC": {"V": 0.875},
		"TGG": {"V": 1.000},
		"CGT": {"V": 0.667},
		"CGC": {"V": 0.667},
		"CGA": {"V": 0.500},
		"CGG": {"V": 0.556},
		"AGT": {"V": 0.889},
		"AGC": {"V": 0.889},
		"AGA": {"V": 0.750},
		"AGG": {"V": 0.778},
		"GGT": {"V": 0.667},
		"GGC": {"V": 0.667},
		"GGA": {"V": 0.625},
		"GGG": {"V": 0.667},
	}

	if aa, ok := AA[codonRes]; ok {
		// AA[codon]["Lname"], AA[codon]["SName"]
		volty = aa["V"]
		// fmt.Println(Lname, Sname)

	}
	return volty
}
