package amino

import (
	// "fmt"
	"strings"
)

func Codon2AA(codon string) (string, string) {
	var Lname, Sname, codon_res string
	codon_res = strings.ToUpper(codon)
	AA := map[string]map[string]string{
		"TCA": map[string]string{
			"LName": "Ser",
			"SName": "S",
		},
		"TCC": map[string]string{
			"LName": "Ser",
			"SName": "S",
		},
		"TCG": map[string]string{
			"LName": "Ser",
			"SName": "S",
		},
		"TCT": map[string]string{
			"LName": "Ser",
			"SName": "S",
		},
		"TTC": map[string]string{
			"LName": "Phe",
			"SName": "F",
		},
		"TTT": map[string]string{
			"LName": "Phe",
			"SName": "F",
		},
		"TTA": map[string]string{
			"LName": "Leu",
			"SName": "L",
		},
		"TTG": map[string]string{
			"LName": "Leu",
			"SName": "L",
		},
		"TAC": map[string]string{
			"LName": "Tyr",
			"SName": "Y",
		},
		"TAT": map[string]string{
			"LName": "Tyr",
			"SName": "Y",
		},
		"TAA": map[string]string{
			"LName": "X",
			"SName": "X",
		},
		"TAG": map[string]string{
			"LName": "X",
			"SName": "X",
		},
		"TGC": map[string]string{
			"LName": "Cys",
			"SName": "C",
		},
		"TGT": map[string]string{
			"LName": "Cys",
			"SName": "C",
		},
		"TGA": map[string]string{
			"LName": "X",
			"SName": "X",
		},
		"TGG": map[string]string{
			"LName": "Trp",
			"SName": "W",
		},
		"CTA": map[string]string{
			"LName": "Leu",
			"SName": "L",
		},
		"CTC": map[string]string{
			"LName": "Leu",
			"SName": "L",
		},
		"CTG": map[string]string{
			"LName": "Leu",
			"SName": "L",
		},
		"CTT": map[string]string{
			"LName": "Leu",
			"SName": "L",
		},
		"CCA": map[string]string{
			"LName": "Pro",
			"SName": "P",
		},
		"CCC": map[string]string{
			"LName": "Pro",
			"SName": "P",
		},
		"CCG": map[string]string{
			"LName": "Pro",
			"SName": "P",
		},
		"CCT": map[string]string{
			"LName": "Pro",
			"SName": "P",
		},

		"CAC": map[string]string{
			"LName": "His",
			"SName": "H",
		},
		"CAT": map[string]string{
			"LName": "His",
			"SName": "H",
		},
		"CAA": map[string]string{
			"LName": "Gln",
			"SName": "Q",
		},
		"CAG": map[string]string{
			"LName": "Gln",
			"SName": "Q",
		},
		"CGA": map[string]string{
			"LName": "Arg",
			"SName": "R",
		},
		"CGC": map[string]string{
			"LName": "Arg",
			"SName": "R",
		},
		"CGG": map[string]string{
			"LName": "Arg",
			"SName": "R",
		},
		"CGT": map[string]string{
			"LName": "Arg",
			"SName": "R",
		},
		"ATA": map[string]string{
			"LName": "Ile",
			"SName": "I",
		},
		"ATC": map[string]string{
			"LName": "Ile",
			"SName": "I",
		},
		"ATT": map[string]string{
			"LName": "Ile",
			"SName": "I",
		},
		"ATG": map[string]string{
			"LName": "Met",
			"SName": "M",
		},
		"ACA": map[string]string{
			"LName": "Thr",
			"SName": "T",
		},
		"ACC": map[string]string{
			"LName": "Thr",
			"SName": "T",
		},
		"ACG": map[string]string{
			"LName": "Thr",
			"SName": "T",
		},
		"ACT": map[string]string{
			"LName": "Thr",
			"SName": "T",
		},
		"AAC": map[string]string{
			"LName": "Asn",
			"SName": "N",
		},
		"AAT": map[string]string{
			"LName": "Asn",
			"SName": "N",
		},
		"AAA": map[string]string{
			"LName": "Lys",
			"SName": "K",
		},
		"AAG": map[string]string{
			"LName": "Lys",
			"SName": "K",
		},
		"AGC": map[string]string{
			"LName": "Ser",
			"SName": "S",
		},
		"AGT": map[string]string{
			"LName": "Ser",
			"SName": "S",
		},
		"AGA": map[string]string{
			"LName": "Arg",
			"SName": "R",
		},
		"AGG": map[string]string{
			"LName": "Arg",
			"SName": "R",
		},
		"GTA": map[string]string{
			"LName": "Val",
			"SName": "V",
		},
		"GTC": map[string]string{
			"LName": "Val",
			"SName": "V",
		},
		"GTG": map[string]string{
			"LName": "Val",
			"SName": "V",
		},
		"GTT": map[string]string{
			"LName": "Val",
			"SName": "V",
		},
		"GCA": map[string]string{
			"LName": "Ala",
			"SName": "A",
		},
		"GCC": map[string]string{
			"LName": "Ala",
			"SName": "A",
		},
		"GCG": map[string]string{
			"LName": "Ala",
			"SName": "A",
		},
		"GCT": map[string]string{
			"LName": "Ala",
			"SName": "A",
		},
		"GAC": map[string]string{
			"LName": "Asp",
			"SName": "D",
		},
		"GAT": map[string]string{
			"LName": "Asp",
			"SName": "D",
		},
		"GAA": map[string]string{
			"LName": "Glu",
			"SName": "E",
		},
		"GAG": map[string]string{
			"LName": "Glu",
			"SName": "E",
		},
		"GGA": map[string]string{
			"LName": "Gly",
			"SName": "G",
		},
		"GGC": map[string]string{
			"LName": "Gly",
			"SName": "G",
		},
		"GGG": map[string]string{
			"LName": "Gly",
			"SName": "G",
		},
		"GGT": map[string]string{
			"LName": "Gly",
			"SName": "G",
		},
	}

	if aa, ok := AA[codon_res]; ok {
		// AA[codon]["Lname"], AA[codon]["SName"]
		Lname = aa["LName"]
		Sname = aa["SName"]
		// fmt.Println(Lname, Sname)

	}
	return Lname, Sname
}
