package amino

import (
	// "fmt"
	"bytes"
	"fmt"
	"strings"
)

var aa = map[string]map[string]string{
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

var tangIdx = map[string]map[string]float64{
	"ST": {"Tang": 2.490},
	"QP": {"Tang": 1.377},
	"EA": {"Tang": 0.901},
	"LH": {"Tang": 0.560},
	"RL": {"Tang": 0.414},
	"VI": {"Tang": 2.415},
	"SG": {"Tang": 1.360},
	"SC": {"Tang": 0.852},
	"KM": {"Tang": 0.559},
	"GC": {"Tang": 0.414},
	"SA": {"Tang": 2.380},
	"AS": {"Tang": 2.380},
	"QH": {"Tang": 1.351},
	"RS": {"Tang": 0.850},
	"RP": {"Tang": 0.559},
	"PL": {"Tang": 0.388},
	"NS": {"Tang": 2.053},
	"VL": {"Tang": 1.329},
	"RT": {"Tang": 0.827},
	"EG": {"Tang": 0.553},
	"RC": {"Tang": 0.382},
	"DE": {"Tang": 2.033},
	"RH": {"Tang": 1.317},
	"IM": {"Tang": 0.827},
	"VF": {"Tang": 0.548},
	"NY": {"Tang": 0.378},
	"IL": {"Tang": 1.726},
	"AP": {"Tang": 1.288},
	"QL": {"Tang": 0.805},
	"EK": {"Tang": 0.548},
	"SW": {"Tang": 0.375},
	"NT": {"Tang": 1.695},
	"KN": {"Tang": 1.075},
	"LW": {"Tang": 0.793},
	"DG": {"Tang": 0.548},
	"SF": {"Tang": 0.365},
	"YF": {"Tang": 1.649},
	"RQ": {"Tang": 1.045},
	"PH": {"Tang": 0.784},
	"IF": {"Tang": 0.545},
	"DV": {"Tang": 0.361},
	"EQ": {"Tang": 1.634},
	"SP": {"Tang": 1.039},
	"TI": {"Tang": 0.750},
	"SI": {"Tang": 0.540},
	"CF": {"Tang": 0.321},
	"LM": {"Tang": 1.601},
	"AV": {"Tang": 1.017},
	"LF": {"Tang": 0.732},
	"GV": {"Tang": 0.539},
	"NI": {"Tang": 0.321},
	"TA": {"Tang": 1.587},
	"DN": {"Tang": 1.015},
	"SL": {"Tang": 0.725},
	"RG": {"Tang": 0.534},
	"CW": {"Tang": 0.271},
	"RK": {"Tang": 1.583},
	"TM": {"Tang": 1.007},
	"KI": {"Tang": 0.688},
	"EV": {"Tang": 0.506},
	"CY": {"Tang": 0.268},
	"KQ": {"Tang": 1.466},
	"TP": {"Tang": 1.001},
	"HY": {"Tang": 0.665},
	"SY": {"Tang": 0.503},
	"RW": {"Tang": 0.263},
	"NH": {"Tang": 1.382},
	"KT": {"Tang": 0.989},
	"DA": {"Tang": 0.657},
	"RI": {"Tang": 0.490},
	"GW": {"Tang": 0.242},
	"GA": {"Tang": 1.379},
	"VM": {"Tang": 0.986},
	"DH": {"Tang": 0.560},
	"RM": {"Tang": 0.470},
	"DY": {"Tang": 0.241},
	"TS": {"Tang": 2.490},
	"PQ": {"Tang": 1.377},
	"AE": {"Tang": 0.906},
	"HL": {"Tang": 0.560},
	"LR": {"Tang": 0.414},
	"IV": {"Tang": 2.415},
	"GS": {"Tang": 1.360},
	"CS": {"Tang": 0.852},
	"MK": {"Tang": 0.559},
	"CG": {"Tang": 0.414},
	"HQ": {"Tang": 1.351},
	"SR": {"Tang": 0.850},
	"PR": {"Tang": 0.559},
	"LP": {"Tang": 0.388},
	"SN": {"Tang": 2.053},
	"LV": {"Tang": 1.329},
	"TR": {"Tang": 0.827},
	"GE": {"Tang": 0.553},
	"CR": {"Tang": 0.382},
	"ED": {"Tang": 2.033},
	"HR": {"Tang": 1.317},
	"MI": {"Tang": 0.827},
	"FV": {"Tang": 0.548},
	"YN": {"Tang": 0.378},
	"LI": {"Tang": 1.726},
	"PA": {"Tang": 1.288},
	"LQ": {"Tang": 0.805},
	"KE": {"Tang": 0.548},
	"WS": {"Tang": 0.375},
	"TN": {"Tang": 1.695},
	"NK": {"Tang": 1.075},
	"WL": {"Tang": 0.793},
	"GD": {"Tang": 0.548},
	"FS": {"Tang": 0.365},
	"FY": {"Tang": 1.649},
	"QR": {"Tang": 1.045},
	"HP": {"Tang": 0.784},
	"FI": {"Tang": 0.545},
	"VD": {"Tang": 0.361},
	"QE": {"Tang": 1.634},
	"PS": {"Tang": 1.039},
	"IT": {"Tang": 0.750},
	"IS": {"Tang": 0.540},
	"FC": {"Tang": 0.321},
	"ML": {"Tang": 1.601},
	"VA": {"Tang": 1.017},
	"FL": {"Tang": 0.732},
	"VG": {"Tang": 0.539},
	"IN": {"Tang": 0.321},
	"AT": {"Tang": 1.587},
	"ND": {"Tang": 1.015},
	"AD": {"Tang": 0.156},
	"LS": {"Tang": 0.725},
	"GR": {"Tang": 0.534},
	"WC": {"Tang": 0.271},
	"KR": {"Tang": 1.583},
	"MT": {"Tang": 1.007},
	"IK": {"Tang": 0.688},
	"VE": {"Tang": 0.506},
	"YC": {"Tang": 0.268},
	"QK": {"Tang": 1.466},
	"PT": {"Tang": 1.001},
	"YH": {"Tang": 0.665},
	"YS": {"Tang": 0.503},
	"WR": {"Tang": 0.263},
	"HN": {"Tang": 1.382},
	"TK": {"Tang": 0.989},
	"AA": {"Tang": 0.657},
	"IR": {"Tang": 0.490},
	"WG": {"Tang": 0.242},
	"AG": {"Tang": 1.379},
	"MV": {"Tang": 0.986},
	"HD": {"Tang": 0.560},
	"MR": {"Tang": 0.470},
	"YD": {"Tang": 0.241},
	"FH": {"Tang": 0},
	"HE": {"Tang": 0},
	"QY": {"Tang": 0},
	"QI": {"Tang": 0},
	"RA": {"Tang": 0},
	"HF": {"Tang": 0},
	"EH": {"Tang": 0},
	"YQ": {"Tang": 0},
	"IQ": {"Tang": 0},
	"AR": {"Tang": 0},
}

var cdaAA = map[string]map[string]float32{
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

func Codon2AA(codon string) (string, string) {
	/*


	 */
	var Lname, Sname, codonRes string
	codonRes = strings.ToUpper(codon)

	if aa, ok := aa[codonRes]; ok {
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

	if aa, ok := cdaAA[codonRes]; ok {
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

// Применение данной формулы рекомендовано автора-
// ми для набора данных с более чем 20 000 аминокислотных замен или,
// по меньшей мере, более чем 2500 замен, что и является основным ограни-
// чением широкого использования метода.
// Таким образом, метод Х. Танга может быть охарактеризован как эм-
// пирический, основанный на анализе эволюционных изменений кодонов,
// применимый для близкородственных видов, и универсальный (для раз-
// личных таксономических групп организмов).

func GetTangInx(refAA string, altAA string) float64 {
	var buffer bytes.Buffer
	var res float64
	buffer.WriteString(refAA)
	buffer.WriteString(altAA)
	// Tang := make(map[string]float64)

	// Tang["ST"] = 2.490
	// Tang["QP"] = 1.377
	// Tang["EA"] = 0.906
	// Tang["LH"] = 0.560
	// Tang["RL"] = 0.414
	// Tang["VI"] = 2.415
	// Tang["SG"] = 1.360
	// Tang["SC"] = 0.852
	// Tang["KM"] = 0.559
	// Tang["GC"] = 0.414
	// Tang["SA"] = 2.380
	// Tang["AS"] = 2.380
	// Tang["QH"] = 1.351
	// Tang["RS"] = 0.850
	// Tang["RP"] = 0.559
	// Tang["PL"] = 0.388
	// Tang["NS"] = 2.053
	// Tang["VL"] = 1.329
	// Tang["RT"] = 0.827
	// Tang["EG"] = 0.553
	// Tang["RC"] = 0.382
	// Tang["DE"] = 2.033
	// Tang["RH"] = 1.317
	// Tang["IM"] = 0.827
	// Tang["VF"] = 0.548
	// Tang["NY"] = 0.378
	// Tang["IL"] = 1.726
	// Tang["AP"] = 1.288
	// Tang["QL"] = 0.805
	// Tang["EK"] = 0.548
	// Tang["SW"] = 0.375
	// Tang["NT"] = 1.695
	// Tang["KN"] = 1.075
	// Tang["LW"] = 0.793
	// Tang["DG"] = 0.548
	// Tang["SF"] = 0.365
	// Tang["YF"] = 1.649
	// Tang["RQ"] = 1.045
	// Tang["PH"] = 0.784
	// Tang["IF"] = 0.545
	// Tang["DV"] = 0.361
	// Tang["EQ"] = 1.634
	// Tang["SP"] = 1.039
	// Tang["TI"] = 0.750
	// Tang["SI"] = 0.540
	// Tang["CF"] = 0.321
	// Tang["LM"] = 1.601
	// Tang["AV"] = 1.017
	// Tang["LF"] = 0.732
	// Tang["GV"] = 0.539
	// Tang["NI"] = 0.321
	// Tang["TA"] = 1.587
	// Tang["DN"] = 1.015
	// Tang["SL"] = 0.725
	// Tang["RG"] = 0.534
	// Tang["CW"] = 0.271
	// Tang["RK"] = 1.583
	// Tang["TM"] = 1.007
	// Tang["KI"] = 0.688
	// Tang["EV"] = 0.506
	// Tang["CY"] = 0.268
	// Tang["KQ"] = 1.466
	// Tang["TP"] = 1.001
	// Tang["HY"] = 0.665
	// Tang["SY"] = 0.503
	// Tang["RW"] = 0.263
	// Tang["NH"] = 1.382
	// Tang["KT"] = 0.989
	// Tang["DA"] = 0.657
	// Tang["RI"] = 0.490
	// Tang["GW"] = 0.242
	// Tang["GA"] = 1.379
	// Tang["VM"] = 0.986
	// Tang["DH"] = 0.560
	// Tang["RM"] = 0.470
	// Tang["DY"] = 0.241
	// Tang["TS"] = 2.490
	// Tang["PQ"] = 1.377
	// Tang["AE"] = 0.906
	// Tang["HL"] = 0.560
	// Tang["LR"] = 0.414
	// Tang["IV"] = 2.415
	// Tang["GS"] = 1.360
	// Tang["CS"] = 0.852
	// Tang["MK"] = 0.559
	// Tang["CG"] = 0.414
	// Tang["AS"] = 2.380
	// Tang["SA"] = 2.380
	// Tang["HQ"] = 1.351
	// Tang["SR"] = 0.850
	// Tang["PR"] = 0.559
	// Tang["LP"] = 0.388
	// Tang["SN"] = 2.053
	// Tang["LV"] = 1.329
	// Tang["TR"] = 0.827
	// Tang["GE"] = 0.553
	// Tang["CR"] = 0.382
	// Tang["ED"] = 2.033
	// Tang["HR"] = 1.317
	// Tang["MI"] = 0.827
	// Tang["FV"] = 0.548
	// Tang["YN"] = 0.378
	// Tang["LI"] = 1.726
	// Tang["PA"] = 1.288
	// Tang["LQ"] = 0.805
	// Tang["KE"] = 0.548
	// Tang["WS"] = 0.375
	// Tang["TN"] = 1.695
	// Tang["NK"] = 1.075
	// Tang["WL"] = 0.793
	// Tang["GD"] = 0.548
	// Tang["FS"] = 0.365
	// Tang["FY"] = 1.649
	// Tang["QR"] = 1.045
	// Tang["HP"] = 0.784
	// Tang["FI"] = 0.545
	// Tang["VD"] = 0.361
	// Tang["QE"] = 1.634
	// Tang["PS"] = 1.039
	// Tang["IT"] = 0.750
	// Tang["IS"] = 0.540
	// Tang["FC"] = 0.321
	// Tang["ML"] = 1.601
	// Tang["VA"] = 1.017
	// Tang["FL"] = 0.732
	// Tang["VG"] = 0.539
	// Tang["IN"] = 0.321
	// Tang["AT"] = 1.587
	// Tang["ND"] = 1.015
	// Tang["AD"] = 0.156
	// Tang["LS"] = 0.725
	// Tang["GR"] = 0.534
	// Tang["WC"] = 0.271
	// Tang["KR"] = 1.583
	// Tang["MT"] = 1.007
	// Tang["IK"] = 0.688
	// Tang["VE"] = 0.506
	// Tang["YC"] = 0.268
	// Tang["QK"] = 1.466
	// Tang["PT"] = 1.001
	// Tang["YH"] = 0.665
	// Tang["YS"] = 0.503
	// Tang["WR"] = 0.263
	// Tang["HN"] = 1.382
	// Tang["TK"] = 0.989
	// Tang["AA"] = 0.657
	// Tang["IR"] = 0.490
	// Tang["WG"] = 0.242
	// Tang["AG"] = 1.379
	// Tang["MV"] = 0.986
	// Tang["HD"] = 0.560
	// Tang["MR"] = 0.470
	// Tang["YD"] = 0.241
	// Tang["FH"] = 0
	// Tang["HE"] = 0
	// Tang["QY"] = 0
	// Tang["QI"] = 0
	// Tang["RA"] = 0
	// Tang["HF"] = 0
	// Tang["EH"] = 0
	// Tang["YQ"] = 0
	// Tang["IQ"] = 0
	// Tang["AR"] = 0
	// // Tang["W"] = "-"
	// Tang["E"] = "-"
	// Tang["G"] = "-"
	// Tang["Q"] = "-"
	// Tang["Y"] = "-"

	// fmt.Printf("%v %v %v\n", refAA, altAA, Tang[buffer.String()])
	if altAA != "" && refAA != "" && altAA != "X" {
		res = tangIdx[buffer.String()]["Tang"]
	} else {
		res = 0
	}
	return res
}

// При использовании модифи-
// цированных физико-химических дистанций замены аминокислот счита-
// ются консервативными при значении GDM выше среднего (более 52,3 для
// замен в целом, более 57,9 для одношаговых замен), в обратном случае —
// радикальными.

func GetGHInx(refAA string, altAA string) int {
	var buffer bytes.Buffer
	var res int
	buffer.WriteString(refAA)
	buffer.WriteString(altAA)
	GH := make(map[string]int)

	GH["AR"] = 112
	GH["RA"] = 112
	GH["AN"] = 111
	GH["NA"] = 111
	GH["AD"] = 126
	GH["DA"] = 126
	GH["AC"] = 195
	GH["CA"] = 195
	GH["AQ"] = 91
	GH["QA"] = 91
	GH["AE"] = 107
	GH["EA"] = 107
	GH["AG"] = 60
	GH["GA"] = 60
	GH["AH"] = 86
	GH["HA"] = 86
	GH["AI"] = 94
	GH["IA"] = 94
	GH["AL"] = 96
	GH["LA"] = 96
	GH["AK"] = 106
	GH["KA"] = 106
	GH["AM"] = 84
	GH["MA"] = 84
	GH["AF"] = 113
	GH["FA"] = 113
	GH["AP"] = 27
	GH["PA"] = 27
	GH["AS"] = 99
	GH["SA"] = 99
	GH["AT"] = 58
	GH["TA"] = 58
	GH["AW"] = 148
	GH["WA"] = 148
	GH["AY"] = 112
	GH["YA"] = 112
	GH["RA"] = 112
	GH["AR"] = 112
	GH["RN"] = 86
	GH["NR"] = 86
	GH["RD"] = 96
	GH["DR"] = 96
	GH["RC"] = 180
	GH["CR"] = 180
	GH["RQ"] = 43
	GH["QR"] = 43
	GH["RE"] = 54
	GH["ER"] = 54
	GH["RG"] = 125
	GH["GR"] = 125
	GH["RH"] = 29
	GH["HR"] = 29
	GH["RI"] = 97
	GH["IR"] = 97
	GH["RL"] = 102
	GH["LR"] = 102
	GH["RK"] = 26
	GH["KR"] = 26
	GH["RM"] = 91
	GH["MR"] = 91
	GH["RF"] = 97
	GH["FR"] = 97
	GH["RP"] = 103
	GH["PR"] = 103
	GH["RS"] = 110
	GH["SR"] = 110
	GH["RT"] = 71
	GH["TR"] = 71
	GH["RW"] = 101
	GH["WR"] = 101
	GH["RY"] = 77
	GH["YR"] = 77
	GH["NA"] = 111
	GH["AN"] = 111
	GH["NR"] = 86
	GH["RN"] = 86
	GH["ND"] = 23
	GH["DN"] = 23
	GH["NC"] = 139
	GH["CN"] = 139
	GH["NQ"] = 46
	GH["QN"] = 46
	GH["NE"] = 42
	GH["EN"] = 42
	GH["NG"] = 80
	GH["GN"] = 80
	GH["NH"] = 68
	GH["HN"] = 68
	GH["NI"] = 149
	GH["IN"] = 149
	GH["NL"] = 153
	GH["LN"] = 153
	GH["NK"] = 94
	GH["KN"] = 94
	GH["NM"] = 142
	GH["MN"] = 142
	GH["NF"] = 158
	GH["FN"] = 158
	GH["NP"] = 91
	GH["PN"] = 91
	GH["NS"] = 46
	GH["SN"] = 46
	GH["NT"] = 65
	GH["TN"] = 65
	GH["NW"] = 174
	GH["WN"] = 174
	GH["NY"] = 143
	GH["YN"] = 143
	GH["DA"] = 126
	GH["AD"] = 126
	GH["DR"] = 96
	GH["RD"] = 96
	GH["DN"] = 23
	GH["ND"] = 23
	GH["DC"] = 154
	GH["CD"] = 154
	GH["DQ"] = 61
	GH["QD"] = 61
	GH["DE"] = 45
	GH["ED"] = 45
	GH["DG"] = 94
	GH["GD"] = 94
	GH["DH"] = 81
	GH["HD"] = 81
	GH["DI"] = 168
	GH["ID"] = 168
	GH["DL"] = 172
	GH["LD"] = 172
	GH["DK"] = 101
	GH["KD"] = 101
	GH["DM"] = 160
	GH["MD"] = 160
	GH["DF"] = 177
	GH["FD"] = 177
	GH["DP"] = 108
	GH["PD"] = 108
	GH["DS"] = 65
	GH["SD"] = 65
	GH["DT"] = 85
	GH["TD"] = 85
	GH["DW"] = 181
	GH["WD"] = 181
	GH["DY"] = 160
	GH["YD"] = 160
	GH["CA"] = 195
	GH["AC"] = 195
	GH["CR"] = 180
	GH["RC"] = 180
	GH["CN"] = 139
	GH["NC"] = 139
	GH["CD"] = 154
	GH["DC"] = 154
	GH["CQ"] = 154
	GH["QC"] = 154
	GH["CE"] = 170
	GH["EC"] = 170
	GH["CG"] = 159
	GH["GC"] = 159
	GH["CH"] = 174
	GH["HC"] = 174
	GH["CI"] = 198
	GH["IC"] = 198
	GH["CL"] = 198
	GH["LC"] = 198
	GH["CK"] = 202
	GH["KC"] = 202
	GH["CM"] = 196
	GH["MC"] = 196
	GH["CF"] = 205
	GH["FC"] = 205
	GH["CP"] = 169
	GH["PC"] = 169
	GH["CS"] = 112
	GH["SC"] = 112
	GH["CT"] = 149
	GH["TC"] = 149
	GH["CW"] = 215
	GH["WC"] = 215
	GH["CY"] = 194
	GH["YC"] = 194
	GH["QA"] = 91
	GH["AQ"] = 91
	GH["QR"] = 43
	GH["RQ"] = 43
	GH["QN"] = 46
	GH["NQ"] = 46
	GH["QD"] = 61
	GH["DQ"] = 61
	GH["QC"] = 154
	GH["CQ"] = 154
	GH["QE"] = 29
	GH["EQ"] = 29
	GH["QG"] = 87
	GH["GQ"] = 87
	GH["QH"] = 24
	GH["HQ"] = 24
	GH["QI"] = 109
	GH["IQ"] = 109
	GH["QL"] = 113
	GH["LQ"] = 113
	GH["QK"] = 53
	GH["KQ"] = 53
	GH["QM"] = 101
	GH["MQ"] = 101
	GH["QF"] = 116
	GH["FQ"] = 116
	GH["QP"] = 76
	GH["PQ"] = 76
	GH["QS"] = 68
	GH["SQ"] = 68
	GH["QT"] = 42
	GH["TQ"] = 42
	GH["QW"] = 130
	GH["WQ"] = 130
	GH["QY"] = 99
	GH["YQ"] = 99
	GH["EA"] = 107
	GH["AE"] = 107
	GH["ER"] = 54
	GH["RE"] = 54
	GH["EN"] = 42
	GH["NE"] = 42
	GH["ED"] = 45
	GH["DE"] = 45
	GH["EC"] = 170
	GH["CE"] = 170
	GH["EQ"] = 29
	GH["QE"] = 29
	GH["EG"] = 98
	GH["GE"] = 98
	GH["EH"] = 40
	GH["HE"] = 40
	GH["EI"] = 134
	GH["IE"] = 134
	GH["EL"] = 138
	GH["LE"] = 138
	GH["EK"] = 56
	GH["KE"] = 56
	GH["EM"] = 126
	GH["ME"] = 126
	GH["EF"] = 140
	GH["FE"] = 140
	GH["EP"] = 93
	GH["PE"] = 93
	GH["ES"] = 80
	GH["SE"] = 80
	GH["ET"] = 65
	GH["TE"] = 65
	GH["EW"] = 152
	GH["WE"] = 152
	GH["EY"] = 122
	GH["YE"] = 122
	GH["GA"] = 60
	GH["AG"] = 60
	GH["GR"] = 125
	GH["RG"] = 125
	GH["GN"] = 80
	GH["NG"] = 80
	GH["GD"] = 94
	GH["DG"] = 94
	GH["GC"] = 159
	GH["CG"] = 159
	GH["GQ"] = 87
	GH["QG"] = 87
	GH["GE"] = 98
	GH["EG"] = 98
	GH["GH"] = 98
	GH["HG"] = 98
	GH["GI"] = 135
	GH["IG"] = 135
	GH["GL"] = 138
	GH["LG"] = 138
	GH["GK"] = 127
	GH["KG"] = 127
	GH["GM"] = 127
	GH["MG"] = 127
	GH["GF"] = 153
	GH["FG"] = 153
	GH["GP"] = 42
	GH["PG"] = 42
	GH["GS"] = 56
	GH["SG"] = 56
	GH["GT"] = 59
	GH["TG"] = 59
	GH["GW"] = 184
	GH["WG"] = 184
	GH["GY"] = 147
	GH["YG"] = 147
	GH["HA"] = 86
	GH["AH"] = 86
	GH["HR"] = 29
	GH["RH"] = 29
	GH["HN"] = 68
	GH["NH"] = 68
	GH["HD"] = 81
	GH["DH"] = 81
	GH["HC"] = 174
	GH["CH"] = 174
	GH["HQ"] = 24
	GH["QH"] = 24
	GH["HE"] = 40
	GH["EH"] = 40
	GH["HG"] = 98
	GH["GH"] = 98
	GH["HI"] = 94
	GH["IH"] = 94
	GH["HL"] = 99
	GH["LH"] = 99
	GH["HK"] = 32
	GH["KH"] = 32
	GH["HM"] = 87
	GH["MH"] = 87
	GH["HF"] = 100
	GH["FH"] = 100
	GH["HP"] = 77
	GH["PH"] = 77
	GH["HS"] = 89
	GH["SH"] = 89
	GH["HT"] = 47
	GH["TH"] = 47
	GH["HW"] = 115
	GH["WH"] = 115
	GH["HY"] = 83
	GH["YH"] = 83
	GH["IA"] = 94
	GH["AI"] = 94
	GH["IR"] = 97
	GH["RI"] = 97
	GH["IN"] = 149
	GH["NI"] = 149
	GH["ID"] = 168
	GH["DI"] = 168
	GH["IC"] = 198
	GH["CI"] = 198
	GH["IQ"] = 109
	GH["QI"] = 109
	GH["IE"] = 134
	GH["EI"] = 134
	GH["IG"] = 135
	GH["GI"] = 135
	GH["IH"] = 94
	GH["HI"] = 94
	GH["IL"] = 5
	GH["LI"] = 5
	GH["IK"] = 102
	GH["KI"] = 102
	GH["IM"] = 10
	GH["MI"] = 10
	GH["IF"] = 21
	GH["FI"] = 21
	GH["IP"] = 95
	GH["PI"] = 95
	GH["IS"] = 142
	GH["SI"] = 142
	GH["IT"] = 89
	GH["TI"] = 89
	GH["IW"] = 61
	GH["WI"] = 61
	GH["IY"] = 33
	GH["YI"] = 33
	GH["LA"] = 96
	GH["AL"] = 96
	GH["LR"] = 102
	GH["RL"] = 102
	GH["LN"] = 153
	GH["NL"] = 153
	GH["LD"] = 172
	GH["DL"] = 172
	GH["LC"] = 198
	GH["CL"] = 198
	GH["LQ"] = 113
	GH["QL"] = 113
	GH["LE"] = 138
	GH["EL"] = 138
	GH["LG"] = 138
	GH["GL"] = 138
	GH["LH"] = 99
	GH["HL"] = 99
	GH["LI"] = 5
	GH["IL"] = 5
	GH["LK"] = 107
	GH["KL"] = 107
	GH["LM"] = 15
	GH["ML"] = 15
	GH["LF"] = 22
	GH["FL"] = 22
	GH["LP"] = 98
	GH["PL"] = 98
	GH["LS"] = 145
	GH["SL"] = 145
	GH["LT"] = 92
	GH["TL"] = 92
	GH["LW"] = 61
	GH["WL"] = 61
	GH["LY"] = 36
	GH["YL"] = 36
	GH["KA"] = 106
	GH["AK"] = 106
	GH["KR"] = 26
	GH["RK"] = 26
	GH["KN"] = 94
	GH["NK"] = 94
	GH["KD"] = 101
	GH["DK"] = 101
	GH["KC"] = 202
	GH["CK"] = 202
	GH["KQ"] = 53
	GH["QK"] = 53
	GH["KE"] = 56
	GH["EK"] = 56
	GH["KG"] = 127
	GH["GK"] = 127
	GH["KH"] = 32
	GH["HK"] = 32
	GH["KI"] = 102
	GH["IK"] = 102
	GH["KL"] = 107
	GH["LK"] = 107
	GH["KM"] = 95
	GH["MK"] = 95
	GH["KF"] = 102
	GH["FK"] = 102
	GH["KP"] = 103
	GH["PK"] = 103
	GH["KS"] = 121
	GH["SK"] = 121
	GH["KT"] = 78
	GH["TK"] = 78
	GH["KW"] = 110
	GH["WK"] = 110
	GH["KY"] = 85
	GH["YK"] = 85
	GH["MA"] = 84
	GH["AM"] = 84
	GH["MR"] = 91
	GH["RM"] = 91
	GH["MN"] = 142
	GH["NM"] = 142
	GH["MD"] = 160
	GH["DM"] = 160
	GH["MC"] = 196
	GH["CM"] = 196
	GH["MQ"] = 101
	GH["QM"] = 101
	GH["ME"] = 126
	GH["EM"] = 126
	GH["MG"] = 127
	GH["GM"] = 127
	GH["MH"] = 87
	GH["HM"] = 87
	GH["MI"] = 10
	GH["IM"] = 10
	GH["ML"] = 15
	GH["LM"] = 15
	GH["MK"] = 95
	GH["KM"] = 95
	GH["MF"] = 28
	GH["FM"] = 28
	GH["MP"] = 87
	GH["PM"] = 87
	GH["MS"] = 135
	GH["SM"] = 135
	GH["MT"] = 81
	GH["TM"] = 81
	GH["MW"] = 67
	GH["WM"] = 67
	GH["MY"] = 36
	GH["YM"] = 36
	GH["FA"] = 113
	GH["AF"] = 113
	GH["FR"] = 97
	GH["RF"] = 97
	GH["FN"] = 158
	GH["NF"] = 158
	GH["FD"] = 177
	GH["DF"] = 177
	GH["FC"] = 205
	GH["CF"] = 205
	GH["FQ"] = 116
	GH["QF"] = 116
	GH["FE"] = 140
	GH["EF"] = 140
	GH["FG"] = 153
	GH["GF"] = 153
	GH["FH"] = 100
	GH["HF"] = 100
	GH["FI"] = 21
	GH["IF"] = 21
	GH["FL"] = 22
	GH["LF"] = 22
	GH["FK"] = 102
	GH["KF"] = 102
	GH["FM"] = 28
	GH["MF"] = 28
	GH["FP"] = 114
	GH["PF"] = 114
	GH["FS"] = 155
	GH["SF"] = 155
	GH["FT"] = 103
	GH["TF"] = 103
	GH["FW"] = 40
	GH["WF"] = 40
	GH["FY"] = 22
	GH["YF"] = 22
	GH["PA"] = 27
	GH["AP"] = 27
	GH["PR"] = 103
	GH["RP"] = 103
	GH["PN"] = 91
	GH["NP"] = 91
	GH["PD"] = 108
	GH["DP"] = 108
	GH["PC"] = 169
	GH["CP"] = 169
	GH["PQ"] = 76
	GH["QP"] = 76
	GH["PE"] = 93
	GH["EP"] = 93
	GH["PG"] = 42
	GH["GP"] = 42
	GH["PH"] = 77
	GH["HP"] = 77
	GH["PI"] = 95
	GH["IP"] = 95
	GH["PL"] = 98
	GH["LP"] = 98
	GH["PK"] = 103
	GH["KP"] = 103
	GH["PM"] = 87
	GH["MP"] = 87
	GH["PF"] = 114
	GH["FP"] = 114
	GH["PS"] = 74
	GH["SP"] = 74
	GH["PT"] = 38
	GH["TP"] = 38
	GH["PW"] = 147
	GH["WP"] = 147
	GH["PY"] = 110
	GH["YP"] = 110
	GH["SA"] = 99
	GH["AS"] = 99
	GH["SR"] = 110
	GH["RS"] = 110
	GH["SN"] = 46
	GH["NS"] = 46
	GH["SD"] = 65
	GH["DS"] = 65
	GH["SC"] = 112
	GH["CS"] = 112
	GH["SQ"] = 68
	GH["QS"] = 68
	GH["SE"] = 80
	GH["ES"] = 80
	GH["SG"] = 56
	GH["GS"] = 56
	GH["SH"] = 89
	GH["HS"] = 89
	GH["SI"] = 142
	GH["IS"] = 142
	GH["SL"] = 145
	GH["LS"] = 145
	GH["SK"] = 121
	GH["KS"] = 121
	GH["SM"] = 135
	GH["MS"] = 135
	GH["SF"] = 155
	GH["FS"] = 155
	GH["SP"] = 74
	GH["PS"] = 74
	GH["ST"] = 58
	GH["TS"] = 58
	GH["SW"] = 177
	GH["WS"] = 177
	GH["SY"] = 144
	GH["YS"] = 144
	GH["TA"] = 58
	GH["AT"] = 58
	GH["TR"] = 71
	GH["RT"] = 71
	GH["TN"] = 65
	GH["NT"] = 65
	GH["TD"] = 85
	GH["DT"] = 85
	GH["TC"] = 149
	GH["CT"] = 149
	GH["TQ"] = 42
	GH["QT"] = 42
	GH["TE"] = 65
	GH["ET"] = 65
	GH["TG"] = 59
	GH["GT"] = 59
	GH["TH"] = 47
	GH["HT"] = 47
	GH["TI"] = 89
	GH["IT"] = 89
	GH["TL"] = 92
	GH["LT"] = 92
	GH["TK"] = 78
	GH["KT"] = 78
	GH["TM"] = 81
	GH["MT"] = 81
	GH["TF"] = 103
	GH["FT"] = 103
	GH["TP"] = 38
	GH["PT"] = 38
	GH["TS"] = 58
	GH["ST"] = 58
	GH["TW"] = 128
	GH["WT"] = 128
	GH["TY"] = 92
	GH["YT"] = 92
	GH["WA"] = 148
	GH["AW"] = 148
	GH["WR"] = 101
	GH["RW"] = 101
	GH["WN"] = 174
	GH["NW"] = 174
	GH["WD"] = 181
	GH["DW"] = 181
	GH["WC"] = 215
	GH["CW"] = 215
	GH["WQ"] = 130
	GH["QW"] = 130
	GH["WE"] = 152
	GH["EW"] = 152
	GH["WG"] = 184
	GH["GW"] = 184
	GH["WH"] = 115
	GH["HW"] = 115
	GH["WI"] = 61
	GH["IW"] = 61
	GH["WL"] = 61
	GH["LW"] = 61
	GH["WK"] = 110
	GH["KW"] = 110
	GH["WM"] = 67
	GH["MW"] = 67
	GH["WF"] = 40
	GH["FW"] = 40
	GH["WP"] = 147
	GH["PW"] = 147
	GH["WS"] = 177
	GH["SW"] = 177
	GH["WT"] = 128
	GH["TW"] = 128
	GH["WY"] = 37
	GH["YW"] = 37
	GH["YA"] = 112
	GH["AY"] = 112
	GH["YR"] = 77
	GH["RY"] = 77
	GH["YN"] = 143
	GH["NY"] = 143
	GH["YD"] = 160
	GH["DY"] = 160
	GH["YC"] = 194
	GH["CY"] = 194
	GH["YQ"] = 99
	GH["QY"] = 99
	GH["YE"] = 122
	GH["EY"] = 122
	GH["YG"] = 147
	GH["GY"] = 147
	GH["YH"] = 83
	GH["HY"] = 83
	GH["YI"] = 33
	GH["IY"] = 33
	GH["YL"] = 36
	GH["LY"] = 36
	GH["YK"] = 85
	GH["KY"] = 85
	GH["YM"] = 36
	GH["MY"] = 36
	GH["YF"] = 22
	GH["FY"] = 22
	GH["YP"] = 110
	GH["PY"] = 110
	GH["YS"] = 144
	GH["SY"] = 144
	GH["YT"] = 92
	GH["TY"] = 92
	GH["YW"] = 37
	GH["WY"] = 37

	// fmt.Printf("%v %v %v\n", refAA, altAA, Tang[buffer.String()])
	if altAA != "" && refAA != "" && altAA != "X" {
		res = GH[buffer.String()]
	} else {
		res = 0
	}
	return res
}

func GetBLOSSInx(refAA string, altAA string) int {
	var buffer bytes.Buffer
	var res int
	buffer.WriteString(refAA)
	buffer.WriteString(altAA)
	BLOSS := make(map[string]int)

	BLOSS["AR"] = -1
	BLOSS["RA"] = -1
	BLOSS["AN"] = -2
	BLOSS["NA"] = -2
	BLOSS["AD"] = -2
	BLOSS["DA"] = -2
	BLOSS["AC"] = 0
	BLOSS["CA"] = 0
	BLOSS["AQ"] = -1
	BLOSS["QA"] = -1
	BLOSS["AE"] = -1
	BLOSS["EA"] = -1
	BLOSS["AG"] = 0
	BLOSS["GA"] = 0
	BLOSS["AH"] = -2
	BLOSS["HA"] = -2
	BLOSS["AI"] = -1
	BLOSS["IA"] = -1
	BLOSS["AL"] = -1
	BLOSS["LA"] = -1
	BLOSS["AK"] = -1
	BLOSS["KA"] = -1
	BLOSS["AM"] = -1
	BLOSS["MA"] = -1
	BLOSS["AF"] = -2
	BLOSS["FA"] = -2
	BLOSS["AP"] = -1
	BLOSS["PA"] = -1
	BLOSS["AS"] = 1
	BLOSS["SA"] = 1
	BLOSS["AT"] = 0
	BLOSS["TA"] = 0
	BLOSS["AW"] = -3
	BLOSS["WA"] = -3
	BLOSS["AY"] = -2
	BLOSS["YA"] = -2
	BLOSS["AV"] = 0
	BLOSS["VA"] = 0
	BLOSS["AB"] = -2
	BLOSS["BA"] = -2
	BLOSS["AZ"] = -1
	BLOSS["ZA"] = -1
	BLOSS["AX"] = 0
	BLOSS["XA"] = 0
	BLOSS["RA"] = -1
	BLOSS["AR"] = -1
	BLOSS["RN"] = 0
	BLOSS["NR"] = 0
	BLOSS["RD"] = -2
	BLOSS["DR"] = -2
	BLOSS["RC"] = -3
	BLOSS["CR"] = -3
	BLOSS["RQ"] = 1
	BLOSS["QR"] = 1
	BLOSS["RE"] = 0
	BLOSS["ER"] = 0
	BLOSS["RG"] = -2
	BLOSS["GR"] = -2
	BLOSS["RH"] = 0
	BLOSS["HR"] = 0
	BLOSS["RI"] = -3
	BLOSS["IR"] = -3
	BLOSS["RL"] = -2
	BLOSS["LR"] = -2
	BLOSS["RK"] = 2
	BLOSS["KR"] = 2
	BLOSS["RM"] = -1
	BLOSS["MR"] = -1
	BLOSS["RF"] = -3
	BLOSS["FR"] = -3
	BLOSS["RP"] = -2
	BLOSS["PR"] = -2
	BLOSS["RS"] = -1
	BLOSS["SR"] = -1
	BLOSS["RT"] = -1
	BLOSS["TR"] = -1
	BLOSS["RW"] = -3
	BLOSS["WR"] = -3
	BLOSS["RY"] = -2
	BLOSS["YR"] = -2
	BLOSS["RV"] = -3
	BLOSS["VR"] = -3
	BLOSS["RB"] = -1
	BLOSS["BR"] = -1
	BLOSS["RZ"] = 0
	BLOSS["ZR"] = 0
	BLOSS["RX"] = -1
	BLOSS["XR"] = -1
	BLOSS["NA"] = -2
	BLOSS["AN"] = -2
	BLOSS["NR"] = 0
	BLOSS["RN"] = 0
	BLOSS["ND"] = 1
	BLOSS["DN"] = 1
	BLOSS["NC"] = -3
	BLOSS["CN"] = -3
	BLOSS["NQ"] = 0
	BLOSS["QN"] = 0
	BLOSS["NE"] = 0
	BLOSS["EN"] = 0
	BLOSS["NG"] = 0
	BLOSS["GN"] = 0
	BLOSS["NH"] = 1
	BLOSS["HN"] = 1
	BLOSS["NI"] = -3
	BLOSS["IN"] = -3
	BLOSS["NL"] = -3
	BLOSS["LN"] = -3
	BLOSS["NK"] = 0
	BLOSS["KN"] = 0
	BLOSS["NM"] = -2
	BLOSS["MN"] = -2
	BLOSS["NF"] = -3
	BLOSS["FN"] = -3
	BLOSS["NP"] = -2
	BLOSS["PN"] = -2
	BLOSS["NS"] = 1
	BLOSS["SN"] = 1
	BLOSS["NT"] = 0
	BLOSS["TN"] = 0
	BLOSS["NW"] = -4
	BLOSS["WN"] = -4
	BLOSS["NY"] = -2
	BLOSS["YN"] = -2
	BLOSS["NV"] = -3
	BLOSS["VN"] = -3
	BLOSS["NB"] = 3
	BLOSS["BN"] = 3
	BLOSS["NZ"] = 0
	BLOSS["ZN"] = 0
	BLOSS["NX"] = -1
	BLOSS["XN"] = -1
	BLOSS["DA"] = -2
	BLOSS["AD"] = -2
	BLOSS["DR"] = -2
	BLOSS["RD"] = -2
	BLOSS["DN"] = 1
	BLOSS["ND"] = 1
	BLOSS["DC"] = -3
	BLOSS["CD"] = -3
	BLOSS["DQ"] = 0
	BLOSS["QD"] = 0
	BLOSS["DE"] = 2
	BLOSS["ED"] = 2
	BLOSS["DG"] = -1
	BLOSS["GD"] = -1
	BLOSS["DH"] = -1
	BLOSS["HD"] = -1
	BLOSS["DI"] = -3
	BLOSS["ID"] = -3
	BLOSS["DL"] = -4
	BLOSS["LD"] = -4
	BLOSS["DK"] = -1
	BLOSS["KD"] = -1
	BLOSS["DM"] = -3
	BLOSS["MD"] = -3
	BLOSS["DF"] = -3
	BLOSS["FD"] = -3
	BLOSS["DP"] = -1
	BLOSS["PD"] = -1
	BLOSS["DS"] = 0
	BLOSS["SD"] = 0
	BLOSS["DT"] = -1
	BLOSS["TD"] = -1
	BLOSS["DW"] = -4
	BLOSS["WD"] = -4
	BLOSS["DY"] = -3
	BLOSS["YD"] = -3
	BLOSS["DV"] = -3
	BLOSS["VD"] = -3
	BLOSS["DB"] = 4
	BLOSS["BD"] = 4
	BLOSS["DZ"] = 1
	BLOSS["ZD"] = 1
	BLOSS["DX"] = -1
	BLOSS["XD"] = -1
	BLOSS["CA"] = 0
	BLOSS["AC"] = 0
	BLOSS["CR"] = -3
	BLOSS["RC"] = -3
	BLOSS["CN"] = -3
	BLOSS["NC"] = -3
	BLOSS["CD"] = -3
	BLOSS["DC"] = -3
	BLOSS["CQ"] = -3
	BLOSS["QC"] = -3
	BLOSS["CE"] = -4
	BLOSS["EC"] = -4
	BLOSS["CG"] = -3
	BLOSS["GC"] = -3
	BLOSS["CH"] = -3
	BLOSS["HC"] = -3
	BLOSS["CI"] = -1
	BLOSS["IC"] = -1
	BLOSS["CL"] = -1
	BLOSS["LC"] = -1
	BLOSS["CK"] = -3
	BLOSS["KC"] = -3
	BLOSS["CM"] = -1
	BLOSS["MC"] = -1
	BLOSS["CF"] = -2
	BLOSS["FC"] = -2
	BLOSS["CP"] = -3
	BLOSS["PC"] = -3
	BLOSS["CS"] = -1
	BLOSS["SC"] = -1
	BLOSS["CT"] = -1
	BLOSS["TC"] = -1
	BLOSS["CW"] = -2
	BLOSS["WC"] = -2
	BLOSS["CY"] = -2
	BLOSS["YC"] = -2
	BLOSS["CV"] = -1
	BLOSS["VC"] = -1
	BLOSS["CB"] = -3
	BLOSS["BC"] = -3
	BLOSS["CZ"] = -3
	BLOSS["ZC"] = -3
	BLOSS["CX"] = -2
	BLOSS["XC"] = -2
	BLOSS["QA"] = -1
	BLOSS["AQ"] = -1
	BLOSS["QR"] = 1
	BLOSS["RQ"] = 1
	BLOSS["QN"] = 0
	BLOSS["NQ"] = 0
	BLOSS["QD"] = 0
	BLOSS["DQ"] = 0
	BLOSS["QC"] = -3
	BLOSS["CQ"] = -3
	BLOSS["QE"] = 2
	BLOSS["EQ"] = 2
	BLOSS["QG"] = -2
	BLOSS["GQ"] = -2
	BLOSS["QH"] = 0
	BLOSS["HQ"] = 0
	BLOSS["QI"] = -3
	BLOSS["IQ"] = -3
	BLOSS["QL"] = -2
	BLOSS["LQ"] = -2
	BLOSS["QK"] = 1
	BLOSS["KQ"] = 1
	BLOSS["QM"] = 0
	BLOSS["MQ"] = 0
	BLOSS["QF"] = -3
	BLOSS["FQ"] = -3
	BLOSS["QP"] = -1
	BLOSS["PQ"] = -1
	BLOSS["QS"] = 0
	BLOSS["SQ"] = 0
	BLOSS["QT"] = -1
	BLOSS["TQ"] = -1
	BLOSS["QW"] = -2
	BLOSS["WQ"] = -2
	BLOSS["QY"] = -1
	BLOSS["YQ"] = -1
	BLOSS["QV"] = -2
	BLOSS["VQ"] = -2
	BLOSS["QB"] = 0
	BLOSS["BQ"] = 0
	BLOSS["QZ"] = 3
	BLOSS["ZQ"] = 3
	BLOSS["QX"] = -1
	BLOSS["XQ"] = -1
	BLOSS["EA"] = -1
	BLOSS["AE"] = -1
	BLOSS["ER"] = 0
	BLOSS["RE"] = 0
	BLOSS["EN"] = 0
	BLOSS["NE"] = 0
	BLOSS["ED"] = 2
	BLOSS["DE"] = 2
	BLOSS["EC"] = -4
	BLOSS["CE"] = -4
	BLOSS["EQ"] = 2
	BLOSS["QE"] = 2
	BLOSS["EG"] = -2
	BLOSS["GE"] = -2
	BLOSS["EH"] = 0
	BLOSS["HE"] = 0
	BLOSS["EI"] = -3
	BLOSS["IE"] = -3
	BLOSS["EL"] = -3
	BLOSS["LE"] = -3
	BLOSS["EK"] = 1
	BLOSS["KE"] = 1
	BLOSS["EM"] = -2
	BLOSS["ME"] = -2
	BLOSS["EF"] = -3
	BLOSS["FE"] = -3
	BLOSS["EP"] = -1
	BLOSS["PE"] = -1
	BLOSS["ES"] = 0
	BLOSS["SE"] = 0
	BLOSS["ET"] = -1
	BLOSS["TE"] = -1
	BLOSS["EW"] = -3
	BLOSS["WE"] = -3
	BLOSS["EY"] = -2
	BLOSS["YE"] = -2
	BLOSS["EV"] = -2
	BLOSS["VE"] = -2
	BLOSS["EB"] = 1
	BLOSS["BE"] = 1
	BLOSS["EZ"] = 4
	BLOSS["ZE"] = 4
	BLOSS["EX"] = -1
	BLOSS["XE"] = -1
	BLOSS["GA"] = 0
	BLOSS["AG"] = 0
	BLOSS["GR"] = -2
	BLOSS["RG"] = -2
	BLOSS["GN"] = 0
	BLOSS["NG"] = 0
	BLOSS["GD"] = -1
	BLOSS["DG"] = -1
	BLOSS["GC"] = -3
	BLOSS["CG"] = -3
	BLOSS["GQ"] = -2
	BLOSS["QG"] = -2
	BLOSS["GE"] = -2
	BLOSS["EG"] = -2
	BLOSS["GH"] = -2
	BLOSS["HG"] = -2
	BLOSS["GI"] = -4
	BLOSS["IG"] = -4
	BLOSS["GL"] = -4
	BLOSS["LG"] = -4
	BLOSS["GK"] = -2
	BLOSS["KG"] = -2
	BLOSS["GM"] = -3
	BLOSS["MG"] = -3
	BLOSS["GF"] = -3
	BLOSS["FG"] = -3
	BLOSS["GP"] = -2
	BLOSS["PG"] = -2
	BLOSS["GS"] = 0
	BLOSS["SG"] = 0
	BLOSS["GT"] = -2
	BLOSS["TG"] = -2
	BLOSS["GW"] = -2
	BLOSS["WG"] = -2
	BLOSS["GY"] = -3
	BLOSS["YG"] = -3
	BLOSS["GV"] = -3
	BLOSS["VG"] = -3
	BLOSS["GB"] = -1
	BLOSS["BG"] = -1
	BLOSS["GZ"] = -2
	BLOSS["ZG"] = -2
	BLOSS["GX"] = -1
	BLOSS["XG"] = -1
	BLOSS["HA"] = -2
	BLOSS["AH"] = -2
	BLOSS["HR"] = 0
	BLOSS["RH"] = 0
	BLOSS["HN"] = 1
	BLOSS["NH"] = 1
	BLOSS["HD"] = -1
	BLOSS["DH"] = -1
	BLOSS["HC"] = -3
	BLOSS["CH"] = -3
	BLOSS["HQ"] = 0
	BLOSS["QH"] = 0
	BLOSS["HE"] = 0
	BLOSS["EH"] = 0
	BLOSS["HG"] = -2
	BLOSS["GH"] = -2
	BLOSS["HI"] = -3
	BLOSS["IH"] = -3
	BLOSS["HL"] = -3
	BLOSS["LH"] = -3
	BLOSS["HK"] = -1
	BLOSS["KH"] = -1
	BLOSS["HM"] = -2
	BLOSS["MH"] = -2
	BLOSS["HF"] = -1
	BLOSS["FH"] = -1
	BLOSS["HP"] = -2
	BLOSS["PH"] = -2
	BLOSS["HS"] = -1
	BLOSS["SH"] = -1
	BLOSS["HT"] = -2
	BLOSS["TH"] = -2
	BLOSS["HW"] = -2
	BLOSS["WH"] = -2
	BLOSS["HY"] = 2
	BLOSS["YH"] = 2
	BLOSS["HV"] = -3
	BLOSS["VH"] = -3
	BLOSS["HB"] = 0
	BLOSS["BH"] = 0
	BLOSS["HZ"] = 0
	BLOSS["ZH"] = 0
	BLOSS["HX"] = -1
	BLOSS["XH"] = -1
	BLOSS["IA"] = -1
	BLOSS["AI"] = -1
	BLOSS["IR"] = -3
	BLOSS["RI"] = -3
	BLOSS["IN"] = -3
	BLOSS["NI"] = -3
	BLOSS["ID"] = -3
	BLOSS["DI"] = -3
	BLOSS["IC"] = -1
	BLOSS["CI"] = -1
	BLOSS["IQ"] = -3
	BLOSS["QI"] = -3
	BLOSS["IE"] = -3
	BLOSS["EI"] = -3
	BLOSS["IG"] = -4
	BLOSS["GI"] = -4
	BLOSS["IH"] = -3
	BLOSS["HI"] = -3
	BLOSS["IL"] = 2
	BLOSS["LI"] = 2
	BLOSS["IK"] = -3
	BLOSS["KI"] = -3
	BLOSS["IM"] = 1
	BLOSS["MI"] = 1
	BLOSS["IF"] = 0
	BLOSS["FI"] = 0
	BLOSS["IP"] = -3
	BLOSS["PI"] = -3
	BLOSS["IS"] = -2
	BLOSS["SI"] = -2
	BLOSS["IT"] = -1
	BLOSS["TI"] = -1
	BLOSS["IW"] = -3
	BLOSS["WI"] = -3
	BLOSS["IY"] = -1
	BLOSS["YI"] = -1
	BLOSS["IV"] = 3
	BLOSS["VI"] = 3
	BLOSS["IB"] = -3
	BLOSS["BI"] = -3
	BLOSS["IZ"] = -3
	BLOSS["ZI"] = -3
	BLOSS["IX"] = -1
	BLOSS["XI"] = -1
	BLOSS["LA"] = -1
	BLOSS["AL"] = -1
	BLOSS["LR"] = -2
	BLOSS["RL"] = -2
	BLOSS["LN"] = -3
	BLOSS["NL"] = -3
	BLOSS["LD"] = -4
	BLOSS["DL"] = -4
	BLOSS["LC"] = -1
	BLOSS["CL"] = -1
	BLOSS["LQ"] = -2
	BLOSS["QL"] = -2
	BLOSS["LE"] = -3
	BLOSS["EL"] = -3
	BLOSS["LG"] = -4
	BLOSS["GL"] = -4
	BLOSS["LH"] = -3
	BLOSS["HL"] = -3
	BLOSS["LI"] = 2
	BLOSS["IL"] = 2
	BLOSS["LK"] = -2
	BLOSS["KL"] = -2
	BLOSS["LM"] = 2
	BLOSS["ML"] = 2
	BLOSS["LF"] = 0
	BLOSS["FL"] = 0
	BLOSS["LP"] = -3
	BLOSS["PL"] = -3
	BLOSS["LS"] = -2
	BLOSS["SL"] = -2
	BLOSS["LT"] = -1
	BLOSS["TL"] = -1
	BLOSS["LW"] = -2
	BLOSS["WL"] = -2
	BLOSS["LY"] = -1
	BLOSS["YL"] = -1
	BLOSS["LV"] = 1
	BLOSS["VL"] = 1
	BLOSS["LB"] = -4
	BLOSS["BL"] = -4
	BLOSS["LZ"] = -3
	BLOSS["ZL"] = -3
	BLOSS["LX"] = -1
	BLOSS["XL"] = -1
	BLOSS["KA"] = -1
	BLOSS["AK"] = -1
	BLOSS["KR"] = 2
	BLOSS["RK"] = 2
	BLOSS["KN"] = 0
	BLOSS["NK"] = 0
	BLOSS["KD"] = -1
	BLOSS["DK"] = -1
	BLOSS["KC"] = -3
	BLOSS["CK"] = -3
	BLOSS["KQ"] = 1
	BLOSS["QK"] = 1
	BLOSS["KE"] = 1
	BLOSS["EK"] = 1
	BLOSS["KG"] = -2
	BLOSS["GK"] = -2
	BLOSS["KH"] = -1
	BLOSS["HK"] = -1
	BLOSS["KI"] = -3
	BLOSS["IK"] = -3
	BLOSS["KL"] = -2
	BLOSS["LK"] = -2
	BLOSS["KM"] = -1
	BLOSS["MK"] = -1
	BLOSS["KF"] = -3
	BLOSS["FK"] = -3
	BLOSS["KP"] = -1
	BLOSS["PK"] = -1
	BLOSS["KS"] = 0
	BLOSS["SK"] = 0
	BLOSS["KT"] = -1
	BLOSS["TK"] = -1
	BLOSS["KW"] = -3
	BLOSS["WK"] = -3
	BLOSS["KY"] = -2
	BLOSS["YK"] = -2
	BLOSS["KV"] = -2
	BLOSS["VK"] = -2
	BLOSS["KB"] = 0
	BLOSS["BK"] = 0
	BLOSS["KZ"] = 1
	BLOSS["ZK"] = 1
	BLOSS["KX"] = -1
	BLOSS["XK"] = -1
	BLOSS["MA"] = -1
	BLOSS["AM"] = -1
	BLOSS["MR"] = -1
	BLOSS["RM"] = -1
	BLOSS["MN"] = -2
	BLOSS["NM"] = -2
	BLOSS["MD"] = -3
	BLOSS["DM"] = -3
	BLOSS["MC"] = -1
	BLOSS["CM"] = -1
	BLOSS["MQ"] = 0
	BLOSS["QM"] = 0
	BLOSS["ME"] = -2
	BLOSS["EM"] = -2
	BLOSS["MG"] = -3
	BLOSS["GM"] = -3
	BLOSS["MH"] = -2
	BLOSS["HM"] = -2
	BLOSS["MI"] = 1
	BLOSS["IM"] = 1
	BLOSS["ML"] = 2
	BLOSS["LM"] = 2
	BLOSS["MK"] = -1
	BLOSS["KM"] = -1
	BLOSS["MF"] = 0
	BLOSS["FM"] = 0
	BLOSS["MP"] = -2
	BLOSS["PM"] = -2
	BLOSS["MS"] = -1
	BLOSS["SM"] = -1
	BLOSS["MT"] = -1
	BLOSS["TM"] = -1
	BLOSS["MW"] = -1
	BLOSS["WM"] = -1
	BLOSS["MY"] = -1
	BLOSS["YM"] = -1
	BLOSS["MV"] = 1
	BLOSS["VM"] = 1
	BLOSS["MB"] = -3
	BLOSS["BM"] = -3
	BLOSS["MZ"] = -1
	BLOSS["ZM"] = -1
	BLOSS["MX"] = -1
	BLOSS["XM"] = -1
	BLOSS["FA"] = -2
	BLOSS["AF"] = -2
	BLOSS["FR"] = -3
	BLOSS["RF"] = -3
	BLOSS["FN"] = -3
	BLOSS["NF"] = -3
	BLOSS["FD"] = -3
	BLOSS["DF"] = -3
	BLOSS["FC"] = -2
	BLOSS["CF"] = -2
	BLOSS["FQ"] = -3
	BLOSS["QF"] = -3
	BLOSS["FE"] = -3
	BLOSS["EF"] = -3
	BLOSS["FG"] = -3
	BLOSS["GF"] = -3
	BLOSS["FH"] = -1
	BLOSS["HF"] = -1
	BLOSS["FI"] = 0
	BLOSS["IF"] = 0
	BLOSS["FL"] = 0
	BLOSS["LF"] = 0
	BLOSS["FK"] = -3
	BLOSS["KF"] = -3
	BLOSS["FM"] = 0
	BLOSS["MF"] = 0
	BLOSS["FP"] = -4
	BLOSS["PF"] = -4
	BLOSS["FS"] = -2
	BLOSS["SF"] = -2
	BLOSS["FT"] = -2
	BLOSS["TF"] = -2
	BLOSS["FW"] = 1
	BLOSS["WF"] = 1
	BLOSS["FY"] = 3
	BLOSS["YF"] = 3
	BLOSS["FV"] = -1
	BLOSS["VF"] = -1
	BLOSS["FB"] = -3
	BLOSS["BF"] = -3
	BLOSS["FZ"] = -3
	BLOSS["ZF"] = -3
	BLOSS["FX"] = -1
	BLOSS["XF"] = -1
	BLOSS["PA"] = -1
	BLOSS["AP"] = -1
	BLOSS["PR"] = -2
	BLOSS["RP"] = -2
	BLOSS["PN"] = -2
	BLOSS["NP"] = -2
	BLOSS["PD"] = -1
	BLOSS["DP"] = -1
	BLOSS["PC"] = -3
	BLOSS["CP"] = -3
	BLOSS["PQ"] = -1
	BLOSS["QP"] = -1
	BLOSS["PE"] = -1
	BLOSS["EP"] = -1
	BLOSS["PG"] = -2
	BLOSS["GP"] = -2
	BLOSS["PH"] = -2
	BLOSS["HP"] = -2
	BLOSS["PI"] = -3
	BLOSS["IP"] = -3
	BLOSS["PL"] = -3
	BLOSS["LP"] = -3
	BLOSS["PK"] = -1
	BLOSS["KP"] = -1
	BLOSS["PM"] = -2
	BLOSS["MP"] = -2
	BLOSS["PF"] = -4
	BLOSS["FP"] = -4
	BLOSS["PS"] = -1
	BLOSS["SP"] = -1
	BLOSS["PT"] = -1
	BLOSS["TP"] = -1
	BLOSS["PW"] = -4
	BLOSS["WP"] = -4
	BLOSS["PY"] = -3
	BLOSS["YP"] = -3
	BLOSS["PV"] = -2
	BLOSS["VP"] = -2
	BLOSS["PB"] = -2
	BLOSS["BP"] = -2
	BLOSS["PZ"] = -1
	BLOSS["ZP"] = -1
	BLOSS["PX"] = -2
	BLOSS["XP"] = -2
	BLOSS["SA"] = 1
	BLOSS["AS"] = 1
	BLOSS["SR"] = -1
	BLOSS["RS"] = -1
	BLOSS["SN"] = 1
	BLOSS["NS"] = 1
	BLOSS["SD"] = 0
	BLOSS["DS"] = 0
	BLOSS["SC"] = -1
	BLOSS["CS"] = -1
	BLOSS["SQ"] = 0
	BLOSS["QS"] = 0
	BLOSS["SE"] = 0
	BLOSS["ES"] = 0
	BLOSS["SG"] = 0
	BLOSS["GS"] = 0
	BLOSS["SH"] = -1
	BLOSS["HS"] = -1
	BLOSS["SI"] = -2
	BLOSS["IS"] = -2
	BLOSS["SL"] = -2
	BLOSS["LS"] = -2
	BLOSS["SK"] = 0
	BLOSS["KS"] = 0
	BLOSS["SM"] = -1
	BLOSS["MS"] = -1
	BLOSS["SF"] = -2
	BLOSS["FS"] = -2
	BLOSS["SP"] = -1
	BLOSS["PS"] = -1
	BLOSS["ST"] = 1
	BLOSS["TS"] = 1
	BLOSS["SW"] = -3
	BLOSS["WS"] = -3
	BLOSS["SY"] = -2
	BLOSS["YS"] = -2
	BLOSS["SV"] = -2
	BLOSS["VS"] = -2
	BLOSS["SB"] = 0
	BLOSS["BS"] = 0
	BLOSS["SZ"] = 0
	BLOSS["ZS"] = 0
	BLOSS["SX"] = 0
	BLOSS["XS"] = 0
	BLOSS["TA"] = 0
	BLOSS["AT"] = 0
	BLOSS["TR"] = -1
	BLOSS["RT"] = -1
	BLOSS["TN"] = 0
	BLOSS["NT"] = 0
	BLOSS["TD"] = -1
	BLOSS["DT"] = -1
	BLOSS["TC"] = -1
	BLOSS["CT"] = -1
	BLOSS["TQ"] = -1
	BLOSS["QT"] = -1
	BLOSS["TE"] = -1
	BLOSS["ET"] = -1
	BLOSS["TG"] = -2
	BLOSS["GT"] = -2
	BLOSS["TH"] = -2
	BLOSS["HT"] = -2
	BLOSS["TI"] = -1
	BLOSS["IT"] = -1
	BLOSS["TL"] = -1
	BLOSS["LT"] = -1
	BLOSS["TK"] = -1
	BLOSS["KT"] = -1
	BLOSS["TM"] = -1
	BLOSS["MT"] = -1
	BLOSS["TF"] = -2
	BLOSS["FT"] = -2
	BLOSS["TP"] = -1
	BLOSS["PT"] = -1
	BLOSS["TS"] = 1
	BLOSS["ST"] = 1
	BLOSS["TW"] = -2
	BLOSS["WT"] = -2
	BLOSS["TY"] = -2
	BLOSS["YT"] = -2
	BLOSS["TV"] = 0
	BLOSS["VT"] = 0
	BLOSS["TB"] = -1
	BLOSS["BT"] = -1
	BLOSS["TZ"] = -1
	BLOSS["ZT"] = -1
	BLOSS["TX"] = 0
	BLOSS["XT"] = 0
	BLOSS["WA"] = -3
	BLOSS["AW"] = -3
	BLOSS["WR"] = -3
	BLOSS["RW"] = -3
	BLOSS["WN"] = -4
	BLOSS["NW"] = -4
	BLOSS["WD"] = -4
	BLOSS["DW"] = -4
	BLOSS["WC"] = -2
	BLOSS["CW"] = -2
	BLOSS["WQ"] = -2
	BLOSS["QW"] = -2
	BLOSS["WE"] = -3
	BLOSS["EW"] = -3
	BLOSS["WG"] = -2
	BLOSS["GW"] = -2
	BLOSS["WH"] = -2
	BLOSS["HW"] = -2
	BLOSS["WI"] = -3
	BLOSS["IW"] = -3
	BLOSS["WL"] = -2
	BLOSS["LW"] = -2
	BLOSS["WK"] = -3
	BLOSS["KW"] = -3
	BLOSS["WM"] = -1
	BLOSS["MW"] = -1
	BLOSS["WF"] = 1
	BLOSS["FW"] = 1
	BLOSS["WP"] = -4
	BLOSS["PW"] = -4
	BLOSS["WS"] = -3
	BLOSS["SW"] = -3
	BLOSS["WT"] = -2
	BLOSS["TW"] = -2
	BLOSS["WY"] = 2
	BLOSS["YW"] = 2
	BLOSS["WV"] = -3
	BLOSS["VW"] = -3
	BLOSS["WB"] = -4
	BLOSS["BW"] = -4
	BLOSS["WZ"] = -3
	BLOSS["ZW"] = -3
	BLOSS["WX"] = -2
	BLOSS["XW"] = -2
	BLOSS["YA"] = -2
	BLOSS["AY"] = -2
	BLOSS["YR"] = -2
	BLOSS["RY"] = -2
	BLOSS["YN"] = -2
	BLOSS["NY"] = -2
	BLOSS["YD"] = -3
	BLOSS["DY"] = -3
	BLOSS["YC"] = -2
	BLOSS["CY"] = -2
	BLOSS["YQ"] = -1
	BLOSS["QY"] = -1
	BLOSS["YE"] = -2
	BLOSS["EY"] = -2
	BLOSS["YG"] = -3
	BLOSS["GY"] = -3
	BLOSS["YH"] = 2
	BLOSS["HY"] = 2
	BLOSS["YI"] = -1
	BLOSS["IY"] = -1
	BLOSS["YL"] = -1
	BLOSS["LY"] = -1
	BLOSS["YK"] = -2
	BLOSS["KY"] = -2
	BLOSS["YM"] = -1
	BLOSS["MY"] = -1
	BLOSS["YF"] = 3
	BLOSS["FY"] = 3
	BLOSS["YP"] = -3
	BLOSS["PY"] = -3
	BLOSS["YS"] = -2
	BLOSS["SY"] = -2
	BLOSS["YT"] = -2
	BLOSS["TY"] = -2
	BLOSS["YW"] = 2
	BLOSS["WY"] = 2
	BLOSS["YV"] = -1
	BLOSS["VY"] = -1
	BLOSS["YB"] = -3
	BLOSS["BY"] = -3
	BLOSS["YZ"] = -2
	BLOSS["ZY"] = -2
	BLOSS["YX"] = -1
	BLOSS["XY"] = -1
	BLOSS["VA"] = 0
	BLOSS["AV"] = 0
	BLOSS["VR"] = -3
	BLOSS["RV"] = -3
	BLOSS["VN"] = -3
	BLOSS["NV"] = -3
	BLOSS["VD"] = -3
	BLOSS["DV"] = -3
	BLOSS["VC"] = -1
	BLOSS["CV"] = -1
	BLOSS["VQ"] = -2
	BLOSS["QV"] = -2
	BLOSS["VE"] = -2
	BLOSS["EV"] = -2
	BLOSS["VG"] = -3
	BLOSS["GV"] = -3
	BLOSS["VH"] = -3
	BLOSS["HV"] = -3
	BLOSS["VI"] = 3
	BLOSS["IV"] = 3
	BLOSS["VL"] = 1
	BLOSS["LV"] = 1
	BLOSS["VK"] = -2
	BLOSS["KV"] = -2
	BLOSS["VM"] = 1
	BLOSS["MV"] = 1
	BLOSS["VF"] = -1
	BLOSS["FV"] = -1
	BLOSS["VP"] = -2
	BLOSS["PV"] = -2
	BLOSS["VS"] = -2
	BLOSS["SV"] = -2
	BLOSS["VT"] = 0
	BLOSS["TV"] = 0
	BLOSS["VW"] = -3
	BLOSS["WV"] = -3
	BLOSS["VY"] = -1
	BLOSS["YV"] = -1
	BLOSS["VB"] = -3
	BLOSS["BV"] = -3
	BLOSS["VZ"] = -2
	BLOSS["ZV"] = -2
	BLOSS["VX"] = -1
	BLOSS["XV"] = -1
	BLOSS["BA"] = -2
	BLOSS["AB"] = -2
	BLOSS["BR"] = -1
	BLOSS["RB"] = -1
	BLOSS["BN"] = 3
	BLOSS["NB"] = 3
	BLOSS["BD"] = 4
	BLOSS["DB"] = 4
	BLOSS["BC"] = -3
	BLOSS["CB"] = -3
	BLOSS["BQ"] = 0
	BLOSS["QB"] = 0
	BLOSS["BE"] = 1
	BLOSS["EB"] = 1
	BLOSS["BG"] = -1
	BLOSS["GB"] = -1
	BLOSS["BH"] = 0
	BLOSS["HB"] = 0
	BLOSS["BI"] = -3
	BLOSS["IB"] = -3
	BLOSS["BL"] = -4
	BLOSS["LB"] = -4
	BLOSS["BK"] = 0
	BLOSS["KB"] = 0
	BLOSS["BM"] = -3
	BLOSS["MB"] = -3
	BLOSS["BF"] = -3
	BLOSS["FB"] = -3
	BLOSS["BP"] = -2
	BLOSS["PB"] = -2
	BLOSS["BS"] = 0
	BLOSS["SB"] = 0
	BLOSS["BT"] = -1
	BLOSS["TB"] = -1
	BLOSS["BW"] = -4
	BLOSS["WB"] = -4
	BLOSS["BY"] = -3
	BLOSS["YB"] = -3
	BLOSS["BV"] = -3
	BLOSS["VB"] = -3
	BLOSS["BZ"] = 1
	BLOSS["ZB"] = 1
	BLOSS["BX"] = -1
	BLOSS["XB"] = -1
	BLOSS["ZA"] = -1
	BLOSS["AZ"] = -1
	BLOSS["ZR"] = 0
	BLOSS["RZ"] = 0
	BLOSS["ZN"] = 0
	BLOSS["NZ"] = 0
	BLOSS["ZD"] = 1
	BLOSS["DZ"] = 1
	BLOSS["ZC"] = -3
	BLOSS["CZ"] = -3
	BLOSS["ZQ"] = 3
	BLOSS["QZ"] = 3
	BLOSS["ZE"] = 4
	BLOSS["EZ"] = 4
	BLOSS["ZG"] = -2
	BLOSS["GZ"] = -2
	BLOSS["ZH"] = 0
	BLOSS["HZ"] = 0
	BLOSS["ZI"] = -3
	BLOSS["IZ"] = -3
	BLOSS["ZL"] = -3
	BLOSS["LZ"] = -3
	BLOSS["ZK"] = 1
	BLOSS["KZ"] = 1
	BLOSS["ZM"] = -1
	BLOSS["MZ"] = -1
	BLOSS["ZF"] = -3
	BLOSS["FZ"] = -3
	BLOSS["ZP"] = -1
	BLOSS["PZ"] = -1
	BLOSS["ZS"] = 0
	BLOSS["SZ"] = 0
	BLOSS["ZT"] = -1
	BLOSS["TZ"] = -1
	BLOSS["ZW"] = -3
	BLOSS["WZ"] = -3
	BLOSS["ZY"] = -2
	BLOSS["YZ"] = -2
	BLOSS["ZV"] = -2
	BLOSS["VZ"] = -2
	BLOSS["ZB"] = 1
	BLOSS["BZ"] = 1
	BLOSS["ZX"] = -1
	BLOSS["XZ"] = -1
	BLOSS["XA"] = 0
	BLOSS["AX"] = 0
	BLOSS["XR"] = -1
	BLOSS["RX"] = -1
	BLOSS["XN"] = -1
	BLOSS["NX"] = -1
	BLOSS["XD"] = -1
	BLOSS["DX"] = -1
	BLOSS["XC"] = -2
	BLOSS["CX"] = -2
	BLOSS["XQ"] = -1
	BLOSS["QX"] = -1
	BLOSS["XE"] = -1
	BLOSS["EX"] = -1
	BLOSS["XG"] = -1
	BLOSS["GX"] = -1
	BLOSS["XH"] = -1
	BLOSS["HX"] = -1
	BLOSS["XI"] = -1
	BLOSS["IX"] = -1
	BLOSS["XL"] = -1
	BLOSS["LX"] = -1
	BLOSS["XK"] = -1
	BLOSS["KX"] = -1
	BLOSS["XM"] = -1
	BLOSS["MX"] = -1
	BLOSS["XF"] = -1
	BLOSS["FX"] = -1
	BLOSS["XP"] = -2
	BLOSS["PX"] = -2
	BLOSS["XS"] = 0
	BLOSS["SX"] = 0
	BLOSS["XT"] = 0
	BLOSS["TX"] = 0
	BLOSS["XW"] = -2
	BLOSS["WX"] = -2
	BLOSS["XY"] = -1
	BLOSS["YX"] = -1
	BLOSS["XV"] = -1
	BLOSS["VX"] = -1
	BLOSS["XB"] = -1
	BLOSS["BX"] = -1
	BLOSS["XZ"] = -1
	BLOSS["ZX"] = -1

	// fmt.Printf("%v %v %v\n", refAA, altAA, Tang[buffer.String()])
	if altAA != "" && refAA != "" && altAA != "X" {
		res = BLOSS[buffer.String()]
	} else {
		res = 0
	}
	return res
}

func GetPAM30Inx(refAA string, altAA string) int {
	var buffer bytes.Buffer
	var res int
	buffer.WriteString(refAA)
	buffer.WriteString(altAA)
	PAM30 := make(map[string]int)

	PAM30["AR"] = -7
	PAM30["RA"] = -7
	PAM30["AN"] = -4
	PAM30["NA"] = -4
	PAM30["AD"] = -3
	PAM30["DA"] = -3
	PAM30["AC"] = -6
	PAM30["CA"] = -6
	PAM30["AQ"] = -4
	PAM30["QA"] = -4
	PAM30["AE"] = -2
	PAM30["EA"] = -2
	PAM30["AG"] = -2
	PAM30["GA"] = -2
	PAM30["AH"] = -7
	PAM30["HA"] = -7
	PAM30["AI"] = -5
	PAM30["IA"] = -5
	PAM30["AL"] = -6
	PAM30["LA"] = -6
	PAM30["AK"] = -7
	PAM30["KA"] = -7
	PAM30["AM"] = -5
	PAM30["MA"] = -5
	PAM30["AF"] = -8
	PAM30["FA"] = -8
	PAM30["AP"] = -2
	PAM30["PA"] = -2
	PAM30["AS"] = 0
	PAM30["SA"] = 0
	PAM30["AT"] = -1
	PAM30["TA"] = -1
	PAM30["AW"] = -13
	PAM30["WA"] = -13
	PAM30["AY"] = -8
	PAM30["YA"] = -8
	PAM30["AV"] = -2
	PAM30["VA"] = -2
	PAM30["AB"] = -3
	PAM30["BA"] = -3
	PAM30["AZ"] = -3
	PAM30["ZA"] = -3
	PAM30["AX"] = -3
	PAM30["XA"] = -3
	PAM30["RA"] = -7
	PAM30["AR"] = -7
	PAM30["RN"] = -6
	PAM30["NR"] = -6
	PAM30["RD"] = -10
	PAM30["DR"] = -10
	PAM30["RC"] = -8
	PAM30["CR"] = -8
	PAM30["RQ"] = -2
	PAM30["QR"] = -2
	PAM30["RE"] = -9
	PAM30["ER"] = -9
	PAM30["RG"] = -9
	PAM30["GR"] = -9
	PAM30["RH"] = -2
	PAM30["HR"] = -2
	PAM30["RI"] = -5
	PAM30["IR"] = -5
	PAM30["RL"] = -8
	PAM30["LR"] = -8
	PAM30["RK"] = 0
	PAM30["KR"] = 0
	PAM30["RM"] = -4
	PAM30["MR"] = -4
	PAM30["RF"] = -9
	PAM30["FR"] = -9
	PAM30["RP"] = -4
	PAM30["PR"] = -4
	PAM30["RS"] = -3
	PAM30["SR"] = -3
	PAM30["RT"] = -6
	PAM30["TR"] = -6
	PAM30["RW"] = -2
	PAM30["WR"] = -2
	PAM30["RY"] = -10
	PAM30["YR"] = -10
	PAM30["RV"] = -8
	PAM30["VR"] = -8
	PAM30["RB"] = -7
	PAM30["BR"] = -7
	PAM30["RZ"] = -4
	PAM30["ZR"] = -4
	PAM30["RX"] = -6
	PAM30["XR"] = -6
	PAM30["NA"] = -4
	PAM30["AN"] = -4
	PAM30["NR"] = -6
	PAM30["RN"] = -6
	PAM30["ND"] = 2
	PAM30["DN"] = 2
	PAM30["NC"] = -11
	PAM30["CN"] = -11
	PAM30["NQ"] = -3
	PAM30["QN"] = -3
	PAM30["NE"] = -2
	PAM30["EN"] = -2
	PAM30["NG"] = -3
	PAM30["GN"] = -3
	PAM30["NH"] = 0
	PAM30["HN"] = 0
	PAM30["NI"] = -5
	PAM30["IN"] = -5
	PAM30["NL"] = -7
	PAM30["LN"] = -7
	PAM30["NK"] = -1
	PAM30["KN"] = -1
	PAM30["NM"] = -9
	PAM30["MN"] = -9
	PAM30["NF"] = -9
	PAM30["FN"] = -9
	PAM30["NP"] = -6
	PAM30["PN"] = -6
	PAM30["NS"] = 0
	PAM30["SN"] = 0
	PAM30["NT"] = -2
	PAM30["TN"] = -2
	PAM30["NW"] = -8
	PAM30["WN"] = -8
	PAM30["NY"] = -4
	PAM30["YN"] = -4
	PAM30["NV"] = -8
	PAM30["VN"] = -8
	PAM30["NB"] = 6
	PAM30["BN"] = 6
	PAM30["NZ"] = -3
	PAM30["ZN"] = -3
	PAM30["NX"] = -3
	PAM30["XN"] = -3
	PAM30["DA"] = -3
	PAM30["AD"] = -3
	PAM30["DR"] = -10
	PAM30["RD"] = -10
	PAM30["DN"] = 2
	PAM30["ND"] = 2
	PAM30["DC"] = -14
	PAM30["CD"] = -14
	PAM30["DQ"] = -2
	PAM30["QD"] = -2
	PAM30["DE"] = 2
	PAM30["ED"] = 2
	PAM30["DG"] = -3
	PAM30["GD"] = -3
	PAM30["DH"] = -4
	PAM30["HD"] = -4
	PAM30["DI"] = -7
	PAM30["ID"] = -7
	PAM30["DL"] = -12
	PAM30["LD"] = -12
	PAM30["DK"] = -4
	PAM30["KD"] = -4
	PAM30["DM"] = -11
	PAM30["MD"] = -11
	PAM30["DF"] = -15
	PAM30["FD"] = -15
	PAM30["DP"] = -8
	PAM30["PD"] = -8
	PAM30["DS"] = -4
	PAM30["SD"] = -4
	PAM30["DT"] = -5
	PAM30["TD"] = -5
	PAM30["DW"] = -15
	PAM30["WD"] = -15
	PAM30["DY"] = -11
	PAM30["YD"] = -11
	PAM30["DV"] = -8
	PAM30["VD"] = -8
	PAM30["DB"] = 6
	PAM30["BD"] = 6
	PAM30["DZ"] = 1
	PAM30["ZD"] = 1
	PAM30["DX"] = -5
	PAM30["XD"] = -5
	PAM30["CA"] = -6
	PAM30["AC"] = -6
	PAM30["CR"] = -8
	PAM30["RC"] = -8
	PAM30["CN"] = -11
	PAM30["NC"] = -11
	PAM30["CD"] = -14
	PAM30["DC"] = -14
	PAM30["CQ"] = -14
	PAM30["QC"] = -14
	PAM30["CE"] = -14
	PAM30["EC"] = -14
	PAM30["CG"] = -9
	PAM30["GC"] = -9
	PAM30["CH"] = -7
	PAM30["HC"] = -7
	PAM30["CI"] = -6
	PAM30["IC"] = -6
	PAM30["CL"] = -15
	PAM30["LC"] = -15
	PAM30["CK"] = -14
	PAM30["KC"] = -14
	PAM30["CM"] = -13
	PAM30["MC"] = -13
	PAM30["CF"] = -13
	PAM30["FC"] = -13
	PAM30["CP"] = -8
	PAM30["PC"] = -8
	PAM30["CS"] = -3
	PAM30["SC"] = -3
	PAM30["CT"] = -8
	PAM30["TC"] = -8
	PAM30["CW"] = -15
	PAM30["WC"] = -15
	PAM30["CY"] = -4
	PAM30["YC"] = -4
	PAM30["CV"] = -6
	PAM30["VC"] = -6
	PAM30["CB"] = -12
	PAM30["BC"] = -12
	PAM30["CZ"] = -14
	PAM30["ZC"] = -14
	PAM30["CX"] = -9
	PAM30["XC"] = -9
	PAM30["QA"] = -4
	PAM30["AQ"] = -4
	PAM30["QR"] = -2
	PAM30["RQ"] = -2
	PAM30["QN"] = -3
	PAM30["NQ"] = -3
	PAM30["QD"] = -2
	PAM30["DQ"] = -2
	PAM30["QC"] = -14
	PAM30["CQ"] = -14
	PAM30["QE"] = 1
	PAM30["EQ"] = 1
	PAM30["QG"] = -7
	PAM30["GQ"] = -7
	PAM30["QH"] = 1
	PAM30["HQ"] = 1
	PAM30["QI"] = -8
	PAM30["IQ"] = -8
	PAM30["QL"] = -5
	PAM30["LQ"] = -5
	PAM30["QK"] = -3
	PAM30["KQ"] = -3
	PAM30["QM"] = -4
	PAM30["MQ"] = -4
	PAM30["QF"] = -13
	PAM30["FQ"] = -13
	PAM30["QP"] = -3
	PAM30["PQ"] = -3
	PAM30["QS"] = -5
	PAM30["SQ"] = -5
	PAM30["QT"] = -5
	PAM30["TQ"] = -5
	PAM30["QW"] = -13
	PAM30["WQ"] = -13
	PAM30["QY"] = -12
	PAM30["YQ"] = -12
	PAM30["QV"] = -7
	PAM30["VQ"] = -7
	PAM30["QB"] = -3
	PAM30["BQ"] = -3
	PAM30["QZ"] = 6
	PAM30["ZQ"] = 6
	PAM30["QX"] = -5
	PAM30["XQ"] = -5
	PAM30["EA"] = -2
	PAM30["AE"] = -2
	PAM30["ER"] = -9
	PAM30["RE"] = -9
	PAM30["EN"] = -2
	PAM30["NE"] = -2
	PAM30["ED"] = 2
	PAM30["DE"] = 2
	PAM30["EC"] = -14
	PAM30["CE"] = -14
	PAM30["EQ"] = 1
	PAM30["QE"] = 1
	PAM30["EG"] = -4
	PAM30["GE"] = -4
	PAM30["EH"] = -5
	PAM30["HE"] = -5
	PAM30["EI"] = -5
	PAM30["IE"] = -5
	PAM30["EL"] = -9
	PAM30["LE"] = -9
	PAM30["EK"] = -4
	PAM30["KE"] = -4
	PAM30["EM"] = -7
	PAM30["ME"] = -7
	PAM30["EF"] = -14
	PAM30["FE"] = -14
	PAM30["EP"] = -5
	PAM30["PE"] = -5
	PAM30["ES"] = -4
	PAM30["SE"] = -4
	PAM30["ET"] = -6
	PAM30["TE"] = -6
	PAM30["EW"] = -17
	PAM30["WE"] = -17
	PAM30["EY"] = -8
	PAM30["YE"] = -8
	PAM30["EV"] = -6
	PAM30["VE"] = -6
	PAM30["EB"] = 1
	PAM30["BE"] = 1
	PAM30["EZ"] = 6
	PAM30["ZE"] = 6
	PAM30["EX"] = -5
	PAM30["XE"] = -5
	PAM30["GA"] = -2
	PAM30["AG"] = -2
	PAM30["GR"] = -9
	PAM30["RG"] = -9
	PAM30["GN"] = -3
	PAM30["NG"] = -3
	PAM30["GD"] = -3
	PAM30["DG"] = -3
	PAM30["GC"] = -9
	PAM30["CG"] = -9
	PAM30["GQ"] = -7
	PAM30["QG"] = -7
	PAM30["GE"] = -4
	PAM30["EG"] = -4
	PAM30["GH"] = -9
	PAM30["HG"] = -9
	PAM30["GI"] = -11
	PAM30["IG"] = -11
	PAM30["GL"] = -10
	PAM30["LG"] = -10
	PAM30["GK"] = -7
	PAM30["KG"] = -7
	PAM30["GM"] = -8
	PAM30["MG"] = -8
	PAM30["GF"] = -9
	PAM30["FG"] = -9
	PAM30["GP"] = -6
	PAM30["PG"] = -6
	PAM30["GS"] = -2
	PAM30["SG"] = -2
	PAM30["GT"] = -6
	PAM30["TG"] = -6
	PAM30["GW"] = -15
	PAM30["WG"] = -15
	PAM30["GY"] = -14
	PAM30["YG"] = -14
	PAM30["GV"] = -5
	PAM30["VG"] = -5
	PAM30["GB"] = -3
	PAM30["BG"] = -3
	PAM30["GZ"] = -5
	PAM30["ZG"] = -5
	PAM30["GX"] = -5
	PAM30["XG"] = -5
	PAM30["HA"] = -7
	PAM30["AH"] = -7
	PAM30["HR"] = -2
	PAM30["RH"] = -2
	PAM30["HN"] = 0
	PAM30["NH"] = 0
	PAM30["HD"] = -4
	PAM30["DH"] = -4
	PAM30["HC"] = -7
	PAM30["CH"] = -7
	PAM30["HQ"] = 1
	PAM30["QH"] = 1
	PAM30["HE"] = -5
	PAM30["EH"] = -5
	PAM30["HG"] = -9
	PAM30["GH"] = -9
	PAM30["HI"] = -9
	PAM30["IH"] = -9
	PAM30["HL"] = -6
	PAM30["LH"] = -6
	PAM30["HK"] = -6
	PAM30["KH"] = -6
	PAM30["HM"] = -10
	PAM30["MH"] = -10
	PAM30["HF"] = -6
	PAM30["FH"] = -6
	PAM30["HP"] = -4
	PAM30["PH"] = -4
	PAM30["HS"] = -6
	PAM30["SH"] = -6
	PAM30["HT"] = -7
	PAM30["TH"] = -7
	PAM30["HW"] = -7
	PAM30["WH"] = -7
	PAM30["HY"] = -3
	PAM30["YH"] = -3
	PAM30["HV"] = -6
	PAM30["VH"] = -6
	PAM30["HB"] = -1
	PAM30["BH"] = -1
	PAM30["HZ"] = -1
	PAM30["ZH"] = -1
	PAM30["HX"] = -5
	PAM30["XH"] = -5
	PAM30["IA"] = -5
	PAM30["AI"] = -5
	PAM30["IR"] = -5
	PAM30["RI"] = -5
	PAM30["IN"] = -5
	PAM30["NI"] = -5
	PAM30["ID"] = -7
	PAM30["DI"] = -7
	PAM30["IC"] = -6
	PAM30["CI"] = -6
	PAM30["IQ"] = -8
	PAM30["QI"] = -8
	PAM30["IE"] = -5
	PAM30["EI"] = -5
	PAM30["IG"] = -11
	PAM30["GI"] = -11
	PAM30["IH"] = -9
	PAM30["HI"] = -9
	PAM30["IL"] = -1
	PAM30["LI"] = -1
	PAM30["IK"] = -6
	PAM30["KI"] = -6
	PAM30["IM"] = -1
	PAM30["MI"] = -1
	PAM30["IF"] = -2
	PAM30["FI"] = -2
	PAM30["IP"] = -8
	PAM30["PI"] = -8
	PAM30["IS"] = -7
	PAM30["SI"] = -7
	PAM30["IT"] = -2
	PAM30["TI"] = -2
	PAM30["IW"] = -14
	PAM30["WI"] = -14
	PAM30["IY"] = -6
	PAM30["YI"] = -6
	PAM30["IV"] = 2
	PAM30["VI"] = 2
	PAM30["IB"] = -6
	PAM30["BI"] = -6
	PAM30["IZ"] = -6
	PAM30["ZI"] = -6
	PAM30["IX"] = -5
	PAM30["XI"] = -5
	PAM30["LA"] = -6
	PAM30["AL"] = -6
	PAM30["LR"] = -8
	PAM30["RL"] = -8
	PAM30["LN"] = -7
	PAM30["NL"] = -7
	PAM30["LD"] = -12
	PAM30["DL"] = -12
	PAM30["LC"] = -15
	PAM30["CL"] = -15
	PAM30["LQ"] = -5
	PAM30["QL"] = -5
	PAM30["LE"] = -9
	PAM30["EL"] = -9
	PAM30["LG"] = -10
	PAM30["GL"] = -10
	PAM30["LH"] = -6
	PAM30["HL"] = -6
	PAM30["LI"] = -1
	PAM30["IL"] = -1
	PAM30["LK"] = -8
	PAM30["KL"] = -8
	PAM30["LM"] = 1
	PAM30["ML"] = 1
	PAM30["LF"] = -3
	PAM30["FL"] = -3
	PAM30["LP"] = -7
	PAM30["PL"] = -7
	PAM30["LS"] = -8
	PAM30["SL"] = -8
	PAM30["LT"] = -7
	PAM30["TL"] = -7
	PAM30["LW"] = -6
	PAM30["WL"] = -6
	PAM30["LY"] = -7
	PAM30["YL"] = -7
	PAM30["LV"] = -2
	PAM30["VL"] = -2
	PAM30["LB"] = -9
	PAM30["BL"] = -9
	PAM30["LZ"] = -7
	PAM30["ZL"] = -7
	PAM30["LX"] = -6
	PAM30["XL"] = -6
	PAM30["KA"] = -7
	PAM30["AK"] = -7
	PAM30["KR"] = 0
	PAM30["RK"] = 0
	PAM30["KN"] = -1
	PAM30["NK"] = -1
	PAM30["KD"] = -4
	PAM30["DK"] = -4
	PAM30["KC"] = -14
	PAM30["CK"] = -14
	PAM30["KQ"] = -3
	PAM30["QK"] = -3
	PAM30["KE"] = -4
	PAM30["EK"] = -4
	PAM30["KG"] = -7
	PAM30["GK"] = -7
	PAM30["KH"] = -6
	PAM30["HK"] = -6
	PAM30["KI"] = -6
	PAM30["IK"] = -6
	PAM30["KL"] = -8
	PAM30["LK"] = -8
	PAM30["KM"] = -2
	PAM30["MK"] = -2
	PAM30["KF"] = -14
	PAM30["FK"] = -14
	PAM30["KP"] = -6
	PAM30["PK"] = -6
	PAM30["KS"] = -4
	PAM30["SK"] = -4
	PAM30["KT"] = -3
	PAM30["TK"] = -3
	PAM30["KW"] = -12
	PAM30["WK"] = -12
	PAM30["KY"] = -9
	PAM30["YK"] = -9
	PAM30["KV"] = -9
	PAM30["VK"] = -9
	PAM30["KB"] = -2
	PAM30["BK"] = -2
	PAM30["KZ"] = -4
	PAM30["ZK"] = -4
	PAM30["KX"] = -5
	PAM30["XK"] = -5
	PAM30["MA"] = -5
	PAM30["AM"] = -5
	PAM30["MR"] = -4
	PAM30["RM"] = -4
	PAM30["MN"] = -9
	PAM30["NM"] = -9
	PAM30["MD"] = -11
	PAM30["DM"] = -11
	PAM30["MC"] = -13
	PAM30["CM"] = -13
	PAM30["MQ"] = -4
	PAM30["QM"] = -4
	PAM30["ME"] = -7
	PAM30["EM"] = -7
	PAM30["MG"] = -8
	PAM30["GM"] = -8
	PAM30["MH"] = -10
	PAM30["HM"] = -10
	PAM30["MI"] = -1
	PAM30["IM"] = -1
	PAM30["ML"] = 1
	PAM30["LM"] = 1
	PAM30["MK"] = -2
	PAM30["KM"] = -2
	PAM30["MF"] = -4
	PAM30["FM"] = -4
	PAM30["MP"] = -8
	PAM30["PM"] = -8
	PAM30["MS"] = -5
	PAM30["SM"] = -5
	PAM30["MT"] = -4
	PAM30["TM"] = -4
	PAM30["MW"] = -13
	PAM30["WM"] = -13
	PAM30["MY"] = -11
	PAM30["YM"] = -11
	PAM30["MV"] = -1
	PAM30["VM"] = -1
	PAM30["MB"] = -10
	PAM30["BM"] = -10
	PAM30["MZ"] = -5
	PAM30["ZM"] = -5
	PAM30["MX"] = -5
	PAM30["XM"] = -5
	PAM30["FA"] = -8
	PAM30["AF"] = -8
	PAM30["FR"] = -9
	PAM30["RF"] = -9
	PAM30["FN"] = -9
	PAM30["NF"] = -9
	PAM30["FD"] = -15
	PAM30["DF"] = -15
	PAM30["FC"] = -13
	PAM30["CF"] = -13
	PAM30["FQ"] = -13
	PAM30["QF"] = -13
	PAM30["FE"] = -14
	PAM30["EF"] = -14
	PAM30["FG"] = -9
	PAM30["GF"] = -9
	PAM30["FH"] = -6
	PAM30["HF"] = -6
	PAM30["FI"] = -2
	PAM30["IF"] = -2
	PAM30["FL"] = -3
	PAM30["LF"] = -3
	PAM30["FK"] = -14
	PAM30["KF"] = -14
	PAM30["FM"] = -4
	PAM30["MF"] = -4
	PAM30["FP"] = -10
	PAM30["PF"] = -10
	PAM30["FS"] = -6
	PAM30["SF"] = -6
	PAM30["FT"] = -9
	PAM30["TF"] = -9
	PAM30["FW"] = -4
	PAM30["WF"] = -4
	PAM30["FY"] = 2
	PAM30["YF"] = 2
	PAM30["FV"] = -8
	PAM30["VF"] = -8
	PAM30["FB"] = -10
	PAM30["BF"] = -10
	PAM30["FZ"] = -13
	PAM30["ZF"] = -13
	PAM30["FX"] = -8
	PAM30["XF"] = -8
	PAM30["PA"] = -2
	PAM30["AP"] = -2
	PAM30["PR"] = -4
	PAM30["RP"] = -4
	PAM30["PN"] = -6
	PAM30["NP"] = -6
	PAM30["PD"] = -8
	PAM30["DP"] = -8
	PAM30["PC"] = -8
	PAM30["CP"] = -8
	PAM30["PQ"] = -3
	PAM30["QP"] = -3
	PAM30["PE"] = -5
	PAM30["EP"] = -5
	PAM30["PG"] = -6
	PAM30["GP"] = -6
	PAM30["PH"] = -4
	PAM30["HP"] = -4
	PAM30["PI"] = -8
	PAM30["IP"] = -8
	PAM30["PL"] = -7
	PAM30["LP"] = -7
	PAM30["PK"] = -6
	PAM30["KP"] = -6
	PAM30["PM"] = -8
	PAM30["MP"] = -8
	PAM30["PF"] = -10
	PAM30["FP"] = -10
	PAM30["PS"] = -2
	PAM30["SP"] = -2
	PAM30["PT"] = -4
	PAM30["TP"] = -4
	PAM30["PW"] = -14
	PAM30["WP"] = -14
	PAM30["PY"] = -13
	PAM30["YP"] = -13
	PAM30["PV"] = -6
	PAM30["VP"] = -6
	PAM30["PB"] = -7
	PAM30["BP"] = -7
	PAM30["PZ"] = -4
	PAM30["ZP"] = -4
	PAM30["PX"] = -5
	PAM30["XP"] = -5
	PAM30["SA"] = 0
	PAM30["AS"] = 0
	PAM30["SR"] = -3
	PAM30["RS"] = -3
	PAM30["SN"] = 0
	PAM30["NS"] = 0
	PAM30["SD"] = -4
	PAM30["DS"] = -4
	PAM30["SC"] = -3
	PAM30["CS"] = -3
	PAM30["SQ"] = -5
	PAM30["QS"] = -5
	PAM30["SE"] = -4
	PAM30["ES"] = -4
	PAM30["SG"] = -2
	PAM30["GS"] = -2
	PAM30["SH"] = -6
	PAM30["HS"] = -6
	PAM30["SI"] = -7
	PAM30["IS"] = -7
	PAM30["SL"] = -8
	PAM30["LS"] = -8
	PAM30["SK"] = -4
	PAM30["KS"] = -4
	PAM30["SM"] = -5
	PAM30["MS"] = -5
	PAM30["SF"] = -6
	PAM30["FS"] = -6
	PAM30["SP"] = -2
	PAM30["PS"] = -2
	PAM30["ST"] = 0
	PAM30["TS"] = 0
	PAM30["SW"] = -5
	PAM30["WS"] = -5
	PAM30["SY"] = -7
	PAM30["YS"] = -7
	PAM30["SV"] = -6
	PAM30["VS"] = -6
	PAM30["SB"] = -1
	PAM30["BS"] = -1
	PAM30["SZ"] = -5
	PAM30["ZS"] = -5
	PAM30["SX"] = -3
	PAM30["XS"] = -3
	PAM30["TA"] = -1
	PAM30["AT"] = -1
	PAM30["TR"] = -6
	PAM30["RT"] = -6
	PAM30["TN"] = -2
	PAM30["NT"] = -2
	PAM30["TD"] = -5
	PAM30["DT"] = -5
	PAM30["TC"] = -8
	PAM30["CT"] = -8
	PAM30["TQ"] = -5
	PAM30["QT"] = -5
	PAM30["TE"] = -6
	PAM30["ET"] = -6
	PAM30["TG"] = -6
	PAM30["GT"] = -6
	PAM30["TH"] = -7
	PAM30["HT"] = -7
	PAM30["TI"] = -2
	PAM30["IT"] = -2
	PAM30["TL"] = -7
	PAM30["LT"] = -7
	PAM30["TK"] = -3
	PAM30["KT"] = -3
	PAM30["TM"] = -4
	PAM30["MT"] = -4
	PAM30["TF"] = -9
	PAM30["FT"] = -9
	PAM30["TP"] = -4
	PAM30["PT"] = -4
	PAM30["TS"] = 0
	PAM30["ST"] = 0
	PAM30["TW"] = -13
	PAM30["WT"] = -13
	PAM30["TY"] = -6
	PAM30["YT"] = -6
	PAM30["TV"] = -3
	PAM30["VT"] = -3
	PAM30["TB"] = -3
	PAM30["BT"] = -3
	PAM30["TZ"] = -6
	PAM30["ZT"] = -6
	PAM30["TX"] = -4
	PAM30["XT"] = -4
	PAM30["WA"] = -13
	PAM30["AW"] = -13
	PAM30["WR"] = -2
	PAM30["RW"] = -2
	PAM30["WN"] = -8
	PAM30["NW"] = -8
	PAM30["WD"] = -15
	PAM30["DW"] = -15
	PAM30["WC"] = -15
	PAM30["CW"] = -15
	PAM30["WQ"] = -13
	PAM30["QW"] = -13
	PAM30["WE"] = -17
	PAM30["EW"] = -17
	PAM30["WG"] = -15
	PAM30["GW"] = -15
	PAM30["WH"] = -7
	PAM30["HW"] = -7
	PAM30["WI"] = -14
	PAM30["IW"] = -14
	PAM30["WL"] = -6
	PAM30["LW"] = -6
	PAM30["WK"] = -12
	PAM30["KW"] = -12
	PAM30["WM"] = -13
	PAM30["MW"] = -13
	PAM30["WF"] = -4
	PAM30["FW"] = -4
	PAM30["WP"] = -14
	PAM30["PW"] = -14
	PAM30["WS"] = -5
	PAM30["SW"] = -5
	PAM30["WT"] = -13
	PAM30["TW"] = -13
	PAM30["WY"] = -5
	PAM30["YW"] = -5
	PAM30["WV"] = -15
	PAM30["VW"] = -15
	PAM30["WB"] = -10
	PAM30["BW"] = -10
	PAM30["WZ"] = -14
	PAM30["ZW"] = -14
	PAM30["WX"] = -11
	PAM30["XW"] = -11
	PAM30["YA"] = -8
	PAM30["AY"] = -8
	PAM30["YR"] = -10
	PAM30["RY"] = -10
	PAM30["YN"] = -4
	PAM30["NY"] = -4
	PAM30["YD"] = -11
	PAM30["DY"] = -11
	PAM30["YC"] = -4
	PAM30["CY"] = -4
	PAM30["YQ"] = -12
	PAM30["QY"] = -12
	PAM30["YE"] = -8
	PAM30["EY"] = -8
	PAM30["YG"] = -14
	PAM30["GY"] = -14
	PAM30["YH"] = -3
	PAM30["HY"] = -3
	PAM30["YI"] = -6
	PAM30["IY"] = -6
	PAM30["YL"] = -7
	PAM30["LY"] = -7
	PAM30["YK"] = -9
	PAM30["KY"] = -9
	PAM30["YM"] = -11
	PAM30["MY"] = -11
	PAM30["YF"] = 2
	PAM30["FY"] = 2
	PAM30["YP"] = -13
	PAM30["PY"] = -13
	PAM30["YS"] = -7
	PAM30["SY"] = -7
	PAM30["YT"] = -6
	PAM30["TY"] = -6
	PAM30["YW"] = -5
	PAM30["WY"] = -5
	PAM30["YV"] = -7
	PAM30["VY"] = -7
	PAM30["YB"] = -6
	PAM30["BY"] = -6
	PAM30["YZ"] = -9
	PAM30["ZY"] = -9
	PAM30["YX"] = -7
	PAM30["XY"] = -7
	PAM30["VA"] = -2
	PAM30["AV"] = -2
	PAM30["VR"] = -8
	PAM30["RV"] = -8
	PAM30["VN"] = -8
	PAM30["NV"] = -8
	PAM30["VD"] = -8
	PAM30["DV"] = -8
	PAM30["VC"] = -6
	PAM30["CV"] = -6
	PAM30["VQ"] = -7
	PAM30["QV"] = -7
	PAM30["VE"] = -6
	PAM30["EV"] = -6
	PAM30["VG"] = -5
	PAM30["GV"] = -5
	PAM30["VH"] = -6
	PAM30["HV"] = -6
	PAM30["VI"] = 2
	PAM30["IV"] = 2
	PAM30["VL"] = -2
	PAM30["LV"] = -2
	PAM30["VK"] = -9
	PAM30["KV"] = -9
	PAM30["VM"] = -1
	PAM30["MV"] = -1
	PAM30["VF"] = -8
	PAM30["FV"] = -8
	PAM30["VP"] = -6
	PAM30["PV"] = -6
	PAM30["VS"] = -6
	PAM30["SV"] = -6
	PAM30["VT"] = -3
	PAM30["TV"] = -3
	PAM30["VW"] = -15
	PAM30["WV"] = -15
	PAM30["VY"] = -7
	PAM30["YV"] = -7
	PAM30["VB"] = -8
	PAM30["BV"] = -8
	PAM30["VZ"] = -6
	PAM30["ZV"] = -6
	PAM30["VX"] = -5
	PAM30["XV"] = -5
	PAM30["BA"] = -3
	PAM30["AB"] = -3
	PAM30["BR"] = -7
	PAM30["RB"] = -7
	PAM30["BN"] = 6
	PAM30["NB"] = 6
	PAM30["BD"] = 6
	PAM30["DB"] = 6
	PAM30["BC"] = -12
	PAM30["CB"] = -12
	PAM30["BQ"] = -3
	PAM30["QB"] = -3
	PAM30["BE"] = 1
	PAM30["EB"] = 1
	PAM30["BG"] = -3
	PAM30["GB"] = -3
	PAM30["BH"] = -1
	PAM30["HB"] = -1
	PAM30["BI"] = -6
	PAM30["IB"] = -6
	PAM30["BL"] = -9
	PAM30["LB"] = -9
	PAM30["BK"] = -2
	PAM30["KB"] = -2
	PAM30["BM"] = -10
	PAM30["MB"] = -10
	PAM30["BF"] = -10
	PAM30["FB"] = -10
	PAM30["BP"] = -7
	PAM30["PB"] = -7
	PAM30["BS"] = -1
	PAM30["SB"] = -1
	PAM30["BT"] = -3
	PAM30["TB"] = -3
	PAM30["BW"] = -10
	PAM30["WB"] = -10
	PAM30["BY"] = -6
	PAM30["YB"] = -6
	PAM30["BV"] = -8
	PAM30["VB"] = -8
	PAM30["BZ"] = 0
	PAM30["ZB"] = 0
	PAM30["BX"] = -5
	PAM30["XB"] = -5
	PAM30["ZA"] = -3
	PAM30["AZ"] = -3
	PAM30["ZR"] = -4
	PAM30["RZ"] = -4
	PAM30["ZN"] = -3
	PAM30["NZ"] = -3
	PAM30["ZD"] = 1
	PAM30["DZ"] = 1
	PAM30["ZC"] = -14
	PAM30["CZ"] = -14
	PAM30["ZQ"] = 6
	PAM30["QZ"] = 6
	PAM30["ZE"] = 6
	PAM30["EZ"] = 6
	PAM30["ZG"] = -5
	PAM30["GZ"] = -5
	PAM30["ZH"] = -1
	PAM30["HZ"] = -1
	PAM30["ZI"] = -6
	PAM30["IZ"] = -6
	PAM30["ZL"] = -7
	PAM30["LZ"] = -7
	PAM30["ZK"] = -4
	PAM30["KZ"] = -4
	PAM30["ZM"] = -5
	PAM30["MZ"] = -5
	PAM30["ZF"] = -13
	PAM30["FZ"] = -13
	PAM30["ZP"] = -4
	PAM30["PZ"] = -4
	PAM30["ZS"] = -5
	PAM30["SZ"] = -5
	PAM30["ZT"] = -6
	PAM30["TZ"] = -6
	PAM30["ZW"] = -14
	PAM30["WZ"] = -14
	PAM30["ZY"] = -9
	PAM30["YZ"] = -9
	PAM30["ZV"] = -6
	PAM30["VZ"] = -6
	PAM30["ZB"] = 0
	PAM30["BZ"] = 0
	PAM30["ZX"] = -5
	PAM30["XZ"] = -5
	PAM30["XA"] = -3
	PAM30["AX"] = -3
	PAM30["XR"] = -6
	PAM30["RX"] = -6
	PAM30["XN"] = -3
	PAM30["NX"] = -3
	PAM30["XD"] = -5
	PAM30["DX"] = -5
	PAM30["XC"] = -9
	PAM30["CX"] = -9
	PAM30["XQ"] = -5
	PAM30["QX"] = -5
	PAM30["XE"] = -5
	PAM30["EX"] = -5
	PAM30["XG"] = -5
	PAM30["GX"] = -5
	PAM30["XH"] = -5
	PAM30["HX"] = -5
	PAM30["XI"] = -5
	PAM30["IX"] = -5
	PAM30["XL"] = -6
	PAM30["LX"] = -6
	PAM30["XK"] = -5
	PAM30["KX"] = -5
	PAM30["XM"] = -5
	PAM30["MX"] = -5
	PAM30["XF"] = -8
	PAM30["FX"] = -8
	PAM30["XP"] = -5
	PAM30["PX"] = -5
	PAM30["XS"] = -3
	PAM30["SX"] = -3
	PAM30["XT"] = -4
	PAM30["TX"] = -4
	PAM30["XW"] = -11
	PAM30["WX"] = -11
	PAM30["XY"] = -7
	PAM30["YX"] = -7
	PAM30["XV"] = -5
	PAM30["VX"] = -5
	PAM30["XB"] = -5
	PAM30["BX"] = -5
	PAM30["XZ"] = -5
	PAM30["ZX"] = -5

	// fmt.Printf("%v %v %v\n", refAA, altAA, Tang[buffer.String()])
	if altAA != "" && refAA != "" && altAA != "X" {
		res = PAM30[buffer.String()]
	} else {
		res = 0
	}
	return res
}

// При использовании данного метода замена считается консерватив-
// ной, если коэффициент Снита больше 0,416 (для одношаговых замен),
// в обратном случае замена считается радикальной

func GetSneathInx(refAA string, altAA string) float64 {
	var buffer bytes.Buffer
	var res float64
	buffer.WriteString(refAA)
	buffer.WriteString(altAA)
	Sneath := make(map[string]float64)

	Sneath["LI"] = 0.889
	Sneath["LV"] = 0.785
	Sneath["LG"] = 0.380
	Sneath["LA"] = 0.643
	Sneath["LP"] = 0.432
	Sneath["LQ"] = 0.501
	Sneath["LN"] = 0.506
	Sneath["LM"] = 0.515
	Sneath["LT"] = 0.432
	Sneath["LS"] = 0.411
	Sneath["LC"] = 0.398
	Sneath["LE"] = 0.333
	Sneath["LD"] = 0.387
	Sneath["LK"] = 0.492
	Sneath["LR"] = 0.360
	Sneath["LY"] = 0.347
	Sneath["LF"] = 0.570
	Sneath["LW"] = 0.368
	Sneath["LH"] = 0.450
	Sneath["IV"] = 0.843
	Sneath["IG"] = 0.371
	Sneath["IA"] = 0.588
	Sneath["IP"] = 0.419
	Sneath["IQ"] = 0.453
	Sneath["IN"] = 0.456
	Sneath["IM"] = 0.494
	Sneath["IT"] = 0.493
	Sneath["IS"] = 0.360
	Sneath["IC"] = 0.348
	Sneath["IE"] = 0.366
	Sneath["ID"] = 0.338
	Sneath["IK"] = 0.477
	Sneath["IR"] = 0.342
	Sneath["IY"] = 0.266
	Sneath["IF"] = 0.487
	Sneath["IW"] = 0.287
	Sneath["IH"] = 0.368
	Sneath["VG"] = 0.437
	Sneath["VA"] = 0.675
	Sneath["VP"] = 0.473
	Sneath["VQ"] = 0.416
	Sneath["VN"] = 0.395
	Sneath["VM"] = 0.465
	Sneath["VT"] = 0.551
	Sneath["VS"] = 0.439
	Sneath["VC"] = 0.430
	Sneath["VE"] = 0.239
	Sneath["VD"] = 0.279
	Sneath["VK"] = 0.419
	Sneath["VR"] = 0.307
	Sneath["VY"] = 0.199
	Sneath["VF"] = 0.380
	Sneath["VW"] = 0.195
	Sneath["VH"] = 0.3
	Sneath["GA"] = 0.659
	Sneath["GP"] = 0.499
	Sneath["GQ"] = 0.163
	Sneath["GN"] = 0.19
	Sneath["GM"] = 0.149
	Sneath["GT"] = 0.396
	Sneath["GS"] = 0.323
	Sneath["GC"] = 0.29
	Sneath["GE"] = 0.049
	Sneath["GD"] = 0.015
	Sneath["GK"] = 0.309
	Sneath["GR"] = 0.149
	Sneath["GY"] = 0.163
	Sneath["GF"] = 0.259
	Sneath["GW"] = 0.138
	Sneath["GH"] = 0.183
	Sneath["AP"] = 0.533
	Sneath["AQ"] = 0.356
	Sneath["AN"] = 0.28
	Sneath["AM"] = 0.421
	Sneath["AT"] = 0.417
	Sneath["AS"] = 0.477
	Sneath["AC"] = 0.578
	Sneath["AE"] = 0.159
	Sneath["AD"] = 0.156
	Sneath["AK"] = 0.426
	Sneath["AR"] = 0.315
	Sneath["AY"] = 0.214
	Sneath["AF"] = 0.356
	Sneath["AW"] = 0.224
	Sneath["AH"] = 0.32
	Sneath["PQ"] = 0.168
	Sneath["PN"] = 0.172
	Sneath["PM"] = 0.265
	Sneath["PT"] = 0.330
	Sneath["PS"] = 0.321
	Sneath["PC"] = 0.318
	Sneath["PE"] = 0.003
	Sneath["PD"] = 0.015
	Sneath["PK"] = 0.295
	Sneath["PR"] = 0.155
	Sneath["PY"] = 0.179
	Sneath["PF"] = 0.282
	Sneath["PW"] = 0.211
	Sneath["PH"] = 0.172
	Sneath["QN"] = 0.589
	Sneath["QM"] = 0.699
	Sneath["QT"] = 0.340
	Sneath["QS"] = 0.501
	Sneath["QC"] = 0.482
	Sneath["QE"] = 0.685
	Sneath["QD"] = 0.492
	Sneath["QK"] = 0.545
	Sneath["QR"] = 0.561
	Sneath["QY"] = 0.368
	Sneath["QF"] = 0.459
	Sneath["QW"] = 0.353
	Sneath["QH"] = 0.406
	Sneath["NM"] = 0.518
	Sneath["NT"] = 0.488
	Sneath["NS"] = 0.581
	Sneath["NC"] = 0.485
	Sneath["NE"] = 0.578
	Sneath["ND"] = 0.637
	Sneath["NK"] = 0.401
	Sneath["NR"] = 0.427
	Sneath["NY"] = 0.391
	Sneath["NF"] = 0.34
	Sneath["NW"] = 0.316
	Sneath["NH"] = 0.459
	Sneath["MT"] = 0.409
	Sneath["MS"] = 0.48
	Sneath["MC"] = 0.612
	Sneath["ME"] = 0.402
	Sneath["MD"] = 0.292
	Sneath["MK"] = 0.482
	Sneath["MR"] = 0.522
	Sneath["MY"] = 0.307
	Sneath["MF"] = 0.465
	Sneath["MW"] = 0.355
	Sneath["MH"] = 0.345
	Sneath["TS"] = 0.668
	Sneath["TC"] = 0.485
	Sneath["TE"] = 0.218
	Sneath["TD"] = 0.254
	Sneath["TK"] = 0.224
	Sneath["TR"] = 0.258
	Sneath["TY"] = 0.285
	Sneath["TF"] = 0.254
	Sneath["TW"] = 0.176
	Sneath["TH"] = 0.208
	Sneath["SC"] = 0.613
	Sneath["SE"] = 0.312
	Sneath["SD"] = 0.330
	Sneath["SK"] = 0.285
	Sneath["SR"] = 0.317
	Sneath["SY"] = 0.354
	Sneath["SF"] = 0.380
	Sneath["SW"] = 0.243
	Sneath["SH"] = 0.342
	Sneath["CE"] = 0.221
	Sneath["CD"] = 0.243
	Sneath["CK"] = 0.269
	Sneath["CR"] = 0.324
	Sneath["CY"] = 0.223
	Sneath["CF"] = 0.289
	Sneath["CW"] = 0.188
	Sneath["CH"] = 0.288
	Sneath["ED"] = 0.84
	Sneath["EK"] = 0.435
	Sneath["ER"] = 0.382
	Sneath["EY"] = 0.261
	Sneath["EF"] = 0.219
	Sneath["EW"] = 0.086
	Sneath["EH"] = 0.201
	Sneath["DK"] = 0.248
	Sneath["DR"] = 0.236
	Sneath["DY"] = 0.287
	Sneath["DF"] = 0.172
	Sneath["DW"] = 0.028
	Sneath["DH"] = 0.2
	Sneath["KR"] = 0.733
	Sneath["KY"] = 0.285
	Sneath["KF"] = 0.381
	Sneath["KW"] = 0.297
	Sneath["KH"] = 0.421
	Sneath["RY"] = 0.407
	Sneath["RF"] = 0.339
	Sneath["RW"] = 0.288
	Sneath["RH"] = 0.396
	Sneath["YF"] = 0.729
	Sneath["YW"] = 0.565
	Sneath["YH"] = 0.504
	Sneath["FW"] = 0.741
	Sneath["FH"] = 0.605
	Sneath["WH"] = 0.484
	// Tang["W"] = "-"
	// Tang["E"] = "-"
	// Tang["G"] = "-"
	// Tang["Q"] = "-"
	// Tang["Y"] = "-"

	// fmt.Printf("%v %v %v\n", refAA, altAA, Tang[buffer.String()])
	if altAA != "" && refAA != "" && altAA != "X" {
		res = Sneath[buffer.String()]
	} else {
		res = 0
	}
	return res
}

// При произвольных заменах аминокислот среднее изменение гидро-
// фобности составляет 1,28 ккал/моль. Замена одной аминокислоты на дру-
// гую считается консервативной в том случае, если разность их гидрофоб-
// ностей ∆Н < 1,28 ккал/моль (для замен в целом) и ∆Н < 1,22 ккал/моль
// (для одношаговых замен — обусловлены заменой одного нуклеотида
// в кодоне ДНК), в противном случае она считается радикальной.

func GetDeltaHInx(refAA string, altAA string) float64 {
	var buffer bytes.Buffer
	var res float64
	buffer.WriteString(refAA)
	buffer.WriteString(altAA)
	deltaH := make(map[string]float64)

	deltaH["PS"] = 2.56
	deltaH["SY"] = 2.83
	deltaH["TI"] = 2.53
	deltaH["ED"] = 0.01
	deltaH["IN"] = 2.96
	deltaH["PT"] = 2.16
	deltaH["ST"] = 0.40
	deltaH["TA"] = 0.19
	deltaH["HR"] = 0.67
	deltaH["VF"] = 0.96
	deltaH["PA"] = 1.97
	deltaH["SI"] = 2.93
	deltaH["ZQ"] = 2.32
	deltaH["HD"] = 0.86
	deltaH["VA"] = 1.06
	deltaH["PL"] = 0.18
	deltaH["SF"] = 2.61
	deltaH["LH"] = 1.02
	deltaH["HN"] = 1.39
	deltaH["VL"] = 0.73
	deltaH["PQ"] = 2.50
	deltaH["SA"] = 0.59
	deltaH["LM"] = 1.12
	deltaH["MR"] = 0.57
	deltaH["TK"] = 1.06
	deltaH["PH"] = 1.20
	deltaH["SL"] = 2.38
	deltaH["LR"] = 1.69
	deltaH["RC"] = 0.08
	deltaH["TM"] = 0.86
	deltaH["PR"] = 1.87
	deltaH["SR"] = 0.69
	deltaH["LW"] = 0.58
	deltaH["RW"] = 2.27
	deltaH["TR"] = 0.29
	deltaH["GS"] = 0.04
	deltaH["SN"] = 0.03
	deltaH["KQ"] = 1.40
	deltaH["DN"] = 0.53
	deltaH["TN"] = 0.44
	deltaH["GV"] = 1.69
	deltaH["SC"] = 0.61
	deltaH["KE"] = 0.95
	deltaH["CW"] = 2.35
	deltaH["AE"] = 0.08
	deltaH["GA"] = 0.63
	deltaH["SW"] = 2.96
	deltaH["KM"] = 0.20
	deltaH["IV"] = 1.28
	deltaH["AD"] = 0.09
	deltaH["GE"] = 0.55
	deltaH["YF"] = 0.22
	deltaH["KR"] = 0.77
	deltaH["IF"] = 0.32
	deltaH["VE"] = 1.14
	deltaH["GR"] = 0.73
	deltaH["YH"] = 1.47
	deltaH["KN"] = 1.49
	deltaH["IL"] = 0.55
	deltaH["VM"] = 0.39
	deltaH["GD"] = 0.54
	deltaH["YD"] = 2.33
	deltaH["QE"] = 0.45
	deltaH["IK"] = 1.47
	deltaH["VD"] = 1.15
	deltaH["GC"] = 0.65
	deltaH["YN"] = 2.86
	deltaH["QH"] = 1.30
	deltaH["IM"] = 1.67
	deltaH["FL"] = 0.23
	deltaH["GW"] = 3.00
	deltaH["YC"] = 2.22
	deltaH["QR"] = 0.63
	deltaH["YR"] = 2.24
	deltaH["FC"] = 2.00

	// fmt.Printf("%v %v %v\n", refAA, altAA, Tang[buffer.String()])
	if altAA != "" && refAA != "" && altAA != "X" {
		res = deltaH[buffer.String()]
	} else {
		res = 0
	}
	return res
}

//GetComplexIndex is....
func GetComplexIndex(refAA string, altAA string, verbose bool) (string, int) {
	var buffer bytes.Buffer
	var Cindex int
	idx := []byte{48, 48, 48, 48} //45
	tang := GetTangInx(refAA, altAA)
	gh := GetGHInx(refAA, altAA)
	// sneath := GetSneathInx(refAA, altAA)
	// deltah := GetDeltaHInx(refAA, altAA)
	pam := GetPAM30Inx(refAA, altAA)
	bloss := GetBLOSSInx(refAA, altAA)

	if tang <= 0.4 {
		idx[0] = 49
		Cindex++
	}
	if gh > 100 {
		idx[1] = 49
		Cindex++
	}
	if pam < 0 {
		idx[2] = 49
		Cindex++
	}
	if bloss < 0 {
		idx[3] = 49
		Cindex++
	}
	if verbose == true {
		buffer.WriteString(fmt.Sprintf("[t:%v,g:%v,p:%v,b:%v]", tang, gh, pam, bloss))
		return fmt.Sprintf("%v%v", bytes.NewBuffer(idx).String(), buffer.String()), Cindex
	} else {
		return bytes.NewBuffer(idx).String(), Cindex
	}

}
