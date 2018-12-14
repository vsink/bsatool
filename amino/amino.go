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
	},
	"TCC": {
		"LName": "Ser",
		"SName": "S",
		"CDA":   "-0.5",
	},
	"TCG": {
		"LName": "Ser",
		"SName": "S",
	},
	"TCT": {
		"LName": "Ser",
		"SName": "S",
	},
	"TTC": {
		"LName": "Phe",
		"SName": "F",
		"CDA":   "0.5",
	},
	"TTT": {
		"LName": "Phe",
		"SName": "F",
	},
	"TTA": {
		"LName": "Leu",
		"SName": "L",
	},
	"TTG": {
		"LName": "Leu",
		"SName": "L",
	},
	"TAC": {
		"LName": "Tyr",
		"SName": "Y",
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

var pam30Idx = map[string]map[string]int{
	"AR": {"PAM30": -7},
	"RA": {"PAM30": -7},
	"AN": {"PAM30": -4},
	"NA": {"PAM30": -4},
	"AD": {"PAM30": -3},
	"DA": {"PAM30": -3},
	"AC": {"PAM30": -6},
	"CA": {"PAM30": -6},
	"AQ": {"PAM30": -4},
	"QA": {"PAM30": -4},
	"AE": {"PAM30": -2},
	"EA": {"PAM30": -2},
	"AG": {"PAM30": -2},
	"GA": {"PAM30": -2},
	"AH": {"PAM30": -7},
	"HA": {"PAM30": -7},
	"AI": {"PAM30": -5},
	"IA": {"PAM30": -5},
	"AL": {"PAM30": -6},
	"LA": {"PAM30": -6},
	"AK": {"PAM30": -7},
	"KA": {"PAM30": -7},
	"AM": {"PAM30": -5},
	"MA": {"PAM30": -5},
	"AF": {"PAM30": -8},
	"FA": {"PAM30": -8},
	"AP": {"PAM30": -2},
	"PA": {"PAM30": -2},
	"AS": {"PAM30": 0},
	"SA": {"PAM30": 0},
	"AT": {"PAM30": -1},
	"TA": {"PAM30": -1},
	"AW": {"PAM30": -13},
	"WA": {"PAM30": -13},
	"AY": {"PAM30": -8},
	"YA": {"PAM30": -8},
	"AV": {"PAM30": -2},
	"VA": {"PAM30": -2},
	"AB": {"PAM30": -3},
	"BA": {"PAM30": -3},
	"AZ": {"PAM30": -3},
	"ZA": {"PAM30": -3},
	"AX": {"PAM30": -3},
	"XA": {"PAM30": -3},
	"RN": {"PAM30": -6},
	"NR": {"PAM30": -6},
	"RD": {"PAM30": -10},
	"DR": {"PAM30": -10},
	"RC": {"PAM30": -8},
	"CR": {"PAM30": -8},
	"RQ": {"PAM30": -2},
	"QR": {"PAM30": -2},
	"RE": {"PAM30": -9},
	"ER": {"PAM30": -9},
	"RG": {"PAM30": -9},
	"GR": {"PAM30": -9},
	"RH": {"PAM30": -2},
	"HR": {"PAM30": -2},
	"RI": {"PAM30": -5},
	"IR": {"PAM30": -5},
	"RL": {"PAM30": -8},
	"LR": {"PAM30": -8},
	"RK": {"PAM30": 0},
	"KR": {"PAM30": 0},
	"RM": {"PAM30": -4},
	"MR": {"PAM30": -4},
	"RF": {"PAM30": -9},
	"FR": {"PAM30": -9},
	"RP": {"PAM30": -4},
	"PR": {"PAM30": -4},
	"RS": {"PAM30": -3},
	"SR": {"PAM30": -3},
	"RT": {"PAM30": -6},
	"TR": {"PAM30": -6},
	"RW": {"PAM30": -2},
	"WR": {"PAM30": -2},
	"RY": {"PAM30": -10},
	"YR": {"PAM30": -10},
	"RV": {"PAM30": -8},
	"VR": {"PAM30": -8},
	"RB": {"PAM30": -7},
	"BR": {"PAM30": -7},
	"RZ": {"PAM30": -4},
	"ZR": {"PAM30": -4},
	"RX": {"PAM30": -6},
	"XR": {"PAM30": -6},
	"ND": {"PAM30": 2},
	"QN": {"PAM30": -3},
	"NE": {"PAM30": -2},
	"EN": {"PAM30": -2},
	"NG": {"PAM30": -3},
	"GN": {"PAM30": -3},
	"NH": {"PAM30": 0},
	"HN": {"PAM30": 0},
	"NI": {"PAM30": -5},
	"IN": {"PAM30": -5},
	"NL": {"PAM30": -7},
	"LN": {"PAM30": -7},
	"NK": {"PAM30": -1},
	"KN": {"PAM30": -1},
	"NM": {"PAM30": -9},
	"MN": {"PAM30": -9},
	"NF": {"PAM30": -9},
	"FN": {"PAM30": -9},
	"NP": {"PAM30": -6},
	"PN": {"PAM30": -6},
	"NS": {"PAM30": 0},
	"SN": {"PAM30": 0},
	"NT": {"PAM30": -2},
	"TN": {"PAM30": -2},
	"NW": {"PAM30": -8},
	"WN": {"PAM30": -8},
	"NY": {"PAM30": -4},
	"YN": {"PAM30": -4},
	"NV": {"PAM30": -8},
	"VN": {"PAM30": -8},
	"NB": {"PAM30": 6},
	"BN": {"PAM30": 6},
	"NZ": {"PAM30": -3},
	"ZN": {"PAM30": -3},
	"NX": {"PAM30": -3},
	"XN": {"PAM30": -3},
	"DC": {"PAM30": -14},
	"CD": {"PAM30": -14},
	"DQ": {"PAM30": -2},
	"QD": {"PAM30": -2},
	"DE": {"PAM30": 2},
	"ED": {"PAM30": 2},
	"DG": {"PAM30": -3},
	"GD": {"PAM30": -3},
	"DH": {"PAM30": -4},
	"HD": {"PAM30": -4},
	"DI": {"PAM30": -7},
	"ID": {"PAM30": -7},
	"DL": {"PAM30": -12},
	"LD": {"PAM30": -12},
	"DK": {"PAM30": -4},
	"KD": {"PAM30": -4},
	"DM": {"PAM30": -11},
	"MD": {"PAM30": -11},
	"DF": {"PAM30": -15},
	"FD": {"PAM30": -15},
	"DP": {"PAM30": -8},
	"PD": {"PAM30": -8},
	"DS": {"PAM30": -4},
	"SD": {"PAM30": -4},
	"DT": {"PAM30": -5},
	"TD": {"PAM30": -5},
	"DW": {"PAM30": -15},
	"WD": {"PAM30": -15},
	"DY": {"PAM30": -11},
	"YD": {"PAM30": -11},
	"DV": {"PAM30": -8},
	"VD": {"PAM30": -8},
	"DB": {"PAM30": 6},
	"BD": {"PAM30": 6},
	"DZ": {"PAM30": 1},
	"ZD": {"PAM30": 1},
	"DX": {"PAM30": -5},
	"XD": {"PAM30": -5},
	"CN": {"PAM30": -11},
	"NC": {"PAM30": -11},
	"CQ": {"PAM30": -14},
	"QC": {"PAM30": -14},
	"CE": {"PAM30": -14},
	"EC": {"PAM30": -14},
	"CG": {"PAM30": -9},
	"GC": {"PAM30": -9},
	"CH": {"PAM30": -7},
	"HC": {"PAM30": -7},
	"CI": {"PAM30": -6},
	"IC": {"PAM30": -6},
	"CL": {"PAM30": -15},
	"LC": {"PAM30": -15},
	"CK": {"PAM30": -14},
	"KC": {"PAM30": -14},
	"CM": {"PAM30": -13},
	"MC": {"PAM30": -13},
	"CF": {"PAM30": -13},
	"FC": {"PAM30": -13},
	"CP": {"PAM30": -8},
	"PC": {"PAM30": -8},
	"CS": {"PAM30": -3},
	"SC": {"PAM30": -3},
	"CT": {"PAM30": -8},
	"TC": {"PAM30": -8},
	"CW": {"PAM30": -15},
	"WC": {"PAM30": -15},
	"CY": {"PAM30": -4},
	"YC": {"PAM30": -4},
	"CV": {"PAM30": -6},
	"VC": {"PAM30": -6},
	"CB": {"PAM30": -12},
	"BC": {"PAM30": -12},
	"CZ": {"PAM30": -14},
	"ZC": {"PAM30": -14},
	"CX": {"PAM30": -9},
	"XC": {"PAM30": -9},
	"NQ": {"PAM30": -3},
	"QE": {"PAM30": 1},
	"EQ": {"PAM30": 1},
	"QG": {"PAM30": -7},
	"GQ": {"PAM30": -7},
	"QH": {"PAM30": 1},
	"HQ": {"PAM30": 1},
	"QI": {"PAM30": -8},
	"IQ": {"PAM30": -8},
	"QL": {"PAM30": -5},
	"LQ": {"PAM30": -5},
	"QK": {"PAM30": -3},
	"KQ": {"PAM30": -3},
	"QM": {"PAM30": -4},
	"MQ": {"PAM30": -4},
	"QF": {"PAM30": -13},
	"FQ": {"PAM30": -13},
	"QP": {"PAM30": -3},
	"PQ": {"PAM30": -3},
	"QS": {"PAM30": -5},
	"SQ": {"PAM30": -5},
	"QT": {"PAM30": -5},
	"TQ": {"PAM30": -5},
	"QW": {"PAM30": -13},
	"WQ": {"PAM30": -13},
	"QY": {"PAM30": -12},
	"YQ": {"PAM30": -12},
	"QV": {"PAM30": -7},
	"VQ": {"PAM30": -7},
	"QB": {"PAM30": -3},
	"BQ": {"PAM30": -3},
	"QZ": {"PAM30": 6},
	"ZQ": {"PAM30": 6},
	"QX": {"PAM30": -5},
	"XQ": {"PAM30": -5},
	"EG": {"PAM30": -4},
	"GE": {"PAM30": -4},
	"EH": {"PAM30": -5},
	"HE": {"PAM30": -5},
	"EI": {"PAM30": -5},
	"IE": {"PAM30": -5},
	"EL": {"PAM30": -9},
	"LE": {"PAM30": -9},
	"EK": {"PAM30": -4},
	"KE": {"PAM30": -4},
	"EM": {"PAM30": -7},
	"ME": {"PAM30": -7},
	"EF": {"PAM30": -14},
	"FE": {"PAM30": -14},
	"EP": {"PAM30": -5},
	"PE": {"PAM30": -5},
	"ES": {"PAM30": -4},
	"SE": {"PAM30": -4},
	"ET": {"PAM30": -6},
	"TE": {"PAM30": -6},
	"EW": {"PAM30": -17},
	"WE": {"PAM30": -17},
	"EY": {"PAM30": -8},
	"YE": {"PAM30": -8},
	"EV": {"PAM30": -6},
	"VE": {"PAM30": -6},
	"EB": {"PAM30": 1},
	"BE": {"PAM30": 1},
	"EZ": {"PAM30": 6},
	"ZE": {"PAM30": 6},
	"EX": {"PAM30": -5},
	"XE": {"PAM30": -5},
	"GH": {"PAM30": -9},
	"HG": {"PAM30": -9},
	"GI": {"PAM30": -11},
	"IG": {"PAM30": -11},
	"GL": {"PAM30": -10},
	"LG": {"PAM30": -10},
	"GK": {"PAM30": -7},
	"KG": {"PAM30": -7},
	"GM": {"PAM30": -8},
	"MG": {"PAM30": -8},
	"GF": {"PAM30": -9},
	"FG": {"PAM30": -9},
	"GP": {"PAM30": -6},
	"PG": {"PAM30": -6},
	"GS": {"PAM30": -2},
	"SG": {"PAM30": -2},
	"GT": {"PAM30": -6},
	"TG": {"PAM30": -6},
	"GW": {"PAM30": -15},
	"WG": {"PAM30": -15},
	"GY": {"PAM30": -14},
	"YG": {"PAM30": -14},
	"GV": {"PAM30": -5},
	"VG": {"PAM30": -5},
	"GB": {"PAM30": -3},
	"BG": {"PAM30": -3},
	"GZ": {"PAM30": -5},
	"ZG": {"PAM30": -5},
	"GX": {"PAM30": -5},
	"XG": {"PAM30": -5},
	"HI": {"PAM30": -9},
	"IH": {"PAM30": -9},
	"HL": {"PAM30": -6},
	"LH": {"PAM30": -6},
	"HK": {"PAM30": -6},
	"KH": {"PAM30": -6},
	"HM": {"PAM30": -10},
	"MH": {"PAM30": -10},
	"HF": {"PAM30": -6},
	"FH": {"PAM30": -6},
	"HP": {"PAM30": -4},
	"PH": {"PAM30": -4},
	"HS": {"PAM30": -6},
	"SH": {"PAM30": -6},
	"HT": {"PAM30": -7},
	"TH": {"PAM30": -7},
	"HW": {"PAM30": -7},
	"WH": {"PAM30": -7},
	"HY": {"PAM30": -3},
	"YH": {"PAM30": -3},
	"HV": {"PAM30": -6},
	"VH": {"PAM30": -6},
	"HB": {"PAM30": -1},
	"BH": {"PAM30": -1},
	"HZ": {"PAM30": -1},
	"ZH": {"PAM30": -1},
	"HX": {"PAM30": -5},
	"XH": {"PAM30": -5},
	"IL": {"PAM30": -1},
	"LI": {"PAM30": -1},
	"IK": {"PAM30": -6},
	"KI": {"PAM30": -6},
	"IM": {"PAM30": -1},
	"MI": {"PAM30": -1},
	"IF": {"PAM30": -2},
	"FI": {"PAM30": -2},
	"IP": {"PAM30": -8},
	"PI": {"PAM30": -8},
	"IS": {"PAM30": -7},
	"SI": {"PAM30": -7},
	"IT": {"PAM30": -2},
	"TI": {"PAM30": -2},
	"IW": {"PAM30": -14},
	"WI": {"PAM30": -14},
	"IY": {"PAM30": -6},
	"YI": {"PAM30": -6},
	"IV": {"PAM30": 2},
	"VI": {"PAM30": 2},
	"IB": {"PAM30": -6},
	"BI": {"PAM30": -6},
	"IZ": {"PAM30": -6},
	"ZI": {"PAM30": -6},
	"IX": {"PAM30": -5},
	"XI": {"PAM30": -5},
	"LK": {"PAM30": -8},
	"KL": {"PAM30": -8},
	"LM": {"PAM30": 1},
	"ML": {"PAM30": 1},
	"LF": {"PAM30": -3},
	"FL": {"PAM30": -3},
	"LP": {"PAM30": -7},
	"PL": {"PAM30": -7},
	"LS": {"PAM30": -8},
	"SL": {"PAM30": -8},
	"LT": {"PAM30": -7},
	"TL": {"PAM30": -7},
	"LW": {"PAM30": -6},
	"WL": {"PAM30": -6},
	"LY": {"PAM30": -7},
	"YL": {"PAM30": -7},
	"LV": {"PAM30": -2},
	"VL": {"PAM30": -2},
	"LB": {"PAM30": -9},
	"BL": {"PAM30": -9},
	"LZ": {"PAM30": -7},
	"ZL": {"PAM30": -7},
	"LX": {"PAM30": -6},
	"XL": {"PAM30": -6},
	"KM": {"PAM30": -2},
	"MK": {"PAM30": -2},
	"KF": {"PAM30": -14},
	"FK": {"PAM30": -14},
	"KP": {"PAM30": -6},
	"PK": {"PAM30": -6},
	"KS": {"PAM30": -4},
	"SK": {"PAM30": -4},
	"KT": {"PAM30": -3},
	"TK": {"PAM30": -3},
	"KW": {"PAM30": -12},
	"WK": {"PAM30": -12},
	"KY": {"PAM30": -9},
	"YK": {"PAM30": -9},
	"KV": {"PAM30": -9},
	"VK": {"PAM30": -9},
	"KB": {"PAM30": -2},
	"BK": {"PAM30": -2},
	"KZ": {"PAM30": -4},
	"ZK": {"PAM30": -4},
	"KX": {"PAM30": -5},
	"XK": {"PAM30": -5},
	"MF": {"PAM30": -4},
	"FM": {"PAM30": -4},
	"MP": {"PAM30": -8},
	"PM": {"PAM30": -8},
	"MS": {"PAM30": -5},
	"SM": {"PAM30": -5},
	"MT": {"PAM30": -4},
	"TM": {"PAM30": -4},
	"MW": {"PAM30": -13},
	"WM": {"PAM30": -13},
	"MY": {"PAM30": -11},
	"YM": {"PAM30": -11},
	"MV": {"PAM30": -1},
	"VM": {"PAM30": -1},
	"MB": {"PAM30": -10},
	"BM": {"PAM30": -10},
	"MZ": {"PAM30": -5},
	"ZM": {"PAM30": -5},
	"MX": {"PAM30": -5},
	"XM": {"PAM30": -5},
	"FP": {"PAM30": -10},
	"PF": {"PAM30": -10},
	"FS": {"PAM30": -6},
	"SF": {"PAM30": -6},
	"FT": {"PAM30": -9},
	"TF": {"PAM30": -9},
	"FW": {"PAM30": -4},
	"WF": {"PAM30": -4},
	"FY": {"PAM30": 2},
	"YF": {"PAM30": 2},
	"FV": {"PAM30": -8},
	"VF": {"PAM30": -8},
	"FB": {"PAM30": -10},
	"BF": {"PAM30": -10},
	"FZ": {"PAM30": -13},
	"ZF": {"PAM30": -13},
	"FX": {"PAM30": -8},
	"XF": {"PAM30": -8},
	"PS": {"PAM30": -2},
	"SP": {"PAM30": -2},
	"PT": {"PAM30": -4},
	"TP": {"PAM30": -4},
	"PW": {"PAM30": -14},
	"WP": {"PAM30": -14},
	"PY": {"PAM30": -13},
	"YP": {"PAM30": -13},
	"PV": {"PAM30": -6},
	"VP": {"PAM30": -6},
	"PB": {"PAM30": -7},
	"BP": {"PAM30": -7},
	"PZ": {"PAM30": -4},
	"ZP": {"PAM30": -4},
	"PX": {"PAM30": -5},
	"XP": {"PAM30": -5},
	"ST": {"PAM30": 0},
	"TS": {"PAM30": 0},
	"SW": {"PAM30": -5},
	"WS": {"PAM30": -5},
	"SY": {"PAM30": -7},
	"YS": {"PAM30": -7},
	"SV": {"PAM30": -6},
	"VS": {"PAM30": -6},
	"SB": {"PAM30": -1},
	"BS": {"PAM30": -1},
	"SZ": {"PAM30": -5},
	"ZS": {"PAM30": -5},
	"SX": {"PAM30": -3},
	"XS": {"PAM30": -3},
	"TW": {"PAM30": -13},
	"WT": {"PAM30": -13},
	"TY": {"PAM30": -6},
	"YT": {"PAM30": -6},
	"TV": {"PAM30": -3},
	"VT": {"PAM30": -3},
	"TB": {"PAM30": -3},
	"BT": {"PAM30": -3},
	"TZ": {"PAM30": -6},
	"ZT": {"PAM30": -6},
	"TX": {"PAM30": -4},
	"XT": {"PAM30": -4},
	"WY": {"PAM30": -5},
	"YW": {"PAM30": -5},
	"WV": {"PAM30": -15},
	"VW": {"PAM30": -15},
	"WB": {"PAM30": -10},
	"BW": {"PAM30": -10},
	"WZ": {"PAM30": -14},
	"ZW": {"PAM30": -14},
	"WX": {"PAM30": -11},
	"XW": {"PAM30": -11},
	"YV": {"PAM30": -7},
	"VY": {"PAM30": -7},
	"YB": {"PAM30": -6},
	"BY": {"PAM30": -6},
	"YZ": {"PAM30": -9},
	"ZY": {"PAM30": -9},
	"YX": {"PAM30": -7},
	"XY": {"PAM30": -7},
	"VB": {"PAM30": -8},
	"BV": {"PAM30": -8},
	"VZ": {"PAM30": -6},
	"ZV": {"PAM30": -6},
	"VX": {"PAM30": -5},
	"XV": {"PAM30": -5},
	"BZ": {"PAM30": 0},
	"ZB": {"PAM30": 0},
	"BX": {"PAM30": -5},
	"XB": {"PAM30": -5},
	"ZX": {"PAM30": -5},
	"XZ": {"PAM30": -5},
}

var blossIdx = map[string]map[string]int{
	"AR": {"BLOSS": -1},
	"RA": {"BLOSS": -1},
	"AN": {"BLOSS": -2},
	"NA": {"BLOSS": -2},
	"AD": {"BLOSS": -2},
	"DA": {"BLOSS": -2},
	"AC": {"BLOSS": 0},
	"CA": {"BLOSS": 0},
	"AQ": {"BLOSS": -1},
	"QA": {"BLOSS": -1},
	"AE": {"BLOSS": -1},
	"EA": {"BLOSS": -1},
	"AG": {"BLOSS": 0},
	"GA": {"BLOSS": 0},
	"AH": {"BLOSS": -2},
	"HA": {"BLOSS": -2},
	"AI": {"BLOSS": -1},
	"IA": {"BLOSS": -1},
	"AL": {"BLOSS": -1},
	"LA": {"BLOSS": -1},
	"AK": {"BLOSS": -1},
	"KA": {"BLOSS": -1},
	"AM": {"BLOSS": -1},
	"MA": {"BLOSS": -1},
	"AF": {"BLOSS": -2},
	"FA": {"BLOSS": -2},
	"AP": {"BLOSS": -1},
	"PA": {"BLOSS": -1},
	"AS": {"BLOSS": 1},
	"SA": {"BLOSS": 1},
	"AT": {"BLOSS": 0},
	"TA": {"BLOSS": 0},
	"AW": {"BLOSS": -3},
	"WA": {"BLOSS": -3},
	"AY": {"BLOSS": -2},
	"YA": {"BLOSS": -2},
	"AV": {"BLOSS": 0},
	"VA": {"BLOSS": 0},
	"AB": {"BLOSS": -2},
	"BA": {"BLOSS": -2},
	"AZ": {"BLOSS": -1},
	"ZA": {"BLOSS": -1},
	"AX": {"BLOSS": 0},
	"XA": {"BLOSS": 0},

	"RN": {"BLOSS": 0},
	"NR": {"BLOSS": 0},
	"RD": {"BLOSS": -2},
	"DR": {"BLOSS": -2},
	"RC": {"BLOSS": -3},
	"CR": {"BLOSS": -3},
	"RQ": {"BLOSS": 1},
	"QR": {"BLOSS": 1},
	"RE": {"BLOSS": 0},
	"ER": {"BLOSS": 0},
	"RG": {"BLOSS": -2},
	"GR": {"BLOSS": -2},
	"RH": {"BLOSS": 0},
	"HR": {"BLOSS": 0},
	"RI": {"BLOSS": -3},
	"IR": {"BLOSS": -3},
	"RL": {"BLOSS": -2},
	"LR": {"BLOSS": -2},
	"RK": {"BLOSS": 2},
	"KR": {"BLOSS": 2},
	"RM": {"BLOSS": -1},
	"MR": {"BLOSS": -1},
	"RF": {"BLOSS": -3},
	"FR": {"BLOSS": -3},
	"RP": {"BLOSS": -2},
	"PR": {"BLOSS": -2},
	"RS": {"BLOSS": -1},
	"SR": {"BLOSS": -1},
	"RT": {"BLOSS": -1},
	"TR": {"BLOSS": -1},
	"RW": {"BLOSS": -3},
	"WR": {"BLOSS": -3},
	"RY": {"BLOSS": -2},
	"YR": {"BLOSS": -2},
	"RV": {"BLOSS": -3},
	"VR": {"BLOSS": -3},
	"RB": {"BLOSS": -1},
	"BR": {"BLOSS": -1},
	"RZ": {"BLOSS": 0},
	"ZR": {"BLOSS": 0},
	"RX": {"BLOSS": -1},
	"XR": {"BLOSS": -1},

	"ND": {"BLOSS": 1},
	"DN": {"BLOSS": 1},
	"NC": {"BLOSS": -3},
	"CN": {"BLOSS": -3},
	"NQ": {"BLOSS": 0},
	"QN": {"BLOSS": 0},
	"NE": {"BLOSS": 0},
	"EN": {"BLOSS": 0},
	"NG": {"BLOSS": 0},
	"GN": {"BLOSS": 0},
	"NH": {"BLOSS": 1},
	"HN": {"BLOSS": 1},
	"NI": {"BLOSS": -3},
	"IN": {"BLOSS": -3},
	"NL": {"BLOSS": -3},
	"LN": {"BLOSS": -3},
	"NK": {"BLOSS": 0},
	"KN": {"BLOSS": 0},
	"NM": {"BLOSS": -2},
	"MN": {"BLOSS": -2},
	"NF": {"BLOSS": -3},
	"FN": {"BLOSS": -3},
	"NP": {"BLOSS": -2},
	"PN": {"BLOSS": -2},
	"NS": {"BLOSS": 1},
	"SN": {"BLOSS": 1},
	"NT": {"BLOSS": 0},
	"TN": {"BLOSS": 0},
	"NW": {"BLOSS": -4},
	"WN": {"BLOSS": -4},
	"NY": {"BLOSS": -2},
	"YN": {"BLOSS": -2},
	"NV": {"BLOSS": -3},
	"VN": {"BLOSS": -3},
	"NB": {"BLOSS": 3},
	"BN": {"BLOSS": 3},
	"NZ": {"BLOSS": 0},
	"ZN": {"BLOSS": 0},
	"NX": {"BLOSS": -1},
	"XN": {"BLOSS": -1},

	"DC": {"BLOSS": -3},
	"CD": {"BLOSS": -3},
	"DQ": {"BLOSS": 0},
	"QD": {"BLOSS": 0},
	"DE": {"BLOSS": 2},
	"ED": {"BLOSS": 2},
	"DG": {"BLOSS": -1},
	"GD": {"BLOSS": -1},
	"DH": {"BLOSS": -1},
	"HD": {"BLOSS": -1},
	"DI": {"BLOSS": -3},
	"ID": {"BLOSS": -3},
	"DL": {"BLOSS": -4},
	"LD": {"BLOSS": -4},
	"DK": {"BLOSS": -1},
	"KD": {"BLOSS": -1},
	"DM": {"BLOSS": -3},
	"MD": {"BLOSS": -3},
	"DF": {"BLOSS": -3},
	"FD": {"BLOSS": -3},
	"DP": {"BLOSS": -1},
	"PD": {"BLOSS": -1},
	"DS": {"BLOSS": 0},
	"SD": {"BLOSS": 0},
	"DT": {"BLOSS": -1},
	"TD": {"BLOSS": -1},
	"DW": {"BLOSS": -4},
	"WD": {"BLOSS": -4},
	"DY": {"BLOSS": -3},
	"YD": {"BLOSS": -3},
	"DV": {"BLOSS": -3},
	"VD": {"BLOSS": -3},
	"DB": {"BLOSS": 4},
	"BD": {"BLOSS": 4},
	"DZ": {"BLOSS": 1},
	"ZD": {"BLOSS": 1},
	"DX": {"BLOSS": -1},
	"XD": {"BLOSS": -1},

	"CQ": {"BLOSS": -3},
	"QC": {"BLOSS": -3},
	"CE": {"BLOSS": -4},
	"EC": {"BLOSS": -4},
	"CG": {"BLOSS": -3},
	"GC": {"BLOSS": -3},
	"CH": {"BLOSS": -3},
	"HC": {"BLOSS": -3},
	"CI": {"BLOSS": -1},
	"IC": {"BLOSS": -1},
	"CL": {"BLOSS": -1},
	"LC": {"BLOSS": -1},
	"CK": {"BLOSS": -3},
	"KC": {"BLOSS": -3},
	"CM": {"BLOSS": -1},
	"MC": {"BLOSS": -1},
	"CF": {"BLOSS": -2},
	"FC": {"BLOSS": -2},
	"CP": {"BLOSS": -3},
	"PC": {"BLOSS": -3},
	"CS": {"BLOSS": -1},
	"SC": {"BLOSS": -1},
	"CT": {"BLOSS": -1},
	"TC": {"BLOSS": -1},
	"CW": {"BLOSS": -2},
	"WC": {"BLOSS": -2},
	"CY": {"BLOSS": -2},
	"YC": {"BLOSS": -2},
	"CV": {"BLOSS": -1},
	"VC": {"BLOSS": -1},
	"CB": {"BLOSS": -3},
	"BC": {"BLOSS": -3},
	"CZ": {"BLOSS": -3},
	"ZC": {"BLOSS": -3},
	"CX": {"BLOSS": -2},
	"XC": {"BLOSS": -2},

	"QE": {"BLOSS": 2},
	"EQ": {"BLOSS": 2},
	"QG": {"BLOSS": -2},
	"GQ": {"BLOSS": -2},
	"QH": {"BLOSS": 0},
	"HQ": {"BLOSS": 0},
	"QI": {"BLOSS": -3},
	"IQ": {"BLOSS": -3},
	"QL": {"BLOSS": -2},
	"LQ": {"BLOSS": -2},
	"QK": {"BLOSS": 1},
	"KQ": {"BLOSS": 1},
	"QM": {"BLOSS": 0},
	"MQ": {"BLOSS": 0},
	"QF": {"BLOSS": -3},
	"FQ": {"BLOSS": -3},
	"QP": {"BLOSS": -1},
	"PQ": {"BLOSS": -1},
	"QS": {"BLOSS": 0},
	"SQ": {"BLOSS": 0},
	"QT": {"BLOSS": -1},
	"TQ": {"BLOSS": -1},
	"QW": {"BLOSS": -2},
	"WQ": {"BLOSS": -2},
	"QY": {"BLOSS": -1},
	"YQ": {"BLOSS": -1},
	"QV": {"BLOSS": -2},
	"VQ": {"BLOSS": -2},
	"QB": {"BLOSS": 0},
	"BQ": {"BLOSS": 0},
	"QZ": {"BLOSS": 3},
	"ZQ": {"BLOSS": 3},
	"QX": {"BLOSS": -1},
	"XQ": {"BLOSS": -1},

	"EG": {"BLOSS": -2},
	"GE": {"BLOSS": -2},
	"EH": {"BLOSS": 0},
	"HE": {"BLOSS": 0},
	"EI": {"BLOSS": -3},
	"IE": {"BLOSS": -3},
	"EL": {"BLOSS": -3},
	"LE": {"BLOSS": -3},
	"EK": {"BLOSS": 1},
	"KE": {"BLOSS": 1},
	"EM": {"BLOSS": -2},
	"ME": {"BLOSS": -2},
	"EF": {"BLOSS": -3},
	"FE": {"BLOSS": -3},
	"EP": {"BLOSS": -1},
	"PE": {"BLOSS": -1},
	"ES": {"BLOSS": 0},
	"SE": {"BLOSS": 0},
	"ET": {"BLOSS": -1},
	"TE": {"BLOSS": -1},
	"EW": {"BLOSS": -3},
	"WE": {"BLOSS": -3},
	"EY": {"BLOSS": -2},
	"YE": {"BLOSS": -2},
	"EV": {"BLOSS": -2},
	"VE": {"BLOSS": -2},
	"EB": {"BLOSS": 1},
	"BE": {"BLOSS": 1},
	"EZ": {"BLOSS": 4},
	"ZE": {"BLOSS": 4},
	"EX": {"BLOSS": -1},
	"XE": {"BLOSS": -1},

	"GH": {"BLOSS": -2},
	"HG": {"BLOSS": -2},
	"GI": {"BLOSS": -4},
	"IG": {"BLOSS": -4},
	"GL": {"BLOSS": -4},
	"LG": {"BLOSS": -4},
	"GK": {"BLOSS": -2},
	"KG": {"BLOSS": -2},
	"GM": {"BLOSS": -3},
	"MG": {"BLOSS": -3},
	"GF": {"BLOSS": -3},
	"FG": {"BLOSS": -3},
	"GP": {"BLOSS": -2},
	"PG": {"BLOSS": -2},
	"GS": {"BLOSS": 0},
	"SG": {"BLOSS": 0},
	"GT": {"BLOSS": -2},
	"TG": {"BLOSS": -2},
	"GW": {"BLOSS": -2},
	"WG": {"BLOSS": -2},
	"GY": {"BLOSS": -3},
	"YG": {"BLOSS": -3},
	"GV": {"BLOSS": -3},
	"VG": {"BLOSS": -3},
	"GB": {"BLOSS": -1},
	"BG": {"BLOSS": -1},
	"GZ": {"BLOSS": -2},
	"ZG": {"BLOSS": -2},
	"GX": {"BLOSS": -1},
	"XG": {"BLOSS": -1},

	"HI": {"BLOSS": -3},
	"IH": {"BLOSS": -3},
	"HL": {"BLOSS": -3},
	"LH": {"BLOSS": -3},
	"HK": {"BLOSS": -1},
	"KH": {"BLOSS": -1},
	"HM": {"BLOSS": -2},
	"MH": {"BLOSS": -2},
	"HF": {"BLOSS": -1},
	"FH": {"BLOSS": -1},
	"HP": {"BLOSS": -2},
	"PH": {"BLOSS": -2},
	"HS": {"BLOSS": -1},
	"SH": {"BLOSS": -1},
	"HT": {"BLOSS": -2},
	"TH": {"BLOSS": -2},
	"HW": {"BLOSS": -2},
	"WH": {"BLOSS": -2},
	"HY": {"BLOSS": 2},
	"YH": {"BLOSS": 2},
	"HV": {"BLOSS": -3},
	"VH": {"BLOSS": -3},
	"HB": {"BLOSS": 0},
	"BH": {"BLOSS": 0},
	"HZ": {"BLOSS": 0},
	"ZH": {"BLOSS": 0},
	"HX": {"BLOSS": -1},
	"XH": {"BLOSS": -1},

	"IL": {"BLOSS": 2},
	"LI": {"BLOSS": 2},
	"IK": {"BLOSS": -3},
	"KI": {"BLOSS": -3},
	"IM": {"BLOSS": 1},
	"MI": {"BLOSS": 1},
	"IF": {"BLOSS": 0},
	"FI": {"BLOSS": 0},
	"IP": {"BLOSS": -3},
	"PI": {"BLOSS": -3},
	"IS": {"BLOSS": -2},
	"SI": {"BLOSS": -2},
	"IT": {"BLOSS": -1},
	"TI": {"BLOSS": -1},
	"IW": {"BLOSS": -3},
	"WI": {"BLOSS": -3},
	"IY": {"BLOSS": -1},
	"YI": {"BLOSS": -1},
	"IV": {"BLOSS": 3},
	"VI": {"BLOSS": 3},
	"IB": {"BLOSS": -3},
	"BI": {"BLOSS": -3},
	"IZ": {"BLOSS": -3},
	"ZI": {"BLOSS": -3},
	"IX": {"BLOSS": -1},
	"XI": {"BLOSS": -1},

	"LK": {"BLOSS": -2},
	"KL": {"BLOSS": -2},
	"LM": {"BLOSS": 2},
	"ML": {"BLOSS": 2},
	"LF": {"BLOSS": 0},
	"FL": {"BLOSS": 0},
	"LP": {"BLOSS": -3},
	"PL": {"BLOSS": -3},
	"LS": {"BLOSS": -2},
	"SL": {"BLOSS": -2},
	"LT": {"BLOSS": -1},
	"TL": {"BLOSS": -1},
	"LW": {"BLOSS": -2},
	"WL": {"BLOSS": -2},
	"LY": {"BLOSS": -1},
	"YL": {"BLOSS": -1},
	"LV": {"BLOSS": 1},
	"VL": {"BLOSS": 1},
	"LB": {"BLOSS": -4},
	"BL": {"BLOSS": -4},
	"LZ": {"BLOSS": -3},
	"ZL": {"BLOSS": -3},
	"LX": {"BLOSS": -1},
	"XL": {"BLOSS": -1},

	"KM": {"BLOSS": -1},
	"MK": {"BLOSS": -1},
	"KF": {"BLOSS": -3},
	"FK": {"BLOSS": -3},
	"KP": {"BLOSS": -1},
	"PK": {"BLOSS": -1},
	"KS": {"BLOSS": 0},
	"SK": {"BLOSS": 0},
	"KT": {"BLOSS": -1},
	"TK": {"BLOSS": -1},
	"KW": {"BLOSS": -3},
	"WK": {"BLOSS": -3},
	"KY": {"BLOSS": -2},
	"YK": {"BLOSS": -2},
	"KV": {"BLOSS": -2},
	"VK": {"BLOSS": -2},
	"KB": {"BLOSS": 0},
	"BK": {"BLOSS": 0},
	"KZ": {"BLOSS": 1},
	"ZK": {"BLOSS": 1},
	"KX": {"BLOSS": -1},
	"XK": {"BLOSS": -1},
	"MF": {"BLOSS": 0},
	"FM": {"BLOSS": 0},
	"MP": {"BLOSS": -2},
	"PM": {"BLOSS": -2},
	"MS": {"BLOSS": -1},
	"SM": {"BLOSS": -1},
	"MT": {"BLOSS": -1},
	"TM": {"BLOSS": -1},
	"MW": {"BLOSS": -1},
	"WM": {"BLOSS": -1},
	"MY": {"BLOSS": -1},
	"YM": {"BLOSS": -1},
	"MV": {"BLOSS": 1},
	"VM": {"BLOSS": 1},
	"MB": {"BLOSS": -3},
	"BM": {"BLOSS": -3},
	"MZ": {"BLOSS": -1},
	"ZM": {"BLOSS": -1},
	"MX": {"BLOSS": -1},
	"XM": {"BLOSS": -1},
	"FP": {"BLOSS": -4},
	"PF": {"BLOSS": -4},
	"FS": {"BLOSS": -2},
	"SF": {"BLOSS": -2},
	"FT": {"BLOSS": -2},
	"TF": {"BLOSS": -2},
	"FW": {"BLOSS": 1},
	"WF": {"BLOSS": 1},
	"FY": {"BLOSS": 3},
	"YF": {"BLOSS": 3},
	"FV": {"BLOSS": -1},
	"VF": {"BLOSS": -1},
	"FB": {"BLOSS": -3},
	"BF": {"BLOSS": -3},
	"FZ": {"BLOSS": -3},
	"ZF": {"BLOSS": -3},
	"FX": {"BLOSS": -1},
	"XF": {"BLOSS": -1},
	"PS": {"BLOSS": -1},
	"SP": {"BLOSS": -1},
	"PT": {"BLOSS": -1},
	"TP": {"BLOSS": -1},
	"PW": {"BLOSS": -4},
	"WP": {"BLOSS": -4},
	"PY": {"BLOSS": -3},
	"YP": {"BLOSS": -3},
	"PV": {"BLOSS": -2},
	"VP": {"BLOSS": -2},
	"PB": {"BLOSS": -2},
	"BP": {"BLOSS": -2},
	"PZ": {"BLOSS": -1},
	"ZP": {"BLOSS": -1},
	"PX": {"BLOSS": -2},
	"XP": {"BLOSS": -2},
	"ST": {"BLOSS": 1},
	"TS": {"BLOSS": 1},
	"SW": {"BLOSS": -3},
	"WS": {"BLOSS": -3},
	"SY": {"BLOSS": -2},
	"YS": {"BLOSS": -2},
	"SV": {"BLOSS": -2},
	"VS": {"BLOSS": -2},
	"SB": {"BLOSS": 0},
	"BS": {"BLOSS": 0},
	"SZ": {"BLOSS": 0},
	"ZS": {"BLOSS": 0},
	"SX": {"BLOSS": 0},
	"XS": {"BLOSS": 0},
	"TW": {"BLOSS": -2},
	"WT": {"BLOSS": -2},
	"TY": {"BLOSS": -2},
	"YT": {"BLOSS": -2},
	"TV": {"BLOSS": 0},
	"VT": {"BLOSS": 0},
	"TB": {"BLOSS": -1},
	"BT": {"BLOSS": -1},
	"TZ": {"BLOSS": -1},
	"ZT": {"BLOSS": -1},
	"TX": {"BLOSS": 0},
	"XT": {"BLOSS": 0},
	"WY": {"BLOSS": 2},
	"YW": {"BLOSS": 2},
	"WV": {"BLOSS": -3},
	"VW": {"BLOSS": -3},
	"WB": {"BLOSS": -4},
	"BW": {"BLOSS": -4},
	"WZ": {"BLOSS": -3},
	"ZW": {"BLOSS": -3},
	"WX": {"BLOSS": -2},
	"XW": {"BLOSS": -2},
	"YV": {"BLOSS": -1},
	"VY": {"BLOSS": -1},
	"YB": {"BLOSS": -3},
	"BY": {"BLOSS": -3},
	"YZ": {"BLOSS": -2},
	"ZY": {"BLOSS": -2},
	"YX": {"BLOSS": -1},
	"XY": {"BLOSS": -1},
	"VB": {"BLOSS": -3},
	"BV": {"BLOSS": -3},
	"VZ": {"BLOSS": -2},
	"ZV": {"BLOSS": -2},
	"VX": {"BLOSS": -1},
	"XV": {"BLOSS": -1},
	"BZ": {"BLOSS": 1},
	"ZB": {"BLOSS": 1},
	"BX": {"BLOSS": -1},
	"XB": {"BLOSS": -1},
	"ZX": {"BLOSS": -1},
	"XZ": {"BLOSS": -1},
}

var ghIdx = map[string]map[string]int{

	"AR": {"GH": 112},
	"RA": {"GH": 112},
	"AN": {"GH": 111},
	"NA": {"GH": 111},
	"AD": {"GH": 126},
	"DA": {"GH": 126},
	"AC": {"GH": 195},
	"CA": {"GH": 195},
	"AQ": {"GH": 91},
	"QA": {"GH": 91},
	"AE": {"GH": 107},
	"EA": {"GH": 107},
	"AG": {"GH": 60},
	"GA": {"GH": 60},
	"AH": {"GH": 86},
	"HA": {"GH": 86},
	"AI": {"GH": 94},
	"IA": {"GH": 94},
	"AL": {"GH": 96},
	"LA": {"GH": 96},
	"AK": {"GH": 106},
	"KA": {"GH": 106},
	"AM": {"GH": 84},
	"MA": {"GH": 84},
	"AF": {"GH": 113},
	"FA": {"GH": 113},
	"AP": {"GH": 27},
	"PA": {"GH": 27},
	"AS": {"GH": 99},
	"SA": {"GH": 99},
	"AT": {"GH": 58},
	"TA": {"GH": 58},
	"AW": {"GH": 148},
	"WA": {"GH": 148},
	"AY": {"GH": 112},
	"YA": {"GH": 112},

	"RN": {"GH": 86},
	"NR": {"GH": 86},
	"RD": {"GH": 96},
	"DR": {"GH": 96},
	"RC": {"GH": 180},
	"CR": {"GH": 180},
	"RQ": {"GH": 43},
	"QR": {"GH": 43},
	"RE": {"GH": 54},
	"ER": {"GH": 54},
	"RG": {"GH": 125},
	"GR": {"GH": 125},
	"RH": {"GH": 29},
	"HR": {"GH": 29},
	"RI": {"GH": 97},
	"IR": {"GH": 97},
	"RL": {"GH": 102},
	"LR": {"GH": 102},
	"RK": {"GH": 26},
	"KR": {"GH": 26},
	"RM": {"GH": 91},
	"MR": {"GH": 91},
	"RF": {"GH": 97},
	"FR": {"GH": 97},
	"RP": {"GH": 103},
	"PR": {"GH": 103},
	"RS": {"GH": 110},
	"SR": {"GH": 110},
	"RT": {"GH": 71},
	"TR": {"GH": 71},
	"RW": {"GH": 101},
	"WR": {"GH": 101},
	"RY": {"GH": 77},
	"YR": {"GH": 77},

	"ND": {"GH": 23},
	"DN": {"GH": 23},
	"NC": {"GH": 139},
	"CN": {"GH": 139},
	"NQ": {"GH": 46},
	"QN": {"GH": 46},
	"NE": {"GH": 42},
	"EN": {"GH": 42},
	"NG": {"GH": 80},
	"GN": {"GH": 80},
	"NH": {"GH": 68},
	"HN": {"GH": 68},
	"NI": {"GH": 149},
	"IN": {"GH": 149},
	"NL": {"GH": 153},
	"LN": {"GH": 153},
	"NK": {"GH": 94},
	"KN": {"GH": 94},
	"NM": {"GH": 142},
	"MN": {"GH": 142},
	"NF": {"GH": 158},
	"FN": {"GH": 158},
	"NP": {"GH": 91},
	"PN": {"GH": 91},
	"NS": {"GH": 46},
	"SN": {"GH": 46},
	"NT": {"GH": 65},
	"TN": {"GH": 65},
	"NW": {"GH": 174},
	"WN": {"GH": 174},
	"NY": {"GH": 143},
	"YN": {"GH": 143},

	"DC": {"GH": 154},
	"CD": {"GH": 154},
	"DQ": {"GH": 61},
	"QD": {"GH": 61},
	"DE": {"GH": 45},
	"ED": {"GH": 45},
	"DG": {"GH": 94},
	"GD": {"GH": 94},
	"DH": {"GH": 81},
	"HD": {"GH": 81},
	"DI": {"GH": 168},
	"ID": {"GH": 168},
	"DL": {"GH": 172},
	"LD": {"GH": 172},
	"DK": {"GH": 101},
	"KD": {"GH": 101},
	"DM": {"GH": 160},
	"MD": {"GH": 160},
	"DF": {"GH": 177},
	"FD": {"GH": 177},
	"DP": {"GH": 108},
	"PD": {"GH": 108},
	"DS": {"GH": 65},
	"SD": {"GH": 65},
	"DT": {"GH": 85},
	"TD": {"GH": 85},
	"DW": {"GH": 181},
	"WD": {"GH": 181},
	"DY": {"GH": 160},
	"YD": {"GH": 160},

	"CQ": {"GH": 154},
	"QC": {"GH": 154},
	"CE": {"GH": 170},
	"EC": {"GH": 170},
	"CG": {"GH": 159},
	"GC": {"GH": 159},
	"CH": {"GH": 174},
	"HC": {"GH": 174},
	"CI": {"GH": 198},
	"IC": {"GH": 198},
	"CL": {"GH": 198},
	"LC": {"GH": 198},
	"CK": {"GH": 202},
	"KC": {"GH": 202},
	"CM": {"GH": 196},
	"MC": {"GH": 196},
	"CF": {"GH": 205},
	"FC": {"GH": 205},
	"CP": {"GH": 169},
	"PC": {"GH": 169},
	"CS": {"GH": 112},
	"SC": {"GH": 112},
	"CT": {"GH": 149},
	"TC": {"GH": 149},
	"CW": {"GH": 215},
	"WC": {"GH": 215},
	"CY": {"GH": 194},
	"YC": {"GH": 194},

	"QE": {"GH": 29},
	"EQ": {"GH": 29},
	"QG": {"GH": 87},
	"GQ": {"GH": 87},
	"QH": {"GH": 24},
	"HQ": {"GH": 24},
	"QI": {"GH": 109},
	"IQ": {"GH": 109},
	"QL": {"GH": 113},
	"LQ": {"GH": 113},
	"QK": {"GH": 53},
	"KQ": {"GH": 53},
	"QM": {"GH": 101},
	"MQ": {"GH": 101},
	"QF": {"GH": 116},
	"FQ": {"GH": 116},
	"QP": {"GH": 76},
	"PQ": {"GH": 76},
	"QS": {"GH": 68},
	"SQ": {"GH": 68},
	"QT": {"GH": 42},
	"TQ": {"GH": 42},
	"QW": {"GH": 130},
	"WQ": {"GH": 130},
	"QY": {"GH": 99},
	"YQ": {"GH": 99},

	"EG": {"GH": 98},
	"GE": {"GH": 98},
	"EH": {"GH": 40},
	"HE": {"GH": 40},
	"EI": {"GH": 134},
	"IE": {"GH": 134},
	"EL": {"GH": 138},
	"LE": {"GH": 138},
	"EK": {"GH": 56},
	"KE": {"GH": 56},
	"EM": {"GH": 126},
	"ME": {"GH": 126},
	"EF": {"GH": 140},
	"FE": {"GH": 140},
	"EP": {"GH": 93},
	"PE": {"GH": 93},
	"ES": {"GH": 80},
	"SE": {"GH": 80},
	"ET": {"GH": 65},
	"TE": {"GH": 65},
	"EW": {"GH": 152},
	"WE": {"GH": 152},
	"EY": {"GH": 122},
	"YE": {"GH": 122},

	"GH": {"GH": 98},
	"HG": {"GH": 98},
	"GI": {"GH": 135},
	"IG": {"GH": 135},
	"GL": {"GH": 138},
	"LG": {"GH": 138},
	"GK": {"GH": 127},
	"KG": {"GH": 127},
	"GM": {"GH": 127},
	"MG": {"GH": 127},
	"GF": {"GH": 153},
	"FG": {"GH": 153},
	"GP": {"GH": 42},
	"PG": {"GH": 42},
	"GS": {"GH": 56},
	"SG": {"GH": 56},
	"GT": {"GH": 59},
	"TG": {"GH": 59},
	"GW": {"GH": 184},
	"WG": {"GH": 184},
	"GY": {"GH": 147},
	"YG": {"GH": 147},

	"HI": {"GH": 94},
	"IH": {"GH": 94},
	"HL": {"GH": 99},
	"LH": {"GH": 99},
	"HK": {"GH": 32},
	"KH": {"GH": 32},
	"HM": {"GH": 87},
	"MH": {"GH": 87},
	"HF": {"GH": 100},
	"FH": {"GH": 100},
	"HP": {"GH": 77},
	"PH": {"GH": 77},
	"HS": {"GH": 89},
	"SH": {"GH": 89},
	"HT": {"GH": 47},
	"TH": {"GH": 47},
	"HW": {"GH": 115},
	"WH": {"GH": 115},
	"HY": {"GH": 83},
	"YH": {"GH": 83},

	"IL": {"GH": 5},
	"LI": {"GH": 5},
	"IK": {"GH": 102},
	"KI": {"GH": 102},
	"IM": {"GH": 10},
	"MI": {"GH": 10},
	"IF": {"GH": 21},
	"FI": {"GH": 21},
	"IP": {"GH": 95},
	"PI": {"GH": 95},
	"IS": {"GH": 142},
	"SI": {"GH": 142},
	"IT": {"GH": 89},
	"TI": {"GH": 89},
	"IW": {"GH": 61},
	"WI": {"GH": 61},
	"IY": {"GH": 33},
	"YI": {"GH": 33},

	"LK": {"GH": 107},
	"KL": {"GH": 107},
	"LM": {"GH": 15},
	"ML": {"GH": 15},
	"LF": {"GH": 22},
	"FL": {"GH": 22},
	"LP": {"GH": 98},
	"PL": {"GH": 98},
	"LS": {"GH": 145},
	"SL": {"GH": 145},
	"LT": {"GH": 92},
	"TL": {"GH": 92},
	"LW": {"GH": 61},
	"WL": {"GH": 61},
	"LY": {"GH": 36},
	"YL": {"GH": 36},

	"KM": {"GH": 95},
	"MK": {"GH": 95},
	"KF": {"GH": 102},
	"FK": {"GH": 102},
	"KP": {"GH": 103},
	"PK": {"GH": 103},
	"KS": {"GH": 121},
	"SK": {"GH": 121},
	"KT": {"GH": 78},
	"TK": {"GH": 78},
	"KW": {"GH": 110},
	"WK": {"GH": 110},
	"KY": {"GH": 85},
	"YK": {"GH": 85},

	"MF": {"GH": 28},
	"FM": {"GH": 28},
	"MP": {"GH": 87},
	"PM": {"GH": 87},
	"MS": {"GH": 135},
	"SM": {"GH": 135},
	"MT": {"GH": 81},
	"TM": {"GH": 81},
	"MW": {"GH": 67},
	"WM": {"GH": 67},
	"MY": {"GH": 36},
	"YM": {"GH": 36},

	"FP": {"GH": 114},
	"PF": {"GH": 114},
	"FS": {"GH": 155},
	"SF": {"GH": 155},
	"FT": {"GH": 103},
	"TF": {"GH": 103},
	"FW": {"GH": 40},
	"WF": {"GH": 40},
	"FY": {"GH": 22},
	"YF": {"GH": 22},

	"PS": {"GH": 74},
	"SP": {"GH": 74},
	"PT": {"GH": 38},
	"TP": {"GH": 38},
	"PW": {"GH": 147},
	"WP": {"GH": 147},
	"PY": {"GH": 110},
	"YP": {"GH": 110},

	"ST": {"GH": 58},
	"TS": {"GH": 58},
	"SW": {"GH": 177},
	"WS": {"GH": 177},
	"SY": {"GH": 144},
	"YS": {"GH": 144},

	"TW": {"GH": 128},
	"WT": {"GH": 128},
	"TY": {"GH": 92},
	"YT": {"GH": 92},

	"YW": {"GH": 37},
	"WY": {"GH": 37},
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
	// GH := make(map[string]int)

	// GH["AR"] = 112
	// GH["RA"] = 112
	// GH["AN"] = 111
	// GH["NA"] = 111
	// GH["AD"] = 126
	// GH["DA"] = 126
	// GH["AC"] = 195
	// GH["CA"] = 195
	// GH["AQ"] = 91
	// GH["QA"] = 91
	// GH["AE"] = 107
	// GH["EA"] = 107
	// GH["AG"] = 60
	// GH["GA"] = 60
	// GH["AH"] = 86
	// GH["HA"] = 86
	// GH["AI"] = 94
	// GH["IA"] = 94
	// GH["AL"] = 96
	// GH["LA"] = 96
	// GH["AK"] = 106
	// GH["KA"] = 106
	// GH["AM"] = 84
	// GH["MA"] = 84
	// GH["AF"] = 113
	// GH["FA"] = 113
	// GH["AP"] = 27
	// GH["PA"] = 27
	// GH["AS"] = 99
	// GH["SA"] = 99
	// GH["AT"] = 58
	// GH["TA"] = 58
	// GH["AW"] = 148
	// GH["WA"] = 148
	// GH["AY"] = 112
	// GH["YA"] = 112
	// GH["RA"] = 112
	// GH["AR"] = 112
	// GH["RN"] = 86
	// GH["NR"] = 86
	// GH["RD"] = 96
	// GH["DR"] = 96
	// GH["RC"] = 180
	// GH["CR"] = 180
	// GH["RQ"] = 43
	// GH["QR"] = 43
	// GH["RE"] = 54
	// GH["ER"] = 54
	// GH["RG"] = 125
	// GH["GR"] = 125
	// GH["RH"] = 29
	// GH["HR"] = 29
	// GH["RI"] = 97
	// GH["IR"] = 97
	// GH["RL"] = 102
	// GH["LR"] = 102
	// GH["RK"] = 26
	// GH["KR"] = 26
	// GH["RM"] = 91
	// GH["MR"] = 91
	// GH["RF"] = 97
	// GH["FR"] = 97
	// GH["RP"] = 103
	// GH["PR"] = 103
	// GH["RS"] = 110
	// GH["SR"] = 110
	// GH["RT"] = 71
	// GH["TR"] = 71
	// GH["RW"] = 101
	// GH["WR"] = 101
	// GH["RY"] = 77
	// GH["YR"] = 77
	// GH["NA"] = 111
	// GH["AN"] = 111
	// GH["NR"] = 86
	// GH["RN"] = 86
	// GH["ND"] = 23
	// GH["DN"] = 23
	// GH["NC"] = 139
	// GH["CN"] = 139
	// GH["NQ"] = 46
	// GH["QN"] = 46
	// GH["NE"] = 42
	// GH["EN"] = 42
	// GH["NG"] = 80
	// GH["GN"] = 80
	// GH["NH"] = 68
	// GH["HN"] = 68
	// GH["NI"] = 149
	// GH["IN"] = 149
	// GH["NL"] = 153
	// GH["LN"] = 153
	// GH["NK"] = 94
	// GH["KN"] = 94
	// GH["NM"] = 142
	// GH["MN"] = 142
	// GH["NF"] = 158
	// GH["FN"] = 158
	// GH["NP"] = 91
	// GH["PN"] = 91
	// GH["NS"] = 46
	// GH["SN"] = 46
	// GH["NT"] = 65
	// GH["TN"] = 65
	// GH["NW"] = 174
	// GH["WN"] = 174
	// GH["NY"] = 143
	// GH["YN"] = 143
	// GH["DA"] = 126
	// GH["AD"] = 126
	// GH["DR"] = 96
	// GH["RD"] = 96
	// GH["DN"] = 23
	// GH["ND"] = 23
	// GH["DC"] = 154
	// GH["CD"] = 154
	// GH["DQ"] = 61
	// GH["QD"] = 61
	// GH["DE"] = 45
	// GH["ED"] = 45
	// GH["DG"] = 94
	// GH["GD"] = 94
	// GH["DH"] = 81
	// GH["HD"] = 81
	// GH["DI"] = 168
	// GH["ID"] = 168
	// GH["DL"] = 172
	// GH["LD"] = 172
	// GH["DK"] = 101
	// GH["KD"] = 101
	// GH["DM"] = 160
	// GH["MD"] = 160
	// GH["DF"] = 177
	// GH["FD"] = 177
	// GH["DP"] = 108
	// GH["PD"] = 108
	// GH["DS"] = 65
	// GH["SD"] = 65
	// GH["DT"] = 85
	// GH["TD"] = 85
	// GH["DW"] = 181
	// GH["WD"] = 181
	// GH["DY"] = 160
	// GH["YD"] = 160
	// GH["CA"] = 195
	// GH["AC"] = 195
	// GH["CR"] = 180
	// GH["RC"] = 180
	// GH["CN"] = 139
	// GH["NC"] = 139
	// GH["CD"] = 154
	// GH["DC"] = 154
	// GH["CQ"] = 154
	// GH["QC"] = 154
	// GH["CE"] = 170
	// GH["EC"] = 170
	// GH["CG"] = 159
	// GH["GC"] = 159
	// GH["CH"] = 174
	// GH["HC"] = 174
	// GH["CI"] = 198
	// GH["IC"] = 198
	// GH["CL"] = 198
	// GH["LC"] = 198
	// GH["CK"] = 202
	// GH["KC"] = 202
	// GH["CM"] = 196
	// GH["MC"] = 196
	// GH["CF"] = 205
	// GH["FC"] = 205
	// GH["CP"] = 169
	// GH["PC"] = 169
	// GH["CS"] = 112
	// GH["SC"] = 112
	// GH["CT"] = 149
	// GH["TC"] = 149
	// GH["CW"] = 215
	// GH["WC"] = 215
	// GH["CY"] = 194
	// GH["YC"] = 194
	// GH["QA"] = 91
	// GH["AQ"] = 91
	// GH["QR"] = 43
	// GH["RQ"] = 43
	// GH["QN"] = 46
	// GH["NQ"] = 46
	// GH["QD"] = 61
	// GH["DQ"] = 61
	// GH["QC"] = 154
	// GH["CQ"] = 154
	// GH["QE"] = 29
	// GH["EQ"] = 29
	// GH["QG"] = 87
	// GH["GQ"] = 87
	// GH["QH"] = 24
	// GH["HQ"] = 24
	// GH["QI"] = 109
	// GH["IQ"] = 109
	// GH["QL"] = 113
	// GH["LQ"] = 113
	// GH["QK"] = 53
	// GH["KQ"] = 53
	// GH["QM"] = 101
	// GH["MQ"] = 101
	// GH["QF"] = 116
	// GH["FQ"] = 116
	// GH["QP"] = 76
	// GH["PQ"] = 76
	// GH["QS"] = 68
	// GH["SQ"] = 68
	// GH["QT"] = 42
	// GH["TQ"] = 42
	// GH["QW"] = 130
	// GH["WQ"] = 130
	// GH["QY"] = 99
	// GH["YQ"] = 99
	// GH["EA"] = 107
	// GH["AE"] = 107
	// GH["ER"] = 54
	// GH["RE"] = 54
	// GH["EN"] = 42
	// GH["NE"] = 42
	// GH["ED"] = 45
	// GH["DE"] = 45
	// GH["EC"] = 170
	// GH["CE"] = 170
	// GH["EQ"] = 29
	// GH["QE"] = 29
	// GH["EG"] = 98
	// GH["GE"] = 98
	// GH["EH"] = 40
	// GH["HE"] = 40
	// GH["EI"] = 134
	// GH["IE"] = 134
	// GH["EL"] = 138
	// GH["LE"] = 138
	// GH["EK"] = 56
	// GH["KE"] = 56
	// GH["EM"] = 126
	// GH["ME"] = 126
	// GH["EF"] = 140
	// GH["FE"] = 140
	// GH["EP"] = 93
	// GH["PE"] = 93
	// GH["ES"] = 80
	// GH["SE"] = 80
	// GH["ET"] = 65
	// GH["TE"] = 65
	// GH["EW"] = 152
	// GH["WE"] = 152
	// GH["EY"] = 122
	// GH["YE"] = 122
	// GH["GA"] = 60
	// GH["AG"] = 60
	// GH["GR"] = 125
	// GH["RG"] = 125
	// GH["GN"] = 80
	// GH["NG"] = 80
	// GH["GD"] = 94
	// GH["DG"] = 94
	// GH["GC"] = 159
	// GH["CG"] = 159
	// GH["GQ"] = 87
	// GH["QG"] = 87
	// GH["GE"] = 98
	// GH["EG"] = 98
	// GH["GH"] = 98
	// GH["HG"] = 98
	// GH["GI"] = 135
	// GH["IG"] = 135
	// GH["GL"] = 138
	// GH["LG"] = 138
	// GH["GK"] = 127
	// GH["KG"] = 127
	// GH["GM"] = 127
	// GH["MG"] = 127
	// GH["GF"] = 153
	// GH["FG"] = 153
	// GH["GP"] = 42
	// GH["PG"] = 42
	// GH["GS"] = 56
	// GH["SG"] = 56
	// GH["GT"] = 59
	// GH["TG"] = 59
	// GH["GW"] = 184
	// GH["WG"] = 184
	// GH["GY"] = 147
	// GH["YG"] = 147
	// GH["HA"] = 86
	// GH["AH"] = 86
	// GH["HR"] = 29
	// GH["RH"] = 29
	// GH["HN"] = 68
	// GH["NH"] = 68
	// GH["HD"] = 81
	// GH["DH"] = 81
	// GH["HC"] = 174
	// GH["CH"] = 174
	// GH["HQ"] = 24
	// GH["QH"] = 24
	// GH["HE"] = 40
	// GH["EH"] = 40
	// GH["HG"] = 98
	// GH["GH"] = 98
	// GH["HI"] = 94
	// GH["IH"] = 94
	// GH["HL"] = 99
	// GH["LH"] = 99
	// GH["HK"] = 32
	// GH["KH"] = 32
	// GH["HM"] = 87
	// GH["MH"] = 87
	// GH["HF"] = 100
	// GH["FH"] = 100
	// GH["HP"] = 77
	// GH["PH"] = 77
	// GH["HS"] = 89
	// GH["SH"] = 89
	// GH["HT"] = 47
	// GH["TH"] = 47
	// GH["HW"] = 115
	// GH["WH"] = 115
	// GH["HY"] = 83
	// GH["YH"] = 83
	// GH["IA"] = 94
	// GH["AI"] = 94
	// GH["IR"] = 97
	// GH["RI"] = 97
	// GH["IN"] = 149
	// GH["NI"] = 149
	// GH["ID"] = 168
	// GH["DI"] = 168
	// GH["IC"] = 198
	// GH["CI"] = 198
	// GH["IQ"] = 109
	// GH["QI"] = 109
	// GH["IE"] = 134
	// GH["EI"] = 134
	// GH["IG"] = 135
	// GH["GI"] = 135
	// GH["IH"] = 94
	// GH["HI"] = 94
	// GH["IL"] = 5
	// GH["LI"] = 5
	// GH["IK"] = 102
	// GH["KI"] = 102
	// GH["IM"] = 10
	// GH["MI"] = 10
	// GH["IF"] = 21
	// GH["FI"] = 21
	// GH["IP"] = 95
	// GH["PI"] = 95
	// GH["IS"] = 142
	// GH["SI"] = 142
	// GH["IT"] = 89
	// GH["TI"] = 89
	// GH["IW"] = 61
	// GH["WI"] = 61
	// GH["IY"] = 33
	// GH["YI"] = 33
	// GH["LA"] = 96
	// GH["AL"] = 96
	// GH["LR"] = 102
	// GH["RL"] = 102
	// GH["LN"] = 153
	// GH["NL"] = 153
	// GH["LD"] = 172
	// GH["DL"] = 172
	// GH["LC"] = 198
	// GH["CL"] = 198
	// GH["LQ"] = 113
	// GH["QL"] = 113
	// GH["LE"] = 138
	// GH["EL"] = 138
	// GH["LG"] = 138
	// GH["GL"] = 138
	// GH["LH"] = 99
	// GH["HL"] = 99
	// GH["LI"] = 5
	// GH["IL"] = 5
	// GH["LK"] = 107
	// GH["KL"] = 107
	// GH["LM"] = 15
	// GH["ML"] = 15
	// GH["LF"] = 22
	// GH["FL"] = 22
	// GH["LP"] = 98
	// GH["PL"] = 98
	// GH["LS"] = 145
	// GH["SL"] = 145
	// GH["LT"] = 92
	// GH["TL"] = 92
	// GH["LW"] = 61
	// GH["WL"] = 61
	// GH["LY"] = 36
	// GH["YL"] = 36
	// GH["KA"] = 106
	// GH["AK"] = 106
	// GH["KR"] = 26
	// GH["RK"] = 26
	// GH["KN"] = 94
	// GH["NK"] = 94
	// GH["KD"] = 101
	// GH["DK"] = 101
	// GH["KC"] = 202
	// GH["CK"] = 202
	// GH["KQ"] = 53
	// GH["QK"] = 53
	// GH["KE"] = 56
	// GH["EK"] = 56
	// GH["KG"] = 127
	// GH["GK"] = 127
	// GH["KH"] = 32
	// GH["HK"] = 32
	// GH["KI"] = 102
	// GH["IK"] = 102
	// GH["KL"] = 107
	// GH["LK"] = 107
	// GH["KM"] = 95
	// GH["MK"] = 95
	// GH["KF"] = 102
	// GH["FK"] = 102
	// GH["KP"] = 103
	// GH["PK"] = 103
	// GH["KS"] = 121
	// GH["SK"] = 121
	// GH["KT"] = 78
	// GH["TK"] = 78
	// GH["KW"] = 110
	// GH["WK"] = 110
	// GH["KY"] = 85
	// GH["YK"] = 85
	// GH["MA"] = 84
	// GH["AM"] = 84
	// GH["MR"] = 91
	// GH["RM"] = 91
	// GH["MN"] = 142
	// GH["NM"] = 142
	// GH["MD"] = 160
	// GH["DM"] = 160
	// GH["MC"] = 196
	// GH["CM"] = 196
	// GH["MQ"] = 101
	// GH["QM"] = 101
	// GH["ME"] = 126
	// GH["EM"] = 126
	// GH["MG"] = 127
	// GH["GM"] = 127
	// GH["MH"] = 87
	// GH["HM"] = 87
	// GH["MI"] = 10
	// GH["IM"] = 10
	// GH["ML"] = 15
	// GH["LM"] = 15
	// GH["MK"] = 95
	// GH["KM"] = 95
	// GH["MF"] = 28
	// GH["FM"] = 28
	// GH["MP"] = 87
	// GH["PM"] = 87
	// GH["MS"] = 135
	// GH["SM"] = 135
	// GH["MT"] = 81
	// GH["TM"] = 81
	// GH["MW"] = 67
	// GH["WM"] = 67
	// GH["MY"] = 36
	// GH["YM"] = 36
	// GH["FA"] = 113
	// GH["AF"] = 113
	// GH["FR"] = 97
	// GH["RF"] = 97
	// GH["FN"] = 158
	// GH["NF"] = 158
	// GH["FD"] = 177
	// GH["DF"] = 177
	// GH["FC"] = 205
	// GH["CF"] = 205
	// GH["FQ"] = 116
	// GH["QF"] = 116
	// GH["FE"] = 140
	// GH["EF"] = 140
	// GH["FG"] = 153
	// GH["GF"] = 153
	// GH["FH"] = 100
	// GH["HF"] = 100
	// GH["FI"] = 21
	// GH["IF"] = 21
	// GH["FL"] = 22
	// GH["LF"] = 22
	// GH["FK"] = 102
	// GH["KF"] = 102
	// GH["FM"] = 28
	// GH["MF"] = 28
	// GH["FP"] = 114
	// GH["PF"] = 114
	// GH["FS"] = 155
	// GH["SF"] = 155
	// GH["FT"] = 103
	// GH["TF"] = 103
	// GH["FW"] = 40
	// GH["WF"] = 40
	// GH["FY"] = 22
	// GH["YF"] = 22
	// GH["PA"] = 27
	// GH["AP"] = 27
	// GH["PR"] = 103
	// GH["RP"] = 103
	// GH["PN"] = 91
	// GH["NP"] = 91
	// GH["PD"] = 108
	// GH["DP"] = 108
	// GH["PC"] = 169
	// GH["CP"] = 169
	// GH["PQ"] = 76
	// GH["QP"] = 76
	// GH["PE"] = 93
	// GH["EP"] = 93
	// GH["PG"] = 42
	// GH["GP"] = 42
	// GH["PH"] = 77
	// GH["HP"] = 77
	// GH["PI"] = 95
	// GH["IP"] = 95
	// GH["PL"] = 98
	// GH["LP"] = 98
	// GH["PK"] = 103
	// GH["KP"] = 103
	// GH["PM"] = 87
	// GH["MP"] = 87
	// GH["PF"] = 114
	// GH["FP"] = 114
	// GH["PS"] = 74
	// GH["SP"] = 74
	// GH["PT"] = 38
	// GH["TP"] = 38
	// GH["PW"] = 147
	// GH["WP"] = 147
	// GH["PY"] = 110
	// GH["YP"] = 110
	// GH["SA"] = 99
	// GH["AS"] = 99
	// GH["SR"] = 110
	// GH["RS"] = 110
	// GH["SN"] = 46
	// GH["NS"] = 46
	// GH["SD"] = 65
	// GH["DS"] = 65
	// GH["SC"] = 112
	// GH["CS"] = 112
	// GH["SQ"] = 68
	// GH["QS"] = 68
	// GH["SE"] = 80
	// GH["ES"] = 80
	// GH["SG"] = 56
	// GH["GS"] = 56
	// GH["SH"] = 89
	// GH["HS"] = 89
	// GH["SI"] = 142
	// GH["IS"] = 142
	// GH["SL"] = 145
	// GH["LS"] = 145
	// GH["SK"] = 121
	// GH["KS"] = 121
	// GH["SM"] = 135
	// GH["MS"] = 135
	// GH["SF"] = 155
	// GH["FS"] = 155
	// GH["SP"] = 74
	// GH["PS"] = 74
	// GH["ST"] = 58
	// GH["TS"] = 58
	// GH["SW"] = 177
	// GH["WS"] = 177
	// GH["SY"] = 144
	// GH["YS"] = 144
	// GH["TA"] = 58
	// GH["AT"] = 58
	// GH["TR"] = 71
	// GH["RT"] = 71
	// GH["TN"] = 65
	// GH["NT"] = 65
	// GH["TD"] = 85
	// GH["DT"] = 85
	// GH["TC"] = 149
	// GH["CT"] = 149
	// GH["TQ"] = 42
	// GH["QT"] = 42
	// GH["TE"] = 65
	// GH["ET"] = 65
	// GH["TG"] = 59
	// GH["GT"] = 59
	// GH["TH"] = 47
	// GH["HT"] = 47
	// GH["TI"] = 89
	// GH["IT"] = 89
	// GH["TL"] = 92
	// GH["LT"] = 92
	// GH["TK"] = 78
	// GH["KT"] = 78
	// GH["TM"] = 81
	// GH["MT"] = 81
	// GH["TF"] = 103
	// GH["FT"] = 103
	// GH["TP"] = 38
	// GH["PT"] = 38
	// GH["TS"] = 58
	// GH["ST"] = 58
	// GH["TW"] = 128
	// GH["WT"] = 128
	// GH["TY"] = 92
	// GH["YT"] = 92
	// GH["WA"] = 148
	// GH["AW"] = 148
	// GH["WR"] = 101
	// GH["RW"] = 101
	// GH["WN"] = 174
	// GH["NW"] = 174
	// GH["WD"] = 181
	// GH["DW"] = 181
	// GH["WC"] = 215
	// GH["CW"] = 215
	// GH["WQ"] = 130
	// GH["QW"] = 130
	// GH["WE"] = 152
	// GH["EW"] = 152
	// GH["WG"] = 184
	// GH["GW"] = 184
	// GH["WH"] = 115
	// GH["HW"] = 115
	// GH["WI"] = 61
	// GH["IW"] = 61
	// GH["WL"] = 61
	// GH["LW"] = 61
	// GH["WK"] = 110
	// GH["KW"] = 110
	// GH["WM"] = 67
	// GH["MW"] = 67
	// GH["WF"] = 40
	// GH["FW"] = 40
	// GH["WP"] = 147
	// GH["PW"] = 147
	// GH["WS"] = 177
	// GH["SW"] = 177
	// GH["WT"] = 128
	// GH["TW"] = 128
	// GH["WY"] = 37
	// GH["YW"] = 37
	// GH["YA"] = 112
	// GH["AY"] = 112
	// GH["YR"] = 77
	// GH["RY"] = 77
	// GH["YN"] = 143
	// GH["NY"] = 143
	// GH["YD"] = 160
	// GH["DY"] = 160
	// GH["YC"] = 194
	// GH["CY"] = 194
	// GH["YQ"] = 99
	// GH["QY"] = 99
	// GH["YE"] = 122
	// GH["EY"] = 122
	// GH["YG"] = 147
	// GH["GY"] = 147
	// GH["YH"] = 83
	// GH["HY"] = 83
	// GH["YI"] = 33
	// GH["IY"] = 33
	// GH["YL"] = 36
	// GH["LY"] = 36
	// GH["YK"] = 85
	// GH["KY"] = 85
	// GH["YM"] = 36
	// GH["MY"] = 36
	// GH["YF"] = 22
	// GH["FY"] = 22
	// GH["YP"] = 110
	// GH["PY"] = 110
	// GH["YS"] = 144
	// GH["SY"] = 144
	// GH["YT"] = 92
	// GH["TY"] = 92
	// GH["YW"] = 37
	// GH["WY"] = 37

	// fmt.Printf("%v %v %v\n", refAA, altAA, Tang[buffer.String()])
	if altAA != "" && refAA != "" && altAA != "X" {
		res = ghIdx[buffer.String()]["GH"]
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
	// BLOSS := make(map[string]int)

	// BLOSS["AR"] = -1
	// BLOSS["RA"] = -1
	// BLOSS["AN"] = -2
	// BLOSS["NA"] = -2
	// BLOSS["AD"] = -2
	// BLOSS["DA"] = -2
	// BLOSS["AC"] = 0
	// BLOSS["CA"] = 0
	// BLOSS["AQ"] = -1
	// BLOSS["QA"] = -1
	// BLOSS["AE"] = -1
	// BLOSS["EA"] = -1
	// BLOSS["AG"] = 0
	// BLOSS["GA"] = 0
	// BLOSS["AH"] = -2
	// BLOSS["HA"] = -2
	// BLOSS["AI"] = -1
	// BLOSS["IA"] = -1
	// BLOSS["AL"] = -1
	// BLOSS["LA"] = -1
	// BLOSS["AK"] = -1
	// BLOSS["KA"] = -1
	// BLOSS["AM"] = -1
	// BLOSS["MA"] = -1
	// BLOSS["AF"] = -2
	// BLOSS["FA"] = -2
	// BLOSS["AP"] = -1
	// BLOSS["PA"] = -1
	// BLOSS["AS"] = 1
	// BLOSS["SA"] = 1
	// BLOSS["AT"] = 0
	// BLOSS["TA"] = 0
	// BLOSS["AW"] = -3
	// BLOSS["WA"] = -3
	// BLOSS["AY"] = -2
	// BLOSS["YA"] = -2
	// BLOSS["AV"] = 0
	// BLOSS["VA"] = 0
	// BLOSS["AB"] = -2
	// BLOSS["BA"] = -2
	// BLOSS["AZ"] = -1
	// BLOSS["ZA"] = -1
	// BLOSS["AX"] = 0
	// BLOSS["XA"] = 0
	// BLOSS["RA"] = -1
	// BLOSS["AR"] = -1
	// BLOSS["RN"] = 0
	// BLOSS["NR"] = 0
	// BLOSS["RD"] = -2
	// BLOSS["DR"] = -2
	// BLOSS["RC"] = -3
	// BLOSS["CR"] = -3
	// BLOSS["RQ"] = 1
	// BLOSS["QR"] = 1
	// BLOSS["RE"] = 0
	// BLOSS["ER"] = 0
	// BLOSS["RG"] = -2
	// BLOSS["GR"] = -2
	// BLOSS["RH"] = 0
	// BLOSS["HR"] = 0
	// BLOSS["RI"] = -3
	// BLOSS["IR"] = -3
	// BLOSS["RL"] = -2
	// BLOSS["LR"] = -2
	// BLOSS["RK"] = 2
	// BLOSS["KR"] = 2
	// BLOSS["RM"] = -1
	// BLOSS["MR"] = -1
	// BLOSS["RF"] = -3
	// BLOSS["FR"] = -3
	// BLOSS["RP"] = -2
	// BLOSS["PR"] = -2
	// BLOSS["RS"] = -1
	// BLOSS["SR"] = -1
	// BLOSS["RT"] = -1
	// BLOSS["TR"] = -1
	// BLOSS["RW"] = -3
	// BLOSS["WR"] = -3
	// BLOSS["RY"] = -2
	// BLOSS["YR"] = -2
	// BLOSS["RV"] = -3
	// BLOSS["VR"] = -3
	// BLOSS["RB"] = -1
	// BLOSS["BR"] = -1
	// BLOSS["RZ"] = 0
	// BLOSS["ZR"] = 0
	// BLOSS["RX"] = -1
	// BLOSS["XR"] = -1
	// BLOSS["NA"] = -2
	// BLOSS["AN"] = -2
	// BLOSS["NR"] = 0
	// BLOSS["RN"] = 0
	// BLOSS["ND"] = 1
	// BLOSS["DN"] = 1
	// BLOSS["NC"] = -3
	// BLOSS["CN"] = -3
	// BLOSS["NQ"] = 0
	// BLOSS["QN"] = 0
	// BLOSS["NE"] = 0
	// BLOSS["EN"] = 0
	// BLOSS["NG"] = 0
	// BLOSS["GN"] = 0
	// BLOSS["NH"] = 1
	// BLOSS["HN"] = 1
	// BLOSS["NI"] = -3
	// BLOSS["IN"] = -3
	// BLOSS["NL"] = -3
	// BLOSS["LN"] = -3
	// BLOSS["NK"] = 0
	// BLOSS["KN"] = 0
	// BLOSS["NM"] = -2
	// BLOSS["MN"] = -2
	// BLOSS["NF"] = -3
	// BLOSS["FN"] = -3
	// BLOSS["NP"] = -2
	// BLOSS["PN"] = -2
	// BLOSS["NS"] = 1
	// BLOSS["SN"] = 1
	// BLOSS["NT"] = 0
	// BLOSS["TN"] = 0
	// BLOSS["NW"] = -4
	// BLOSS["WN"] = -4
	// BLOSS["NY"] = -2
	// BLOSS["YN"] = -2
	// BLOSS["NV"] = -3
	// BLOSS["VN"] = -3
	// BLOSS["NB"] = 3
	// BLOSS["BN"] = 3
	// BLOSS["NZ"] = 0
	// BLOSS["ZN"] = 0
	// BLOSS["NX"] = -1
	// BLOSS["XN"] = -1
	// BLOSS["DA"] = -2
	// BLOSS["AD"] = -2
	// BLOSS["DR"] = -2
	// BLOSS["RD"] = -2
	// BLOSS["DN"] = 1
	// BLOSS["ND"] = 1
	// BLOSS["DC"] = -3
	// BLOSS["CD"] = -3
	// BLOSS["DQ"] = 0
	// BLOSS["QD"] = 0
	// BLOSS["DE"] = 2
	// BLOSS["ED"] = 2
	// BLOSS["DG"] = -1
	// BLOSS["GD"] = -1
	// BLOSS["DH"] = -1
	// BLOSS["HD"] = -1
	// BLOSS["DI"] = -3
	// BLOSS["ID"] = -3
	// BLOSS["DL"] = -4
	// BLOSS["LD"] = -4
	// BLOSS["DK"] = -1
	// BLOSS["KD"] = -1
	// BLOSS["DM"] = -3
	// BLOSS["MD"] = -3
	// BLOSS["DF"] = -3
	// BLOSS["FD"] = -3
	// BLOSS["DP"] = -1
	// BLOSS["PD"] = -1
	// BLOSS["DS"] = 0
	// BLOSS["SD"] = 0
	// BLOSS["DT"] = -1
	// BLOSS["TD"] = -1
	// BLOSS["DW"] = -4
	// BLOSS["WD"] = -4
	// BLOSS["DY"] = -3
	// BLOSS["YD"] = -3
	// BLOSS["DV"] = -3
	// BLOSS["VD"] = -3
	// BLOSS["DB"] = 4
	// BLOSS["BD"] = 4
	// BLOSS["DZ"] = 1
	// BLOSS["ZD"] = 1
	// BLOSS["DX"] = -1
	// BLOSS["XD"] = -1
	// BLOSS["CA"] = 0
	// BLOSS["AC"] = 0
	// BLOSS["CR"] = -3
	// BLOSS["RC"] = -3
	// BLOSS["CN"] = -3
	// BLOSS["NC"] = -3
	// BLOSS["CD"] = -3
	// BLOSS["DC"] = -3
	// BLOSS["CQ"] = -3
	// BLOSS["QC"] = -3
	// BLOSS["CE"] = -4
	// BLOSS["EC"] = -4
	// BLOSS["CG"] = -3
	// BLOSS["GC"] = -3
	// BLOSS["CH"] = -3
	// BLOSS["HC"] = -3
	// BLOSS["CI"] = -1
	// BLOSS["IC"] = -1
	// BLOSS["CL"] = -1
	// BLOSS["LC"] = -1
	// BLOSS["CK"] = -3
	// BLOSS["KC"] = -3
	// BLOSS["CM"] = -1
	// BLOSS["MC"] = -1
	// BLOSS["CF"] = -2
	// BLOSS["FC"] = -2
	// BLOSS["CP"] = -3
	// BLOSS["PC"] = -3
	// BLOSS["CS"] = -1
	// BLOSS["SC"] = -1
	// BLOSS["CT"] = -1
	// BLOSS["TC"] = -1
	// BLOSS["CW"] = -2
	// BLOSS["WC"] = -2
	// BLOSS["CY"] = -2
	// BLOSS["YC"] = -2
	// BLOSS["CV"] = -1
	// BLOSS["VC"] = -1
	// BLOSS["CB"] = -3
	// BLOSS["BC"] = -3
	// BLOSS["CZ"] = -3
	// BLOSS["ZC"] = -3
	// BLOSS["CX"] = -2
	// BLOSS["XC"] = -2
	// BLOSS["QA"] = -1
	// BLOSS["AQ"] = -1
	// BLOSS["QR"] = 1
	// BLOSS["RQ"] = 1
	// BLOSS["QN"] = 0
	// BLOSS["NQ"] = 0
	// BLOSS["QD"] = 0
	// BLOSS["DQ"] = 0
	// BLOSS["QC"] = -3
	// BLOSS["CQ"] = -3
	// BLOSS["QE"] = 2
	// BLOSS["EQ"] = 2
	// BLOSS["QG"] = -2
	// BLOSS["GQ"] = -2
	// BLOSS["QH"] = 0
	// BLOSS["HQ"] = 0
	// BLOSS["QI"] = -3
	// BLOSS["IQ"] = -3
	// BLOSS["QL"] = -2
	// BLOSS["LQ"] = -2
	// BLOSS["QK"] = 1
	// BLOSS["KQ"] = 1
	// BLOSS["QM"] = 0
	// BLOSS["MQ"] = 0
	// BLOSS["QF"] = -3
	// BLOSS["FQ"] = -3
	// BLOSS["QP"] = -1
	// BLOSS["PQ"] = -1
	// BLOSS["QS"] = 0
	// BLOSS["SQ"] = 0
	// BLOSS["QT"] = -1
	// BLOSS["TQ"] = -1
	// BLOSS["QW"] = -2
	// BLOSS["WQ"] = -2
	// BLOSS["QY"] = -1
	// BLOSS["YQ"] = -1
	// BLOSS["QV"] = -2
	// BLOSS["VQ"] = -2
	// BLOSS["QB"] = 0
	// BLOSS["BQ"] = 0
	// BLOSS["QZ"] = 3
	// BLOSS["ZQ"] = 3
	// BLOSS["QX"] = -1
	// BLOSS["XQ"] = -1
	// BLOSS["EA"] = -1
	// BLOSS["AE"] = -1
	// BLOSS["ER"] = 0
	// BLOSS["RE"] = 0
	// BLOSS["EN"] = 0
	// BLOSS["NE"] = 0
	// BLOSS["ED"] = 2
	// BLOSS["DE"] = 2
	// BLOSS["EC"] = -4
	// BLOSS["CE"] = -4
	// BLOSS["EQ"] = 2
	// BLOSS["QE"] = 2
	// BLOSS["EG"] = -2
	// BLOSS["GE"] = -2
	// BLOSS["EH"] = 0
	// BLOSS["HE"] = 0
	// BLOSS["EI"] = -3
	// BLOSS["IE"] = -3
	// BLOSS["EL"] = -3
	// BLOSS["LE"] = -3
	// BLOSS["EK"] = 1
	// BLOSS["KE"] = 1
	// BLOSS["EM"] = -2
	// BLOSS["ME"] = -2
	// BLOSS["EF"] = -3
	// BLOSS["FE"] = -3
	// BLOSS["EP"] = -1
	// BLOSS["PE"] = -1
	// BLOSS["ES"] = 0
	// BLOSS["SE"] = 0
	// BLOSS["ET"] = -1
	// BLOSS["TE"] = -1
	// BLOSS["EW"] = -3
	// BLOSS["WE"] = -3
	// BLOSS["EY"] = -2
	// BLOSS["YE"] = -2
	// BLOSS["EV"] = -2
	// BLOSS["VE"] = -2
	// BLOSS["EB"] = 1
	// BLOSS["BE"] = 1
	// BLOSS["EZ"] = 4
	// BLOSS["ZE"] = 4
	// BLOSS["EX"] = -1
	// BLOSS["XE"] = -1
	// BLOSS["GA"] = 0
	// BLOSS["AG"] = 0
	// BLOSS["GR"] = -2
	// BLOSS["RG"] = -2
	// BLOSS["GN"] = 0
	// BLOSS["NG"] = 0
	// BLOSS["GD"] = -1
	// BLOSS["DG"] = -1
	// BLOSS["GC"] = -3
	// BLOSS["CG"] = -3
	// BLOSS["GQ"] = -2
	// BLOSS["QG"] = -2
	// BLOSS["GE"] = -2
	// BLOSS["EG"] = -2
	// BLOSS["GH"] = -2
	// BLOSS["HG"] = -2
	// BLOSS["GI"] = -4
	// BLOSS["IG"] = -4
	// BLOSS["GL"] = -4
	// BLOSS["LG"] = -4
	// BLOSS["GK"] = -2
	// BLOSS["KG"] = -2
	// BLOSS["GM"] = -3
	// BLOSS["MG"] = -3
	// BLOSS["GF"] = -3
	// BLOSS["FG"] = -3
	// BLOSS["GP"] = -2
	// BLOSS["PG"] = -2
	// BLOSS["GS"] = 0
	// BLOSS["SG"] = 0
	// BLOSS["GT"] = -2
	// BLOSS["TG"] = -2
	// BLOSS["GW"] = -2
	// BLOSS["WG"] = -2
	// BLOSS["GY"] = -3
	// BLOSS["YG"] = -3
	// BLOSS["GV"] = -3
	// BLOSS["VG"] = -3
	// BLOSS["GB"] = -1
	// BLOSS["BG"] = -1
	// BLOSS["GZ"] = -2
	// BLOSS["ZG"] = -2
	// BLOSS["GX"] = -1
	// BLOSS["XG"] = -1
	// BLOSS["HA"] = -2
	// BLOSS["AH"] = -2
	// BLOSS["HR"] = 0
	// BLOSS["RH"] = 0
	// BLOSS["HN"] = 1
	// BLOSS["NH"] = 1
	// BLOSS["HD"] = -1
	// BLOSS["DH"] = -1
	// BLOSS["HC"] = -3
	// BLOSS["CH"] = -3
	// BLOSS["HQ"] = 0
	// BLOSS["QH"] = 0
	// BLOSS["HE"] = 0
	// BLOSS["EH"] = 0
	// BLOSS["HG"] = -2
	// BLOSS["GH"] = -2
	// BLOSS["HI"] = -3
	// BLOSS["IH"] = -3
	// BLOSS["HL"] = -3
	// BLOSS["LH"] = -3
	// BLOSS["HK"] = -1
	// BLOSS["KH"] = -1
	// BLOSS["HM"] = -2
	// BLOSS["MH"] = -2
	// BLOSS["HF"] = -1
	// BLOSS["FH"] = -1
	// BLOSS["HP"] = -2
	// BLOSS["PH"] = -2
	// BLOSS["HS"] = -1
	// BLOSS["SH"] = -1
	// BLOSS["HT"] = -2
	// BLOSS["TH"] = -2
	// BLOSS["HW"] = -2
	// BLOSS["WH"] = -2
	// BLOSS["HY"] = 2
	// BLOSS["YH"] = 2
	// BLOSS["HV"] = -3
	// BLOSS["VH"] = -3
	// BLOSS["HB"] = 0
	// BLOSS["BH"] = 0
	// BLOSS["HZ"] = 0
	// BLOSS["ZH"] = 0
	// BLOSS["HX"] = -1
	// BLOSS["XH"] = -1
	// BLOSS["IA"] = -1
	// BLOSS["AI"] = -1
	// BLOSS["IR"] = -3
	// BLOSS["RI"] = -3
	// BLOSS["IN"] = -3
	// BLOSS["NI"] = -3
	// BLOSS["ID"] = -3
	// BLOSS["DI"] = -3
	// BLOSS["IC"] = -1
	// BLOSS["CI"] = -1
	// BLOSS["IQ"] = -3
	// BLOSS["QI"] = -3
	// BLOSS["IE"] = -3
	// BLOSS["EI"] = -3
	// BLOSS["IG"] = -4
	// BLOSS["GI"] = -4
	// BLOSS["IH"] = -3
	// BLOSS["HI"] = -3
	// BLOSS["IL"] = 2
	// BLOSS["LI"] = 2
	// BLOSS["IK"] = -3
	// BLOSS["KI"] = -3
	// BLOSS["IM"] = 1
	// BLOSS["MI"] = 1
	// BLOSS["IF"] = 0
	// BLOSS["FI"] = 0
	// BLOSS["IP"] = -3
	// BLOSS["PI"] = -3
	// BLOSS["IS"] = -2
	// BLOSS["SI"] = -2
	// BLOSS["IT"] = -1
	// BLOSS["TI"] = -1
	// BLOSS["IW"] = -3
	// BLOSS["WI"] = -3
	// BLOSS["IY"] = -1
	// BLOSS["YI"] = -1
	// BLOSS["IV"] = 3
	// BLOSS["VI"] = 3
	// BLOSS["IB"] = -3
	// BLOSS["BI"] = -3
	// BLOSS["IZ"] = -3
	// BLOSS["ZI"] = -3
	// BLOSS["IX"] = -1
	// BLOSS["XI"] = -1
	// BLOSS["LA"] = -1
	// BLOSS["AL"] = -1
	// BLOSS["LR"] = -2
	// BLOSS["RL"] = -2
	// BLOSS["LN"] = -3
	// BLOSS["NL"] = -3
	// BLOSS["LD"] = -4
	// BLOSS["DL"] = -4
	// BLOSS["LC"] = -1
	// BLOSS["CL"] = -1
	// BLOSS["LQ"] = -2
	// BLOSS["QL"] = -2
	// BLOSS["LE"] = -3
	// BLOSS["EL"] = -3
	// BLOSS["LG"] = -4
	// BLOSS["GL"] = -4
	// BLOSS["LH"] = -3
	// BLOSS["HL"] = -3
	// BLOSS["LI"] = 2
	// BLOSS["IL"] = 2
	// BLOSS["LK"] = -2
	// BLOSS["KL"] = -2
	// BLOSS["LM"] = 2
	// BLOSS["ML"] = 2
	// BLOSS["LF"] = 0
	// BLOSS["FL"] = 0
	// BLOSS["LP"] = -3
	// BLOSS["PL"] = -3
	// BLOSS["LS"] = -2
	// BLOSS["SL"] = -2
	// BLOSS["LT"] = -1
	// BLOSS["TL"] = -1
	// BLOSS["LW"] = -2
	// BLOSS["WL"] = -2
	// BLOSS["LY"] = -1
	// BLOSS["YL"] = -1
	// BLOSS["LV"] = 1
	// BLOSS["VL"] = 1
	// BLOSS["LB"] = -4
	// BLOSS["BL"] = -4
	// BLOSS["LZ"] = -3
	// BLOSS["ZL"] = -3
	// BLOSS["LX"] = -1
	// BLOSS["XL"] = -1
	// BLOSS["KA"] = -1
	// BLOSS["AK"] = -1
	// BLOSS["KR"] = 2
	// BLOSS["RK"] = 2
	// BLOSS["KN"] = 0
	// BLOSS["NK"] = 0
	// BLOSS["KD"] = -1
	// BLOSS["DK"] = -1
	// BLOSS["KC"] = -3
	// BLOSS["CK"] = -3
	// BLOSS["KQ"] = 1
	// BLOSS["QK"] = 1
	// BLOSS["KE"] = 1
	// BLOSS["EK"] = 1
	// BLOSS["KG"] = -2
	// BLOSS["GK"] = -2
	// BLOSS["KH"] = -1
	// BLOSS["HK"] = -1
	// BLOSS["KI"] = -3
	// BLOSS["IK"] = -3
	// BLOSS["KL"] = -2
	// BLOSS["LK"] = -2
	// BLOSS["KM"] = -1
	// BLOSS["MK"] = -1
	// BLOSS["KF"] = -3
	// BLOSS["FK"] = -3
	// BLOSS["KP"] = -1
	// BLOSS["PK"] = -1
	// BLOSS["KS"] = 0
	// BLOSS["SK"] = 0
	// BLOSS["KT"] = -1
	// BLOSS["TK"] = -1
	// BLOSS["KW"] = -3
	// BLOSS["WK"] = -3
	// BLOSS["KY"] = -2
	// BLOSS["YK"] = -2
	// BLOSS["KV"] = -2
	// BLOSS["VK"] = -2
	// BLOSS["KB"] = 0
	// BLOSS["BK"] = 0
	// BLOSS["KZ"] = 1
	// BLOSS["ZK"] = 1
	// BLOSS["KX"] = -1
	// BLOSS["XK"] = -1
	// BLOSS["MA"] = -1
	// BLOSS["AM"] = -1
	// BLOSS["MR"] = -1
	// BLOSS["RM"] = -1
	// BLOSS["MN"] = -2
	// BLOSS["NM"] = -2
	// BLOSS["MD"] = -3
	// BLOSS["DM"] = -3
	// BLOSS["MC"] = -1
	// BLOSS["CM"] = -1
	// BLOSS["MQ"] = 0
	// BLOSS["QM"] = 0
	// BLOSS["ME"] = -2
	// BLOSS["EM"] = -2
	// BLOSS["MG"] = -3
	// BLOSS["GM"] = -3
	// BLOSS["MH"] = -2
	// BLOSS["HM"] = -2
	// BLOSS["MI"] = 1
	// BLOSS["IM"] = 1
	// BLOSS["ML"] = 2
	// BLOSS["LM"] = 2
	// BLOSS["MK"] = -1
	// BLOSS["KM"] = -1
	// BLOSS["MF"] = 0
	// BLOSS["FM"] = 0
	// BLOSS["MP"] = -2
	// BLOSS["PM"] = -2
	// BLOSS["MS"] = -1
	// BLOSS["SM"] = -1
	// BLOSS["MT"] = -1
	// BLOSS["TM"] = -1
	// BLOSS["MW"] = -1
	// BLOSS["WM"] = -1
	// BLOSS["MY"] = -1
	// BLOSS["YM"] = -1
	// BLOSS["MV"] = 1
	// BLOSS["VM"] = 1
	// BLOSS["MB"] = -3
	// BLOSS["BM"] = -3
	// BLOSS["MZ"] = -1
	// BLOSS["ZM"] = -1
	// BLOSS["MX"] = -1
	// BLOSS["XM"] = -1
	// BLOSS["FA"] = -2
	// BLOSS["AF"] = -2
	// BLOSS["FR"] = -3
	// BLOSS["RF"] = -3
	// BLOSS["FN"] = -3
	// BLOSS["NF"] = -3
	// BLOSS["FD"] = -3
	// BLOSS["DF"] = -3
	// BLOSS["FC"] = -2
	// BLOSS["CF"] = -2
	// BLOSS["FQ"] = -3
	// BLOSS["QF"] = -3
	// BLOSS["FE"] = -3
	// BLOSS["EF"] = -3
	// BLOSS["FG"] = -3
	// BLOSS["GF"] = -3
	// BLOSS["FH"] = -1
	// BLOSS["HF"] = -1
	// BLOSS["FI"] = 0
	// BLOSS["IF"] = 0
	// BLOSS["FL"] = 0
	// BLOSS["LF"] = 0
	// BLOSS["FK"] = -3
	// BLOSS["KF"] = -3
	// BLOSS["FM"] = 0
	// BLOSS["MF"] = 0
	// BLOSS["FP"] = -4
	// BLOSS["PF"] = -4
	// BLOSS["FS"] = -2
	// BLOSS["SF"] = -2
	// BLOSS["FT"] = -2
	// BLOSS["TF"] = -2
	// BLOSS["FW"] = 1
	// BLOSS["WF"] = 1
	// BLOSS["FY"] = 3
	// BLOSS["YF"] = 3
	// BLOSS["FV"] = -1
	// BLOSS["VF"] = -1
	// BLOSS["FB"] = -3
	// BLOSS["BF"] = -3
	// BLOSS["FZ"] = -3
	// BLOSS["ZF"] = -3
	// BLOSS["FX"] = -1
	// BLOSS["XF"] = -1
	// BLOSS["PA"] = -1
	// BLOSS["AP"] = -1
	// BLOSS["PR"] = -2
	// BLOSS["RP"] = -2
	// BLOSS["PN"] = -2
	// BLOSS["NP"] = -2
	// BLOSS["PD"] = -1
	// BLOSS["DP"] = -1
	// BLOSS["PC"] = -3
	// BLOSS["CP"] = -3
	// BLOSS["PQ"] = -1
	// BLOSS["QP"] = -1
	// BLOSS["PE"] = -1
	// BLOSS["EP"] = -1
	// BLOSS["PG"] = -2
	// BLOSS["GP"] = -2
	// BLOSS["PH"] = -2
	// BLOSS["HP"] = -2
	// BLOSS["PI"] = -3
	// BLOSS["IP"] = -3
	// BLOSS["PL"] = -3
	// BLOSS["LP"] = -3
	// BLOSS["PK"] = -1
	// BLOSS["KP"] = -1
	// BLOSS["PM"] = -2
	// BLOSS["MP"] = -2
	// BLOSS["PF"] = -4
	// BLOSS["FP"] = -4
	// BLOSS["PS"] = -1
	// BLOSS["SP"] = -1
	// BLOSS["PT"] = -1
	// BLOSS["TP"] = -1
	// BLOSS["PW"] = -4
	// BLOSS["WP"] = -4
	// BLOSS["PY"] = -3
	// BLOSS["YP"] = -3
	// BLOSS["PV"] = -2
	// BLOSS["VP"] = -2
	// BLOSS["PB"] = -2
	// BLOSS["BP"] = -2
	// BLOSS["PZ"] = -1
	// BLOSS["ZP"] = -1
	// BLOSS["PX"] = -2
	// BLOSS["XP"] = -2
	// BLOSS["SA"] = 1
	// BLOSS["AS"] = 1
	// BLOSS["SR"] = -1
	// BLOSS["RS"] = -1
	// BLOSS["SN"] = 1
	// BLOSS["NS"] = 1
	// BLOSS["SD"] = 0
	// BLOSS["DS"] = 0
	// BLOSS["SC"] = -1
	// BLOSS["CS"] = -1
	// BLOSS["SQ"] = 0
	// BLOSS["QS"] = 0
	// BLOSS["SE"] = 0
	// BLOSS["ES"] = 0
	// BLOSS["SG"] = 0
	// BLOSS["GS"] = 0
	// BLOSS["SH"] = -1
	// BLOSS["HS"] = -1
	// BLOSS["SI"] = -2
	// BLOSS["IS"] = -2
	// BLOSS["SL"] = -2
	// BLOSS["LS"] = -2
	// BLOSS["SK"] = 0
	// BLOSS["KS"] = 0
	// BLOSS["SM"] = -1
	// BLOSS["MS"] = -1
	// BLOSS["SF"] = -2
	// BLOSS["FS"] = -2
	// BLOSS["SP"] = -1
	// BLOSS["PS"] = -1
	// BLOSS["ST"] = 1
	// BLOSS["TS"] = 1
	// BLOSS["SW"] = -3
	// BLOSS["WS"] = -3
	// BLOSS["SY"] = -2
	// BLOSS["YS"] = -2
	// BLOSS["SV"] = -2
	// BLOSS["VS"] = -2
	// BLOSS["SB"] = 0
	// BLOSS["BS"] = 0
	// BLOSS["SZ"] = 0
	// BLOSS["ZS"] = 0
	// BLOSS["SX"] = 0
	// BLOSS["XS"] = 0
	// BLOSS["TA"] = 0
	// BLOSS["AT"] = 0
	// BLOSS["TR"] = -1
	// BLOSS["RT"] = -1
	// BLOSS["TN"] = 0
	// BLOSS["NT"] = 0
	// BLOSS["TD"] = -1
	// BLOSS["DT"] = -1
	// BLOSS["TC"] = -1
	// BLOSS["CT"] = -1
	// BLOSS["TQ"] = -1
	// BLOSS["QT"] = -1
	// BLOSS["TE"] = -1
	// BLOSS["ET"] = -1
	// BLOSS["TG"] = -2
	// BLOSS["GT"] = -2
	// BLOSS["TH"] = -2
	// BLOSS["HT"] = -2
	// BLOSS["TI"] = -1
	// BLOSS["IT"] = -1
	// BLOSS["TL"] = -1
	// BLOSS["LT"] = -1
	// BLOSS["TK"] = -1
	// BLOSS["KT"] = -1
	// BLOSS["TM"] = -1
	// BLOSS["MT"] = -1
	// BLOSS["TF"] = -2
	// BLOSS["FT"] = -2
	// BLOSS["TP"] = -1
	// BLOSS["PT"] = -1
	// BLOSS["TS"] = 1
	// BLOSS["ST"] = 1
	// BLOSS["TW"] = -2
	// BLOSS["WT"] = -2
	// BLOSS["TY"] = -2
	// BLOSS["YT"] = -2
	// BLOSS["TV"] = 0
	// BLOSS["VT"] = 0
	// BLOSS["TB"] = -1
	// BLOSS["BT"] = -1
	// BLOSS["TZ"] = -1
	// BLOSS["ZT"] = -1
	// BLOSS["TX"] = 0
	// BLOSS["XT"] = 0
	// BLOSS["WA"] = -3
	// BLOSS["AW"] = -3
	// BLOSS["WR"] = -3
	// BLOSS["RW"] = -3
	// BLOSS["WN"] = -4
	// BLOSS["NW"] = -4
	// BLOSS["WD"] = -4
	// BLOSS["DW"] = -4
	// BLOSS["WC"] = -2
	// BLOSS["CW"] = -2
	// BLOSS["WQ"] = -2
	// BLOSS["QW"] = -2
	// BLOSS["WE"] = -3
	// BLOSS["EW"] = -3
	// BLOSS["WG"] = -2
	// BLOSS["GW"] = -2
	// BLOSS["WH"] = -2
	// BLOSS["HW"] = -2
	// BLOSS["WI"] = -3
	// BLOSS["IW"] = -3
	// BLOSS["WL"] = -2
	// BLOSS["LW"] = -2
	// BLOSS["WK"] = -3
	// BLOSS["KW"] = -3
	// BLOSS["WM"] = -1
	// BLOSS["MW"] = -1
	// BLOSS["WF"] = 1
	// BLOSS["FW"] = 1
	// BLOSS["WP"] = -4
	// BLOSS["PW"] = -4
	// BLOSS["WS"] = -3
	// BLOSS["SW"] = -3
	// BLOSS["WT"] = -2
	// BLOSS["TW"] = -2
	// BLOSS["WY"] = 2
	// BLOSS["YW"] = 2
	// BLOSS["WV"] = -3
	// BLOSS["VW"] = -3
	// BLOSS["WB"] = -4
	// BLOSS["BW"] = -4
	// BLOSS["WZ"] = -3
	// BLOSS["ZW"] = -3
	// BLOSS["WX"] = -2
	// BLOSS["XW"] = -2
	// BLOSS["YA"] = -2
	// BLOSS["AY"] = -2
	// BLOSS["YR"] = -2
	// BLOSS["RY"] = -2
	// BLOSS["YN"] = -2
	// BLOSS["NY"] = -2
	// BLOSS["YD"] = -3
	// BLOSS["DY"] = -3
	// BLOSS["YC"] = -2
	// BLOSS["CY"] = -2
	// BLOSS["YQ"] = -1
	// BLOSS["QY"] = -1
	// BLOSS["YE"] = -2
	// BLOSS["EY"] = -2
	// BLOSS["YG"] = -3
	// BLOSS["GY"] = -3
	// BLOSS["YH"] = 2
	// BLOSS["HY"] = 2
	// BLOSS["YI"] = -1
	// BLOSS["IY"] = -1
	// BLOSS["YL"] = -1
	// BLOSS["LY"] = -1
	// BLOSS["YK"] = -2
	// BLOSS["KY"] = -2
	// BLOSS["YM"] = -1
	// BLOSS["MY"] = -1
	// BLOSS["YF"] = 3
	// BLOSS["FY"] = 3
	// BLOSS["YP"] = -3
	// BLOSS["PY"] = -3
	// BLOSS["YS"] = -2
	// BLOSS["SY"] = -2
	// BLOSS["YT"] = -2
	// BLOSS["TY"] = -2
	// BLOSS["YW"] = 2
	// BLOSS["WY"] = 2
	// BLOSS["YV"] = -1
	// BLOSS["VY"] = -1
	// BLOSS["YB"] = -3
	// BLOSS["BY"] = -3
	// BLOSS["YZ"] = -2
	// BLOSS["ZY"] = -2
	// BLOSS["YX"] = -1
	// BLOSS["XY"] = -1
	// BLOSS["VA"] = 0
	// BLOSS["AV"] = 0
	// BLOSS["VR"] = -3
	// BLOSS["RV"] = -3
	// BLOSS["VN"] = -3
	// BLOSS["NV"] = -3
	// BLOSS["VD"] = -3
	// BLOSS["DV"] = -3
	// BLOSS["VC"] = -1
	// BLOSS["CV"] = -1
	// BLOSS["VQ"] = -2
	// BLOSS["QV"] = -2
	// BLOSS["VE"] = -2
	// BLOSS["EV"] = -2
	// BLOSS["VG"] = -3
	// BLOSS["GV"] = -3
	// BLOSS["VH"] = -3
	// BLOSS["HV"] = -3
	// BLOSS["VI"] = 3
	// BLOSS["IV"] = 3
	// BLOSS["VL"] = 1
	// BLOSS["LV"] = 1
	// BLOSS["VK"] = -2
	// BLOSS["KV"] = -2
	// BLOSS["VM"] = 1
	// BLOSS["MV"] = 1
	// BLOSS["VF"] = -1
	// BLOSS["FV"] = -1
	// BLOSS["VP"] = -2
	// BLOSS["PV"] = -2
	// BLOSS["VS"] = -2
	// BLOSS["SV"] = -2
	// BLOSS["VT"] = 0
	// BLOSS["TV"] = 0
	// BLOSS["VW"] = -3
	// BLOSS["WV"] = -3
	// BLOSS["VY"] = -1
	// BLOSS["YV"] = -1
	// BLOSS["VB"] = -3
	// BLOSS["BV"] = -3
	// BLOSS["VZ"] = -2
	// BLOSS["ZV"] = -2
	// BLOSS["VX"] = -1
	// BLOSS["XV"] = -1
	// BLOSS["BA"] = -2
	// BLOSS["AB"] = -2
	// BLOSS["BR"] = -1
	// BLOSS["RB"] = -1
	// BLOSS["BN"] = 3
	// BLOSS["NB"] = 3
	// BLOSS["BD"] = 4
	// BLOSS["DB"] = 4
	// BLOSS["BC"] = -3
	// BLOSS["CB"] = -3
	// BLOSS["BQ"] = 0
	// BLOSS["QB"] = 0
	// BLOSS["BE"] = 1
	// BLOSS["EB"] = 1
	// BLOSS["BG"] = -1
	// BLOSS["GB"] = -1
	// BLOSS["BH"] = 0
	// BLOSS["HB"] = 0
	// BLOSS["BI"] = -3
	// BLOSS["IB"] = -3
	// BLOSS["BL"] = -4
	// BLOSS["LB"] = -4
	// BLOSS["BK"] = 0
	// BLOSS["KB"] = 0
	// BLOSS["BM"] = -3
	// BLOSS["MB"] = -3
	// BLOSS["BF"] = -3
	// BLOSS["FB"] = -3
	// BLOSS["BP"] = -2
	// BLOSS["PB"] = -2
	// BLOSS["BS"] = 0
	// BLOSS["SB"] = 0
	// BLOSS["BT"] = -1
	// BLOSS["TB"] = -1
	// BLOSS["BW"] = -4
	// BLOSS["WB"] = -4
	// BLOSS["BY"] = -3
	// BLOSS["YB"] = -3
	// BLOSS["BV"] = -3
	// BLOSS["VB"] = -3
	// BLOSS["BZ"] = 1
	// BLOSS["ZB"] = 1
	// BLOSS["BX"] = -1
	// BLOSS["XB"] = -1
	// BLOSS["ZA"] = -1
	// BLOSS["AZ"] = -1
	// BLOSS["ZR"] = 0
	// BLOSS["RZ"] = 0
	// BLOSS["ZN"] = 0
	// BLOSS["NZ"] = 0
	// BLOSS["ZD"] = 1
	// BLOSS["DZ"] = 1
	// BLOSS["ZC"] = -3
	// BLOSS["CZ"] = -3
	// BLOSS["ZQ"] = 3
	// BLOSS["QZ"] = 3
	// BLOSS["ZE"] = 4
	// BLOSS["EZ"] = 4
	// BLOSS["ZG"] = -2
	// BLOSS["GZ"] = -2
	// BLOSS["ZH"] = 0
	// BLOSS["HZ"] = 0
	// BLOSS["ZI"] = -3
	// BLOSS["IZ"] = -3
	// BLOSS["ZL"] = -3
	// BLOSS["LZ"] = -3
	// BLOSS["ZK"] = 1
	// BLOSS["KZ"] = 1
	// BLOSS["ZM"] = -1
	// BLOSS["MZ"] = -1
	// BLOSS["ZF"] = -3
	// BLOSS["FZ"] = -3
	// BLOSS["ZP"] = -1
	// BLOSS["PZ"] = -1
	// BLOSS["ZS"] = 0
	// BLOSS["SZ"] = 0
	// BLOSS["ZT"] = -1
	// BLOSS["TZ"] = -1
	// BLOSS["ZW"] = -3
	// BLOSS["WZ"] = -3
	// BLOSS["ZY"] = -2
	// BLOSS["YZ"] = -2
	// BLOSS["ZV"] = -2
	// BLOSS["VZ"] = -2
	// BLOSS["ZB"] = 1
	// BLOSS["BZ"] = 1
	// BLOSS["ZX"] = -1
	// BLOSS["XZ"] = -1
	// BLOSS["XA"] = 0
	// BLOSS["AX"] = 0
	// BLOSS["XR"] = -1
	// BLOSS["RX"] = -1
	// BLOSS["XN"] = -1
	// BLOSS["NX"] = -1
	// BLOSS["XD"] = -1
	// BLOSS["DX"] = -1
	// BLOSS["XC"] = -2
	// BLOSS["CX"] = -2
	// BLOSS["XQ"] = -1
	// BLOSS["QX"] = -1
	// BLOSS["XE"] = -1
	// BLOSS["EX"] = -1
	// BLOSS["XG"] = -1
	// BLOSS["GX"] = -1
	// BLOSS["XH"] = -1
	// BLOSS["HX"] = -1
	// BLOSS["XI"] = -1
	// BLOSS["IX"] = -1
	// BLOSS["XL"] = -1
	// BLOSS["LX"] = -1
	// BLOSS["XK"] = -1
	// BLOSS["KX"] = -1
	// BLOSS["XM"] = -1
	// BLOSS["MX"] = -1
	// BLOSS["XF"] = -1
	// BLOSS["FX"] = -1
	// BLOSS["XP"] = -2
	// BLOSS["PX"] = -2
	// BLOSS["XS"] = 0
	// BLOSS["SX"] = 0
	// BLOSS["XT"] = 0
	// BLOSS["TX"] = 0
	// BLOSS["XW"] = -2
	// BLOSS["WX"] = -2
	// BLOSS["XY"] = -1
	// BLOSS["YX"] = -1
	// BLOSS["XV"] = -1
	// BLOSS["VX"] = -1
	// BLOSS["XB"] = -1
	// BLOSS["BX"] = -1
	// BLOSS["XZ"] = -1
	// BLOSS["ZX"] = -1

	// fmt.Printf("%v %v %v\n", refAA, altAA, Tang[buffer.String()])
	if altAA != "" && refAA != "" && altAA != "X" {
		res = blossIdx[buffer.String()]["BLOSS"]
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

	// fmt.Printf("%v %v %v\n", refAA, altAA, Tang[buffer.String()])
	if altAA != "" && refAA != "" && altAA != "X" {
		res = pam30Idx[buffer.String()]["PAM30"]
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

var aaLong = map[string]map[string]string{
	"Ala": {
		"Name":  "Alanine",
		"Short": "A",
	},
	"Arg": {
		"Name":  "Arginine",
		"Short": "R",
	},
	"Asn": {
		"Name":  "Asparagine",
		"Short": "N",
	},
	"Asp": {
		"Name":  "Aspartic acid",
		"Short": "D",
	},
	"Cys": {
		"Name":  "Cysteine",
		"Short": "C",
	},
	"Gln": {
		"Name":  "Glutamine",
		"Short": "Q",
	},
	"Glu": {
		"Name":  "Glutamic acid",
		"Short": "E",
	},
	"Gly": {
		"Name":  "Glycine",
		"Short": "G",
	},
	"His": {
		"Name":  "Histidine",
		"Short": "H",
	},
	"Ile": {
		"Name":  "Isoleucine",
		"Short": "I",
	},
	"Leu": {
		"Name":  "Leucine",
		"Short": "L",
	},
	"Lys": {
		"Name":  "Lysisne",
		"Short": "K",
	},
	"Met": {
		"Name":  "Methionine",
		"Short": "M",
	},
	"Phe": {
		"Name":  "Phenylalanine",
		"Short": "F",
	},
	"Pro": {
		"Name":  "Proline",
		"Short": "P",
	},
	"Ser": {
		"Name":  "Serine",
		"Short": "S",
	},
	"Thr": {
		"Name":  "Threonine",
		"Short": "T",
	},
	"Trp": {
		"Name":  "Tryptophan",
		"Short": "W",
	},
	"Tyr": {
		"Name":  "Tyrosine",
		"Short": "Y",
	},
	"Val": {
		"Name":  "Valine",
		"Short": "V",
	},
}

func GetShortNameAA(longName string) (shortName string) {
	if len(longName) == 3 {
		shortName =aaLong[longName]["Short"]
	}
	return shortName
}
