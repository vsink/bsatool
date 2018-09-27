package codon

import (
	"fmt"
	"math"

	// "math"
	"math/big"
	"regexp"
	"strconv"
	"strings"

	"../amino"
)

// //CodonCompInf is.....
// type CodonCompInf struct {

// }

//DnDs struct is....
type DnDs struct {
	N, S, PN, PS, DN, DS, DNDS float64
	ND, NS                     int
}

type DnDsQuery struct {
	OutChan        chan DnDs
	RefSeq, AltSeq string
}

func (q *DnDsQuery) Request() {
	q.OutChan <- CalcDnDs(q.RefSeq, q.AltSeq)

}

// CalcDnDs is....
func CalcDnDs(refSeq string, altSeq string) DnDs {
	// var codonPos int
	var refSeqCodons = make(map[int]string)
	var altSeqCodons = make(map[int]string)
	// var mutType string
	var nd []int
	var sd []int
	// var expN, expS float64
	var N, S, expN, expS float64
	threeFour := big.NewRat(3, 4)
	threeFour2Float, _ := threeFour.Float64()

	re := regexp.MustCompile("\\w{3}")

	refCodons := re.FindAllStringSubmatch(refSeq, -1)
	altCodons := re.FindAllStringSubmatch(altSeq, -1)
	for i, rCodon := range refCodons {
		refSeqCodons[i] = strings.ToUpper(rCodon[0])
		N, S = expectedSites(strings.ToUpper(rCodon[0]))

		expN = expN + N
		expS = expS + S

	}
	for i, aCodon := range altCodons {
		altSeqCodons[i] = strings.ToUpper(aCodon[0])
	}
	for key, val := range refSeqCodons {
		for key1, val1 := range altSeqCodons {
			if key == key1 && val == val1 {

				// fmt.Println("the same:", key+1, val, val1)
			} else if key == key1 && val != val1 {
				_, refAA := amino.Codon2AA(val)
				_, altAA := amino.Codon2AA(val1)
				_, diff := codonCompare(val, val1)
				if refAA == altAA && diff == 1 {
					// mutType = "s"
					sd = append(sd, 1)
				} else if refAA != altAA && diff == 1 {
					// mutType = "n"
					nd = append(nd, 1)
				} else if refAA != altAA && diff > 1 || refAA == altAA && diff > 1 {
					// mutType = "n/s"
					nd = append(nd, 1)
					sd = append(sd, 1)
				}

			}
		}
	}
	Nd := 0
	for _, value := range nd {
		Nd += value
	}
	Ns := 0
	for _, value := range sd {
		Ns += value
	}
	pN := float64(Nd) / expN
	pS := float64(Ns) / expS
	dn := -threeFour2Float * math.Log(float64(1-(4*float64(pN)/3)))
	ds := -threeFour2Float * math.Log(float64(1-(4*float64(pS)/3)))
	dNdS := dn / ds

	newDNDS := DnDs{N: expN, S: expS, ND: Nd, NS: Ns, PN: pN, PS: pS, DN: dn, DS: ds, DNDS: dNdS}
	return newDNDS
}

func codonCompare(codon1, codon2 string) (int, int) {

	codon1Nuc := strings.SplitAfter(codon1, "")
	codon2Nuc := strings.SplitAfter(codon2, "")
	// fmt.Println(res)
	var diffC, exactC int

	for j, C1 := range codon1Nuc {

		for jj, C2 := range codon2Nuc {

			if C1 == C2 && j == jj {
				exactC++
				// fmt.Println(C1, C2, j+1, "pos is exact")
			} else if C1 != C2 && j == jj {
				diffC++
				// fmt.Println(C1, C2, j+1, "pos is diff")
			}
		}
	}

	return exactC, diffC
	// fmt.Printf("diff:%v exact:%v\n", diffC, exactC)
}

// GcCodonCalc is ....
func GcCodonCalc(seq string) (float64, float64, float64, float64) {
	var codonsMap = make(map[int]string)
	var i, gcCount1, gcCount2, gcCount3, gcCount int
	var gc1, gc2, gc3, gc float64
	// var codons []string
	var buffer strings.Builder
	for _, val := range seq {

		// if i < 2 {
		// codons = append(codons, string(val))

		buffer.WriteString(strings.ToUpper(string(val)))
		if buffer.Len() == 3 {
			i++
			// codons = append(codons, buffer.String())
			codonsMap[i] = buffer.String()
			// fmt.Println(buffer.String())
			buffer.Reset()
		}

	}

	for _, val := range codonsMap {

		for pos, nuc := range val {

			if pos == 0 && nuc == 'C' || pos == 0 && nuc == 'G' {
				gcCount1 = gcCount1 + 1

			} else if pos == 1 && nuc == 'C' || pos == 1 && nuc == 'G' {
				gcCount2 = gcCount2 + 1

			} else if pos == 2 && nuc == 'C' || pos == 2 && nuc == 'G' {
				gcCount3 = gcCount3 + 1

			}
			if nuc == 'C' || nuc == 'G' {
				gcCount = gcCount + 1

			}
			// fmt.Println(pos, string(nuc))
		}
	}
	gc = (float64(gcCount) / float64(len(codonsMap)*3)) * 100
	gc1 = (float64(gcCount1) / float64(len(codonsMap)*3)) * 100
	gc2 = (float64(gcCount2) / float64(len(codonsMap)*3)) * 100
	gc3 = (float64(gcCount3) / float64(len(codonsMap)*3)) * 100
	// fmt.Printf("gc:%.2f%% gc1:%.2f%% gc2:%.2f%% gc3:%.2f%%\n", gc, gc1, gc2, gc3)
	// fmt.Println(codonsMap)
	return gc, gc1, gc2, gc3

}

func expectedSites(codon string) (float64, float64) {
	var buffer strings.Builder
	var codons []string
	var n, s float64
	oneThree := big.NewRat(1, 3)
	oneThree2Float, _ := oneThree.Float64()
	_, refAA := amino.Codon2AA(codon)
	codon2Nuc := strings.SplitAfter(codon, "")
	nucs := []string{"C", "T", "G", "A"}

	// fmt.Printf("codon:%v\n", refAA)
	for i, val := range codon2Nuc {
		//pos 1
		if i == 0 {
			for _, nuc := range nucs {
				if nuc != val {
					buffer.WriteString(nuc)
					buffer.WriteString(strings.Join(codon2Nuc[1:3], ""))
					if buffer.Len() == 3 {
						codons = append(codons, buffer.String())
						// fmt.Println(buffer.String())
						buffer.Reset()
					}
					// fmt.Printf("pos:%v codon:%v %v\n", i, nuc, strings.Join(codon2Nuc[1:3], ""))
				}
			}
		} else if i == 1 {
			for _, nuc := range nucs {
				if nuc != val {
					buffer.WriteString(string(codon2Nuc[0]))
					buffer.WriteString(nuc)
					buffer.WriteString(strings.Join(codon2Nuc[2:3], ""))
					if buffer.Len() == 3 {
						codons = append(codons, buffer.String())
						// fmt.Println(buffer.String())
						buffer.Reset()
					}
				}
			}
		} else if i == 2 {
			for _, nuc := range nucs {
				if nuc != val {
					buffer.WriteString(strings.Join(codon2Nuc[0:2], ""))
					buffer.WriteString(nuc)

					if buffer.Len() == 3 {
						codons = append(codons, buffer.String())
						// fmt.Println(buffer.String())
						buffer.Reset()
					}
				}
			}
		}
	}
	for _, val := range codons {
		_, altAA := amino.Codon2AA(val)
		if refAA != altAA {
			n = n + oneThree2Float
			// fmt.Printf("%.5g\n", n)
		} else if refAA == altAA {
			n = n + 0
			// fmt.Printf("%.5g\n", n)
		}

	}
	s = 3 - n
	s, _ = strconv.ParseFloat(fmt.Sprintf("%.5g", s), 64)
	n, _ = strconv.ParseFloat(fmt.Sprintf("%.5g", n), 64)

	return n, s
}
