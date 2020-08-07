package codon

import (
	"fmt"
	"math"

	// "math"
	// "math/big"
	// "regexp"
	"strconv"
	"strings"

	"../amino"
)

// //CodonCompInf is.....
// type CodonCompInf struct {

// }

//DnDs struct is....
type DnDs struct {
	N, S, PN, PS, DN, DS, DNDS, ND, NS float64
	// ND, NS                     int
}

type DnDsQuery struct {
	OutChan        chan DnDs
	RefSeq, AltSeq string
}

func (q *DnDsQuery) Request() {
	q.OutChan <- CalcDnDs(q.RefSeq, q.AltSeq)

}

var (
	nucs = []string{"C", "T", "G", "A"}
)

// CalcDnDs is....
func CalcDnDs(refSeq string, altSeq string) DnDs {
	// var codonPos int
	var refSeqCodons = make(map[int]string)
	var altSeqCodons = make(map[int]string)
	// var mutType string

	// var expNS []float64
	var expN, expS float64

	// var expN, expS float64
	// var Nd, Ns int
	// threeFour := big.NewRat(3, 4)
	// threeFour, _ := threeFour.Float64()
	// threeFour := float64(3) / 4
	// re := regexp.MustCompile("\\w{3}")

	// refCodons := re.FindAllStringSubmatch(refSeq, -1)
	// altCodons := re.FindAllStringSubmatch(altSeq, -1)
	refAlphabet := strings.Split(refSeq, "")
	altAlphabet := strings.Split(altSeq, "")
	var testRefCodons []string
	var testAltCodons []string
	j := 0
	var codon []string
	for i := 0; i < len(refAlphabet); i++ {
		j++
		codon = append(codon, refAlphabet[i])
		if j == 3 {
			j = 0
			// fmt.Println(strings.Join(codon, ""))
			testRefCodons = append(testRefCodons, strings.Join(codon, ""))
			codon = nil
		}

	}
	for i := 0; i < len(altAlphabet); i++ {
		j++
		codon = append(codon, altAlphabet[i])
		if j == 3 {
			j = 0
			// fmt.Println(strings.Join(codon, ""))
			testAltCodons = append(testAltCodons, strings.Join(codon, ""))
			codon = nil
		}

	}
	// fmt.Println(testAltCodons, testRefCodons)
	// refCodons := testRefCodons
	// altCodons := testAltCodons

	for i := 0; i < len(testRefCodons); i++ {
		refSeqCodons[i] = strings.ToUpper(testRefCodons[i])
		expNSReq := make(chan []float64)
		go func() {
			// fmt.Println("Gorutine 1 started")
			// expNSReq <- expectedSites(strings.ToUpper(rCodon[0]))
			expNSReq <- expectedSites(strings.ToUpper(testRefCodons[i]))
		}()
		// expNS = expectedSites(strings.ToUpper(rCodon[0]))
		expNS := <-expNSReq
		close(expNSReq)

		expN = expN + expNS[0] //N
		expS = expS + expNS[1] //S
		// fmt.Println(expNS[0], expNS[1])
	}

	// for i, rCodon := range refCodons {
	// 	refSeqCodons[i] = strings.ToUpper(rCodon[0])
	// 	expNSReq := make(chan []float64)
	// 	go func() {
	// 		// fmt.Println("Gorutine 1 started")
	// 		expNSReq <- expectedSites(strings.ToUpper(rCodon[0]))
	// 	}()
	// 	// expNS = expectedSites(strings.ToUpper(rCodon[0]))
	// 	expNS := <-expNSReq
	// 	close(expNSReq)
	// 	expN = expN + expNS[0] //N
	// 	expS = expS + expNS[1] //S

	// }

	// for i, aCodon := range altCodons {
	// 	altSeqCodons[i] = strings.ToUpper(aCodon[0])
	// }

	for i := 0; i < len(testAltCodons); i++ {
		altSeqCodons[i] = strings.ToUpper(testAltCodons[i])
	}

	ndnsResChan := make(chan []float64)

	go func() {
		// fmt.Println("Gorutine 2 started")
		// ndnsResChan <- getNdNs(refSeqCodons, altSeqCodons)
		ndnsResChan <- getNdNsTest(testRefCodons, testAltCodons)
	}()

	ndnsRes := <-ndnsResChan
	// close(ndnsResChan)

	ndnsRes0, _ := strconv.ParseFloat(fmt.Sprintf("%8f", ndnsRes[0]), 64)
	ndnsRes1, _ := strconv.ParseFloat(fmt.Sprintf("%8f", ndnsRes[1]), 64)
	// Nd, Ns := ndnsRes[0], ndnsRes[1]
	Nd, Ns := ndnsRes0, ndnsRes1

	pN, _ := strconv.ParseFloat(fmt.Sprintf("%8f", Nd/expN), 64)
	pS, _ := strconv.ParseFloat(fmt.Sprintf("%8f", Ns/expS), 64)

	// test1 := fmt.Sprintf("%8f", pN)
	// if s, err := strconv.ParseFloat(test1, 64); err == nil {
	// 	fmt.Printf("%T, %v\n", s, s)
	// } else {
	// 	fmt.Println(err)
	// }
	// test2, _ := strconv.ParseFloat(test1, 64)
	// fmt.Printf("%v", test1)
	dnLong := -float64(3) / 4 * (math.Log(float64(1 - (4 * float64(pN) / 3))))
	dsLong := -float64(3) / 4 * (math.Log(float64(1 - (4 * float64(pS) / 3))))

	dn, _ := strconv.ParseFloat(fmt.Sprintf("%8f", dnLong), 64)
	ds, _ := strconv.ParseFloat(fmt.Sprintf("%8f", dsLong), 64)

	dNdS := dn / ds

	newDNDS := DnDs{N: expN, S: expS, ND: Nd, NS: Ns, PN: pN, PS: pS, DN: dn, DS: ds, DNDS: dNdS}

	// fmt.Printf("N:%v S:%v ND:%v NS:%v PN:%v PS:%v DN:%v DS:%v\n", expN, expS, Nd, Ns, pN, pS, dn, ds)
	// fmt.Println(dNdS, newDNDS)
	return newDNDS
}

func getNdNsTest(refSeqCodons []string, altSeqCodons []string) (ndnsRes []float64) {
	var (
		// nd, sd []int
		Nd = 0.0
		Ns = 0.0
	)

	for i := 0; i < len(refSeqCodons); i++ {

		_, refAA := amino.Codon2AA(refSeqCodons[i])
		_, altAA := amino.Codon2AA(altSeqCodons[i])
		_, diff := codonCompare(refSeqCodons[i], altSeqCodons[i])

		if refAA == altAA && diff == 1 {
			// mutType = "s"
			// sd = append(sd, 1)
			Ns++
			// fmt.Println(refAA, altAA, diff, refSeqCodons[i], altSeqCodons[i], "Syn")
		} else if refAA != altAA && diff == 1 {
			// mutType = "n"
			// nd = append(nd, 1)
			Nd++
			// fmt.Println(refAA, altAA, diff, refSeqCodons[i], altSeqCodons[i], "Nonsyn")
		} else if diff == 2 {
			// mutType = "n/s"
			// nd = append(nd, 1)
			// sd = append(sd, 1)
			Nd++
			Ns++
			// fmt.Println(refAA, altAA, diff, refSeqCodons[i], altSeqCodons[i], "Both")
		} else if diff == 3 {
			// mutType = "n/s"
			// nd = append(nd, 1)
			// sd = append(sd, 1)
			Nd = Nd + 2.333
			Ns = Ns + 0.666
			// fmt.Println(refAA, altAA, diff, refSeqCodons[i], altSeqCodons[i], "Total")
		}
	}

	// }

	// for key, val := range refSeqCodons {
	// 	for key1, val1 := range altSeqCodons {
	// 		if key == key1 && val != val1 {
	// 			_, refAA := amino.Codon2AA(val)
	// 			_, altAA := amino.Codon2AA(val1)
	// 			_, diff := codonCompare(val, val1)
	// 			// fmt.Println(refAA, altAA)

	// 			if refAA == altAA && diff == 1 {
	// 				// mutType = "s"
	// 				// sd = append(sd, 1)
	// 				Ns++
	// 			} else if refAA != altAA && diff == 1 {
	// 				// mutType = "n"
	// 				// nd = append(nd, 1)
	// 				Nd++
	// 			} else if refAA != altAA && diff > 1 || refAA == altAA && diff > 1 {
	// 				// mutType = "n/s"
	// 				// nd = append(nd, 1)
	// 				// sd = append(sd, 1)
	// 				Nd++
	// 				Ns++
	// 			}

	// 		}
	// 	}
	// }

	// for i := 0; i < len(nd); i++ {
	// 	Nd += nd[i]
	// }
	// for _, value := range nd {
	// 	Nd += value
	// }

	// for _, value := range sd {
	// 	Ns += value
	// }
	ndnsRes = append(ndnsRes, Nd) //Nd
	ndnsRes = append(ndnsRes, Ns) //Ns
	return ndnsRes

}

func getNdNs(refSeqCodons map[int]string, altSeqCodons map[int]string) (ndnsRes []int) {
	var (
		// nd, sd []int
		Nd = 0
		Ns = 0
	)

	for key, val := range refSeqCodons {
		for key1, val1 := range altSeqCodons {
			if key == key1 && val != val1 {
				_, refAA := amino.Codon2AA(val)
				_, altAA := amino.Codon2AA(val1)
				_, diff := codonCompare(val, val1)
				// fmt.Println(refAA, altAA)

				if refAA == altAA && diff == 1 {
					// mutType = "s"
					// sd = append(sd, 1)
					Ns++
				} else if refAA != altAA && diff == 1 {
					// mutType = "n"
					// nd = append(nd, 1)
					Nd++
				} else if refAA != altAA && diff > 1 || refAA == altAA && diff > 1 {
					// mutType = "n/s"
					// nd = append(nd, 1)
					// sd = append(sd, 1)
					Nd++
					Ns++
				}

			}
		}
	}

	// for i := 0; i < len(nd); i++ {
	// 	Nd += nd[i]
	// }
	// for _, value := range nd {
	// 	Nd += value
	// }

	// for _, value := range sd {
	// 	Ns += value
	// }
	ndnsRes = append(ndnsRes, Nd) //Nd
	ndnsRes = append(ndnsRes, Ns) //Ns
	return ndnsRes

}

func codonCompare(refCodon, altCodon string) (int, int) {

	codon1Nuc := strings.SplitAfter(refCodon, "")
	codon2Nuc := strings.SplitAfter(altCodon, "")
	// fmt.Println(res)
	var diffC, exactC int
	if altCodon == "TAG" || altCodon == "TAA" || altCodon == "TGA" {
		exactC = 0
		diffC = 0
	} else {
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

func expectedSites(codon string) (expNS []float64) {
	var buffer strings.Builder
	var codons []string
	var n, s float64
	// oneThree := big.NewRat(1, 3)
	// oneThree2Float, _ := oneThree.Float64()
	_, refAA := amino.Codon2AA(codon)
	codon2Nuc := strings.SplitAfter(codon, "")
	// fmt.Println(codon2Nuc)
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
			n = n + float64(1)/3
			// fmt.Printf("%.5g\n", n)
		} else if refAA == altAA {
			n = n + 0
			// fmt.Printf("%.5g\n", n)
		}

	}
	s = 3 - n
	s, _ = strconv.ParseFloat(fmt.Sprintf("%4f", s), 64)
	n, _ = strconv.ParseFloat(fmt.Sprintf("%4f", n), 64)
	expNS = append(expNS, n)
	expNS = append(expNS, s)
	return expNS
}
