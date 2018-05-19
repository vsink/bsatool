package gene

// "bufio"
// "bytes"
// "fmt"
// "io"
// "os"
// "regexp"
// "strings"

// "fmt"
// "strings"

//AllGenesVal массив информации об генах, загруженный из файла базы данных
// var AllGenesVal []Gene

// type SNP interface {
// 	ReadGB()
// }
// var AllGenesVal []Gene

// var GenomeSeqSlice []string // информация об генах, загруженная из базы

//Gene is structrue of imported from GB file fields
type Gene struct {
	Locus, Name, Product, Direction, GeneID, ProteinID, Note, GOA, TypeOf string
	Start, End                                                            int
}

// SNPinfo is structure of annotatted SNP
type SNPinfo struct {
	/*
				APos абсолютная позиция в геноме
					PosInGene позиция в гене
		PosInCodonG позиция в буквы в кодоне (0-первая, 1-средняя, 2-последняя)

	*/
	APos, PosInGene, PosInCodonG, CodonNbrInG, GeneLen, Start, End int
	RefCodon, AltCodon, RefAA, AltAA, Locus,
	Direction, NucInPos, Product, Name,
	RefAAShort, AltAAShort, Mutation, Tang, Alt, Note, ReportType, ProteinID, GeneID, GOA, TiTv, TypeOf string
}

type INDELinfo struct {
	/*
	 */
	Apos     int
	Ref, Alt string
}

type SNPcheck struct {
	Locus, PosInGene, CodonNbrInG, Ref, Alt, Name, TypeOf, APos, AASref, AASalt, AALref, AALalt, Raw string
}

type GeneInfo struct {
	GeneLen   int
	Direction string
}
