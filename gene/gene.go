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
var AllGenesVal []Gene

// var GenomeSeqSlice []string // информация об генах, загруженная из базы

//Gene is structrue of imported from GB file fields
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
	GOA       string `json:"goa"`
}

// SNPinfo is structure of annotatted SNP
type SNPinfo struct {
	/*
				APos абсолютная позиция в геноме
					PosInGene позиция в гене
		PosInCodonG позиция в буквы в кодоне (0-первая, 1-средняя, 2-последняя)

	*/
	APos, PosInGene, PosInCodonG, CodonNbrInG, GeneLen int
	RefCodon, AltCodon, RefAA, AltAA, Locus,
	Direction, NucInPos, Product, Name, Start,
	RefAAShort, AltAAShort, End, Mutation, Tang, Alt, Note, ReportType, ProteinID, GeneID, GOA, TiTv string
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
