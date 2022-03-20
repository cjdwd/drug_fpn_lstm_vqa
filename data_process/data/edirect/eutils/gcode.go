// ===========================================================================
//
//                            PUBLIC DOMAIN NOTICE
//            National Center for Biotechnology Information (NCBI)
//
//  This software/database is a "United States Government Work" under the
//  terms of the United States Copyright Act. It was written as part of
//  the author's official duties as a United States Government employee and
//  thus cannot be copyrighted. This software/database is freely available
//  to the public for use. The National Library of Medicine and the U.S.
//  Government do not place any restriction on its use or reproduction.
//  We would, however, appreciate having the NCBI and the author cited in
//  any work or product based on this material.
//
//  Although all reasonable efforts have been taken to ensure the accuracy
//  and reliability of the software and data, the NLM and the U.S.
//  Government do not and cannot warrant the performance or results that
//  may be obtained by using this software or data. The NLM and the U.S.
//  Government disclaim all warranties, express or implied, including
//  warranties of performance, merchantability or fitness for any particular
//  purpose.
//
// ===========================================================================
//
// File Name:  gcode.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"fmt"
	"os"
	"sort"
	"strings"
)

// When a new genetic code is added, update the raw genetic code tables below,
// recompile, run 'transmute -degenerate > gdata.go' to generate new data tables,
// replace the old gdata.go source code file with the new version, and then
// recompile again to bring all executables up to date.

// bases are ncbi4na bit flags: A 0, C 1, G 2, T 4

const (
	//                                      tgca
	baseGap = iota //    0     -    -       0000
	baseA          //    1     A    A       0001
	baseC          //    2     C    C       0010
	baseM          //    3     M    AC      0011
	baseG          //    4     G    G       0100
	baseR          //    5     R    AG      0101
	baseS          //    6     S    CG      0110
	baseV          //    7     V    ACG     0111
	baseT          //    8     T    T       1000
	baseW          //    9     W    AT      1001
	baseY          //    10    Y    CT      1010
	baseH          //    11    H    ACT     1011
	baseK          //    12    K    GT      1100
	baseD          //    13    D    AGT     1101
	baseB          //    14    B    CGT     1110
	baseN          //    15    N    ACGT    1111
)

// Base
//    1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
//    2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
//    3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

var ncbieaaCode = map[int]string{
	1:  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	2:  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
	3:  "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	4:  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	5:  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
	6:  "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	9:  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	10: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	12: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	13: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
	14: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	15: "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	16: "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	21: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	22: "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	23: "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	24: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
	25: "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	26: "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	27: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	28: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	29: "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	30: "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	31: "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	32: "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	33: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
}

// Base
//    1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
//    2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
//    3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

var sncbieaaCode = map[int]string{
	1:  "---M------**--*----M---------------M----------------------------",
	2:  "----------**--------------------MMMM----------**---M------------",
	3:  "----------**----------------------MM---------------M------------",
	4:  "--MM------**-------M------------MMMM---------------M------------",
	5:  "---M------**--------------------MMMM---------------M------------",
	6:  "--------------*--------------------M----------------------------",
	9:  "----------**-----------------------M---------------M------------",
	10: "----------**-----------------------M----------------------------",
	11: "---M------**--*----M------------MMMM---------------M------------",
	12: "----------**--*----M---------------M----------------------------",
	13: "---M------**----------------------MM---------------M------------",
	14: "-----------*-----------------------M----------------------------",
	15: "----------*---*--------------------M----------------------------",
	16: "----------*---*--------------------M----------------------------",
	21: "----------**-----------------------M---------------M------------",
	22: "------*---*---*--------------------M----------------------------",
	23: "--*-------**--*-----------------M--M---------------M------------",
	24: "---M------**-------M---------------M---------------M------------",
	25: "---M------**-----------------------M---------------M------------",
	26: "----------**--*----M---------------M----------------------------",
	27: "--------------*--------------------M----------------------------",
	28: "----------**--*--------------------M----------------------------",
	29: "--------------*--------------------M----------------------------",
	30: "--------------*--------------------M----------------------------",
	31: "----------**-----------------------M----------------------------",
	32: "---M------*---*----M------------MMMM---------------M------------",
	33: "---M-------*-------M---------------M---------------M------------",
}

var genCodeNames = map[int]string{
	1:  "Standard",
	2:  "Vertebrate Mitochondrial",
	3:  "Yeast Mitochondrial",
	4:  "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma",
	5:  "Invertebrate Mitochondrial",
	6:  "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
	9:  "Echinoderm Mitochondrial; Flatworm Mitochondrial",
	10: "Euplotid Nuclear",
	11: "Bacterial, Archaeal and Plant Plastid",
	12: "Alternative Yeast Nuclear",
	13: "Ascidian Mitochondrial",
	14: "Alternative Flatworm Mitochondrial",
	15: "Blepharisma Macronuclear",
	16: "Chlorophycean Mitochondrial",
	21: "Trematode Mitochondrial",
	22: "Scenedesmus obliquus Mitochondrial",
	23: "Thraustochytrium Mitochondrial",
	24: "Rhabdopleuridae Mitochondrial",
	25: "Candidate Division SR1 and Gracilibacteria",
	26: "Pachysolen tannophilus Nuclear",
	27: "Karyorelict Nuclear",
	28: "Condylostoma Nuclear",
	29: "Mesodinium Nuclear",
	30: "Peritrich Nuclear",
	31: "Blastocrithidia Nuclear",
	32: "Balanophoraceae Plastid",
	33: "Cephalodiscidae Mitochondrial",
}

// finite state machine into lookup table for simultaneous two-strand translation

func correctGenCode(genCode int) int {

	switch genCode {
	case 0:
		genCode = 1
	case 7:
		genCode = 4
	case 8:
		genCode = 1
	}

	return genCode
}

// GenCodeName returns full name of genetic code
func GenCodeName(genCode int) string {

	genCode = correctGenCode(genCode)

	return genCodeNames[genCode]
}

// SetCodonState initializes a state from three nucleotide letters
func SetCodonState(ch1, ch2, ch3 int) int {

	if ch1 < 0 || ch1 > 255 || ch2 < 0 || ch2 > 255 || ch3 < 0 || ch3 > 255 {
		return 0
	}

	return 256*baseToIdx[ch1] + 16*baseToIdx[ch2] + baseToIdx[ch3]
}

// GetCodonFromState returns the three-letter codon from the state value
func GetCodonFromState(state int) string {

	if state < 0 || state > 4095 {
		return ""
	}

	upperToBase := "-ACMGRSVTWYHKDBN"

	ch3 := upperToBase[state&15]
	state = state >> 4
	ch2 := upperToBase[state&15]
	state = state >> 4
	ch1 := upperToBase[state&15]

	str := string(ch1) + string(ch2) + string(ch3)

	return str
}

// NextCodonState uses the old state plus the new letter
func NextCodonState(state, ch int) int {

	if state < 0 || state > 4095 {
		return 0
	}

	return nextState[state] + baseToIdx[ch]
}

// RevCompState finds the state for the last three bases on the minus strand
func RevCompState(state int) int {

	if state < 0 || state > 4095 {
		return 0
	}

	return rvCpState[state]
}

// translation codon lookup functions

// GetCodonResidue returns the amino acid for the state representing the current codon
func GetCodonResidue(genCode, state int) int {

	genCode = correctGenCode(genCode)

	if state < 0 || state > 4095 {
		return 'X'
	}

	// look for override in alternative genetic code
	gc, ok := aminoAcidMaps[genCode]
	if ok {
		aa, ok := gc[state]
		if ok {
			return aa
		}
	}

	if genCode != 1 {
		// now fetch residue from standard genetic code
		gc, ok = aminoAcidMaps[1]
		if ok {
			aa, ok := gc[state]
			if ok {
				return aa
			}
		}
	}

	// no entry in either map defaults to unknown residue X
	return 'X'
}

// GetStartResidue returns the amino acid for the first codon
func GetStartResidue(genCode, state int) int {

	genCode = correctGenCode(genCode)

	if state < 0 || state > 4095 {
		return '-'
	}

	gc, ok := orfStartMaps[genCode]
	if ok {
		aa, ok := gc[state]
		if ok {
			return aa
		}
	}

	return '-'
}

// GetStopResidue returns the amino acid for the last codon
func GetStopResidue(genCode, state int) int {

	genCode = correctGenCode(genCode)

	if state < 0 || state > 4095 {
		return '-'
	}

	gc, ok := orfStopMaps[genCode]
	if ok {
		aa, ok := gc[state]
		if ok {
			return aa
		}
	}

	return '-'
}

// IsOrfStart returns true for a start codon
func IsOrfStart(genCode, state int) bool {

	genCode = correctGenCode(genCode)

	if state < 0 || state > 4095 {
		return false
	}

	return GetStartResidue(genCode, state) == 'M'
}

// IsOrfStop returns true for a stop codon
func IsOrfStop(genCode, state int) bool {

	genCode = correctGenCode(genCode)

	if state < 0 || state > 4095 {
		return false
	}

	return GetStopResidue(genCode, state) == '*'
}

// IsATGStart returns true for an ATG start codon
func IsATGStart(genCode, state int) bool {

	genCode = correctGenCode(genCode)

	const ATGState = 388

	if state < 0 || state > 4095 {
		return false
	}

	return IsOrfStart(genCode, state) && state == ATGState
}

// IsAltStart returns true for a non-ATG start codon
func IsAltStart(genCode, state int) bool {

	genCode = correctGenCode(genCode)

	const ATGState = 388

	if state < 0 || state > 4095 {
		return false
	}

	return IsOrfStart(genCode, state) && state != ATGState
}

// TranslateCdRegion converts a single coding region to a protein sequence
func TranslateCdRegion(seq string, genCode, frame int, includeStop, doEveryCodon, removeTrailingX, is5primeComplete, is3primeComplete bool, between string) string {

	genCode = correctGenCode(genCode)

	usableSize := len(seq) - frame
	mod := usableSize % 3
	length := usableSize / 3
	checkStart := is5primeComplete && frame == 0

	aa := 0
	state := 0
	// start_state := 0
	firstTime := true

	// standard and bacterial code use built-in array for fast access
	fast := genCode == 1 || genCode == 11

	pos := frame

	var aminos []string

	// first codon has extra logic, process separately
	if length > 0 {

		// loop through first codon to initialize state
		for k := 0; k < 3; k++ {

			// internal version of NextCodonState without value reality check
			ch := int(seq[pos])
			state = nextState[state] + baseToIdx[ch]
			// state = NextCodonState(state, int(seq[pos]))
			pos++
		}

		// start_state = state

		// save first translated amino acid
		if checkStart {
			aa = GetStartResidue(genCode, state)
		} else {
			aa = GetCodonResidue(genCode, state)
		}

		aminos = append(aminos, string(rune(aa)))

		firstTime = false
	}

	// start at offset 1, continue with rest of sequence
	for i := 1; i < length; i++ {

		// loop through next triplet codon
		for k := 0; k < 3; k++ {

			// internal version of NextCodonState without value reality check
			ch := int(seq[pos])
			state = nextState[state] + baseToIdx[ch]
			// state = NextCodonState(state, int(seq[pos]))
			pos++
		}

		// save translated amino acid
		if fast {
			// use high-efficiency built-in lookup table for standard code
			aa = int(stdGenCode[state])
		} else {
			aa = GetCodonResidue(genCode, state)
		}

		aminos = append(aminos, string(rune(aa)))

		firstTime = false
	}

	if mod > 0 {

		var last []string

		k := 0
		for ; k < mod; k++ {
			ch := int(seq[pos])
			state = NextCodonState(state, ch)
			pos++
			last = append(last, string(rune(ch)))
		}

		for ; k < 3; k++ {
			state = NextCodonState(state, int('N'))
			last = append(last, "N")
		}

		if firstTime {
			// start_state = state
		}

		// save translated amino acid
		ch := GetCodonResidue(genCode, state)

		aa = 0
		if firstTime && checkStart {
			aa = GetStartResidue(genCode, state)
		} else if ch != 'X' {
			aa = GetCodonResidue(genCode, state)
		}

		if aa != 0 {
			aminos = append(aminos, string(rune(aa)))
		}
	}

	txt := strings.Join(aminos, between)

	// check for stop codon that normally encodes an amino acid
	if aa != '*' && includeStop && mod > 0 && len(txt) > 0 && is3primeComplete {

		aa = GetStopResidue(genCode, state)
		if aa == '*' {
			rgt := len(txt) - 1
			txt = txt[:rgt] + "*"
		}
	}

	if doEveryCodon {

		// read through all stop codons

	} else if includeStop {

		idx := strings.Index(txt, "*")
		if idx >= 0 && idx+1 < len(txt) {
			txt = txt[:idx+1]
		}

	} else {

		idx := strings.Index(txt, "*")
		if idx >= 0 {
			txt = txt[:idx]
		}
	}

	if removeTrailingX {

		for strings.HasSuffix(txt, "X") || strings.HasSuffix(txt, "*") {
			txt = strings.TrimRight(txt, "*")
			txt = strings.TrimRight(txt, "X")
		}
	}

	return txt
}

// NucProtCodonReport displays triplet codons above the translated amino acid
func NucProtCodonReport(nuc, prt string, frame int, threeLetter bool) string {

	usableSize := len(nuc) - frame
	mod := usableSize % 3
	length := usableSize / 3

	if length < 1 {
		return ""
	}

	pos := frame
	lag := pos

	var codons []string

	for i := 0; i < length; i++ {
		pos += 3
		codons = append(codons, nuc[lag:pos])
		lag = pos
	}

	if mod > 0 {

		var last []string

		k := 0
		for ; k < mod; k++ {
			ch := int(nuc[pos])
			pos++
			last = append(last, string(rune(ch)))
		}

		lag = pos

		for ; k < 3; k++ {
			last = append(last, "N")
		}

		cdn := strings.Join(last, "")
		codons = append(codons, cdn)
	}

	var buffer strings.Builder
	cumulative := 0

	for len(codons) > 0 || len(prt) > 0 {
		k := 0
		for ; k < 15 && len(codons) > 0; k++ {
			cdn := codons[0]
			codons = codons[1:]
			buffer.WriteString(cdn)
			buffer.WriteString(" ")
		}
		cumulative += k
		for ; k < 15; k++ {
			buffer.WriteString("    ")
		}
		num := fmt.Sprintf("%6d", cumulative)
		buffer.WriteString(num)
		buffer.WriteString("\n")
		k = 0
		for ; k < 15 && len(prt) > 0; k++ {
			res := string(prt[0])
			prt = prt[1:]
			if threeLetter {
				res = aaTo3[res]
			} else {
				res = " " + res + " "
			}
			buffer.WriteString(res)
			buffer.WriteString(" ")
		}
		buffer.WriteString("\n\n")
	}

	txt := buffer.String()

	return txt
}

const gdataHeading = "" +
	`// ===========================================================================
//
//                            PUBLIC DOMAIN NOTICE
//            National Center for Biotechnology Information (NCBI)
//
//  This software/database is a "United States Government Work" under the
//  terms of the United States Copyright Act. It was written as part of
//  the author's official duties as a United States Government employee and
//  thus cannot be copyrighted. This software/database is freely available
//  to the public for use. The National Library of Medicine and the U.S.
//  Government do not place any restriction on its use or reproduction.
//  We would, however, appreciate having the NCBI and the author cited in
//  any work or product based on this material.
//
//  Although all reasonable efforts have been taken to ensure the accuracy
//  and reliability of the software and data, the NLM and the U.S.
//  Government do not and cannot warrant the performance or results that
//  may be obtained by using this software or data. The NLM and the U.S.
//  Government disclaim all warranties, express or implied, including
//  warranties of performance, merchantability or fitness for any particular
//  purpose.
//
// ===========================================================================
//
// File Name:  gdata.go
//
// Author:  transmute -degenerate > gdata.go
//
// ==========================================================================`

// GenerateGeneticCodeMaps regenerates static protein translation maps
func GenerateGeneticCodeMaps() {

	fmt.Fprintf(os.Stdout, "%s\n\n", gdataHeading)

	fmt.Fprintf(os.Stdout, "package eutils\n\n")

	// generate baseToIdx for source code

	var baseToIdx [256]int

	upperToBase := "-ACMGRSVTWYHKDBN"
	lowerToBase := "-acmgrsvtwyhkdbn"

	// illegal characters map to 0
	for i := 0; i < 256; i++ {
		baseToIdx[i] = baseGap
	}

	// map iupacna alphabet to int
	for i := baseGap; i <= baseN; i++ {
		ch := upperToBase[i]
		baseToIdx[int(ch)] = i
		ch = lowerToBase[i]
		baseToIdx[int(ch)] = i
	}
	baseToIdx['U'] = baseT
	baseToIdx['u'] = baseT
	baseToIdx['X'] = baseN
	baseToIdx['x'] = baseN

	// also map ncbi4na alphabet to int
	for i := baseGap; i <= baseN; i++ {
		baseToIdx[i] = i
	}

	spaces := ""
	fmt.Fprintf(os.Stdout, "var baseToIdx = map[int]int{\n")
	for i := 0; i < 256; i++ {
		if i < 10 {
			spaces = "  "
		} else if i < 100 {
			spaces = " "
		} else {
			spaces = ""
		}
		if baseToIdx[i] != 0 {
			fmt.Fprintf(os.Stdout, "\t%d:%s %d,\n", i, spaces, baseToIdx[i])
		}
	}
	fmt.Fprintf(os.Stdout, "}\n\n")

	// populate nextState and rvCpState arrays

	var nextState [4096]int
	var rvCpState [4096]int

	baseToComp := "-TGKCYSBAWRDMHVN"

	// states 0 through 4095 are triple letter states (---, --A, ..., NNT, NNN)
	for i, st := baseGap, 0; i <= baseN; i++ {
		for j, nx := baseGap, 0; j <= baseN; j++ {
			for k := baseGap; k <= baseN; k++ {
				nextState[st] = nx
				p := baseToIdx[int(baseToComp[k])]
				q := baseToIdx[int(baseToComp[j])]
				r := baseToIdx[int(baseToComp[i])]
				rvCpState[st] = int(256*p + 16*q + r)
				st++
				nx += 16
			}
		}
	}

	// NextCodonState indexes through nextState array and adds base* value
	fmt.Fprintf(os.Stdout, "var nextState = [4096]int{\n")
	for i := 0; i < 4096; i += 16 {
		prefix := "\t"
		for j := 0; j < 16; j++ {
			fmt.Fprintf(os.Stdout, "%s%4d,", prefix, nextState[i+j])
			prefix = " "
		}
		fmt.Fprintf(os.Stdout, "\n")
	}
	fmt.Fprintf(os.Stdout, "}\n\n")

	// RevCompState is a direct index to the reverse complement of a state
	fmt.Fprintf(os.Stdout, "var rvCpState = [4096]int{\n")
	for i := 0; i < 4096; i += 16 {
		prefix := "\t"
		for j := 0; j < 16; j++ {
			fmt.Fprintf(os.Stdout, "%s%4d,", prefix, rvCpState[i+j])
			prefix = " "
		}
		fmt.Fprintf(os.Stdout, "\n")
	}
	fmt.Fprintf(os.Stdout, "}\n\n")

	// translation table structure
	type TransTable struct {
		aminoAcid [4096]int
		orfStart  [4096]int
		orfStop   [4096]int
	}

	// create translation table for given genetic code
	newTransTable := func(genCode int) *TransTable {

		tbl := new(TransTable)
		if tbl == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to allocate translation table for genetic code %d\n", genCode)
			os.Exit(1)
		}

		// genetic code number corrections
		if genCode == 7 {
			genCode = 4
		} else if genCode == 8 {
			genCode = 1
		} else if genCode == 0 {
			genCode = 1
		}

		// return if unable to find ncbieaa and sncbieaa strings
		ncbieaa, ok := ncbieaaCode[genCode]
		if !ok {
			fmt.Fprintf(os.Stderr, "\nERROR: Genetic code %d does not exist\n", genCode)
			os.Exit(1)
		}
		sncbieaa, ok := sncbieaaCode[genCode]
		if !ok {
			fmt.Fprintf(os.Stderr, "\nERROR: Genetic code %d does not exist\n", genCode)
			os.Exit(1)
		}

		// also check length of ncbieaa and sncbieaa strings
		if len(ncbieaa) != 64 || len(sncbieaa) != 64 {
			fmt.Fprintf(os.Stderr, "\nERROR: Genetic code %d length mismatch\n", genCode)
			os.Exit(1)
		}

		// ambiguous codons map to unknown amino acid or not start
		for i := 0; i < 4096; i++ {
			tbl.aminoAcid[i] = int('X')
			tbl.orfStart[i] = int('-')
			tbl.orfStop[i] = int('-')
		}

		var expansions = [4]int{baseA, baseC, baseG, baseT}
		// T = 0, C = 1, A = 2, G = 3
		var codonIdx = [9]int{0, 2, 1, 0, 3, 0, 0, 0, 0}

		// lookup amino acid for each codon in genetic code table
		for i, st := baseGap, 0; i <= baseN; i++ {
			for j := baseGap; j <= baseN; j++ {
				for k := baseGap; k <= baseN; k++ {

					aa := 0
					orf := 0
					goOn := true

					// expand ambiguous IJK nucleotide symbols into component bases XYZ
					for p := 0; p < 4 && goOn; p++ {
						x := expansions[p]
						if (x & i) != 0 {
							for q := 0; q < 4 && goOn; q++ {
								y := expansions[q]
								if (y & j) != 0 {
									for r := 0; r < 4 && goOn; r++ {
										z := expansions[r]
										if (z & k) != 0 {

											// calculate offset in genetic code string

											// the T = 0, C = 1, A = 2, G = 3 order is
											// necessary because the genetic code
											// strings are presented in TCAG order
											cd := 16*codonIdx[x] + 4*codonIdx[y] + codonIdx[z]

											// lookup amino acid for codon XYZ
											ch := int(ncbieaa[cd])
											if aa == 0 {
												aa = ch
											} else if aa != ch {
												// allow Asx (Asp or Asn) and Glx (Glu or Gln)
												if (aa == 'B' || aa == 'D' || aa == 'N') && (ch == 'D' || ch == 'N') {
													aa = 'B'
												} else if (aa == 'Z' || aa == 'E' || aa == 'Q') && (ch == 'E' || ch == 'Q') {
													aa = 'Z'
												} else if (aa == 'J' || aa == 'I' || aa == 'L') && (ch == 'I' || ch == 'L') {
													aa = 'J'
												} else {
													aa = 'X'
												}
											}

											// lookup translation start flag
											ch = int(sncbieaa[cd])
											if orf == 0 {
												orf = ch
											} else if orf != ch {
												orf = 'X'
											}

											// drop out of loop as soon as answer is known
											if aa == 'X' && orf == 'X' {
												goOn = false
											}
										}
									}
								}
							}
						}
					}

					// assign amino acid
					if aa != 0 {
						tbl.aminoAcid[st] = aa
					}
					// assign orf start/stop
					if orf == '*' {
						tbl.orfStop[st] = orf
					} else if orf != 0 && orf != '-' && orf != 'X' {
						tbl.orfStart[st] = orf
					}

					st++
				}
			}
		}

		return tbl
	}

	// generate aminoAcidMaps, orfStartMaps, and orfStopMaps for source code

	// sort genetic code keys in numerical order
	var keys []int
	for ky := range sncbieaaCode {
		keys = append(keys, ky)
	}
	sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })

	tbl := newTransTable(1)
	if tbl == nil {
		fmt.Fprintf(os.Stderr, "newTransTable failed\n")
		return
	}

	fmt.Fprintf(os.Stdout, "// standard genetic code [1] shows only non-X residues,\n")
	fmt.Fprintf(os.Stdout, "// remaining tables show differences from the standard\n\n")
	fmt.Fprintf(os.Stdout, "var aminoAcidMaps = map[int]map[int]int{\n")

	spaces = ""
	fmt.Fprintf(os.Stdout, "\t1: {\n")
	for i := 0; i < 4096; i++ {
		aa := tbl.aminoAcid[i]
		if aa != 0 && aa != 'X' {
			if i < 10 {
				spaces = "   "
			} else if i < 100 {
				spaces = "  "
			} else if i < 1000 {
				spaces = " "
			} else {
				spaces = ""
			}
			fmt.Fprintf(os.Stdout, "\t\t%d:%s '%c', // %s\n", i, spaces, aa, GetCodonFromState(i))
		}
	}
	fmt.Fprintf(os.Stdout, "\t},\n")

	for _, id := range keys {

		// skip standard code
		if id == 1 {
			continue
		}

		tb := newTransTable(id)
		if tb == nil {
			fmt.Fprintf(os.Stderr, "newTransTable failed for code %d, skipping\n", id)
			continue
		}

		spaces = ""
		fmt.Fprintf(os.Stdout, "\t%d: {\n", id)
		for i := 0; i < 4096; i++ {
			aa := tbl.aminoAcid[i]
			nw := tb.aminoAcid[i]
			if aa == nw {
				continue
			}
			// allow 'X' to override specific amino acid in standard code
			if aa != 0 {
				if i < 10 {
					spaces = "   "
				} else if i < 100 {
					spaces = "  "
				} else if i < 1000 {
					spaces = " "
				} else {
					spaces = ""
				}
				fmt.Fprintf(os.Stdout, "\t\t%d:%s '%c', // %s\n", i, spaces, nw, GetCodonFromState(i))
			}
		}
		fmt.Fprintf(os.Stdout, "\t},\n")
	}

	fmt.Fprintf(os.Stdout, "}\n\n")

	// continue with start and stop tables, but all are instantiated

	fmt.Fprintf(os.Stdout, "// start and stop tables show all relevant entries for each genetic code\n\n")
	fmt.Fprintf(os.Stdout, "var orfStartMaps = map[int]map[int]int{\n")

	for _, id := range keys {

		tb := newTransTable(id)
		if tb == nil {
			fmt.Fprintf(os.Stderr, "newTransTable failed for code %d, skipping\n", id)
			continue
		}

		spaces = ""
		fmt.Fprintf(os.Stdout, "\t%d: {\n", id)
		for i := 0; i < 4096; i++ {
			aa := tb.orfStart[i]
			if aa != 0 && aa != '*' && aa != '-' && aa != 'X' {
				if i < 10 {
					spaces = "   "
				} else if i < 100 {
					spaces = "  "
				} else if i < 1000 {
					spaces = " "
				} else {
					spaces = ""
				}
				fmt.Fprintf(os.Stdout, "\t\t%d:%s '%c', // %s\n", i, spaces, aa, GetCodonFromState(i))
			}
		}
		fmt.Fprintf(os.Stdout, "\t},\n")
	}

	fmt.Fprintf(os.Stdout, "}\n\n")

	fmt.Fprintf(os.Stdout, "var orfStopMaps = map[int]map[int]int{\n")

	for _, id := range keys {

		tb := newTransTable(id)
		if tb == nil {
			fmt.Fprintf(os.Stderr, "newTransTable failed for code %d, skipping\n", id)
			continue
		}

		spaces = ""
		fmt.Fprintf(os.Stdout, "\t%d: {\n", id)
		for i := 0; i < 4096; i++ {
			aa := tb.orfStop[i]
			if aa == '*' {
				if i < 10 {
					spaces = "   "
				} else if i < 100 {
					spaces = "  "
				} else if i < 1000 {
					spaces = " "
				} else {
					spaces = ""
				}
				fmt.Fprintf(os.Stdout, "\t\t%d:%s '%c', // %s\n", i, spaces, aa, GetCodonFromState(i))
			}
		}
		fmt.Fprintf(os.Stdout, "\t},\n")
	}

	fmt.Fprintf(os.Stdout, "}\n\n")

	// generate static string of standard genetic code array for performance

	fmt.Fprintf(os.Stdout, "// standard genetic code array for performance\n\n")
	fmt.Fprintf(os.Stdout, "const stdGenCode = \"\" +\n")

	for i := 0; i < 4096; i += 64 {
		fmt.Fprintf(os.Stdout, "\t\"")
		for j := 0; j < 64; j++ {
			aa := tbl.aminoAcid[i+j]
			fmt.Fprintf(os.Stdout, "%c", aa)
		}
		fmt.Fprintf(os.Stdout, "\"")
		if i < 4032 {
			fmt.Fprintf(os.Stdout, " +")
		}
		fmt.Fprintf(os.Stdout, "\n")
	}

	/*
		// print codon to state table for reference

		for _, ch1 := range upperToBase {
			for  _, ch2 := range upperToBase {
				for  _, ch3 := range upperToBase {
					st := SetCodonState(int(ch1), int(ch2), int(ch3))
					fmt.Fprintf(os.Stderr, "%c%c%c\t%d\n", ch1, ch2, ch3, st)
				}
			}
		}
	*/
}
