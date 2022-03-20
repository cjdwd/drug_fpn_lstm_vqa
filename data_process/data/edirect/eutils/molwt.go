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
// File Name:  molwt.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"strconv"
	"strings"
)

var numC = map[rune]int{
	'A': 3,
	'B': 4,
	'C': 3,
	'D': 4,
	'E': 5,
	'F': 9,
	'G': 2,
	'H': 6,
	'I': 6,
	'J': 6,
	'K': 6,
	'L': 6,
	'M': 5,
	'N': 4,
	'O': 12,
	'P': 5,
	'Q': 5,
	'R': 6,
	'S': 3,
	'T': 4,
	'U': 3,
	'V': 5,
	'W': 11,
	'X': 0,
	'Y': 9,
	'Z': 5,
}

var numH = map[rune]int{
	'A': 5,
	'B': 5,
	'C': 5,
	'D': 5,
	'E': 7,
	'F': 9,
	'G': 3,
	'H': 7,
	'I': 11,
	'J': 11,
	'K': 12,
	'L': 11,
	'M': 9,
	'N': 6,
	'O': 19,
	'P': 7,
	'Q': 8,
	'R': 12,
	'S': 5,
	'T': 7,
	'U': 5,
	'V': 9,
	'W': 10,
	'X': 0,
	'Y': 9,
	'Z': 7,
}

var numN = map[rune]int{
	'A': 1,
	'B': 1,
	'C': 1,
	'D': 1,
	'E': 1,
	'F': 1,
	'G': 1,
	'H': 3,
	'I': 1,
	'J': 1,
	'K': 2,
	'L': 1,
	'M': 1,
	'N': 2,
	'O': 3,
	'P': 1,
	'Q': 2,
	'R': 4,
	'S': 1,
	'T': 1,
	'U': 1,
	'V': 1,
	'W': 2,
	'X': 0,
	'Y': 1,
	'Z': 1,
}

var numO = map[rune]int{
	'A': 1,
	'B': 3,
	'C': 1,
	'D': 3,
	'E': 3,
	'F': 1,
	'G': 1,
	'H': 1,
	'I': 1,
	'J': 1,
	'K': 1,
	'L': 1,
	'M': 1,
	'N': 2,
	'O': 2,
	'P': 1,
	'Q': 2,
	'R': 1,
	'S': 2,
	'T': 2,
	'U': 1,
	'V': 1,
	'W': 1,
	'X': 0,
	'Y': 2,
	'Z': 3,
}

var numS = map[rune]int{
	'A': 0,
	'B': 0,
	'C': 1,
	'D': 0,
	'E': 0,
	'F': 0,
	'G': 0,
	'H': 0,
	'I': 0,
	'J': 0,
	'K': 0,
	'L': 0,
	'M': 1,
	'N': 0,
	'O': 0,
	'P': 0,
	'Q': 0,
	'R': 0,
	'S': 0,
	'T': 0,
	'U': 0,
	'V': 0,
	'W': 0,
	'X': 0,
	'Y': 0,
	'Z': 0,
}

var numSe = map[rune]int{
	'A': 0,
	'B': 0,
	'C': 0,
	'D': 0,
	'E': 0,
	'F': 0,
	'G': 0,
	'H': 0,
	'I': 0,
	'J': 0,
	'K': 0,
	'L': 0,
	'M': 0,
	'N': 0,
	'O': 0,
	'P': 0,
	'Q': 0,
	'R': 0,
	'S': 0,
	'T': 0,
	'U': 1,
	'V': 0,
	'W': 0,
	'X': 0,
	'Y': 0,
	'Z': 0,
}

// ProteinWeight calculates the molecular weight of a peptide sequence
func ProteinWeight(str string, trimLeadingMet bool) string {

	// Start with water (H2O)
	c := 0
	h := 2
	n := 0
	o := 1
	s := 0
	se := 0

	str = strings.ToUpper(str)

	if trimLeadingMet {
		// leading methionine usually removed by post-translational modification
		str = strings.TrimPrefix(str, "M")
	}

	// add number of carbon, hydrogen, nitrogen, oxygen, sulfur, and selenium atoms per amino acid
	for _, ch := range str {
		c += numC[ch]
		h += numH[ch]
		n += numN[ch]
		o += numO[ch]
		s += numS[ch]
		se += numSe[ch]
	}

	// calculate molecular weight
	wt := 12.01115*float64(c) +
		1.0079*float64(h) +
		14.0067*float64(n) +
		15.9994*float64(o) +
		32.064*float64(s) +
		78.96*float64(se)

	str = strconv.Itoa(int(wt + 0.5))

	return str
}
