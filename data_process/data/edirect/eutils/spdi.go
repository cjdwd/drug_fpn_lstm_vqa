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
// File Name:  spdi.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"fmt"
	"os"
	"strconv"
	"strings"
)

// SequenceRemove removes characters from the beginning, the end, or both
func SequenceRemove(seq, first, last string) string {

	if seq == "" {
		return ""
	}

	ln := len(seq)

	if IsAllDigits(first) {
		val, err := strconv.Atoi(first)
		if err == nil && val > 0 {
			if val <= ln {
				seq = seq[val:]
				ln = len(seq)
			} else {
				fmt.Fprintf(os.Stderr, "%s ERROR: %s -first argument %d is greater than sequence length %d%s\n", INVT, LOUD, val, ln, INIT)
				seq = ""
			}
		}
	} else {
		val := len(first)
		if val > 0 {
			if val <= ln {
				// warn if existing sequence does not match deletion argument
				ext := seq[:val]
				if strings.ToUpper(first) != strings.ToUpper(ext) {
					fmt.Fprintf(os.Stderr, "%s WARNING: %s -first argument %s does not match existing sequence %s%s\n", INVT, LOUD, strings.ToUpper(first), strings.ToUpper(ext), INIT)
				}
				// delete characters
				seq = seq[val:]
				ln = len(seq)
			} else {
				fmt.Fprintf(os.Stderr, "%s ERROR: %s -first argument %d is greater than sequence length %d%s\n", INVT, LOUD, val, ln, INIT)
				seq = ""
			}
		}
	}

	if IsAllDigits(last) {
		val, err := strconv.Atoi(last)
		if err == nil && val > 0 {
			if val <= ln {
				seq = seq[:ln-val]
			} else {
				fmt.Fprintf(os.Stderr, "%s ERROR: %s -last argument %d is greater than remaining sequence length %d%s\n", INVT, LOUD, val, ln, INIT)
				seq = ""
			}
		}
	} else {
		val := len(last)
		if val > 0 {
			if val <= ln {
				// warn if existing sequence does not match deletion argument
				ext := seq[ln-val:]
				if strings.ToUpper(last) != strings.ToUpper(ext) {
					fmt.Fprintf(os.Stderr, "%s WARNING: %s -last argument %s does not match existing sequence %s%s\n", INVT, LOUD, strings.ToUpper(last), strings.ToUpper(ext), INIT)
				}
				// delete characters
				seq = seq[:ln-val]
				ln = len(seq)
			} else {
				fmt.Fprintf(os.Stderr, "%s ERROR: %s -last argument %d is greater than sequence length %d%s\n", INVT, LOUD, val, ln, INIT)
				seq = ""
			}
		}
	}

	return seq
}

// SequenceRetain keeps leading or trailing sequence characters
func SequenceRetain(seq string, lead, trail int) string {

	if seq == "" {
		return ""
	}

	ln := len(seq)

	if lead > 0 && trail > 0 {
		fmt.Fprintf(os.Stderr, "%s ERROR: %s Cannot have both -leading and -trailing arguments%s\n", INVT, LOUD, INIT)
		seq = ""
	} else if lead > 0 {
		if lead <= ln {
			seq = seq[:lead]
		} else {
			fmt.Fprintf(os.Stderr, "%s ERROR: %s -leading argument %d is greater than sequence length %d%s\n", INVT, LOUD, lead, ln, INIT)
			seq = ""
		}
	} else if trail > 0 {
		if trail <= ln {
			seq = seq[ln-trail:]
		} else {
			fmt.Fprintf(os.Stderr, "%s ERROR: %s -trailing argument %d is greater than sequence length %d%s\n", INVT, LOUD, trail, ln, INIT)
			seq = ""
		}
	}

	return seq
}

// SequenceReplace applies SPDI instructions
func SequenceReplace(seq string, pos int, del, ins string) string {

	if seq == "" {
		return ""
	}

	ln := len(seq)

	if del == "" && ins == "" {
		fmt.Fprintf(os.Stderr, "%s ERROR: %s -replace command requires either -delete or -insert%s\n", INVT, LOUD, INIT)
		return ""
	}

	if pos > ln {

		if pos == ln+1 && del == "" && ins != "" {

			// append to end of sequence
			seq = seq[:] + ins

		} else {
			fmt.Fprintf(os.Stderr, "%s ERROR: %s -replace position %d is greater than sequence length %d%s\n", INVT, LOUD, pos, ln, INIT)
			return ""
		}

	} else {

		if IsAllDigits(del) {
			val, err := strconv.Atoi(del)
			if err == nil && val > 0 {
				if val <= ln-pos {
					seq = seq[:pos] + seq[pos+val:]
				} else {
					fmt.Fprintf(os.Stderr, "%s ERROR: %s -replace deletion %d is greater than remaining sequence length %d%s\n", INVT, LOUD, val, ln-pos, INIT)
					return ""
				}
			}
		} else {
			val := len(del)
			if val > 0 {
				if val <= ln-pos {
					// warn if existing sequence does not match deletion argument
					ext := seq[pos : pos+val]
					if strings.ToUpper(del) != strings.ToUpper(ext) {
						fmt.Fprintf(os.Stderr, "%s WARNING: %s -replace deletion %s does not match existing sequence %s%s\n", INVT, LOUD, strings.ToUpper(del), strings.ToUpper(ext), INIT)
					}
					// delete characters
					seq = seq[:pos] + seq[pos+val:]
				} else {
					fmt.Fprintf(os.Stderr, "%s ERROR: %s -replace deletion %d is greater than remaining sequence length %d%s\n", INVT, LOUD, val, ln-pos, INIT)
					return ""
				}
			}
		}

		ln = len(seq)
		if ins != "" {
			if pos <= ln {
				seq = seq[:pos] + ins + seq[pos:]
			} else {
				fmt.Fprintf(os.Stderr, "%s ERROR: %s -replace position %d is greater than remaining sequence length %d%s\n", INVT, LOUD, pos, ln-pos, INIT)
				return ""
			}
		}
	}

	return seq
}

// SequenceExtract returns the sequence under the intervals of a feature location
func SequenceExtract(seq, featLoc string, isOneBased bool) string {

	if seq == "" {
		return ""
	}

	ln := len(seq)

	// split intervals, e.g., "201..224,1550..1920,1986..2085,2317..2404,2466..2629"
	comma := strings.Split(featLoc, ",")

	var buffer strings.Builder

	for _, item := range comma {

		// also allow dash separator, e.g., "201-224,1550-1920"
		item = strings.Replace(item, "-", "..", -1)
		// and colon separator, e.g., "201:224,1550:1920"
		item = strings.Replace(item, ":", "..", -1)

		fr, to := SplitInTwoLeft(item, "..")

		fr = strings.TrimSpace(fr)
		to = strings.TrimSpace(to)

		min, err := strconv.Atoi(fr)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s ERROR: %s Unrecognized number '%s'%s\n", INVT, LOUD, fr, INIT)
			os.Exit(1)
		}
		if min < 1 || min > ln {
			fmt.Fprintf(os.Stderr, "%s ERROR: %s Starting point '%s' out of range%s\n", INVT, LOUD, fr, INIT)
			os.Exit(1)
		}

		max, err := strconv.Atoi(to)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s ERROR: %s Unrecognized number '%s'%s\n", INVT, LOUD, to, INIT)
			os.Exit(1)
		}
		if max < 1 || max > ln {
			fmt.Fprintf(os.Stderr, "%s ERROR: %s Ending point '%s' out of range%s\n", INVT, LOUD, to, INIT)
			os.Exit(1)
		}

		if !isOneBased {
			min++
			max++
		}

		if min < max {
			min--
			sub := seq[min:max]
			buffer.WriteString(sub)
		} else if min > max {
			max--
			sub := seq[max:min]
			sub = ReverseComplement(sub)
			buffer.WriteString(sub)
		} else {
			// need more information to know strand if single point
		}
	}

	return buffer.String()
}

// ReverseComplement returns the reverse complement of a sequence
func ReverseComplement(seq string) string {

	runes := []rune(seq)
	// reverse sequence letters - middle base in odd-length sequence is not touched
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	found := false
	// now complement every base, also handling uracil, leaving case intact
	for i, ch := range runes {
		runes[i], found = revComp[ch]
		if !found {
			runes[i] = 'X'
		}
	}
	seq = string(runes)

	return seq
}

// SequenceReverse reverses a sequence, but does not complement the bases
func SequenceReverse(seq string) string {

	runes := []rune(seq)
	// reverse sequence letters - middle base in odd-length sequence is not touched
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	seq = string(runes)

	return seq
}
