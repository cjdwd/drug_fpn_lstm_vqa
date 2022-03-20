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
// File Name:  search.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"container/list"
	"fmt"
	"os"
	"strings"
	"unicode"
)

// Searching for a single substring uses the Boyer-Moore-Horspool algorithm,
// as described in Niklaus Wirth, Algorithms and Data Structures, Prentice-Hall,
// 1986, p. 69.
// The original published code had a typographical error, where:
//   UNTIL (j < 0) OR (p[j] # s[i])
// should have been:
//   UNTIL (j < 0) OR (p[j] # s[k])
// This was corrected in later editions.

// Simultaneous searching for multiple substrings uses a Finite State Machine,
// as described in Andrew Binstock and John Rex, Practical Algorithms for
// Programmers, Addison-Wesley, 1995, p. 111.

// Boyer-Moore-Horspool structure
type bmhData struct {
	skipTable [256]int
	pattern   string
	alias     string
	patlen    int
	sensitive bool
}

// Finite State Machine structures
type transitEntry struct {
	char rune
	next int
}

type matchesEntry struct {
	match string
	alias string
}

type stateEntry struct {
	transit []transitEntry
	failure int
	matches []matchesEntry
}

type fsmData struct {
	stateArray []stateEntry
	maxpatlen  int
}

// Searcher contains private fields for Boyer-Moore-Horspool or Finite State Machine searching
type Searcher struct {
	relaxed  bool
	circular bool
	useBMH   bool
	bmhData
	fsmData
}

const failState = -1

// RelaxString removes all punctuation from patterns and search text.
func RelaxString(str string) string {

	if str == "" {
		return str
	}

	terms := strings.FieldsFunc(str, func(c rune) bool {
		return (!unicode.IsLetter(c) && !unicode.IsDigit(c)) || c > 127
	})

	str = strings.Join(terms, " ")
	str = CompressRunsOfSpaces(str)
	str = strings.TrimSpace(str)

	return str
}

func initializeBMH(pattern string, caseSensitive, relaxed, isSequence bool) bmhData {

	alias := ""

	if isSequence {
		pattern, alias = SplitInTwoLeft(pattern, ":")
	}

	if relaxed {
		pattern = RelaxString(pattern)
	}

	if alias == "" {
		alias = pattern
		if isSequence {
			alias = strings.ToUpper(alias)
		}
	}

	patlen := len(pattern)

	// position of last character in pattern
	last := patlen - 1

	// displacement table
	var skip [256]int

	// initialize bad character displacement table
	for i := range skip {
		skip[i] = patlen
	}
	for i := 0; i < last; i++ {
		ch := pattern[i]
		if caseSensitive {
			skip[ch] = last - i
		} else {
			// same displacement for upper and lower case character
			lwr := unicode.ToLower(rune(ch))
			upr := unicode.ToUpper(rune(ch))
			skip[lwr] = last - i
			skip[upr] = last - i
		}
	}

	return bmhData{
		skipTable: skip,
		pattern:   pattern,
		alias:     alias,
		patlen:    patlen,
		sensitive: caseSensitive,
	}
}

func initializeFSM(patterns []string, caseSensitive, relaxed, isSequence, bothStrands bool) fsmData {

	highState := 0

	stateTable := make(map[int]map[rune]int)
	stateFails := make(map[int]int)
	stateMatch := make(map[int][]string)
	stateAlias := make(map[string]string)

	maxpatlen := 0

	gotoState := func(state int, ch rune, failZero bool) int {

		st, ok := stateTable[state]
		if ok {
			for lt, nx := range st {
				if ch == lt {
					return nx
				}
			}
		}
		if state == 0 && failZero {
			return 0
		}
		return failState
	}

	addTransition := func(state int, ch rune, next int) {

		st := stateTable[state]
		if st == nil {
			st = make(map[rune]int)
			stateTable[state] = st
		}
		if caseSensitive {
			st[ch] = next
		} else {
			// same transition for upper and lower case character
			lwr := unicode.ToLower(rune(ch))
			upr := unicode.ToUpper(rune(ch))
			st[lwr] = next
			st[upr] = next
		}
	}

	addOutput := func(state int, word string) {

		mt, ok := stateMatch[state]
		if !ok {
			mt = make([]string, 0, 1)
		}
		mt = append(mt, word)
		stateMatch[state] = mt
	}

	enterWord := func(word string) {

		state := 0

		// try to overlay beginning of word onto existing table
		i := 0
		for _, ch := range word {
			nx := gotoState(state, ch, false)
			if nx == failState {
				break
			}
			i++
			state = nx
		}

		// create new states for remaining characters in word
		rem := word[i:]
		for _, ch := range rem {
			highState++
			addTransition(state, ch, highState)
			state = highState
		}

		// at end of word record match information
		addOutput(state, word)
	}

	findFail := func(state, newState int, ch rune) {

		nx := 0

		// traverse existing failure path
		for {
			nx = gotoState(state, ch, true)
			if nx != failState {
				break
			}
			state = stateFails[state]
		}

		// add new failure state
		stateFails[newState] = nx

		// add matches of substring at new state
		mt, ok := stateMatch[nx]
		if ok {
			for _, wrd := range mt {
				addOutput(newState, wrd)
			}
		}
	}

	computeFail := func() {

		// first-in first-out queue for computing failure states
		queue := list.New()

		// queue up states reached directly from state 0 (depth 1)
		stateFails[0] = 0

		st, ok := stateTable[0]
		if ok {
			for _, s := range st {
				stateFails[s] = 0
				queue.PushBack(s)
			}
		}

		for queue.Len() > 0 {

			e := queue.Front()
			// type assertion
			r := e.Value.(int)
			queue.Remove(e)

			// depth 1 states beget depth 2 states, etc.
			st, ok := stateTable[r]
			if ok {
				for ch, s := range st {
					queue.PushBack(s)
					state := stateFails[r]
					findFail(state, s, ch)
				}
			}
		}
	}

	alreadySeen := make(map[string]bool)

	// initialize state table with search patterns
	for _, txt := range patterns {

		alias := ""

		if isSequence {
			txt, alias = SplitInTwoLeft(txt, ":")
		}

		if relaxed {
			txt = RelaxString(txt)
		}

		// skip if duplicate pattern
		if alreadySeen[txt] {
			if isSequence {
				saw := stateAlias[txt]
				// but if sequence alias previously generated for reverse complement of a pattern
				if saw != "" && strings.HasPrefix(saw, "(") && strings.HasSuffix(saw, ")") {
					if alias == "" {
						alias = strings.ToUpper(txt)
					}
					// print top strand match under current entry
					stateAlias[txt] = alias
				}
			}
			continue
		}
		alreadySeen[txt] = true

		if alias == "" {
			alias = txt
			if isSequence {
				alias = strings.ToUpper(alias)
			}
		}
		stateAlias[txt] = alias

		enterWord(txt)

		if bothStrands {
			// also enter non-palindromic reverse complemenet
			rev := ReverseComplement(txt)
			if strings.ToUpper(rev) != strings.ToUpper(txt) {
				alreadySeen[rev] = true
				// search for lacZ ATGACCATGATTACGGATT on E. coli minus strand
				// shows top strand match (AATCCGTAATCATGGTCAT) in parentheses
				stateAlias[rev] = "(" + alias + ")"
				enterWord(rev)
			}
		}

		// track longest pattern, duplicate that much past end if circular molecule
		mx := len(txt)
		if maxpatlen < mx {
			maxpatlen = mx
		}
	}

	// compute states to go to on failure to match a character
	computeFail()

	/*
		for i := 0; i <= highState; i++ {
			fmt.Fprintf(os.Stderr, "State %d\n", i)
			st, ok := stateTable[i]
			if ok {
				for ch, nx := range st {
					fmt.Fprintf(os.Stderr, "  Character %c, Next State %d\n", ch, nx)
				}
			}
			fl, ok := stateFails[i]
			if ok {
				fmt.Fprintf(os.Stderr, "  Fail State %d\n", fl)
			}
			mt, ok := stateMatch[i]
			if ok {
				for _, wrd := range mt {
					fmt.Fprintf(os.Stderr, "  Match %s\n", wrd)
				}
			}
		}
	*/

	// replace maps with arrays of structures for faster performance

	// create master array
	stateArray := make([]stateEntry, highState+1)

	// populate master array from data in maps
	for state := 0; state <= highState; state++ {

		tsit := make([]transitEntry, 0, 1)
		st, ok := stateTable[state]
		if ok {
			for lt, nx := range st {
				te := transitEntry{char: lt, next: nx}
				tsit = append(tsit, te)
			}
		}
		stateArray[state].transit = tsit

		stateArray[state].failure = stateFails[state]

		mtch := make([]matchesEntry, 0, 1)
		for _, mt := range stateMatch[state] {
			al := stateAlias[mt]
			if al == "" {
				al = mt
			}
			me := matchesEntry{match: mt, alias: al}
			mtch = append(mtch, me)
		}
		stateArray[state].matches = mtch
	}

	return fsmData{
		stateArray: stateArray,
		maxpatlen:  maxpatlen,
	}
}

// newSearcher primes tables for searching on one or more strings.
func newSearcher(patterns []string, caseSensitive, relaxed, isSequence, circular, bothStrands bool) *Searcher {

	if patterns == nil {
		return nil
	}

	if isSequence {
		caseSensitive = false
		relaxed = false
	} else {
		circular = false
		bothStrands = false
	}

	// use Boyer-Moore-Horspool for single pattern, Finite State Machine for multiple patterns
	useBMH := true
	if len(patterns) > 1 {
		useBMH = false
	} else if bothStrands {
		// isolate single sequence pattern from optional alias
		pat, _ := SplitInTwoLeft(patterns[0], ":")
		rev := ReverseComplement(pat)
		// sequence search is not case-sensitive
		if strings.ToUpper(rev) != strings.ToUpper(pat) {
			// if non-palindromic sequence, FSM will also scan for reverse complement pattern
			useBMH = false
		}
	}

	if useBMH {

		bmhObj := initializeBMH(patterns[0], caseSensitive, relaxed, isSequence)

		return &Searcher{
			relaxed:  relaxed,
			circular: circular,
			useBMH:   useBMH,
			bmhData:  bmhObj,
		}
	}

	fsmObj := initializeFSM(patterns, caseSensitive, relaxed, isSequence, bothStrands)

	return &Searcher{
		relaxed:  relaxed,
		circular: circular,
		useBMH:   useBMH,
		fsmData:  fsmObj,
	}
}

// StringSearcher primes tables for searching on one or more strings.
func StringSearcher(patterns []string, caseSensitive, removePunctuation bool) *Searcher {

	return newSearcher(patterns, caseSensitive, removePunctuation, false, false, false)
}

// SequenceSearcher primes tables for searching on one or more nucleotide or protein sequences.
func SequenceSearcher(patterns []string, isProtein, isCircular, topStrandOnly bool) *Searcher {

	if patterns == nil {
		return nil
	}

	if isProtein {
		topStrandOnly = true
	}

	// expands nucleotide ambiguity characters into individual instantiated patterns
	expandPattern := func(pat string) ([]string, bool) {

		if pat == "" {
			return nil, false
		}

		var expanded []string

		overflowed := false

		// recursive function definition
		var expandNext func(prev, next string, level int)

		expandNext = func(prev, next string, level int) {

			// limits on recursion depth and number of expanded patterns
			if overflowed {
				return
			} else if level > 256 {
				overflowed = true
			} else if len(expanded) > 256 {
				overflowed = true
			} else if next != "" {
				// take next character
				curr := next[:1]
				next = next[1:]
				// get set of bases if ambiguous
				exp, ok := expandNuc[curr]
				if !ok {
					exp = curr
				}
				// recursively expand for each base at this position
				for _, ch := range exp {
					expandNext(prev+string(ch), next, level+1)
				}
			} else {
				// at end of string, record one unambiguous pattern
				expanded = append(expanded, prev)
			}
		}

		expandNext("", pat, 1)

		return expanded, overflowed
	}

	var arry []string

	for _, pat := range patterns {

		// each pattern can optionally be followed by a colon and an alias
		txt, alias := SplitInTwoLeft(pat, ":")

		if isProtein {
			// bypass nucleotide character expansion

			if alias != "" {
				txt += ":" + alias
			} else {
				txt += ":" + txt
			}
			arry = append(arry, txt)

			// go to next pattern
			continue
		}

		// expand nucleotide ambiguity characters in pattern to instantiate all matching sequences
		expanded, overflowed := expandPattern(txt)

		if overflowed {
			fmt.Fprintf(os.Stderr, "ERROR: Ignoring pattern expansion of '%s' due to overflow\n", pat)
			continue
		}

		for _, str := range expanded {
			if alias == "+" {
				// RCCGGY:+ prints RCCGGY
				str += ":" + txt
			} else if alias == "-" {
				// RCCGGY:- prints RCCGGY-ACCGGC, etc.
				str += ":" + txt + "-" + str
			} else if alias != "" {
				// RCCGGY:BsrFI-v2 prints BsrFI-v2
				str += ":" + alias
			} else {
				// RCCGGY prints ACCGGC, ACCGGT, GCCGGC, or GCCGGT
				str += ":" + str
			}
			arry = append(arry, str)
		}
	}

	return newSearcher(arry, false, false, true, isCircular, !topStrandOnly)
}

// SearchHit contains position and match or alias
type SearchHit struct {
	Point int
	Match string
}

// Search uses precomputed Searcher tables to search a string or sequence.
func (srch *Searcher) Search(text string) []SearchHit {

	if text == "" {
		return nil
	}

	// failure to create searcher, due to overflow when expanding pattern, results in nil method variable
	if srch == nil {
		return nil
	}

	var arry []SearchHit

	if srch.relaxed {
		text = RelaxString(text)
	}

	if srch.useBMH {
		// single pattern, use Boyer-Moore-Horspool

		cutoff := len(text)

		if cutoff < srch.patlen {
			fmt.Fprintf(os.Stderr, "ERROR: Search text is shorter than pattern\n")
			return nil
		}

		if srch.circular {
			// for circular DNA molecule, copy initial characters and add them to the end of the text
			overhang := text[:srch.patlen]
			text += overhang
		}

		txtlen := len(text)
		max := txtlen - srch.patlen
		last := srch.patlen - 1

		i := 0
		for i <= max {

			j := last
			k := i + last
			if srch.sensitive {
				for j >= 0 && text[k] == srch.pattern[j] {
					j--
					k--
				}
			} else {
				for j >= 0 && unicode.ToLower(rune(text[k])) == unicode.ToLower(rune(srch.pattern[j])) {
					j--
					k--
				}
			}
			if j < 0 {
				// send result if not past end of original text
				if i < cutoff {
					hit := SearchHit{Point: i, Match: srch.alias}
					arry = append(arry, hit)
				}
			}
			// find character in text above last character in pattern
			ch := text[i+last]
			// displacement table can shift pattern by one or more positions
			i += srch.skipTable[ch]
		}

	} else {
		// multiple patterns, use Finite State Machine
		// (can coerce use of FSM on single pattern by entering the same pattern twice)

		cutoff := len(text)

		if cutoff < srch.maxpatlen {
			fmt.Fprintf(os.Stderr, "ERROR: Search text is shorter than pattern\n")
			return nil
		}

		if srch.circular {
			// for circular DNA molecule, copy initial characters and add them to the end of the text
			overhang := text[:srch.maxpatlen]
			text += overhang
		}

		// not the same as version used in creation - no failZero argument, and works on final arrays
		gotoState := func(state int, ch rune) int {

			tbl := srch.stateArray[state]
			if tbl.transit != nil {
				for _, te := range tbl.transit {
					if ch == te.char {
						return te.next
					}
				}
			}
			if state == 0 {
				return 0
			}
			return failState
		}

		state := 0
		for pos, ch := range text {

			nx := 0

			for {
				nx = gotoState(state, ch)
				if nx != failState {
					break
				}
				tbl := srch.stateArray[state]
				state = tbl.failure
			}

			// fmt.Fprintf(os.Stdout, "%5d  %c  %5d -> %5d\n", pos, ch, state, nx)

			state = nx

			tbl := srch.stateArray[state]
			if tbl.matches != nil {
				for _, me := range tbl.matches {
					// send result if not past end of original text
					point := pos - len(me.match) + 1
					if point < cutoff {
						hit := SearchHit{Point: point, Match: me.alias}
						arry = append(arry, hit)
					}
				}
			}
		}
	}

	return arry
}
