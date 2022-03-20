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
// File Name:  phrase.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"github.com/surgebase/porter2"
	"io"
	"os"
	"path"
	"regexp"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"unicode"
)

// Master points to a term and to its postings data
type Master struct {
	TermOffset int32
	PostOffset int32
}

// Arrays contains postings lists and word offsets
type Arrays struct {
	Data []int32
	Ofst [][]int16
	Dist int
}

func commonOpenFile(dpath, fname string) (*os.File, int64) {

	fpath := path.Join(dpath, fname)
	if fpath == "" {
		return nil, 0
	}

	inFile, err := os.Open(fpath)
	if err != nil && os.IsNotExist(err) {
		return nil, 0
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil, 0
	}

	fi, err := inFile.Stat()
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil, 0
	}

	size := fi.Size()

	return inFile, size
}

func readMasterIndex(dpath, key, field string) []Master {

	inFile, size := commonOpenFile(dpath, key+"."+field+".mst")
	if inFile == nil {
		return nil
	}

	defer inFile.Close()

	data := make([]Master, size/8)
	if data == nil || len(data) < 1 {
		return nil
	}

	err := binary.Read(inFile, binary.LittleEndian, &data)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	return data
}

func readTermList(dpath, key, field string) []byte {

	inFile, size := commonOpenFile(dpath, key+"."+field+".trm")
	if inFile == nil {
		return nil
	}

	defer inFile.Close()

	data := make([]byte, size)
	if data == nil || len(data) < 1 {
		return nil
	}

	err := binary.Read(inFile, binary.LittleEndian, &data)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	return data
}

func readPostingData(dpath, key, field string, offset int32, size int32) []int32 {

	inFile, _ := commonOpenFile(dpath, key+"."+field+".pst")
	if inFile == nil {
		return nil
	}

	defer inFile.Close()

	data := make([]int32, size/4)
	if data == nil || len(data) < 1 {
		return nil
	}

	_, err := inFile.Seek(int64(offset), io.SeekStart)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	err = binary.Read(inFile, binary.LittleEndian, data)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	return data
}

func readPositionIndex(dpath, key, field string, offset int32, size int32) []int32 {

	inFile, _ := commonOpenFile(dpath, key+"."+field+".uqi")
	if inFile == nil {
		return nil
	}

	defer inFile.Close()

	data := make([]int32, size/4)
	if data == nil || len(data) < 1 {
		return nil
	}

	_, err := inFile.Seek(int64(offset), io.SeekStart)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	err = binary.Read(inFile, binary.LittleEndian, data)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	return data
}

func readOffsetData(dpath, key, field string, offset int32, size int32) []int16 {

	inFile, _ := commonOpenFile(dpath, key+"."+field+".ofs")
	if inFile == nil {
		return nil
	}

	defer inFile.Close()

	data := make([]int16, size/2)
	if data == nil || len(data) < 1 {
		return nil
	}

	_, err := inFile.Seek(int64(offset), io.SeekStart)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	err = binary.Read(inFile, binary.LittleEndian, data)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	return data
}

func readMasterIndexFuture(dpath, key, field string) <-chan []Master {

	out := make(chan []Master, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create master index channel\n")
		os.Exit(1)
	}

	// masterIndexFuture asynchronously gets the master file and sends results through channel
	masterIndexFuture := func(dpath, key, field string, out chan<- []Master) {

		data := readMasterIndex(dpath, key, field)

		out <- data

		close(out)
	}

	// launch single future goroutine
	go masterIndexFuture(dpath, key, field, out)

	return out
}

func readTermListFuture(dpath, key, field string) <-chan []byte {

	out := make(chan []byte, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create term list channel\n")
		os.Exit(1)
	}

	// termListFuture asynchronously gets posting IDs and sends results through channel
	termListFuture := func(dpath, key, field string, out chan<- []byte) {

		data := readTermList(dpath, key, field)

		out <- data

		close(out)
	}

	// launch single future goroutine
	go termListFuture(dpath, key, field, out)

	return out
}

func getPostingIDs(prom, term, field string, simple bool) ([]int32, [][]int16) {

	var (
		arry [516]rune
	)

	dpath, key := PostingPath(prom, field, term, arry)
	if dpath == "" {
		return nil, nil
	}

	// schedule asynchronous fetching
	mi := readMasterIndexFuture(dpath, key, field)

	tl := readTermListFuture(dpath, key, field)

	// fetch master index and term list
	indx := <-mi

	trms := <-tl

	if indx == nil || len(indx) < 1 {
		return nil, nil
	}

	if trms == nil || len(trms) < 1 {
		return nil, nil
	}

	// master index is padded with phantom term and postings position
	numTerms := len(indx) - 1

	strs := make([]string, numTerms)
	if strs == nil || len(strs) < 1 {
		return nil, nil
	}

	retlength := int32(len("\n"))

	// populate array of strings from term list
	for i, j := 0, 1; i < numTerms; i++ {
		from := indx[i].TermOffset
		to := indx[j].TermOffset - retlength
		j++
		txt := string(trms[from:to])
		strs[i] = txt
	}

	// change protecting underscore to space
	term = strings.Replace(term, "_", " ", -1)

	// if term ends with dollar sign, use porter2 stemming, then add asterisk
	if strings.HasSuffix(term, "$") && term != "$" {
		term = strings.TrimSuffix(term, "$")
		term = porter2.Stem(term)
		term += "*"
	}

	isWildCard := false
	if strings.HasSuffix(term, "*") && term != "*" {
		tlen := len(term)
		isWildCard = true
		term = strings.TrimSuffix(term, "*")
		pdlen := len(PostingDir(term))
		if tlen < pdlen {
			fmt.Fprintf(os.Stderr, "Wildcard term '%s' must be at least %d characters long - ignoring this word\n", term, pdlen)
			return nil, nil
		}
	}

	// binary search in term list
	L, R := 0, numTerms-1
	for L < R {
		mid := (L + R) / 2
		if strs[mid] < term {
			L = mid + 1
		} else {
			R = mid
		}
	}

	// wild card search scans term lists, fuses adjacent postings lists
	if isWildCard {
		if R < numTerms && strings.HasPrefix(strs[R], term) {
			offset := indx[R].PostOffset
			for R < numTerms && strings.HasPrefix(strs[R], term) {
				R++
			}
			size := indx[R].PostOffset - offset

			// read relevant postings list section
			data := readPostingData(dpath, key, field, offset, size)
			if data == nil || len(data) < 1 {
				return nil, nil
			}

			if simple {

				merged := make(map[int32]bool)

				// combine all postings in term range
				for _, val := range data {
					merged[val] = true
				}

				fused := make([]int32, len(merged))

				// convert map to slice
				i := 0
				for num := range merged {
					fused[i] = num
					i++
				}

				sort.Slice(fused, func(i, j int) bool { return fused[i] < fused[j] })

				return fused, nil
			}

			// read relevant word position section, includes phantom offset at end
			uqis := readPositionIndex(dpath, key, field, offset, size+4)
			if uqis == nil {
				return nil, nil
			}
			ulen := len(uqis)
			if ulen < 1 {
				return nil, nil
			}

			from := uqis[0]
			to := uqis[ulen-1]

			// read offset section
			ofst := readOffsetData(dpath, key, field, from, to-from)
			if ofst == nil {
				return nil, nil
			}

			combo := make(map[int32][]int16)

			addPositions := func(uid int32, pos int16) {

				arrs, ok := combo[uid]
				if !ok {
					arrs = make([]int16, 0, 1)
				}
				arrs = append(arrs, pos)
				combo[uid] = arrs
			}

			// populate array of positions per UID
			for i, j, k := 0, 1, int32(0); i < ulen-1; i++ {
				uid := data[i]
				num := (uqis[j] - uqis[i]) / 2
				j++
				for q := k; q < k+num; q++ {
					addPositions(uid, ofst[q])
				}
				k += num
			}

			fused := make([]int32, len(combo))

			// convert map to slice
			i := 0
			for num := range combo {
				fused[i] = num
				i++
			}

			sort.Slice(fused, func(i, j int) bool { return fused[i] < fused[j] })

			// make array of int16 arrays, populate for each UID
			arrs := make([][]int16, ulen-1)
			if arrs == nil {
				return nil, nil
			}

			for j, uid := range fused {
				posn := combo[uid]

				if len(posn) > 1 {
					sort.Slice(posn, func(i, j int) bool { return posn[i] < posn[j] })
				}

				arrs[j] = posn
			}

			return fused, arrs
		}

		return nil, nil
	}

	// regular search requires exact match from binary search
	if R < numTerms && strs[R] == term {

		offset := indx[R].PostOffset
		size := indx[R+1].PostOffset - offset

		// read relevant postings list section
		data := readPostingData(dpath, key, field, offset, size)
		if data == nil || len(data) < 1 {
			return nil, nil
		}

		if simple {
			return data, nil
		}

		// read relevant word position section, includes phantom offset at end
		uqis := readPositionIndex(dpath, key, field, offset, size+4)
		if uqis == nil {
			return nil, nil
		}
		ulen := len(uqis)
		if ulen < 1 {
			return nil, nil
		}

		from := uqis[0]
		to := uqis[ulen-1]

		// read offset section
		ofst := readOffsetData(dpath, key, field, from, to-from)
		if ofst == nil {
			return nil, nil
		}

		// make array of int16 arrays, populate for each UID
		arrs := make([][]int16, ulen)
		if arrs == nil || len(arrs) < 1 {
			return nil, nil
		}

		// populate array of positions per UID
		for i, j, k := 0, 1, int32(0); i < ulen-1; i++ {
			num := (uqis[j] - uqis[i]) / 2
			j++
			arrs[i] = ofst[k : k+num]
			k += num
		}

		return data, arrs
	}

	return nil, nil
}

func printTermCount(base, term, field string) int {

	data, _ := getPostingIDs(base, term, field, true)
	size := len(data)
	fmt.Fprintf(os.Stdout, "%d\t%s\n", size, term)

	return size
}

func printTermCounts(base, term, field string) int {

	pdlen := len(PostingDir(term))

	if len(term) < pdlen {
		fmt.Fprintf(os.Stderr, "\nERROR: Term count argument must be at least %d characters\n", pdlen)
		os.Exit(1)
	}

	if strings.Contains(term[:pdlen], "*") {
		fmt.Fprintf(os.Stderr, "\nERROR: Wildcard asterisk must not be in first %d characters\n", pdlen)
		os.Exit(1)
	}

	var arry [516]rune
	dpath, key := PostingPath(base, field, term, arry)
	if dpath == "" {
		return 0
	}

	// schedule asynchronous fetching
	mi := readMasterIndexFuture(dpath, key, field)

	tl := readTermListFuture(dpath, key, field)

	// fetch master index and term list
	indx := <-mi

	trms := <-tl

	if indx == nil || len(indx) < 1 {
		return 0
	}

	if trms == nil || len(trms) < 1 {
		return 0
	}

	// master index is padded with phantom term and postings position
	numTerms := len(indx) - 1

	strs := make([]string, numTerms)
	if strs == nil || len(strs) < 1 {
		return 0
	}

	retlength := int32(len("\n"))

	// populate array of strings from term list
	for i, j := 0, 1; i < numTerms; i++ {
		from := indx[i].TermOffset
		to := indx[j].TermOffset - retlength
		j++
		txt := string(trms[from:to])
		strs[i] = txt
	}

	// change protecting underscore to space
	term = strings.Replace(term, "_", " ", -1)

	// flank pattern with start-of-string and end-of-string symbols
	pat := "^" + term + "$"

	// change asterisk in query to dot + star for regular expression
	pat = strings.Replace(pat, "*", ".*", -1)

	re, err := regexp.Compile(pat)

	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return 0
	}

	count := 0

	for R, str := range strs {
		if re.MatchString(str) {
			offset := indx[R].PostOffset
			size := indx[R+1].PostOffset - offset
			fmt.Fprintf(os.Stdout, "%d\t%s\n", size/4, str)
			count++
		}
	}

	return count
}

func printTermPositions(base, term, field string) int {

	data, ofst := getPostingIDs(base, term, field, false)
	size := len(data)
	fmt.Fprintf(os.Stdout, "\n%d\t%s\n\n", size, term)

	for i := 0; i < len(data); i++ {
		fmt.Fprintf(os.Stdout, "%12d\t", data[i])
		pos := ofst[i]
		sep := ""
		for j := 0; j < len(pos); j++ {
			fmt.Fprintf(os.Stdout, "%s%d", sep, pos[j])
			sep = ","
		}
		fmt.Fprintf(os.Stdout, "\n")
	}

	return size
}

// BOOLEAN OPERATIONS FOR POSTINGS LISTS

func extendPositionalIDs(N []int32, np [][]int16, M []int32, mp [][]int16, delta int, proc func(pn, pm []int16, dlt int16) []int16) ([]int32, [][]int16) {

	if proc == nil {
		return nil, nil
	}

	if N == nil || len(N) < 1 || np == nil || len(np) < 1 {
		return M, mp
	}
	if M == nil || len(M) < 1 || mp == nil || len(mp) < 1 {
		return N, np
	}

	n, m := len(N), len(M)

	// order matters when extending phrase or testing proximity, do not swap lists based on size

	sz := n
	if sz > m {
		sz = m
	}

	if sz < 1 {
		return N, np
	}

	res := make([]int32, sz)
	ofs := make([][]int16, sz)

	if res == nil || len(res) < 1 || ofs == nil || len(ofs) < 1 {
		return nil, nil
	}

	i, j, k := 0, 0, 0

	// use local variables for speed
	en, em := N[i], M[j]

	for {
		// do inequality tests first
		if en < em {
			i++
			if i == n {
				break
			}
			en = N[i]
		} else if en > em {
			j++
			if j == m {
				break
			}
			em = M[j]
		} else {
			// specific callbacks test position arrays to match terms by adjacency or phrases by proximity
			adj := proc(np[i], mp[j], int16(delta))
			if adj != nil && len(adj) > 0 {
				res[k] = en
				ofs[k] = adj
				k++
			}
			i++
			j++
			if i == n || j == m {
				break
			}
			en = N[i]
			em = M[j]
		}
	}

	// truncate output arrays to actual size of intersection
	res = res[:k]
	ofs = ofs[:k]

	return res, ofs
}

func intersectIDs(N, M []int32) []int32 {

	if N == nil {
		return M
	}
	if M == nil {
		return N
	}

	n, m := len(N), len(M)

	// swap to make M the smaller list
	if n < m {
		N, M = M, N
		n, m = m, n
	}

	if m < 1 {
		return N
	}

	res := make([]int32, m)

	i, j, k := 0, 0, 0

	// use local variables for speed
	en, em := N[i], M[j]

	for {
		// do inequality tests first
		if en < em {
			// index to larger list most likely to be advanced
			i++
			if i == n {
				break
			}
			en = N[i]
		} else if en > em {
			j++
			if j == m {
				break
			}
			em = M[j]
		} else {
			// equality (intersection match) least likely
			res[k] = en
			k++
			i++
			j++
			if i == n || j == m {
				break
			}
			en = N[i]
			em = M[j]
		}
	}

	// truncate output array to actual size of intersection
	res = res[:k]

	return res
}

// if m * log(n) < m + n, binary search has fewer comparisons, but processor memory caches make linear algorithm faster
/*
func intersectBinary(N, M []int32) []int32 {

	if N == nil {
		return M
	}
	if M == nil {
		return N
	}

	n, m := len(N), len(M)

	// swap to make M the smaller list
	if n < m {
		N, M = M, N
		n, m = m, n
	}

	if m < 1 {
		return N
	}

	k := 0

	res := make([]int32, m)

	for _, uid := range M {
		// inline binary search is faster than sort.Search
		L, R := 0, n-1
		for L < R {
			mid := (L + R) / 2
			if N[mid] < uid {
				L = mid + 1
			} else {
				R = mid
			}
		}
		// R := sort.Search(len(N), func(i int) bool { return N[i] >= uid })
		if R < n && N[R] == uid {
			res[k] = uid
			k++
			// remove leading part of N for slight speed gain
			N = N[R:]
			n = len(N)
		}
	}

	res = res[:k]

	return res
}
*/

func combineIDs(N, M []int32) []int32 {

	if N == nil {
		return M
	}
	if M == nil {
		return N
	}

	n, m := len(N), len(M)

	// swap to make M the smaller list
	if n < m {
		N, M = M, N
		n, m = m, n
	}

	if m < 1 {
		return N
	}

	i, j, k := 0, 0, 0

	res := make([]int32, n+m)

	for i < n && j < m {
		if N[i] < M[j] {
			res[k] = N[i]
			k++
			i++
		} else if N[i] > M[j] {
			res[k] = M[j]
			k++
			j++
		} else {
			res[k] = N[i]
			k++
			i++
			j++
		}
	}
	for i < n {
		res[k] = N[i]
		k++
		i++
	}
	for j < m {
		res[k] = M[j]
		k++
		j++
	}

	res = res[:k]

	return res
}

func excludeIDs(N, M []int32) []int32 {

	if N == nil {
		return nil
	}
	if M == nil {
		return N
	}

	n, m := len(N), len(M)

	if m < 1 {
		return N
	}

	res := make([]int32, n)

	i, j, k := 0, 0, 0

	// use local variables for speed
	en, em := N[i], M[j]

	for {
		// do inequality tests first
		if en < em {
			// item is not excluded
			res[k] = en
			k++
			i++
			if i == n {
				break
			}
			en = N[i]
		} else if en > em {
			// advance second list pointer
			j++
			if j == m {
				break
			}
			em = M[j]
		} else {
			// exclude
			i++
			j++
			if i == n || j == m {
				break
			}
			en = N[i]
			em = M[j]
		}
	}

	// truncate output array to actual size of result
	res = res[:k]

	return res
}

// QUERY EVALUATION FUNCTION

func postingIDsFuture(base, term, field string, dist int) <-chan Arrays {

	out := make(chan Arrays, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create postings channel\n")
		os.Exit(1)
	}

	// postingFuture asynchronously gets posting IDs and sends results through channel
	postingFuture := func(base, term, field string, dist int, out chan<- Arrays) {

		data, ofst := getPostingIDs(base, term, field, false)

		out <- Arrays{Data: data, Ofst: ofst, Dist: dist}

		close(out)
	}

	// launch single future goroutine
	go postingFuture(base, term, field, dist, out)

	return out
}

func evaluateQuery(base string, clauses []string) int {

	if clauses == nil || clauses[0] == "" {
		return 0
	}

	count := 0

	// flag set if no tildes, indicates no proximity tests in query
	noProx := true
	for _, tkn := range clauses {
		if strings.HasPrefix(tkn, "~") {
			noProx = false
		}
	}

	phrasePositions := func(pn, pm []int16, dlt int16) []int16 {

		var arry []int16

		ln, lm := len(pn), len(pm)

		q, r := 0, 0

		vn, vm := pn[q], pm[r]
		vnd := vn + dlt

		for {
			if vnd > vm {
				r++
				if r == lm {
					break
				}
				vm = pm[r]
			} else if vnd < vm {
				q++
				if q == ln {
					break
				}
				vn = pn[q]
				vnd = vn + dlt
			} else {
				// store position of first word in current growing phrase
				arry = append(arry, vn)
				q++
				r++
				if q == ln || r == lm {
					break
				}
				vn = pn[q]
				vm = pm[r]
				vnd = vn + dlt
			}
		}

		return arry
	}

	proximityPositions := func(pn, pm []int16, dlt int16) []int16 {

		var arry []int16

		ln, lm := len(pn), len(pm)

		q, r := 0, 0

		vn, vm := pn[q], pm[r]
		vnd := vn + dlt

		for {
			if vnd < vm {
				q++
				if q == ln {
					break
				}
				vn = pn[q]
				vnd = vn + dlt
			} else if vn < vm {
				// store position of first word in downstream phrase that passes proximity test
				arry = append(arry, vm)
				q++
				r++
				if q == ln || r == lm {
					break
				}
				vn = pn[q]
				vm = pm[r]
				vnd = vn + dlt
			} else {
				r++
				if r == lm {
					break
				}
				vm = pm[r]
			}
		}

		return arry
	}

	eval := func(str string) ([]int32, [][]int16, int) {

		// extract optional [FIELD] qualifier
		field := "TIAB"

		if strings.HasSuffix(str, "]") {
			pos := strings.Index(str, "[")
			if pos >= 0 {
				field = str[pos:]
				field = strings.TrimPrefix(field, "[")
				field = strings.TrimSuffix(field, "]")
				str = str[:pos]
				str = strings.TrimSpace(str)
			}
			switch field {
			case "NORM":
				field = "TIAB"
			case "STEM", "TIAB", "TITL":
			case "PIPE":
				// esearch -db pubmed -query "complement system proteins [MESH]" -pub clinical |
				// efetch -format uid | phrase-search -query "[PIPE] AND L [THME]"
				var data []int32
				// read UIDs from stdin
				uidq := CreateUIDReader(os.Stdin)
				for ext := range uidq {

					val, err := strconv.Atoi(ext.Text)
					if err != nil {
						fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized UID %s\n", ext.Text)
						os.Exit(1)
					}

					data = append(data, int32(val))
				}
				// sort UIDs before returning
				sort.Slice(data, func(i, j int) bool { return data[i] < data[j] })
				return data, nil, 0
			default:
				str = strings.Replace(str, " ", "_", -1)
			}
		}

		words := strings.Fields(str)

		if words == nil || len(words) < 1 {
			return nil, nil, 0
		}

		// if no tilde proximity tests, and not building up phrase from multiple words,
		// no need to use more expensive position tests when calculating intersection
		if noProx && len(words) == 1 {
			term := words[0]
			if strings.HasPrefix(term, "+") {
				return nil, nil, 0
			}
			term = strings.Replace(term, "_", " ", -1)
			data, _ := getPostingIDs(base, term, field, true)
			count++
			return data, nil, 1
		}

		dist := 0

		var intersect []Arrays

		var futures []<-chan Arrays

		// schedule asynchronous fetching
		for _, term := range words {

			term = strings.Replace(term, "_", " ", -1)

			if strings.HasPrefix(term, "+") {
				dist += strings.Count(term, "+")
				// run of stop words or explicit plus signs skip past one or more words in phrase
				continue
			}

			fetch := postingIDsFuture(base, term, field, dist)

			futures = append(futures, fetch)

			dist++
		}

		runtime.Gosched()

		for _, chn := range futures {

			// fetch postings data
			fut := <-chn

			if len(fut.Data) < 1 {
				// bail if word not present
				return nil, nil, 0
			}

			// append posting and positions
			intersect = append(intersect, fut)

			runtime.Gosched()
		}

		if len(intersect) < 1 {
			return nil, nil, 0
		}

		// start phrase with first word
		data, ofst, dist := intersect[0].Data, intersect[0].Ofst, intersect[0].Dist+1

		if len(intersect) == 1 {
			return data, ofst, dist
		}

		for i := 1; i < len(intersect); i++ {

			// add subsequent words, keep starting positions of phrases that contain all words in proper position
			data, ofst = extendPositionalIDs(data, ofst, intersect[i].Data, intersect[i].Ofst, intersect[i].Dist, phrasePositions)
			if len(data) < 1 {
				// bail if phrase not present
				return nil, nil, 0
			}
			dist = intersect[i].Dist + 1
		}

		count += len(intersect)

		// return UIDs and all positions of current phrase
		return data, ofst, dist
	}

	prevTkn := ""

	nextToken := func() string {

		if len(clauses) < 1 {
			return ""
		}

		// remove next token from slice
		tkn := clauses[0]
		clauses = clauses[1:]

		if tkn == "(" && prevTkn != "" && prevTkn != "&" && prevTkn != "|" && prevTkn != "!" {
			fmt.Fprintf(os.Stderr, "\nERROR: Tokens '%s' and '%s' should be separated by AND, OR, or NOT\n", prevTkn, tkn)
			os.Exit(1)
		}

		if prevTkn == ")" && tkn != "" && tkn != "&" && tkn != "|" && tkn != "!" && tkn != ")" {
			fmt.Fprintf(os.Stderr, "\nERROR: Tokens '%s' and '%s' should be separated by AND, OR, or NOT\n", prevTkn, tkn)
			os.Exit(1)
		}

		prevTkn = tkn

		return tkn
	}

	// recursive definitions
	var fact func() ([]int32, [][]int16, int, string)
	var prox func() ([]int32, string)
	var excl func() ([]int32, string)
	var term func() ([]int32, string)
	var expr func() ([]int32, string)

	fact = func() ([]int32, [][]int16, int, string) {

		var (
			data  []int32
			ofst  [][]int16
			delta int
			tkn   string
		)

		tkn = nextToken()

		if tkn == "(" {
			// recursively process expression in parentheses
			data, tkn = expr()
			if tkn == ")" {
				tkn = nextToken()
			} else {
				fmt.Fprintf(os.Stderr, "\nERROR: Expected ')' but received '%s'\n", tkn)
				os.Exit(1)
			}
		} else if tkn == ")" {
			fmt.Fprintf(os.Stderr, "\nERROR: Unexpected ')' token\n")
			os.Exit(1)
		} else if tkn == "&" || tkn == "|" || tkn == "!" {
			fmt.Fprintf(os.Stderr, "\nERROR: Unexpected operator '%s' in expression\n", tkn)
			os.Exit(1)
		} else if tkn == "" {
			fmt.Fprintf(os.Stderr, "\nERROR: Unexpected end of expression\n")
			os.Exit(1)
		} else {
			// evaluate current phrase
			data, ofst, delta = eval(tkn)
			tkn = nextToken()
		}

		return data, ofst, delta, tkn
	}

	prox = func() ([]int32, string) {

		var (
			next []int32
			noff [][]int16
			ndlt int
		)

		data, ofst, delta, tkn := fact()
		if len(data) < 1 {
			return nil, ""
		}

		for strings.HasPrefix(tkn, "~") {
			dist := strings.Count(tkn, "~")
			next, noff, ndlt, tkn = fact()
			if len(next) < 1 {
				return nil, ""
			}
			// next phrase must be within specified distance after the previous phrase
			data, ofst = extendPositionalIDs(data, ofst, next, noff, delta+dist, proximityPositions)
			if len(data) < 1 {
				return nil, ""
			}
			delta = ndlt
		}

		return data, tkn
	}

	excl = func() ([]int32, string) {

		var next []int32

		data, tkn := prox()
		for tkn == "!" {
			next, tkn = prox()
			data = excludeIDs(data, next)
		}

		return data, tkn
	}

	term = func() ([]int32, string) {

		var next []int32

		data, tkn := excl()
		for tkn == "&" {
			next, tkn = excl()
			data = intersectIDs(data, next)
		}

		return data, tkn
	}

	expr = func() ([]int32, string) {

		var next []int32

		data, tkn := term()
		for tkn == "|" {
			next, tkn = term()
			data = combineIDs(data, next)
		}

		return data, tkn
	}

	// enter recursive descent parser
	result, tkn := expr()

	if tkn != "" {
		fmt.Fprintf(os.Stderr, "\nERROR: Unexpected token '%s' at end of expression\n", tkn)
		os.Exit(1)
	}

	// sort final result
	sort.Slice(result, func(i, j int) bool { return result[i] < result[j] })

	// use buffers to speed up uid printing
	var buffer strings.Builder

	wrtr := bufio.NewWriter(os.Stdout)

	for _, pmid := range result {
		val := strconv.Itoa(int(pmid))
		buffer.WriteString(val[:])
		buffer.WriteString("\n")
	}

	txt := buffer.String()
	if txt != "" {
		// print buffer
		wrtr.WriteString(txt[:])
	}

	wrtr.Flush()

	runtime.Gosched()

	return count
}

// QUERY PARSING FUNCTIONS

func prepareQuery(str string) string {

	if str == "" {
		return ""
	}

	if strings.HasPrefix(str, "[PIPE]") {
		str = "stdin " + str
	}

	str = CleanupQuery(str, false, true)

	str = strings.Replace(str, "~ ~", "~~", -1)
	str = strings.Replace(str, "~ ~", "~~", -1)

	str = strings.TrimSpace(str)

	// temporarily flank with spaces to detect misplaced operators at ends
	str = " " + str + " "

	str = strings.Replace(str, " AND ", " & ", -1)
	str = strings.Replace(str, " OR ", " | ", -1)
	str = strings.Replace(str, " NOT ", " ! ", -1)

	str = strings.Replace(str, "(", " ( ", -1)
	str = strings.Replace(str, ")", " ) ", -1)
	str = strings.Replace(str, "&", " & ", -1)
	str = strings.Replace(str, "|", " | ", -1)
	str = strings.Replace(str, "!", " ! ", -1)

	// ensure that bracketed fields are flanked by spaces
	str = strings.Replace(str, "[", " [", -1)
	str = strings.Replace(str, "]", "] ", -1)

	// remove temporary flanking spaces
	str = strings.TrimSpace(str)

	str = strings.ToLower(str)

	str = strings.Replace(str, "_", " ", -1)

	hasPlusOrMinus := func(str string) bool {

		for _, ch := range str {
			if ch == '-' || ch == '+' {
				return true
			}
		}

		return false
	}

	fixThemeCases := func(str string) string {

		if !strings.Contains(str, "[thme]") && !strings.Contains(str, "[conv]") {
			return str
		}

		var arry []string

		terms := strings.Fields(str)

		for _, item := range terms {

			switch item {
			case "a+":
				arry = append(arry, "ap")
			case "e+":
				arry = append(arry, "ep")
			case "ec+":
				arry = append(arry, "ecp")
			case "eg+":
				arry = append(arry, "egp")
			case "v+":
				arry = append(arry, "vp")
			case "a-":
				arry = append(arry, "am")
			case "e-":
				arry = append(arry, "em")
			case "ec-":
				arry = append(arry, "ecm")
			default:
				arry = append(arry, item)
			}
		}

		// reconstruct string from transformed words
		str = strings.Join(arry, " ")

		return str
	}

	if hasPlusOrMinus(str) {
		str = fixThemeCases(str)
	}

	if HasHyphenOrApostrophe(str) {
		str = FixSpecialCases(str)
	}

	str = strings.Replace(str, "-", " ", -1)

	// break terms at punctuation, and at non-ASCII characters, allowing brackets for field names,
	// along with Boolean control symbols, underscore for protected terms, asterisk to indicate
	// truncation wildcard, tilde for maximum proximity, and plus sign for exactly one wildcard word
	terms := strings.FieldsFunc(str, func(c rune) bool {
		return (!unicode.IsLetter(c) && !unicode.IsDigit(c) &&
			c != '_' && c != '*' && c != '~' && c != '+' &&
			c != '$' && c != '&' && c != '|' && c != '!' &&
			c != '(' && c != ')' && c != '[' && c != ']') || c > 127
	})

	// rejoin into processed sentence
	tmp := strings.Join(terms, " ")

	tmp = CompressRunsOfSpaces(tmp)
	tmp = strings.TrimSpace(tmp)

	return tmp
}

func prepareExact(str, sfx string, deStop bool) string {

	if str == "" {
		return ""
	}

	if str == "[Not Available]." || str == "Health." {
		return ""
	}

	str = CleanupQuery(str, true, true)

	str = strings.Replace(str, "(", " ", -1)
	str = strings.Replace(str, ")", " ", -1)

	str = strings.Replace(str, "_", " ", -1)

	if HasHyphenOrApostrophe(str) {
		str = FixSpecialCases(str)
	}

	str = strings.Replace(str, "-", " ", -1)

	// remove trailing punctuation from each word
	var arry []string

	terms := strings.Fields(str)
	for _, item := range terms {
		max := len(item)
		for max > 1 {
			ch := item[max-1]
			if ch != '.' && ch != ',' && ch != ':' && ch != ';' {
				break
			}
			// trim trailing period, comma, colon, and semicolon
			item = item[:max-1]
			// continue checking for runs of punctuation at end
			max--
		}
		if item == "" {
			continue
		}
		arry = append(arry, item)
	}

	// rejoin into string
	cleaned := strings.Join(arry, " ")

	// break clauses at punctuation other than space or underscore, and at non-ASCII characters
	clauses := strings.FieldsFunc(cleaned, func(c rune) bool {
		return (!unicode.IsLetter(c) && !unicode.IsDigit(c)) && c != ' ' && c != '_' || c > 127
	})

	// space replaces plus sign to separate runs of unpunctuated words
	phrases := strings.Join(clauses, " ")

	var chain []string

	// break phrases into individual words
	words := strings.Fields(phrases)

	for _, item := range words {

		// skip at site of punctuation break
		if item == "+" {
			chain = append(chain, "+")
			continue
		}

		// skip terms that are all digits
		if IsAllDigitsOrPeriod(item) {
			chain = append(chain, "+")
			continue
		}

		// optional stop word removal
		if deStop && IsStopWord(item) {
			chain = append(chain, "+")
			continue
		}

		// index single normalized term
		chain = append(chain, item)
	}

	// rejoin into processed sentence
	tmp := strings.Join(chain, " ")

	tmp = strings.Replace(tmp, "+ +", "++", -1)
	tmp = strings.Replace(tmp, "+ +", "++", -1)

	tmp = CompressRunsOfSpaces(tmp)
	tmp = strings.TrimSpace(tmp)

	if tmp != "" && !strings.HasSuffix(tmp, "]") {
		tmp += " " + sfx
	}

	return tmp
}

func processStopWords(str string, deStop bool) string {

	if str == "" {
		return ""
	}

	var chain []string

	terms := strings.Fields(str)

	nextField := func(terms []string) (string, int) {

		for j, item := range terms {
			if strings.HasPrefix(item, "[") && strings.HasSuffix(item, "]") {
				return strings.ToUpper(item), j + 1
			}
		}

		return "", 0
	}

	// replace unwanted and stop words with plus sign
	for len(terms) > 0 {

		item := terms[0]
		terms = terms[1:]

		fld, j := nextField(terms)

		// with addition of TITL field, switch from TIAB to NORM
		if fld == "[NORM]" {
			fld = "[TIAB]"
		}

		stps := false
		rlxd := false
		if fld == "[TITL]" || fld == "[TIAB]" {
			stps = true
		} else if fld == "[STEM]" {
			stps = true
			rlxd = true
		} else if fld == "" {
			stps = true
		}

		addOneTerm := func(itm string) {

			if stps {
				if IsAllDigitsOrPeriod(itm) {
					// skip terms that are all digits
					chain = append(chain, "+")
				} else if deStop && IsStopWord(itm) {
					// skip if stop word, breaking phrase chain
					chain = append(chain, "+")
				} else if rlxd {
					isWildCard := strings.HasSuffix(itm, "*")
					if isWildCard {
						// temporarily remove trailing asterisk
						itm = strings.TrimSuffix(itm, "*")
					}

					itm = porter2.Stem(itm)
					itm = strings.TrimSpace(itm)

					if isWildCard {
						// do wildcard search in stemmed term list
						itm += "*"
					}
					chain = append(chain, itm)
				} else {
					// record single unmodified term
					chain = append(chain, itm)
				}
			} else {
				// do not treat non-TIAB terms as stop words
				chain = append(chain, itm)
			}
		}

		if j == 0 {
			// index single normalized term
			addOneTerm(item)
			continue
		}

		for j > 0 {

			addOneTerm(item)

			j--
			item = terms[0]
			terms = terms[1:]
		}

		if fld != "" {
			chain = append(chain, fld)
		}
	}

	// rejoin into processed sentence
	tmp := strings.Join(chain, " ")

	tmp = strings.Replace(tmp, "+ +", "++", -1)
	tmp = strings.Replace(tmp, "+ +", "++", -1)

	tmp = strings.Replace(tmp, "~ +", "~+", -1)
	tmp = strings.Replace(tmp, "+ ~", "+~", -1)

	for strings.Contains(tmp, "~+") {
		tmp = strings.Replace(tmp, "~+", "~~", -1)
	}
	for strings.Contains(tmp, "+~") {
		tmp = strings.Replace(tmp, "+~", "~~", -1)
	}

	tmp = CompressRunsOfSpaces(tmp)
	tmp = strings.TrimSpace(tmp)

	return tmp
}

func partitionQuery(str string) []string {

	if str == "" {
		return nil
	}

	str = CompressRunsOfSpaces(str)
	str = strings.TrimSpace(str)

	str = " " + str + " "

	// flank all operators with caret
	str = strings.Replace(str, " ( ", " ^ ( ^ ", -1)
	str = strings.Replace(str, " ) ", " ^ ) ^ ", -1)
	str = strings.Replace(str, " & ", " ^ & ^ ", -1)
	str = strings.Replace(str, " | ", " ^ | ^ ", -1)
	str = strings.Replace(str, " ! ", " ^ ! ^ ", -1)
	str = strings.Replace(str, " ~", " ^ ~", -1)
	str = strings.Replace(str, "~ ", "~ ^ ", -1)

	str = CompressRunsOfSpaces(str)
	str = strings.TrimSpace(str)

	str = strings.Replace(str, "^ ^", "^", -1)

	if strings.HasPrefix(str, "^ ") {
		str = str[2:]
	}
	if strings.HasSuffix(str, " ^") {
		max := len(str)
		str = str[:max-2]
	}

	str = strings.Replace(str, "~ ^ +", "~+", -1)
	str = strings.Replace(str, "+ ^ ~", "+~", -1)

	str = strings.Replace(str, "~ +", "~+", -1)
	str = strings.Replace(str, "+ ~", "+~", -1)

	for strings.Contains(str, "~+") {
		str = strings.Replace(str, "~+", "~~", -1)
	}
	for strings.Contains(str, "+~") {
		str = strings.Replace(str, "+~", "~~", -1)
	}

	// split into non-broken phrase segments or operator symbols
	tmp := strings.Split(str, " ^ ")

	return tmp
}

func setFieldQualifiers(clauses []string, rlxd bool) []string {

	var res []string

	if clauses == nil {
		return nil
	}

	for _, str := range clauses {

		// pass control symbols unchanged
		if str == "(" || str == ")" || str == "&" || str == "|" || str == "!" || strings.HasPrefix(str, "~") {
			res = append(res, str)
			continue
		}

		// pass angle bracket content delimiters (for -phrase, -require, -exclude)
		if str == "<" || str == ">" {
			res = append(res, str)
			continue
		}

		if strings.HasSuffix(str, " [YEAR]") {

			slen := len(str)
			str = str[:slen-7]

			// regular 4-digit year
			if len(str) == 4 && IsAllDigitsOrPeriod(str) {
				res = append(res, str+" [YEAR]")
				continue
			}

			// check for year wildcard
			if len(str) == 4 && str[3] == '*' && IsAllDigitsOrPeriod(str[:3]) {

				fmt.Fprintf(os.Stderr, "\nERROR: Wildcards not supported for years - use ####:#### range instead\n")
				os.Exit(1)
			}

			// check for year range
			if len(str) == 9 && str[4] == ' ' && IsAllDigitsOrPeriod(str[:4]) && IsAllDigitsOrPeriod(str[5:]) {
				start, err := strconv.Atoi(str[:4])
				if err != nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to recognize starting year '%s'\n", str[:4])
					os.Exit(1)
				}
				stop, err := strconv.Atoi(str[5:])
				if err != nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to recognize stopping year '%s'\n", str[5:])
					os.Exit(1)
				}
				if start > stop {
					continue
				}
				// expand year range into individual year-by-year queries
				pfx := "("
				sfx := ")"
				for start <= stop {
					res = append(res, pfx)
					pfx = "|"
					yr := strconv.Itoa(start)
					res = append(res, yr+" [year]")
					start++
				}
				res = append(res, sfx)
				continue
			}

			fmt.Fprintf(os.Stderr, "\nERROR: Unable to recognize year expression '%s'\n", str)
			os.Exit(1)

		} else if strings.HasSuffix(str, " [TREE]") {

			slen := len(str)
			str = str[:slen-7]

			// pad if top-level mesh tree wildcard uses four character trie
			if len(str) == 4 && str[3] == '*' {
				key := str[:2]
				num, ok := TrieLen[key]
				if ok && num > 3 {
					str = str[0:3] + " " + "*"
				}
			}

			str = strings.Replace(str, " ", ".", -1)
			tmp := str
			tmp = strings.TrimSuffix(tmp, "*")
			if len(tmp) > 2 && unicode.IsLower(rune(tmp[0])) && IsAllDigitsOrPeriod(tmp[1:]) {
				str = strings.Replace(str, ".", " ", -1)
				res = append(res, str+" [TREE]")
				continue
			}

			fmt.Fprintf(os.Stderr, "\nERROR: Unable to recognize mesh code expression '%s'\n", str)
			os.Exit(1)
		}

		// remove leading and trailing plus signs and spaces
		for strings.HasPrefix(str, "+") || strings.HasPrefix(str, " ") {
			str = str[1:]
		}
		for strings.HasSuffix(str, "+") || strings.HasSuffix(str, " ") {
			slen := len(str)
			str = str[:slen-1]
		}

		res = append(res, str)
	}

	return res
}

// SEARCH TERM LISTS FOR PHRASES OR NORMALIZED TERMS, OR MATCH BY PATTERN

// ProcessSearch evaluates query, returns list of PMIDs
func ProcessSearch(base, phrase string, xact, titl, rlxd, deStop bool) int {

	if phrase == "" {
		return 0
	}

	if titl {
		phrase = prepareExact(phrase, "[titl]", deStop)
	} else if xact {
		phrase = prepareExact(phrase, "[tiab]", deStop)
	} else {
		phrase = prepareQuery(phrase)
	}

	phrase = processStopWords(phrase, deStop)

	clauses := partitionQuery(phrase)

	clauses = setFieldQualifiers(clauses, rlxd)

	return evaluateQuery(base, clauses)
}

// ProcessMock shows individual steps in processing query for evaluation
func ProcessMock(base, phrase string, xact, titl, rlxd, deStop bool) int {

	if phrase == "" {
		return 0
	}

	fmt.Fprintf(os.Stdout, "processSearch:\n\n%s\n\n", phrase)

	if titl {
		phrase = prepareExact(phrase, "[titl]", deStop)

		fmt.Fprintf(os.Stdout, "prepareExact:\n\n%s\n\n", phrase)
	} else if xact {
		phrase = prepareExact(phrase, "[tiab]", deStop)

		fmt.Fprintf(os.Stdout, "prepareExact:\n\n%s\n\n", phrase)
	} else {
		phrase = prepareQuery(phrase)

		fmt.Fprintf(os.Stdout, "prepareQuery:\n\n%s\n\n", phrase)
	}

	phrase = processStopWords(phrase, deStop)

	fmt.Fprintf(os.Stdout, "processStopWords:\n\n%s\n\n", phrase)

	clauses := partitionQuery(phrase)

	fmt.Fprintf(os.Stdout, "partitionQuery:\n\n")
	for _, tkn := range clauses {
		fmt.Fprintf(os.Stdout, "%s\n", tkn)
	}
	fmt.Fprintf(os.Stdout, "\n")

	clauses = setFieldQualifiers(clauses, rlxd)

	fmt.Fprintf(os.Stdout, "setFieldQualifiers:\n\n")
	for _, tkn := range clauses {
		fmt.Fprintf(os.Stdout, "%s\n", tkn)
	}
	fmt.Fprintf(os.Stdout, "\n")

	return 0
}

// ProcessCount prints document count for each term, also supports terminal wildcards
func ProcessCount(base, phrase string, plrl, psns, rlxd, deStop bool) int {

	if phrase == "" {
		return 0
	}

	phrase = prepareQuery(phrase)

	phrase = processStopWords(phrase, deStop)

	clauses := partitionQuery(phrase)

	clauses = setFieldQualifiers(clauses, rlxd)

	if clauses == nil {
		return 0
	}

	count := 0

	splitIntoWords := func(str string) []string {

		if str == "" {
			return nil
		}

		var arry []string

		parts := strings.Split(str, "+")

		for _, segment := range parts {

			segment = strings.TrimSpace(segment)

			if segment == "" {
				continue
			}

			words := strings.Fields(segment)

			for _, item := range words {
				if strings.HasPrefix(item, "~") {
					continue
				}
				arry = append(arry, item)
			}
		}

		return arry
	}

	parseField := func(str string) (string, string) {

		field := "TIAB"

		if strings.HasSuffix(str, "]") {
			pos := strings.Index(str, "[")
			if pos >= 0 {
				field = str[pos:]
				field = strings.TrimPrefix(field, "[")
				field = strings.TrimSuffix(field, "]")
				str = str[:pos]
				str = strings.TrimSpace(str)
			}
			switch field {
			case "NORM":
				field = "TIAB"
			case "STEM", "TIAB", "TITL":
			case "PIPE":
			default:
				str = strings.Replace(str, " ", "_", -1)
			}
		}

		return field, str
	}

	checkTermCounts := func(txt string) {

		field, str := parseField(txt)

		var words []string

		words = splitIntoWords(str)

		if words == nil || len(words) < 1 {
			return
		}

		for _, term := range words {

			term = strings.Replace(term, "_", " ", -1)

			if psns {
				count += printTermPositions(base, term, field)
			} else if plrl {
				count += printTermCounts(base, term, field)
			} else {
				count += printTermCount(base, term, field)
			}
		}
	}

	for _, item := range clauses {

		// skip control symbols
		if item == "(" || item == ")" || item == "&" || item == "|" || item == "!" {
			continue
		}

		checkTermCounts(item)
	}

	runtime.Gosched()

	return count
}
