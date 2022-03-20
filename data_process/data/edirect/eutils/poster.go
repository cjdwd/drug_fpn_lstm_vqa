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
// File Name:  poster.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"encoding/binary"
	"fmt"
	"github.com/rainycape/unidecode"
	"io"
	"os"
	"path"
	"strconv"
	"strings"
	"sync"
	"unicode"
)

// MakeArchiveTrie allows a short prefix of letters with an optional underscore, and splits the remainder into character pairs
func MakeArchiveTrie(str string, arry [132]rune) string {

	if len(str) > 64 {
		return ""
	}

	if IsAllDigits(str) {

		// pad numeric identifier to 8 characters with leading zeros
		ln := len(str)
		if ln < 8 {
			zeros := "00000000"
			str = zeros[ln:] + str
		}
	}

	if IsAllDigitsOrPeriod(str) {

		// limit trie to first 6 characters
		if len(str) > 6 {
			str = str[:6]
		}
	}

	max := 4
	k := 0
	for _, ch := range str {
		if unicode.IsLetter(ch) {
			k++
			continue
		}
		if ch == '_' {
			k++
			max = 6
		}
		break
	}

	// prefix is up to three letters if followed by digits, or up to four letters if followed by an underscore
	pfx := str[:k]
	if len(pfx) < max {
		str = str[k:]
	} else {
		pfx = ""
	}

	i := 0

	if pfx != "" {
		for _, ch := range pfx {
			arry[i] = ch
			i++
		}
		arry[i] = '/'
		i++
	}

	between := 0
	doSlash := false

	// remainder is divided in character pairs, e.g., NP_/06/00/51 for NP_060051.2
	for _, ch := range str {
		// break at period separating accession from version
		if ch == '.' {
			break
		}
		if doSlash {
			arry[i] = '/'
			i++
			doSlash = false
		}
		if ch == ' ' {
			ch = '_'
		}
		if !unicode.IsLetter(ch) && !unicode.IsDigit(ch) {
			ch = '_'
		}
		arry[i] = ch
		i++
		between++
		if between > 1 {
			doSlash = true
			between = 0
		}
	}

	res := string(arry[:i])

	if !strings.HasSuffix(res, "/") {
		arry[i] = '/'
		i++
		res = string(arry[:i])
	}

	return strings.ToUpper(res)
}

// MakePostingsTrie splits a string into characters, separated by path delimiting slashes
func MakePostingsTrie(str string, arry [516]rune) string {

	if len(str) > 256 {
		return ""
	}

	// expand Greek letters, anglicize characters in other alphabets
	if IsNotASCII(str) {
		if HasGreek(str) {
			str = SpellGreek(str)
			str = CompressRunsOfSpaces(str)
		}
		str = unidecode.Unidecode(str)
		str = strings.TrimSpace(str)
	}

	i := 0
	doSlash := false

	for _, ch := range str {
		if doSlash {
			arry[i] = '/'
			i++
		}
		if ch == ' ' {
			ch = '_'
		}
		if !unicode.IsLetter(ch) && !unicode.IsDigit(ch) {
			ch = '_'
		}
		arry[i] = ch
		i++
		doSlash = true
	}

	return strings.ToLower(string(arry[:i]))
}

// PostingPath constructs a Postings directory subpath for a given term prefix
func PostingPath(prom, field, term string, arry [516]rune) (string, string) {

	// use first few characters of identifier for directory
	dir := IdentifierKey(term)

	trie := MakePostingsTrie(dir, arry)
	if trie == "" {
		return "", ""
	}

	dpath := path.Join(prom, field, trie)

	return dpath, dir
}

// CreatePromoters creates term lists and postings files from merged inverted
// index files. All specified fields (e.g., "TITL TIAB YEAR TREE") are completed
// in a single scan of the inverted files.
//
// Postings consist of three files (.mst, .trm, and .pst) for all terms, plus
// two additional files (.uqi and .ofs) for terms with position data.
//
// Master index files (with .mst suffixes) contain pairs of 32-bit values, in
// little endian form, pointing to an offset into the term list (.trm files,
// saved as lines of text words), and an offset into the postings list (.pst
// files,  containing 32-bit PMIDs).
//
// For position data, the .uqi file is parallel to the .pst file (one entry
// for each PMID associated with a given term), and contains 32-bit offsets
// to 16-bit paragraph position values.
//
// An extra entry at the end, pointing just past the end of data, allows
// the length of a term, or the size of a postings list, or the number of
// positions for a term in a given PMID, to be calculated as the difference
// between two adjacent pointers.
//
// The number of term positions per PMID is the term frequency (TF). The
// number of PMIDs per term is the document frequency (DF). All that remains
// for calculating TF-IDF term weights, which can support ranked retrieval,
// is the total number of live PubMed documents, which could easily be saved
// during indexing.
//
func CreatePromoters(prom, fields string, files []string) <-chan string {

	if files == nil {
		return nil
	}

	out := make(chan string, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create promoter channel\n")
		os.Exit(1)
	}

	flds := strings.Split(fields, " ")

	// xmlPromoter saves records in a single set of term/posting files
	xmlPromoter := func(wg *sync.WaitGroup, fileName string, out chan<- string) {

		defer wg.Done()

		f, err := os.Open(fileName)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
			os.Exit(1)
		}

		// close input file when all records have been processed
		defer f.Close()

		var in io.Reader

		in = f

		// if suffix is ".gz", use decompressor
		iszip := false
		if strings.HasSuffix(fileName, ".gz") {
			iszip = true
		}

		if iszip {
			brd := bufio.NewReader(f)
			if brd == nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create buffered reader on '%s'\n", fileName)
				os.Exit(1)
			}
			zpr, err := gzip.NewReader(brd)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create decompressor on '%s'\n", fileName)
				os.Exit(1)
			}

			// close decompressor when all records have been processed
			defer zpr.Close()

			// use decompressor for reading file
			in = zpr
		}

		rdr := CreateXMLStreamer(in)

		if rdr == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
			os.Exit(1)
		}

		getOnePosting := func(field, text string) (string, []int32, []string) {

			var data []int32
			var atts []string

			term := ""

			doPromote := func(tag, attr, content string) {

				if tag == "InvKey" {

					// term used for postings file name
					term = content

					term = strings.ToLower(term)

				} else if tag == field {

					// convert UID string to integer
					if content == "" {
						fmt.Fprintf(os.Stderr, "\nERROR: Empty UID for term '%s'\n", term)
						return
					}
					value, err := strconv.ParseInt(content, 10, 32)
					if err != nil {
						fmt.Fprintf(os.Stderr, "%s\n", err.Error())
						return
					}
					data = append(data, int32(value))

					if strings.HasPrefix(attr, "pos=\"") {
						attr = attr[5:]
						lgth := len(attr)
						if lgth > 1 && attr[lgth-1] == '"' {
							// "
							attr = attr[:lgth-1]
						}
						atts = append(atts, attr)
					}
				}
			}

			// explore data fields
			StreamValues(text[:], "InvDocument", doPromote)

			if term == "" || len(data) < 1 {
				return "", nil, nil
			}

			return term, data, atts
		}

		var (
			termPos int32
			postPos int32
			ofstPos int32

			indxList bytes.Buffer
			termList bytes.Buffer
			postList bytes.Buffer
			uqidList bytes.Buffer
			ofstList bytes.Buffer
		)

		retlength := len("\n")

		addOnePosting := func(term string, data []int32, atts []string) {

			tlength := len(term)
			dlength := len(data)
			alength := len(atts)

			// write to term list buffer
			termList.WriteString(term[:])
			termList.WriteString("\n")

			// write to postings buffer
			binary.Write(&postList, binary.LittleEndian, data)

			// write to master index buffer
			binary.Write(&indxList, binary.LittleEndian, termPos)
			binary.Write(&indxList, binary.LittleEndian, postPos)

			postPos += int32(dlength * 4)
			termPos += int32(tlength + retlength)

			// return if no position attributes
			if alength < 1 {
				return
			}
			if dlength != alength {
				fmt.Fprintf(os.Stderr, "dlength %d, alength %d\n", dlength, alength)
				return
			}

			// write term offset list for each UID
			for _, attr := range atts {

				binary.Write(&uqidList, binary.LittleEndian, ofstPos)

				atrs := strings.Split(attr, ",")
				atln := len(atrs)
				for _, att := range atrs {
					if att == "" {
						continue
					}
					value, err := strconv.ParseInt(att, 10, 32)
					if err != nil {
						fmt.Fprintf(os.Stderr, "%s\n", err.Error())
						return
					}
					binary.Write(&ofstList, binary.LittleEndian, int16(value))
				}

				ofstPos += int32(atln * 2)
			}
		}

		topOffMaster := func() {

			// phantom term and postings positions eliminates special case calculation at end
			binary.Write(&indxList, binary.LittleEndian, termPos)
			binary.Write(&indxList, binary.LittleEndian, postPos)
			binary.Write(&uqidList, binary.LittleEndian, ofstPos)
		}

		writeFile := func(dpath, fname string, bfr bytes.Buffer) {

			fpath := path.Join(dpath, fname)
			if fpath == "" {
				return
			}

			// overwrites and truncates existing file
			fl, err := os.Create(fpath)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
				return
			}

			data := bfr.Bytes()

			wrtr := bufio.NewWriter(fl)

			_, err = wrtr.Write(data)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			}

			wrtr.Flush()

			// fl.Sync()

			fl.Close()
		}

		writeFiveFiles := func(field, key string) {

			var arry [516]rune
			dpath, key := PostingPath(prom, field, key, arry)
			if dpath == "" {
				return
			}

			// make subdirectories, if necessary
			err := os.MkdirAll(dpath, os.ModePerm)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
				return
			}

			writeFile(dpath, key+"."+field+".trm", termList)

			writeFile(dpath, key+"."+field+".pst", postList)

			writeFile(dpath, key+"."+field+".mst", indxList)

			// do not write position index and offset data files for fields with no position attributes recorded
			if uqidList.Len() > 0 && ofstList.Len() > 0 {

				writeFile(dpath, key+"."+field+".uqi", uqidList)

				writeFile(dpath, key+"."+field+".ofs", ofstList)
			}
		}

		processOneField := func(field string, recs []string) {

			tag := ""

			for _, str := range recs {

				term, data, atts := getOnePosting(field, str)

				if term == "" || data == nil {
					continue
				}

				// use first few characters of identifier
				if tag == "" {
					tag = IdentifierKey(term)
				}

				addOnePosting(term, data, atts)
			}

			if tag != "" {
				topOffMaster()
				writeFiveFiles(field, tag)
			}

			// reset buffers and position counters
			termPos = 0
			postPos = 0
			ofstPos = 0

			indxList.Reset()
			termList.Reset()
			postList.Reset()
			uqidList.Reset()
			ofstList.Reset()
		}

		find := ParseIndex("InvKey")

		currTag := ""
		prevTag := ""

		var arry []string

		// read next array of InvDocument records with same key
		PartitionPattern("InvDocument", "", false, rdr,
			func(str string) {

				id := FindIdentifier(str[:], "InvDocument", find)
				if id == "" {
					return
				}

				// use first few characters of identifier
				currTag = IdentifierKey(id)

				if prevTag != currTag {

					// after IdentifierKey converts space to underscore,
					// okay that xxx_ and xxx0 will be out of alphabetical order

					// records with same identifier key as a unit
					if prevTag != "" {
						for _, fld := range flds {
							processOneField(fld, arry)
						}
						out <- prevTag
					}

					// empty the slice
					arry = nil
				}

				// collect next InvDocument record
				arry = append(arry, str[:])

				prevTag = currTag
			})

		if arry != nil {

			// remaining records with last identifier key
			for _, fld := range flds {
				processOneField(fld, arry)
			}
			out <- prevTag
		}
	}

	var wg sync.WaitGroup

	// launch multiple promoter goroutines
	for _, str := range files {
		wg.Add(1)
		go xmlPromoter(&wg, str, out)
	}

	// launch separate anonymous goroutine to wait until all promoters are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}
