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
// File Name:  invert.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"fmt"
	"github.com/rainycape/unidecode"
	"os"
	"runtime"
	"sort"
	"strings"
	"sync"
	"unicode"
)

// CreateDispensers collects field, uid, positions, for each term
func CreateDispensers(inp <-chan XMLRecord) <-chan []string {

	if inp == nil {
		return nil
	}

	out := make(chan []string, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create dispenser channel\n")
		os.Exit(1)
	}

	var ilock sync.Mutex

	// map for inverted index
	inverted := make(map[string][]string)

	// add single posting
	addPost := func(fld, term, pos, uid string) {

		// protect map with mutex
		ilock.Lock()

		data, ok := inverted[term]
		if !ok {
			data = make([]string, 0, 4)
			// first entry on new slice is term
			data = append(data, term)
		}
		data = append(data, fld)
		data = append(data, uid)
		data = append(data, pos)
		// always need to update inverted, since data may be reallocated
		inverted[term] = data

		// unlock at end to avoid defer overhead
		ilock.Unlock()
	}

	// xmlDispenser prepares UID, term, and position strings for inversion
	xmlDispenser := func(wg *sync.WaitGroup, inp <-chan XMLRecord, out chan<- []string) {

		defer wg.Done()

		currUID := ""

		doDispense := func(tag, attr, content string) {

			if tag == "IdxUid" {
				currUID = content
			} else {

				// expand Greek letters, anglicize characters in other alphabets
				if IsNotASCII(content) {

					if HasGreek(content) {
						content = SpellGreek(content)
						content = CompressRunsOfSpaces(content)
					}
					content = unidecode.Unidecode(content)

					content = DoAccentTransform(content)
					content = UnicodeToASCII(content)

					content = strings.TrimSpace(content)
				}

				content = strings.ToLower(content)

				// remove punctuation from term
				content = strings.Map(func(c rune) rune {
					if !unicode.IsLetter(c) && !unicode.IsDigit(c) && c != ' ' && c != '-' && c != '_' {
						return -1
					}
					return c
				}, content)

				content = strings.Replace(content, "_", " ", -1)
				content = strings.Replace(content, "-", " ", -1)

				content = CompressRunsOfSpaces(content)
				content = strings.TrimSpace(content)

				if content != "" && currUID != "" {
					addPost(tag, content, attr, currUID)
				}
			}
		}

		// read partitioned XML from producer channel
		for ext := range inp {

			StreamValues(ext.Text[:], "IdxDocument", doDispense)
		}
	}

	var wg sync.WaitGroup

	// launch multiple dispenser goroutines
	for i := 0; i < NumServe(); i++ {
		wg.Add(1)
		go xmlDispenser(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all dispensers are done
	go func() {
		wg.Wait()

		// send results to inverters
		for _, data := range inverted {
			out <- data

			runtime.Gosched()
		}

		close(out)
	}()

	return out
}

// CreateInverters sorts UIDs and positions for each term in each field
func CreateInverters(inp <-chan []string) <-chan XMLRecord {

	if inp == nil {
		return nil
	}

	out := make(chan XMLRecord, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create inverter channel\n")
		os.Exit(1)
	}

	// xmlInverter sorts and prints one posting list
	xmlInverter := func(wg *sync.WaitGroup, inp <-chan []string, out chan<- XMLRecord) {

		defer wg.Done()

		var buffer strings.Builder

		printPosting := func(key string, data []string) string {

			fields := make(map[string]map[string]string)

			for len(data) > 1 {
				fld := data[0]
				uid := data[1]
				att := data[2]
				positions, ok := fields[fld]
				if !ok {
					positions = make(map[string]string)
					fields[fld] = positions
				}
				// store position attribute string by uid
				positions[uid] = att
				// skip to next position
				data = data[3:]
			}

			buffer.Reset()

			buffer.WriteString("  <InvDocument>\n")
			buffer.WriteString("    <InvKey>")
			buffer.WriteString(key)
			buffer.WriteString("</InvKey>\n")
			buffer.WriteString("    <InvIDs>\n")

			// sort fields in alphabetical order
			var keys []string
			for ky := range fields {
				keys = append(keys, ky)
			}
			sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })

			for _, fld := range keys {

				positions := fields[fld]

				var arry []string

				for item := range positions {
					arry = append(arry, item)
				}

				if len(arry) > 1 {
					sort.Slice(arry, func(i, j int) bool {
						// numeric sort on strings checks lengths first
						lni := len(arry[i])
						lnj := len(arry[j])
						// shorter string is numerically less, assuming no leading zeros
						if lni < lnj {
							return true
						}
						if lni > lnj {
							return false
						}
						// same length, can now do string comparison on contents
						return arry[i] < arry[j]
					})
				}

				// print list of UIDs, skipping duplicates
				prev := ""
				for _, uid := range arry {
					if uid == prev {
						continue
					}

					buffer.WriteString("      <")
					buffer.WriteString(fld)
					atr := positions[uid]
					if atr != "" {
						buffer.WriteString(" ")
						buffer.WriteString(atr)
					}
					buffer.WriteString(">")
					buffer.WriteString(uid)
					buffer.WriteString("</")
					buffer.WriteString(fld)
					buffer.WriteString(">\n")

					prev = uid
				}
			}

			buffer.WriteString("    </InvIDs>\n")
			buffer.WriteString("  </InvDocument>\n")

			str := buffer.String()

			return str
		}

		for inv := range inp {

			key := inv[0]
			data := inv[1:]

			str := printPosting(key, data)

			out <- XMLRecord{Ident: key, Text: str}

			runtime.Gosched()
		}
	}

	var wg sync.WaitGroup

	// launch multiple inverter goroutines
	for i := 0; i < NumServe(); i++ {
		wg.Add(1)
		go xmlInverter(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all inverters are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

// CreateResolver sorts postings by identifier prefix to prepare for multi-file merge
func CreateResolver(inp <-chan XMLRecord) <-chan string {

	if inp == nil {
		return nil
	}

	out := make(chan string, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create resolver channel\n")
		os.Exit(1)
	}

	// xmlResolver prints inverted postings alphabetized by identifier prefix
	xmlResolver := func(inp <-chan XMLRecord, out chan<- string) {

		// close channel when all records have been processed
		defer close(out)

		// map for inverted index
		inverted := make(map[string]string)

		// drain channel, populate map for alphabetizing
		for curr := range inp {

			inverted[curr.Ident] = curr.Text
		}

		var ordered []string

		for item := range inverted {
			ordered = append(ordered, item)
		}

		if len(ordered) > 1 {
			sort.Slice(ordered, func(i, j int) bool { return ordered[i] < ordered[j] })
		}

		// iterate through alphabetized results
		for _, curr := range ordered {

			txt := inverted[curr]

			// send result to output
			out <- txt

			runtime.Gosched()
		}
	}

	// launch single resolver goroutine
	go xmlResolver(inp, out)

	return out
}
