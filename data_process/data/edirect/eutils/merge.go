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
// File Name:  merge.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"bufio"
	"compress/gzip"
	"container/heap"
	"fmt"
	"io"
	"os"
	"path"
	"runtime"
	"runtime/debug"
	"sort"
	"strings"
	"sync"
	"unicode"
)

// TrieLen directory depth parameters are based on the observed size distribution of PubMed indices
var TrieLen = map[string]int{
	"19": 4,
	"20": 4,
	"a1": 3,
	"ab": 3,
	"ac": 4,
	"ad": 3,
	"af": 4,
	"ag": 3,
	"al": 3,
	"an": 4,
	"ap": 4,
	"ar": 3,
	"as": 4,
	"b0": 3,
	"ba": 4,
	"be": 4,
	"bi": 3,
	"br": 3,
	"c0": 3,
	"c1": 3,
	"ca": 4,
	"ce": 4,
	"ch": 4,
	"cl": 4,
	"co": 4,
	"cr": 3,
	"cy": 3,
	"d0": 4,
	"d1": 4,
	"d2": 3,
	"da": 4,
	"de": 4,
	"di": 4,
	"do": 3,
	"dr": 3,
	"e0": 3,
	"ef": 4,
	"en": 3,
	"ev": 3,
	"ex": 4,
	"fa": 3,
	"fi": 3,
	"fo": 4,
	"fr": 4,
	"fu": 4,
	"g0": 3,
	"ge": 4,
	"gr": 4,
	"he": 4,
	"hi": 4,
	"im": 3,
	"in": 4,
	"la": 3,
	"le": 3,
	"li": 3,
	"lo": 3,
	"ma": 3,
	"me": 4,
	"mi": 3,
	"mo": 4,
	"mu": 3,
	"mz": 3,
	"n0": 3,
	"ne": 3,
	"no": 4,
	"ob": 3,
	"on": 3,
	"oz": 3,
	"pa": 4,
	"pe": 4,
	"ph": 3,
	"pl": 4,
	"po": 4,
	"pr": 4,
	"ra": 3,
	"re": 4,
	"ri": 3,
	"rz": 3,
	"se": 3,
	"si": 4,
	"sp": 4,
	"st": 4,
	"su": 4,
	"sy": 4,
	"te": 3,
	"th": 3,
	"ti": 3,
	"tr": 4,
	"tw": 4,
	"un": 3,
	"va": 3,
	"ve": 3,
	"vi": 3,
	"we": 3,
	"wh": 3,
}

// MergLen directory depth parameters are based on the observed size distribution of PubMed indices
var MergLen = map[string]int{
	"ana": 4,
	"app": 4,
	"ass": 4,
	"can": 4,
	"cas": 4,
	"cha": 4,
	"cli": 4,
	"com": 4,
	"con": 4,
	"d00": 4,
	"d01": 4,
	"d02": 4,
	"d12": 4,
	"dam": 4,
	"dat": 4,
	"dec": 4,
	"ded": 4,
	"del": 4,
	"dem": 4,
	"dep": 4,
	"des": 4,
	"det": 4,
	"dif": 4,
	"dis": 4,
	"eff": 4,
	"exp": 4,
	"for": 4,
	"gen": 4,
	"gro": 4,
	"hea": 4,
	"hig": 4,
	"inc": 4,
	"ind": 4,
	"int": 4,
	"inv": 4,
	"met": 4,
	"mod": 4,
	"pat": 4,
	"per": 4,
	"pre": 4,
	"pro": 4,
	"rel": 4,
	"rep": 4,
	"res": 4,
	"sig": 4,
	"sta": 4,
	"str": 4,
	"stu": 4,
	"tre": 4,
}

// PostingDir returns directory trie (without slashes) for location of indices for a given term
func PostingDir(term string) string {

	if len(term) < 3 {
		return term
	}

	key := term[:2]

	num, ok := TrieLen[key]
	if ok && len(term) >= num {
		return term[:num]
	}

	switch term[0] {
	case 'u', 'v', 'w', 'x', 'y', 'z':
		return term[:2]
	}

	return term[:3]
}

// IdentifierKey cleans up a term then returns the posting directory
func IdentifierKey(term string) string {

	// remove punctuation from term
	key := strings.Map(func(c rune) rune {
		if !unicode.IsLetter(c) && !unicode.IsDigit(c) && c != ' ' && c != '-' && c != '_' {
			return -1
		}
		return c
	}, term)

	key = strings.Replace(key, " ", "_", -1)
	key = strings.Replace(key, "-", "_", -1)

	// use first 2, 3, or 4 characters of identifier for directory
	key = PostingDir(key)

	return key
}

// Plex allows distribution of indexing
type Plex struct {
	Which int
	Ident string
	Text  string
	Index int
	Sibs  []string
}

// PlexHeap methods satisfy heap.Interface
type PlexHeap []Plex

func (h PlexHeap) Len() int {
	return len(h)
}
func (h PlexHeap) Less(i, j int) bool {
	return h[i].Ident < h[j].Ident
}
func (h PlexHeap) Swap(i, j int) {
	h[i], h[j] = h[j], h[i]
}

// Push works on pointer to PlexHeap
func (h *PlexHeap) Push(x interface{}) {
	*h = append(*h, x.(Plex))
}

// Pop works on pointer to PlexHeap
func (h *PlexHeap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

// CreatePresenters creates one channel per input file
func CreatePresenters(files []string) []<-chan Plex {

	if files == nil {
		return nil
	}

	numFiles := len(files)
	if numFiles < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Not enough inverted files to merge\n")
		os.Exit(1)
	}

	chns := make([]<-chan Plex, numFiles)
	if chns == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create presenter channel array\n")
		os.Exit(1)
	}

	// xmlPresenter sends partitioned XML strings through channel
	xmlPresenter := func(fileNum int, fileName string, out chan<- Plex) {

		// close channel when all records have been processed
		defer close(out)

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

		find := ParseIndex("InvKey")

		// partition all input by pattern and send XML substring through channel
		PartitionPattern("InvDocument", "", false, rdr,
			func(str string) {
				id := FindIdentifier(str[:], "InvDocument", find)

				out <- Plex{fileNum, id, str, 0, nil}
			})
	}

	// launch multiple presenter goroutines
	for i, str := range files {

		chn := make(chan Plex, ChanDepth())
		if chn == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create presenter channel\n")
			os.Exit(1)
		}

		go xmlPresenter(i, str, chn)

		chns[i] = chn
	}

	// no need for separate anonymous goroutine to wait until all presenters are done

	return chns
}

// CreateManifold reads from each file, sends merged postings in sorted order
func CreateManifold(inp []<-chan Plex) <-chan Plex {

	if inp == nil {
		return nil
	}

	out := make(chan Plex, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create manifold channel\n")
		os.Exit(1)
	}

	// xmlManifold restores alphabetical order of merged postings
	xmlManifold := func(inp []<-chan Plex, out chan<- Plex) {

		// close channel when all records have been processed
		defer close(out)

		// initialize empty heap
		hp := &PlexHeap{}
		heap.Init(hp)

		// read first object from all input channels in turn
		for _, chn := range inp {
			plx, ok := <-chn
			if ok {
				heap.Push(hp, plx)
			}
		}

		// array to collect strings with same identifier
		var arry []string

		prevIdent := ""
		rec := 0

		// reading from heap returns objects in alphabetical order
		for hp.Len() > 0 {

			// remove lowest item from heap, use interface type assertion
			curr := heap.Pop(hp).(Plex)

			// compare adjacent record identifiers
			if prevIdent == curr.Ident {

				// save next inverted object string in slice
				arry = append(arry, curr.Text)

			} else {

				if len(arry) > 0 {

					rec++
					// send set from previous identifier to output channel
					out <- Plex{0, prevIdent, "", rec, arry}

					// empty the slice
					arry = nil

					runtime.Gosched()
				}

				// remember new identifier
				prevIdent = curr.Ident

				// save first inverted object with this identifier
				arry = append(arry, curr.Text)
			}

			// read next object from channel that just supplied lowest item
			chn := inp[curr.Which]
			plx, ok := <-chn
			if ok {
				heap.Push(hp, plx)
			}
		}

		if len(arry) > 0 {

			rec++
			// send last record
			out <- Plex{0, prevIdent, "", rec, arry}

			arry = nil
		}
	}

	// launch single manifold goroutine
	go xmlManifold(inp, out)

	return out
}

// CreateFusers collects all inverted indices for a given term
func CreateFusers(inp <-chan XMLRecord) <-chan Plex {

	if inp == nil {
		return nil
	}

	out := make(chan Plex, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create fuser channel\n")
		os.Exit(1)
	}

	var flock sync.Mutex

	// map for combining inverted indices
	inverted := make(map[string][]string)

	addInverts := func(id, str string) {

		// protect map with mutex
		flock.Lock()

		data, ok := inverted[id]
		if !ok {
			data = make([]string, 0, 1)
		}

		data = append(data, str)
		// always need to update inverted, since data may be reallocated
		inverted[id] = data

		// unlock at end to avoid defer overhead
		flock.Unlock()
	}

	xmlFuser := func(wg *sync.WaitGroup, inp <-chan XMLRecord, out chan<- Plex) {

		defer wg.Done()

		find := ParseIndex("InvKey")

		// read partitioned XML from producer channel
		for ext := range inp {

			str := ext.Text[:]
			id := FindIdentifier(str, "InvDocument", find)
			addInverts(id, str)
		}
	}

	var wg sync.WaitGroup

	// launch multiple fuser goroutines
	for i := 0; i < NumServe(); i++ {
		wg.Add(1)
		go xmlFuser(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all fusers are done
	go func() {
		wg.Wait()

		// sort id keys in alphabetical order
		var keys []string
		for ky := range inverted {
			keys = append(keys, ky)
		}
		sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })

		rec := 0

		for _, id := range keys {

			arry := inverted[id]

			rec++
			// send array of records with same identifier to output channel
			out <- Plex{0, id, "", rec, arry}

			// empty the slice
			arry = nil
		}

		close(out)
	}()

	return out
}

// CreateMergers combines collected indices for the same term
func CreateMergers(inp <-chan Plex) <-chan XMLRecord {

	if inp == nil {
		return nil
	}

	out := make(chan XMLRecord, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create merger channel\n")
		os.Exit(1)
	}

	// xmlMerger fuses adjacent InvDocument records with the same identifier
	xmlMerger := func(wg *sync.WaitGroup, inp <-chan Plex, out chan<- XMLRecord) {

		defer wg.Done()

		var buffer strings.Builder

		fusePostings := func(key string, data []string) string {

			fields := make(map[string]map[string]string)

			addIdents := func(fld, pos, uid string) {

				// no need for mutex here
				positions, ok := fields[fld]
				if !ok {
					positions = make(map[string]string)
					fields[fld] = positions
				}

				positions[uid] = pos
			}

			addUID := func(tag, attr, content string) {

				if tag != "InvKey" {

					addIdents(tag, attr, content)
				}
			}

			for _, str := range data {
				StreamValues(str[:], "InvDocument", addUID)
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
				last := ""
				for _, uid := range arry {
					// detect duplicate UIDs, now in same list after conversion of one term entry from foreign alphabet
					if uid == last {
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

					last = uid
				}
			}

			buffer.WriteString("    </InvIDs>\n")
			buffer.WriteString("  </InvDocument>\n")

			txt := buffer.String()

			return txt
		}

		for plx := range inp {

			rec := plx.Index
			key := plx.Ident
			data := plx.Sibs

			if len(data) < 1 {
				continue
			}

			str := fusePostings(key, data)

			out <- XMLRecord{Index: rec, Ident: key, Text: str}

			runtime.Gosched()
		}
	}

	var wg sync.WaitGroup

	// launch multiple merger goroutines
	for i := 0; i < NumServe(); i++ {
		wg.Add(1)
		go xmlMerger(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all mergers are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

// CreateSplitter distributes adjacent records with the same identifier prefix
func CreateSplitter(mergePath string, zipp bool, inp <-chan XMLRecord) <-chan string {

	if inp == nil {
		return nil
	}

	out := make(chan string, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create splitter channel\n")
		os.Exit(1)
	}

	openSaver := func(mergePath, key string, zipp bool) (*os.File, *bufio.Writer, *gzip.Writer) {

		var (
			fl   *os.File
			wrtr *bufio.Writer
			zpr  *gzip.Writer
			err  error
		)

		sfx := ".mrg"
		if zipp {
			sfx += ".gz"
		}

		fpath := path.Join(mergePath, key+sfx)
		if fpath == "" {
			return nil, nil, nil
		}

		// overwrites and truncates existing file
		fl, err = os.Create(fpath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return nil, nil, nil
		}

		var out io.Writer

		out = fl

		if zipp {

			zpr, err = gzip.NewWriterLevel(fl, gzip.BestSpeed)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
				return nil, nil, nil
			}

			out = zpr
		}

		// create buffered writer layer
		wrtr = bufio.NewWriter(out)
		if wrtr == nil {
			fmt.Fprintf(os.Stderr, "Unable to create bufio.NewWriter\n")
			return nil, nil, nil
		}

		return fl, wrtr, zpr
	}

	closeSaver := func(fl *os.File, wrtr *bufio.Writer, zpr *gzip.Writer) {

		wrtr.Flush()
		if zpr != nil {
			zpr.Close()
		}
		// fl.Sync()

		err := fl.Close()
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		}
	}

	// xmlSplitter distributes adjacent records with the same identifier prefix
	xmlSplitter := func(inp <-chan XMLRecord, out chan<- string) {

		// close channel when all records have been processed
		defer close(out)

		var (
			fl   *os.File
			wrtr *bufio.Writer
			zpr  *gzip.Writer
		)

		currTag := ""
		prevTag := ""

		// remember previous record
		prev := XMLRecord{}

		for curr := range inp {

			// use first few characters of identifier
			currTag = IdentifierKey(curr.Ident)
			if currTag == "" {
				continue
			}

			// then truncate to 2, 3, or 4 character prefix
			if len(currTag) > 2 {
				key := currTag[:2]
				num, ok := TrieLen[key]
				if ok {
					if num > 3 && len(currTag) > 3 {
						key = currTag[:3]
						num, ok = MergLen[key]
						if ok && num > 3 {
							currTag = currTag[:4]
						} else {
							currTag = currTag[:3]
						}
					} else if num > 2 {
						currTag = currTag[:3]
					} else {
						currTag = currTag[:2]
					}
				} else {
					currTag = currTag[:2]
				}
			}

			if fl == nil {
				// open initial file
				fl, wrtr, zpr = openSaver(mergePath, currTag, zipp)
				if wrtr == nil {
					continue
				}

				// send first opening tag and indent
				wrtr.WriteString("<InvDocumentSet>\n  ")
			}

			// compare keys from adjacent term lists
			if prev.Text != "" && prevTag != currTag {

				// after IdentifierKey converts space to underscore,
				// okay that x_ and x0 will be out of alphabetical order

				// send closing tag
				wrtr.WriteString("</InvDocumentSet>\n")

				closeSaver(fl, wrtr, zpr)

				out <- currTag

				// force garbage collection
				runtime.GC()
				debug.FreeOSMemory()

				runtime.Gosched()

				// open next file
				fl, wrtr, zpr = openSaver(mergePath, currTag, zipp)
				if wrtr == nil {
					continue
				}

				// send opening tag and indent
				wrtr.WriteString("<InvDocumentSet>\n  ")
			}

			// send one InvDocument
			str := strings.TrimSpace(curr.Text)

			wrtr.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				wrtr.WriteString("\n")
			}

			// now remember this record
			prev = curr

			prevTag = currTag
		}

		if prev.Text != "" {

			// send last closing tag
			wrtr.WriteString("</InvDocumentSet>\n")

			closeSaver(fl, wrtr, zpr)

			out <- currTag

			// force garbage collection
			runtime.GC()
			debug.FreeOSMemory()

			runtime.Gosched()
		}
	}

	// launch single splitter goroutine
	go xmlSplitter(inp, out)

	return out
}
