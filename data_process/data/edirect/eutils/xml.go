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
// File Name:  xml.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"bytes"
	"container/heap"
	"fmt"
	"io"
	"os"
)

// XMLBlock is a string that begins with a left angle bracket and is trimmed back to
// end with a right angle bracket. The excluded characters are saved and prepended
// to the next buffer. Providing complete object tags simplifies subsequent parsing.
type XMLBlock string

// CreateXMLStreamer reads XML input into a channel of trimmed strings that are
// then split by PartitionPattern into individual records (which can be processed
// concurrently), or parsed directly into a channel of tokens by CreateTokenizer.
func CreateXMLStreamer(in io.Reader) <-chan XMLBlock {

	if in == nil {
		return nil
	}

	out := make(chan XMLBlock, chanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML block reader channel\n")
		os.Exit(1)
	}

	// xmlReader sends trimmed XML blocks through the output channel.
	xmlReader := func(in io.Reader, out chan<- XMLBlock) {

		// close channel when all blocks have been processed
		defer close(out)

		// 65536 appears to be the maximum number of characters presented to io.Reader
		// when input is piped from stdin. Increasing the buffer size when input is from
		// a file does not improve program performance. An additional 16384 bytes are
		// reserved for copying the previous remainder to the beginning of the buffer
		// before the next read.
		const XMLBUFSIZE = 65536 + 16384

		buffer := make([]byte, XMLBUFSIZE)
		remainder := ""
		position := int64(0)
		delta := 0
		isClosed := false

		// htmlBehind is used in strict mode to trim back further when a lower-case tag
		// is encountered. This may be a formatting decoration, such as <i> or </i> for
		// italics. Processing HTML, which may have embedded mixed content, requires use
		// of mixed mode.
		htmlBehind := func(bufr []byte, pos, txtlen int) bool {

			for pos >= 0 {
				if bufr[pos] == '<' {
					// detect lower-case markup tags, or DispFormula in PubMed
					return HTMLAhead(string(bufr), pos, txtlen) != 0
				}
				pos--
			}

			return false
		}

		// nextBuffer reads one buffer, trims back to the right-most > character, and
		// retains the remainder for prepending in the next call. It also signals if
		// there was no > character, resulting in subsequent calls to nextBuffer to
		// continue reading a large content string.
		nextBuffer := func() ([]byte, bool, bool) {

			if isClosed {
				return nil, false, true
			}

			// prepend previous remainder to beginning of buffer
			m := copy(buffer, remainder)
			remainder = ""
			if m > 16384 {
				// previous remainder is larger than reserved section,
				// write and signal the need to continue reading.
				return buffer[:m], true, false
			}

			// read next block, append behind copied remainder from previous read
			n, err := in.Read(buffer[m:])
			// with data piped through stdin, read function may not always return the
			// same number of bytes each time
			if err != nil {
				if err != io.EOF {
					// real error.
					fmt.Fprintf(os.Stderr, "\n%sERROR: %s%s\n", RED, err.Error(), INIT)
					// ignore bytes - non-conforming implementations of io.Reader may
					// return mangled data on non-EOF errors
					isClosed = true
					return nil, false, true
				}
				// end of file.
				isClosed = true
				if n == 0 {
					// if EOF and no more data, do not send final remainder (not terminated
					// by right angle bracket that is used as a sentinel)
					return nil, false, true
				}
			}
			if n < 0 {
				// reality check - non-conforming implementations of io.Reader may return -1
				fmt.Fprintf(os.Stderr, "\n%sERROR: io.Reader returned negative count %d%s\n", RED, n, INIT)
				// treat as n == 0 in order to update file offset and avoid losing previous remainder
				n = 0
			}

			// keep track of file offset
			position += int64(delta)
			delta = n

			// slice of actual characters read
			bufr := buffer[:n+m]

			// Look for last > character. It is safe to back up on UTF-8 rune array when looking
			// for a 7-bit ASCII character.
			pos := -1
			for pos = len(bufr) - 1; pos >= 0; pos-- {
				if bufr[pos] == '>' {
					if doStrict {
						// optionally skip backwards past embedded i, b, u, sub, and sup
						// HTML open, close, and empty tags, and MathML instructions
						if htmlBehind(bufr, pos, len(bufr)) {
							continue
						}
					}
					// found end of XML tag, break
					break
				}
			}

			// trim back to last > character, save remainder for next buffer
			if pos > -1 {
				pos++
				remainder = string(bufr[pos:])
				return bufr[:pos], false, false
			}

			// no > found, signal need to continue reading long content
			return bufr[:], true, false
		}

		// nextBlock reads buffer, concatenates if necessary to place long element content
		// into a single string. All result strings end in > character that is used as a
		// sentinel in subsequent code.
		nextBlock := func() string {

			// read next buffer
			line, cont, closed := nextBuffer()

			if closed {
				// no sentinel in remainder at end of file
				return ""
			}

			if cont {
				// current line does not end with > character
				var buff bytes.Buffer

				// keep reading long content blocks
				for {
					if len(line) > 0 {
						buff.Write(line)
					}
					if !cont {
						// last buffer ended with sentinel
						break
					}
					line, cont, closed = nextBuffer()
					if closed {
						// no sentinel in multi-block buffer at end of file
						return ""
					}
				}

				// concatenate blocks
				return buff.String()
			}

			return string(line)
		}

		// read XML and send blocks through channel
		for {
			str := nextBlock()

			// trimming spaces here would throw off line tracking

			// optionally compress/cleanup tags/attributes and contents
			if doCleanup {
				if HasBadSpace(str) {
					str = CleanupBadSpaces(str)
				}
				if HasAdjacentSpaces(str) {
					str = CompressRunsOfSpaces(str)
				}
			}

			out <- XMLBlock(str)

			// bail after sending empty string sentinel
			if str == "" {
				return
			}
		}
	}

	// launch single block reader goroutine
	go xmlReader(in, out)

	return out
}

// XMLRecord wraps a numbered XML record or the results of data extraction on
// that record. The Index field stores the record's original position in the
// input stream. The Data field is used for binary compressed PubmedArticle XML.
type XMLRecord struct {
	Index int
	Ident string
	Text  string
	Data  []byte
}

// CreateXMLProducer partitions an XML set and sends records down a channel.
// After processing asynchronously in multiple concurrent go routines, the
// original order can be restored by passage through the XMLUnshuffler.
func CreateXMLProducer(pat, star string, turbo bool, rdr <-chan XMLBlock) <-chan XMLRecord {

	if rdr == nil {
		return nil
	}

	out := make(chan XMLRecord, chanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML producer channel\n")
		os.Exit(1)
	}

	// xmlProducer sends partitioned XML strings through channel.
	xmlProducer := func(pat, star string, turbo bool, rdr <-chan XMLBlock, out chan<- XMLRecord) {

		// close channel when all records have been processed
		defer close(out)

		rec := 0

		// partition all input by pattern and send XML substring to available consumer through channel
		PartitionPattern(pat, star, turbo, rdr,
			func(str string) {
				rec++
				out <- XMLRecord{rec, "", str, nil}
			})
	}

	// launch single producer goroutine
	go xmlProducer(pat, star, turbo, rdr, out)

	return out
}

// xmlRecordHeap collects asynchronous processing results for presentation in the original order.
type xmlRecordHeap []XMLRecord

// methods that satisfy heap.Interface
func (h xmlRecordHeap) Len() int {
	return len(h)
}
func (h xmlRecordHeap) Less(i, j int) bool {
	return h[i].Index < h[j].Index
}
func (h xmlRecordHeap) Swap(i, j int) {
	h[i], h[j] = h[j], h[i]
}
func (h *xmlRecordHeap) Push(x interface{}) {
	*h = append(*h, x.(XMLRecord))
}
func (h *xmlRecordHeap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

// CreateXMLUnshuffler passes the output of multiple concurrent processors to
// a heap, which releases results in the same order as the original records.
func CreateXMLUnshuffler(inp <-chan XMLRecord) <-chan XMLRecord {

	if inp == nil {
		return nil
	}

	out := make(chan XMLRecord, chanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML unshuffler channel\n")
		os.Exit(1)
	}

	// xmlUnshuffler restores original order with heap.
	xmlUnshuffler := func(inp <-chan XMLRecord, out chan<- XMLRecord) {

		// close channel when all records have been processed
		defer close(out)

		// initialize empty heap
		hp := &xmlRecordHeap{}
		heap.Init(hp)

		// index of next desired result
		next := 1

		delay := 0

		for ext := range inp {

			// push result onto heap
			heap.Push(hp, ext)

			// Read several values before checking to see if next record to print has been processed.
			// The default heapSize value has been tuned by experiment for maximum performance.
			if delay < heapSize {
				delay++
				continue
			}

			delay = 0

			for hp.Len() > 0 {

				// remove lowest item from heap, use interface type assertion
				curr := heap.Pop(hp).(XMLRecord)

				if curr.Index > next {

					// record should be printed later, push back onto heap
					heap.Push(hp, curr)
					// and go back to waiting on input channel
					break
				}

				// send even if empty to get all record counts for reordering
				out <- XMLRecord{curr.Index, curr.Ident, curr.Text, curr.Data}

				// prevent ambiguous -limit filter from clogging heap (deprecated)
				if curr.Index == next {
					// increment index for next expected match
					next++
				}

				// continue to check heap to see if next result is already available
			}
		}

		// flush remainder of heap to output
		for hp.Len() > 0 {
			curr := heap.Pop(hp).(XMLRecord)

			out <- XMLRecord{curr.Index, curr.Ident, curr.Text, curr.Data}
		}
	}

	// launch single unshuffler goroutine
	go xmlUnshuffler(inp, out)

	return out
}
