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
// File Name:  fasta.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"fmt"
	"io"
	"os"
	"strings"
)

// FASTARecord contains parsed data from FASTA format
type FASTARecord struct {
	SeqID    string
	Title    string
	Length   int
	Sequence string
}

// FASTAConverter partitions a FASTA set and sends records down a channel.
func FASTAConverter(inp io.Reader, caseSensitive bool) <-chan FASTARecord {

	if inp == nil {
		return nil
	}

	tks := make(chan string, chanDepth)
	if tks == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create FASTA tokenizer channel\n")
		os.Exit(1)
	}

	out := make(chan FASTARecord, chanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create FASTA streamer channel\n")
		os.Exit(1)
	}

	// fastaTokenizer splits FASTA input stream into tokens
	fastaTokenizer := func(inp io.Reader, tks chan<- string) {

		// close channel when all tokens have been sent
		defer close(tks)

		const FASTABUFSIZE = 65536

		buffer := make([]byte, FASTABUFSIZE)
		isClosed := false

		nextBuffer := func() string {

			if isClosed {
				return ""
			}

			n, err := inp.Read(buffer[:])

			if err != nil {
				if err != io.EOF {
					// real error.
					fmt.Fprintf(os.Stderr, "\n%sERROR: %s%s\n", RED, err.Error(), INIT)
					// ignore bytes - non-conforming implementations of io.Reader may
					// return mangled data on non-EOF errors
					isClosed = true
					return ""
				}
				// end of file.
				isClosed = true
				if n == 0 {
					// if EOF and no more data
					return ""
				}
			}
			if n < 0 {
				// reality check - non-conforming implementations of io.Reader may return -1
				fmt.Fprintf(os.Stderr, "\n%sERROR: io.Reader returned negative count %d%s\n", RED, n, INIT)
				// treat as n == 0
				return ""
			}

			// slice of actual characters read
			bufr := buffer[:n]

			return string(bufr[:])
		}

		line := ""

		for {

			if line == "" {
				line = nextBuffer()
				if line == "" {
					break
				}
			}

			// look for start of FASTA defline
			pos := strings.Index(line, ">")

			if pos < 0 {
				// no angle bracket, send sequence buffer
				tks <- line
				line = ""
				continue
			}

			if pos > 0 {
				// send sequence buffer up to angle bracket
				str := line[:pos]
				tks <- str
				line = line[pos:]
				continue
			}

			// look for end of FASTA defline
			pos = strings.Index(line, "\n")

			if pos > 0 {
				// send full defline within buffer
				str := line[:pos]
				tks <- str
				line = line[pos+1:]
				continue
			}

			// defline continues into next buffer
			defln := line

			for {

				// read next buffer
				line = nextBuffer()
				if line == "" {
					// file ends in defline
					tks <- defln
					break
				}

				pos = strings.Index(line, "\n")

				if pos < 0 {
					// add full buffer to defline
					defln += line
					continue
				}

				// send constructed defline
				defln += line[:pos]
				tks <- defln
				line = line[pos+1:]
				break
			}
		}
	}

	// fastaStreamer sends FASTA records down a channel
	fastaStreamer := func(tks <-chan string, out chan<- FASTARecord) {

		// close channel when all records have been processed
		defer close(out)

		const FASTABUFSIZE = 65536

		seqid := ""
		title := ""

		var fasta []string

		sendFasta := func() {

			seq := strings.Join(fasta, "")
			seqlen := len(seq)

			if seqlen > 0 {
				out <- FASTARecord{SeqID: seqid, Title: title, Length: seqlen, Sequence: seq[:]}
			}

			seqid = ""
			title = ""
			// reset sequence accumulator
			fasta = nil
		}

		for line := range tks {

			if strings.HasPrefix(line, ">") {

				// send current record
				sendFasta()

				// process next defline
				line = line[1:]
				seqid, title = SplitInTwoLeft(line, " ")

				continue
			}

			if !caseSensitive {
				// optionally convert FASTA letters to upper case
				line = strings.ToUpper(line)
			}

			// leave only letters, asterisk, or hyphen
			line = strings.Map(func(c rune) rune {
				if c >= 'A' && c <= 'Z' {
					return c
				}
				if c >= 'a' && c <= 'z' {
					return c
				}
				if c == '*' || c == '-' {
					return c
				}
				return -1
			}, line)

			// append current line
			fasta = append(fasta, line)
		}

		// send final record
		sendFasta()
	}

	// launch single fasta tokenizer goroutine
	go fastaTokenizer(inp, tks)

	// launch single fasta streamer goroutine
	go fastaStreamer(tks, out)

	return out
}
