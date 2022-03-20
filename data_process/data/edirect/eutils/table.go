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
// File Name:  table.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"bufio"
	"fmt"
	"html"
	"io"
	"os"
	"runtime"
	"strings"
)

// TableConverter parses tab-delimited or comma-separated values files into XML object stream
func TableConverter(inp io.Reader, delim, set, rec string, skip int, header, lower, upper, indent bool, fields []string) <-chan string {

	if inp == nil {
		return nil
	}

	head := ""
	tail := ""

	hd := ""
	tl := ""

	if set != "" && set != "-" {
		head = "<" + set + ">"
		tail = "</" + set + ">"
	}

	if rec != "" && rec != "-" {
		hd = "<" + rec + ">"
		tl = "</" + rec + ">"
	}

	numFlds := len(fields)

	if numFlds < 1 && !header {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient arguments for table converter\n")
		os.Exit(1)
	}

	out := make(chan string, chanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create table converter channel\n")
		os.Exit(1)
	}

	convertTable := func(inp io.Reader, out chan<- string) {

		// close channel when all records have been sent
		defer close(out)

		okay := false
		row := 0

		var buffer strings.Builder

		scanr := bufio.NewScanner(inp)

		if head != "" {
			buffer.WriteString(head)
			buffer.WriteString("\n")
		}

		if header {

			// uses fields from first row for column names
			for scanr.Scan() {

				line := scanr.Text()

				row++

				if skip > 0 {
					skip--
					continue
				}

				cols := strings.Split(line, delim)

				for _, str := range cols {
					fields = append(fields, str)
					numFlds++
				}
				break
			}

			if numFlds < 1 {
				fmt.Fprintf(os.Stderr, "\nERROR: Line with column names not found\n")
				os.Exit(1)
			}
		}

		for scanr.Scan() {

			line := scanr.Text()

			row++

			if skip > 0 {
				skip--
				continue
			}

			cols := strings.Split(line, delim)

			if len(cols) != numFlds {
				fmt.Fprintf(os.Stderr, "Mismatched columns in row %d - '%s'\n", row, line)
				continue
			}

			if hd != "" {
				if indent {
					buffer.WriteString("  ")
				}
				buffer.WriteString(hd)
				buffer.WriteString("\n")
			}

			for i, fld := range fields {
				val := cols[i]
				if lower {
					val = strings.ToLower(val)
				}
				if upper {
					val = strings.ToUpper(val)
				}
				if fld[0] == '*' {
					fld = fld[1:]
				} else {
					val = html.EscapeString(val)
				}
				val = strings.TrimSpace(val)
				if indent {
					buffer.WriteString("    ")
				}
				buffer.WriteString("<")
				buffer.WriteString(fld)
				buffer.WriteString(">")
				buffer.WriteString(val)
				buffer.WriteString("</")
				buffer.WriteString(fld)
				buffer.WriteString(">")
				buffer.WriteString("\n")
			}

			if tl != "" {
				if indent {
					buffer.WriteString("  ")
				}
				buffer.WriteString(tl)
				buffer.WriteString("\n")
			}

			okay = true
		}

		if tail != "" {
			buffer.WriteString(tail)
			buffer.WriteString("\n")
		}

		if okay {
			txt := buffer.String()
			if txt != "" {
				// send remaining result through output channel
				out <- txt
			}
		}

		buffer.Reset()

		runtime.Gosched()
	}

	go convertTable(inp, out)

	return out
}
