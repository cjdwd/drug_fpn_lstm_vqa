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
// File Name:  asn.go
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

// ASN1Converter parses text ASN.1 records into XML objects
func ASN1Converter(inp io.Reader, set, rec string) <-chan string {

	if inp == nil {
		return nil
	}

	tks := make(chan string, chanDepth)
	if tks == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create ASN1 tokenizer channel\n")
		os.Exit(1)
	}

	out := make(chan string, chanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create ASN1 converter channel\n")
		os.Exit(1)
	}

	// tokenizeASN1 sends ASN1 tokens down a channel
	tokenizeASN1 := func(inp io.Reader, tks chan<- string) {

		// close channel when all tokens have been sent
		defer close(tks)

		var buf strings.Builder

		scanr := bufio.NewScanner(inp)

		row := 0
		idx := 0
		line := ""

		sentinel := string(rune(0))

		nextLine := func() string {

			for scanr.Scan() {
				// read line
				line := scanr.Text()
				row++
				if line == "" {
					// ignore blank lines
					continue
				}
				// add sentinel
				line += sentinel
				return line
			}

			// end of data
			return sentinel
		}

		readRestOfAsnString := func() string {

			// continue reading additional lines of string
			for {

				line = nextLine()
				if line == "" || line == sentinel {
					break
				}
				idx = 0

				ch := line[idx]

				for inAsnString[ch] {
					idx++
					ch = line[idx]
				}
				str := line[:idx]
				line = line[idx:]
				idx = 0

				buf.WriteString(str)

				if strings.HasPrefix(line, "\"") {
					// "
					// skip past closing quote
					line = line[1:]
					idx = 0

					// break out of continuation loop
					return line
				}

				// string continues on additional lines
			}

			return line
		}

		readRestOfAsnBits := func() string {

			// continue reading additional lines of string
			for {

				line = nextLine()
				if line == "" || line == sentinel {
					break
				}
				idx = 0

				ch := line[idx]

				for inAsnBits[ch] {
					idx++
					ch = line[idx]
				}
				str := line[:idx]
				line = line[idx:]
				idx = 0

				buf.WriteString(str)

				if strings.HasPrefix(line, "'") {
					// skip past closing apostrophe
					line = line[1:]
					idx = 0

					// break out of continuation loop
					return line
				}

				// string continues on additional lines
			}

			return line
		}

		for {

			line = nextLine()
			if line == "" || line == sentinel {
				break
			}

			for {

				if line == "" || line == sentinel {
					break
				}
				idx = 0

				// trim leading blanks
				ch := line[idx]
				for inBlank[ch] {
					idx++
					ch = line[idx]
				}
				line = line[idx:]
				idx = 0

				if ch == ',' {
					tks <- string(ch)
					line = line[1:]
					continue
				}

				if ch == '{' {
					// start structure
					tks <- string(ch)
					line = line[1:]
					continue
				}

				if ch == '}' {
					// end structure
					tks <- string(ch)
					line = line[1:]
					continue
				}

				if ch == '"' {
					// "
					// start of string
					buf.Reset()

					// skip past opening quote
					line = line[1:]
					idx = 0
					ch = line[idx]

					// read to closing quote or sentinel at end of line
					for inAsnString[ch] {
						idx++
						ch = line[idx]
					}
					str := line[:idx]
					line = line[idx:]
					idx = 0

					buf.WriteString(str)

					if strings.HasPrefix(line, "\"") {
						// "
						// skip past closing quote
						line = line[1:]
						idx = 0

					} else {

						// continue reading additional lines of string
						line = readRestOfAsnString()
						idx = 0
					}

					tmp := buf.String()
					if tmp == "" {
						// encode empty string
						tmp = "\"\""
					}

					tks <- tmp
					buf.Reset()

					continue
				}

				if ch == '\'' {
					// start of bit string
					buf.Reset()

					// skip past opening apostrophe
					line = line[1:]
					idx = 0
					ch = line[idx]

					// read to closing apostrophe or sentinel at end of line
					for inAsnBits[ch] {
						idx++
						ch = line[idx]
					}
					str := line[:idx]
					line = line[idx:]
					idx = 0

					buf.WriteString(str)

					if strings.HasPrefix(line, "'") {
						// skip past closing apostrophe
						line = line[1:]
						idx = 0

					} else {

						// continue reading additional lines of bit string
						line = readRestOfAsnBits()
						idx = 0
					}

					// if apostrophe is at end of line, read next line
					if line == "" || line == sentinel {
						line = nextLine()
					}

					// then skip past trailing hex or binary indicator
					if strings.HasPrefix(line, "H") || strings.HasPrefix(line, "B") {
						line = line[1:]
						idx = 0
					}

					tks <- buf.String()
					buf.Reset()

					continue
				}

				if ch == ':' && strings.HasPrefix(line, "::=") {
					// start of record contents
					tks <- "::="
					line = line[3:]
					idx = 0
					continue
				}

				if ch == '-' && strings.HasPrefix(line, "--") {
					// skip comments
					break
				}
				if ch == ';' {
					// skip comments
					break
				}

				// read token or unquoted numeric value
				idx = 0
				for inAsnTag[ch] {
					idx++
					ch = line[idx]
				}
				tkn := line[:idx]
				line = line[idx:]
				idx = 0

				tks <- tkn
			}
		}
	}

	// convertASN1 sends XML records down a channel
	convertASN1 := func(inp <-chan string, out chan<- string) {

		// close channel when all tokens have been processed
		defer close(out)

		// ensure that XML tags are legal
		fixTag := func(tag string) string {

			if tag == "" {
				return tag
			}

			okay := true
			for _, ch := range tag {
				if !inElement[ch] {
					okay = false
				}
			}
			if okay {
				return tag
			}

			var temp strings.Builder

			// replace illegal characters with underscore
			for _, ch := range tag {
				if inElement[ch] {
					temp.WriteRune(ch)
				} else {
					temp.WriteRune('_')
				}
			}

			return temp.String()
		}

		nextToken := func() string {

			for {
				tkn, ok := <-inp
				if !ok {
					break
				}
				if tkn == "" {
					// ignore blank tokens
					continue
				}

				return tkn
			}

			// end of data
			return ""
		}

		// collects tags until next brace or comma
		var arry []string

		// builds XML output for current record
		var buffer strings.Builder

		// array to speed up indentation
		indentSpaces := []string{
			"",
			"  ",
			"    ",
			"      ",
			"        ",
			"          ",
			"            ",
			"              ",
			"                ",
			"                  ",
		}

		count := 0

		indent := 0
		if set != "" {
			indent = 1
		}

		// indent a specified number of spaces
		doIndent := func(indt int) {
			i := indt
			for i > 9 {
				buffer.WriteString("                    ")
				i -= 10
			}
			if i < 0 {
				return
			}
			buffer.WriteString(indentSpaces[i])
		}

		printOpeningTag := func(tag string) {

			tag = fixTag(tag)
			doIndent(indent)
			indent++
			buffer.WriteString("<")
			buffer.WriteString(tag)
			buffer.WriteString(">\n")
		}

		printClosingTag := func(tag string) {

			tag = fixTag(tag)
			indent--
			doIndent(indent)
			buffer.WriteString("</")
			buffer.WriteString(tag)
			buffer.WriteString(">\n")
		}

		printContent := func(tag, tkn string) {

			if tkn == "\"\"" {
				return
			}
			tag = fixTag(tag)
			doIndent(indent)
			buffer.WriteString("<")
			buffer.WriteString(tag)
			buffer.WriteString(">")
			tkn = strings.TrimSpace(tkn)
			tkn = html.EscapeString(tkn)
			buffer.WriteString(tkn)
			buffer.WriteString("</")
			buffer.WriteString(tag)
			buffer.WriteString(">\n")
		}

		popFromArry := func() (string, string, string) {

			fst, sec, trd := "", "", ""
			switch len(arry) {
			case 1:
				fst = arry[0]
			case 2:
				fst, sec = arry[0], arry[1]
			case 3:
				fst, sec, trd = arry[0], arry[1], arry[2]
			}
			arry = nil
			return fst, sec, trd
		}

		// recursive function definition
		var parseAsnObject func(prnt string, lvl int)

		parseAsnObject = func(prnt string, lvl int) {

			for {
				tkn := nextToken()
				if tkn == "" {
					return
				}

				switch tkn {
				case "{":
					fst, sec, trd := popFromArry()
					tag := fst
					if fst == "" {
						fst = prnt
					}
					if fst == "" {
						return
					}
					printOpeningTag(fst)
					tag = fst + "_E"
					if sec != "" {
						printOpeningTag(sec)
						tag = sec + "_E"
					}
					if trd != "" {
						printOpeningTag(trd)
						tag = trd
					}
					parseAsnObject(tag, lvl+1)
					if trd != "" {
						printClosingTag(trd)
					}
					if sec != "" {
						printClosingTag(sec)
					}
					printClosingTag(fst)
					if lvl == 0 {
						return
					}
				case ",":
					fst, sec, trd := popFromArry()
					if trd != "" {
						printOpeningTag(fst)
						printContent(sec, trd)
						printClosingTag(fst)
					} else if sec != "" {
						printContent(fst, sec)
					} else if fst != "" {
						printContent(prnt, fst)
					}
				case "}":
					fst, sec, trd := popFromArry()
					if fst == "" {
						return
					}
					if trd != "" {
						printOpeningTag(fst)
						printContent(sec, trd)
						printClosingTag(fst)
					} else if sec != "" {
						printContent(fst, sec)
					} else if fst != "" {
						printContent(prnt, fst)
					}
					return
				case "::=":
					fmt.Fprintf(os.Stderr, "\nERROR: Unexpected ::= token found\n")
					os.Exit(1)
				default:
					arry = append(arry, tkn)
				}

				count++
				if count > 1000 {
					count = 0
					txt := buffer.String()
					if txt != "" {
						// send current result through output channel
						out <- txt
					}
					buffer.Reset()
				}
			}
		}

		if set != "" {
			out <- "<" + set + ">"
		}

		// process stream of catenated top-level ASN1 records
		for {
			arry = nil

			top := nextToken()
			if top == "" {
				break
			}

			if rec != "" {
				top = rec
			}

			arry = append(arry, top)

			tkn := nextToken()
			if tkn == "" {
				fmt.Fprintf(os.Stderr, "\nERROR: Incomplete ASN1 starting with '%s'\n", top)
				os.Exit(1)
			}
			if tkn != "::=" {
				fmt.Fprintf(os.Stderr, "\nERROR: ASN1 message missing expected ::= token, found '%s'\n", tkn)
				os.Exit(1)
			}

			parseAsnObject(top, 0)

			txt := buffer.String()
			if txt != "" {
				// send remaining result through output channel
				out <- txt
			}

			buffer.Reset()

			runtime.Gosched()
		}

		if set != "" {
			out <- "</" + set + ">"
		}
	}

	// launch single tokenizer goroutine
	go tokenizeASN1(inp, tks)

	// launch single converter goroutine
	go convertASN1(tks, out)

	return out
}
