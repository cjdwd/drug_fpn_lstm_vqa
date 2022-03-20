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
// File Name:  json.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"encoding/json"
	"fmt"
	"github.com/gedex/inflector"
	"html"
	"io"
	"os"
	"runtime"
	"strconv"
	"strings"
)

// JSONConverter parses JSON stream into XML object stream
func JSONConverter(inp io.Reader, set, rec, nest string) <-chan string {

	if inp == nil {
		return nil
	}

	tks := make(chan string, chanDepth)
	if tks == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create JSON tokenizer channel\n")
		os.Exit(1)
	}

	out := make(chan string, chanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create JSON converter channel\n")
		os.Exit(1)
	}

	tokenizeJSON := func(inp io.Reader, tks chan<- string) {

		// close channel when all tokens have been sent
		defer close(tks)

		// use token decoder from encoding/json package
		dec := json.NewDecoder(inp)
		if dec == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create JSON Decoder\n")
			os.Exit(1)
		}
		dec.UseNumber()

		for {
			t, err := dec.Token()
			if err == io.EOF {
				return
			}
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to read JSON token '%s'\n", err)
				os.Exit(1)
			}

			// type switch performs sequential type assertions until match is found
			switch v := t.(type) {
			case json.Delim:
				// opening or closing braces (for objects) or brackets (for arrays)
				tks <- string(v)
			case string:
				str := v
				if HasAdjacentSpacesOrNewline(str) {
					str = CompressRunsOfSpaces(str)
				}
				tks <- str
			case json.Number:
				tks <- v.String()
			case float64:
				tks <- strconv.FormatFloat(v, 'f', -1, 64)
			case bool:
				if v {
					tks <- "true"
				} else {
					tks <- "false"
				}
			case nil:
				tks <- "null"
			default:
				tks <- t.(string)
			}
		}
	}

	// opt is used for anonymous top-level objects, anon for anonymous top-level arrays
	opt := "opt"
	anon := "anon"
	if rec != "" {
		// override record delimiter
		opt = rec
		anon = rec
	}

	flat := false
	elem := false
	depth := false
	plural := false

	switch nest {
	case "flat":
		flat = true
	case "element", "elem", "_E":
		elem = true
	case "depth", "deep", "level":
		depth = true
	case "plural", "name":
		plural = true
	case "recurse", "recursive", "same":
	}

	// convertJSON sends XML records down a channel
	convertJSON := func(tks <-chan string, out chan<- string) {

		// close channel when all tokens have been processed
		defer close(out)

		// ensure that XML tags are legal (initial digit allowed by xtract for biological data in JSON)
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

		// closure silently places local variable pointer onto inner function call stack
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

		count := 0

		// recursive function definitions
		var parseObject func(tag string)
		var parseArray func(tag, pfx string, lvl int)

		// recursive descent parser uses mutual recursion
		parseValue := func(tag, pfx, tkn string, lvl int) {

			switch tkn {
			case "{":
				parseObject(tag)
				// no break needed, would use fallthrough to explicitly cause program control to flow to the next case
			case "[":
				if flat {
					parseArray(tag, pfx, lvl+1)
				} else if lvl > 0 {
					// nested JSON arrays create recursive XML objects
					doIndent(indent)
					indent++
					tg := tag
					if plural {
						tg = inflector.Pluralize(tag)
					}
					buffer.WriteString("<")
					buffer.WriteString(tg)
					buffer.WriteString(">\n")
					if depth {
						parseArray(pfx+"_"+strconv.Itoa(lvl), tag, lvl+1)
					} else if elem {
						sfx := ""
						for i := 0; i < lvl; i++ {
							sfx += "_E"
						}
						parseArray(pfx+sfx, tag, lvl+1)
					} else {
						parseArray(tag, pfx, lvl+1)
					}
					indent--
					doIndent(indent)
					buffer.WriteString("</")
					buffer.WriteString(tg)
					buffer.WriteString(">\n")
				} else {
					parseArray(tag, pfx, lvl+1)
				}
			case "}", "]":
				// should not get here, decoder tracks nesting of braces and brackets
			case "":
				// empty value string generates self-closing object
				doIndent(indent)
				buffer.WriteString("<")
				buffer.WriteString(tag)
				buffer.WriteString("/>\n")
			default:
				// write object and contents to string builder
				doIndent(indent)
				tkn = strings.TrimSpace(tkn)
				tkn = html.EscapeString(tkn)
				buffer.WriteString("<")
				buffer.WriteString(tag)
				buffer.WriteString(">")
				buffer.WriteString(tkn)
				buffer.WriteString("</")
				buffer.WriteString(tag)
				buffer.WriteString(">\n")
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

		parseObject = func(tag string) {

			doIndent(indent)
			indent++
			buffer.WriteString("<")
			buffer.WriteString(tag)
			buffer.WriteString(">\n")

			for {
				// shadowing tag variable inside for loop does not step on value of tag argument in outer scope
				tag, ok := <-tks
				if !ok {
					break
				}

				if tag == "}" || tag == "]" {
					break
				}

				tag = fixTag(tag)

				tkn, ok := <-tks
				if !ok {
					break
				}

				if tkn == "}" || tkn == "]" {
					break
				}

				parseValue(tag, tag, tkn, 0)
			}

			indent--
			doIndent(indent)
			buffer.WriteString("</")
			buffer.WriteString(tag)
			buffer.WriteString(">\n")
		}

		parseArray = func(tag, pfx string, lvl int) {

			for {
				tkn, ok := <-tks
				if !ok {
					break
				}

				if tkn == "}" || tkn == "]" {
					break
				}

				parseValue(tag, pfx, tkn, lvl)
			}
		}

		if set != "" {
			out <- "<" + set + ">"
		}

		// process stream of catenated top-level JSON objects or arrays
		for {
			tkn, ok := <-tks
			if !ok {
				break
			}
			if tkn == "{" {
				parseObject(opt)
			} else if tkn == "[" {
				parseArray(anon, anon, 0)
			} else {
				break
			}

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
	go tokenizeJSON(inp, tks)

	// launch single converter goroutine
	go convertJSON(tks, out)

	return out
}
