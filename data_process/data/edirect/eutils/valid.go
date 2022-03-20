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
// File Name:  valid.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"fmt"
	"os"
	"strings"
)

// ValidateXML checks for well-formed XML
func ValidateXML(rdr <-chan XMLBlock, fnd string, html bool) int {

	if rdr == nil {
		return 0
	}

	countLines = true

	tknq := CreateTokenizer(rdr)

	if tknq == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create validator tokenizer\n")
		os.Exit(1)
	}

	find := ParseIndex(fnd)

	// report unexpectedly large maximum nesting depth (undocumented)
	maxDepth := 0
	depthLine := 0
	depthID := ""
	maxLine := 0

	// warn if HTML tags are not well-formed
	unbalancedHTML := func(text string) bool {

		var arry []string

		idx := 0
		txtlen := len(text)

		inTag := false
		start := 0
		stop := -1
		var last byte

		for idx < txtlen {
			ch := text[idx]
			if ch == '<' {
				if inTag {
					return true
				}
				inTag = true
				start = idx
				stop = -1
			} else if ch == ' ' || ch == '\n' {
				// space before attributes marks end of tag
				if stop < 0 {
					stop = idx
				}
			} else if ch == '>' {
				if !inTag {
					return true
				}
				inTag = false
				if stop < 0 {
					stop = idx
				}
				curr := text[start+1 : stop]
				if strings.HasPrefix(curr, "/") {
					curr = curr[1:]
					if len(arry) < 1 {
						return true
					}
					prev := arry[len(arry)-1]
					if curr != prev {
						return true
					}
					arry = arry[:len(arry)-1]
				} else if last == '/' {
					// ignore self-closing tag
				} else {
					arry = append(arry, curr)
				}
			}
			last = ch
			idx++
		}

		if inTag {
			return true
		}

		if len(arry) > 0 {
			return true
		}

		return false
	}

	// warn if HTML tags are encoded
	encodedHTML := func(str string) bool {

		lookAhead := func(txt string, to int) string {

			mx := len(txt)
			if to > mx {
				to = mx
			}
			pos := strings.Index(txt[:to], "gt;")
			if pos > 0 {
				to = pos + 3
			}
			return txt[:to]
		}

		for i, ch := range str {
			if ch == '<' {
				continue
			} else if ch != '&' {
				continue
			} else if strings.HasPrefix(str[i:], "&lt;") {
				sub := lookAhead(str[i:], 14)
				_, ok := HTMLRepair(sub)
				if ok {
					return true
				}
			} else if strings.HasPrefix(str[i:], "&amp;lt;") {
				sub := lookAhead(str[i:], 22)
				_, ok := HTMLRepair(sub)
				if ok {
					return true
				}
			} else if strings.HasPrefix(str[i:], "&amp;amp;") {
				return true
			}
		}

		return false
	}

	currID := ""

	// verifyLevel recursive definition
	var verifyLevel func(string, string, int)

	// verify integrity of XML object nesting (well-formed)
	verifyLevel = func(parent, prev string, level int) {

		status := START
		for tkn := range tknq {

			tag := tkn.Tag
			name := tkn.Name
			line := tkn.Line
			maxLine = line

			if level > maxDepth {
				maxDepth = level
				depthLine = line
				depthID = currID
			}

			switch tag {
			case STARTTAG:
				if status == CHAR && !doMixed {
					fmt.Fprintf(os.Stdout, "%s%8d\t<%s> not expected after contents\n", currID, line, name)
				}
				verifyLevel(name, parent, level+1)
				// returns here after recursion
				status = STOP
			case SELFTAG:
				status = OTHER
			case STOPTAG:
				if parent != name && parent != "" {
					fmt.Fprintf(os.Stdout, "%s%8d\tExpected </%s>, found </%s>\n", currID, line, parent, name)
				}
				if level < 1 {
					fmt.Fprintf(os.Stdout, "%s%8d\tUnexpected </%s> at end of XML\n", currID, line, name)
				}
				// break recursion
				return
			case CONTENTTAG:
				// check for content index match
				if find != nil && find.Index != "" {
					if parent == find.Match || find.Match == "" {
						if find.Parent == "" || prev == find.Parent {
							currID = name + "\t"
						}
					}
				}
				if status != START && !doMixed {
					fmt.Fprintf(os.Stdout, "%s%8d\tContents not expected before </%s>\n", currID, line, parent)
				}
				if allowEmbed {
					if unbalancedHTML(name) {
						fmt.Fprintf(os.Stdout, "%s%8d\tUnbalanced mixed-content tags in <%s>\n", currID, line, parent)
					}
					if html && encodedHTML(name) {
						fmt.Fprintf(os.Stdout, "%s%8d\tEncoded mixed-content markup in <%s>\n", currID, line, parent)
					}
				}
				status = CHAR
			case CDATATAG, COMMENTTAG:
				status = OTHER
			case DOCTYPETAG:
			case NOTAG:
			case ISCLOSED:
				if level > 0 {
					fmt.Fprintf(os.Stdout, "%s%8d\tUnexpected end of data\n", currID, line)
				}
				return
			default:
				status = OTHER
			}
		}
	}

	verifyLevel("", "", 0)

	if maxDepth > 25 {
		fmt.Fprintf(os.Stdout, "%s%8d\tMaximum nesting, %d levels\n", depthID, depthLine, maxDepth)
	}

	return maxLine
}
