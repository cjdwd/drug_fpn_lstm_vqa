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
// File Name:  format.go
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

// FormatArgs contains XML format customization arguments
type FormatArgs struct {
	Format  string
	XML     string
	Doctype string
	Unicode string
	Script  string
	Mathml  string
	Combine bool
	Self    bool
	Comment bool
	Cdata   bool
}

// xmlFormatter reformats a record string or a stream of XML tokens
func xmlFormatter(rcrd, prnt string, inp <-chan XMLToken, offset int, doXML bool, args FormatArgs) <-chan string {

	if rcrd == "" && inp == nil {
		return nil
	}

	out := make(chan string, chanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "Unable to create formatter channel\n")
		os.Exit(1)
	}

	compRecrd := false
	flushLeft := false
	wrapAttrs := false
	ret := "\n"

	switch args.Format {
	case "compact", "compacted", "compress", "compressed", "terse", "*":
		// compress to one record per line
		compRecrd = true
		ret = ""
	case "flush", "flushed", "left":
		// suppress line indentation
		flushLeft = true
	case "expand", "expanded", "extend", "extended", "verbose", "@":
		// each attribute on its own line
		wrapAttrs = true
	case "indent", "indented", "normal", "default", "":
		// default behavior
	default:
		fmt.Fprintf(os.Stderr, "Unrecognized format '%s'\n", args.Format)
		os.Exit(1)
	}

	fuseTopSets := args.Combine
	keepSelfClosing := args.Self

	doComment := args.Comment
	doCdata := args.Cdata

	// formatXML goroutine processes one record
	formatXML := func(rcrd, prnt string, inp <-chan XMLToken, offset int, doXML bool, out chan<- string) {

		// close channel when all chunks have been sent
		defer close(out)

		xml := args.XML

		doctype := ""
		customDoctype := false

		if args.Doctype != "" {
			customDoctype = true
			if !strings.HasPrefix(args.Doctype, "-") {
				doctype = args.Doctype
			}
		}

		var buffer strings.Builder

		count := 0

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

		indent := offset

		// indent a specified number of spaces
		doIndent := func(indt int) {
			if compRecrd || flushLeft {
				return
			}
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

		badAttributeSpacing := func(attr string) bool {

			if len(attr) < 2 {
				return false
			}

			var prev rune

			for _, ch := range attr {
				if ch == '=' && prev == ' ' {
					return true
				}
				if ch == ' ' && prev == '=' {
					return true
				}
				prev = ch
			}

			return false
		}

		// print attributes
		printAttributes := func(attr string) {

			if attr == "" {
				return
			}
			attr = strings.TrimSpace(attr)
			attr = CompressRunsOfSpaces(attr)

			if deAccent {
				if IsNotASCII(attr) {
					attr = DoAccentTransform(attr)
				}
			}
			if doASCII {
				if IsNotASCII(attr) {
					attr = UnicodeToASCII(attr)
				}
			}

			if wrapAttrs {

				start := 0
				idx := 0
				inQuote := false

				attlen := len(attr)

				for idx < attlen {
					ch := attr[idx]
					if ch == '=' && !inQuote {
						inQuote = true
						str := strings.TrimSpace(attr[start:idx])
						buffer.WriteString("\n")
						doIndent(indent)
						buffer.WriteString(" ")
						buffer.WriteString(str)
						// skip past equal sign
						idx++
						ch = attr[idx]
						if ch != '"' && ch != '\'' {
							// "
							// skip past unexpected blanks
							for inBlank[ch] {
								idx++
								ch = attr[idx]
							}
						}
						// skip past leading double quote
						idx++
						start = idx
					} else if ch == '"' || ch == '\'' {
						// "
						inQuote = !inQuote
						str := strings.TrimSpace(attr[start:idx])
						buffer.WriteString("=\"")
						buffer.WriteString(str)
						buffer.WriteString("\"")
						// skip past trailing double quote and (possible) space
						idx += 2
						start = idx
					} else {
						idx++
					}
				}

				buffer.WriteString("\n")
				doIndent(indent)

			} else if badAttributeSpacing(attr) {

				buffer.WriteString(" ")

				start := 0
				idx := 0

				attlen := len(attr)

				for idx < attlen {
					ch := attr[idx]
					if ch == '=' {
						str := strings.TrimSpace(attr[start:idx])
						buffer.WriteString(str)
						// skip past equal sign
						idx++
						ch = attr[idx]
						if ch != '"' && ch != '\'' {
							// "
							// skip past unexpected blanks
							for inBlank[ch] {
								idx++
								ch = attr[idx]
							}
						}
						// skip past leading double quote
						idx++
						start = idx
					} else if ch == '"' || ch == '\'' {
						// "
						str := strings.TrimSpace(attr[start:idx])
						buffer.WriteString("=\"")
						buffer.WriteString(str)
						buffer.WriteString("\"")
						// skip past trailing double quote and (possible) space
						idx += 2
						start = idx
					} else {
						idx++
					}
				}

			} else {

				buffer.WriteString(" ")
				buffer.WriteString(attr)
			}
		}

		parent := ""
		pfx := ""
		skip := 0
		okIndent := true

		printXMLAndDoctype := func(xml, doctype, parent string) {

			// check for xml line explicitly set in argument
			if xml != "" {
				xml = strings.TrimSpace(xml)
				if strings.HasPrefix(xml, "<") {
					xml = xml[1:]
				}
				if strings.HasPrefix(xml, "?") {
					xml = xml[1:]
				}
				if strings.HasPrefix(xml, "xml") {
					xml = xml[3:]
				}
				if strings.HasPrefix(xml, " ") {
					xml = xml[1:]
				}
				if strings.HasSuffix(xml, "?>") {
					xlen := len(xml)
					xml = xml[:xlen-2]
				}
				xml = strings.TrimSpace(xml)

				buffer.WriteString("<?xml ")
				buffer.WriteString(xml)
				buffer.WriteString(" ?>")
			} else {
				buffer.WriteString("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>")
			}

			buffer.WriteString("\n")

			// check for doctype taken from XML file or explicitly set in argument
			if doctype != "" {
				doctype = strings.TrimSpace(doctype)
				if strings.HasPrefix(doctype, "<") {
					doctype = doctype[1:]
				}
				if strings.HasPrefix(doctype, "!") {
					doctype = doctype[1:]
				}
				if strings.HasPrefix(doctype, "DOCTYPE") {
					doctype = doctype[7:]
				}
				if strings.HasPrefix(doctype, " ") {
					doctype = doctype[1:]
				}
				doctype = strings.TrimSuffix(doctype, ">")
				doctype = strings.TrimSpace(doctype)

				buffer.WriteString("<!DOCTYPE ")
				buffer.WriteString(doctype)
				buffer.WriteString(">")
			} else {
				buffer.WriteString("<!DOCTYPE ")
				buffer.WriteString(parent)
				buffer.WriteString(">")
			}

			buffer.WriteString("\n")
		}

		cleanupMixed := func(txt string) string {

			txt = DoTrimFlankingHTML(txt)
			if HasBadSpace(txt) {
				txt = CleanupBadSpaces(txt)
			}
			if HasAdjacentSpaces(txt) {
				txt = CompressRunsOfSpaces(txt)
			}
			if NeedsTightening(txt) {
				txt = TightenParentheses(txt)
			}

			return txt
		}

		cleanToken := func(tkn XMLToken, nxtTag int, nxtName, nxtAttr string, lastContent bool) {

			if skip > 0 {
				skip--
				return
			}

			name := tkn.Name

			switch tkn.Tag {
			case STARTTAG:
				// detect first start tag, print xml and doctype parent
				if indent == 0 && parent == "" {
					parent = name
					if doXML {
						printXMLAndDoctype(xml, doctype, parent)
					}
					// do not fuse <opt> or <anon> top-level objects (converted from JSON)
					if parent == "opt" || parent == "anon" {
						fuseTopSets = false
					}
				}
				// convert start-stop to self-closing tag if attributes are present, otherwise skip
				if nxtTag == STOPTAG && nxtName == name {
					if tkn.Attr != "" || keepSelfClosing {
						buffer.WriteString(pfx)
						doIndent(indent)
						buffer.WriteString("<")
						buffer.WriteString(name)
						printAttributes(tkn.Attr)
						buffer.WriteString("/>")
						pfx = ret
						okIndent = true
					}
					skip++
					return
				}
				buffer.WriteString(pfx)
				if doMixed {
					if !lastContent {
						doIndent(indent)
					}
				} else {
					doIndent(indent)
				}
				indent++
				buffer.WriteString("<")
				buffer.WriteString(name)
				printAttributes(tkn.Attr)
				buffer.WriteString(">")
				pfx = ret
				okIndent = true
				if compRecrd && indent == 1 {
					buffer.WriteString("\n")
				}
			case SELFTAG:
				if tkn.Attr != "" || keepSelfClosing {
					buffer.WriteString(pfx)
					doIndent(indent)
					buffer.WriteString("<")
					buffer.WriteString(name)
					printAttributes(tkn.Attr)
					buffer.WriteString("/>")
					pfx = ret
					okIndent = true
				}
			case STOPTAG:
				// skip internal copies of top-level </parent><parent> tags (to fuse multiple chunks returned by efetch)
				// do not skip if attributes are present, unless pattern ends in "Set" (e.g., <DocumentSummarySet status="OK">)
				if nxtTag == STARTTAG && indent == 1 && fuseTopSets && nxtName == parent {
					if nxtAttr == "" || (strings.HasSuffix(parent, "Set") && len(parent) > 3) {
						skip++
						return
					}
				}
				buffer.WriteString(pfx)
				indent--
				if okIndent {
					doIndent(indent)
				}
				buffer.WriteString("</")
				buffer.WriteString(name)
				buffer.WriteString(">")
				if doMixed && nxtTag == CONTENTTAG && nxtName != "." {
					buffer.WriteString(" ")
				}
				pfx = ret
				okIndent = true
				if compRecrd && indent < 2 {
					buffer.WriteString("\n")
				}
			case CONTENTTAG:
				if nxtTag == STARTTAG || nxtTag == SELFTAG {
					if doStrict {
						fmt.Fprintf(os.Stderr, "%s ERROR: %s UNRECOGNIZED MIXED CONTENT <%s> IN <%s>%s\n", INVT, LOUD, nxtName, name, INIT)
					} else if !doMixed {
						fmt.Fprintf(os.Stderr, "%s ERROR: %s UNEXPECTED MIXED CONTENT <%s> IN <%s>%s\n", INVT, LOUD, nxtName, name, INIT)
					}
				}
				if len(name) > 0 && IsNotJustWhitespace(name) {
					// support for all content processing flags
					if doStrict || doMixed || doCompress {
						ctype := tkn.Cont
						name = CleanupContents(name, (ctype&ASCII) != 0, (ctype&AMPER) != 0, (ctype&MIXED) != 0)
					}
					if doMixed {
						name = cleanupMixed(name)
					}
					buffer.WriteString(name)
				}
				if (doStrict || doMixed) && nxtTag == STARTTAG {
					buffer.WriteString(" ")
				}
				pfx = ""
				okIndent = false
			case CDATATAG:
				if doCdata {
					buffer.WriteString(pfx)
					doIndent(indent)
					buffer.WriteString("<![CDATA[")
					buffer.WriteString(name)
					buffer.WriteString("]]>")
					pfx = ret
					okIndent = true
				}
			case COMMENTTAG:
				if doComment {
					buffer.WriteString(pfx)
					doIndent(indent)
					buffer.WriteString("<!--")
					buffer.WriteString(name)
					buffer.WriteString("-->")
					pfx = ret
					okIndent = true
				}
			case DOCTYPETAG:
				if customDoctype && doctype == "" {
					doctype = name
				}
			case NOTAG:
			case ISCLOSED:
				// now handled at end of calling function
			default:
				buffer.WriteString(pfx)
				pfx = ""
				okIndent = false
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

		var prev XMLToken
		primed := false
		skipDoctype := false
		lastContent := false

		// track adjacent pairs to give look-ahead at next token
		doPair := func(tkn XMLToken) {

			if tkn.Tag == NOTAG {
				return
			}

			if tkn.Tag == DOCTYPETAG {
				if skipDoctype {
					return
				}
				skipDoctype = true
			}

			if primed {
				cleanToken(prev, tkn.Tag, tkn.Name, tkn.Attr, lastContent)
			}

			lastContent = (prev.Tag == CONTENTTAG)
			prev = XMLToken{tkn.Tag, tkn.Cont, tkn.Name, tkn.Attr, tkn.Index, tkn.Line}
			primed = true
		}

		if inp != nil {

			// track adjacent pairs to give look-ahead at next token
			for tkn := range inp {

				doPair(tkn)
			}

		} else {

			parseXML(rcrd, prnt, nil, doPair, nil, nil)
		}

		// isclosed tag
		if primed {
			buffer.WriteString(pfx)
			txt := buffer.String()
			if txt != "" {
				// send remaining result through output channel
				out <- txt
			}
		}
	}

	// launch single formatter goroutine
	go formatXML(rcrd, prnt, inp, offset, doXML, out)

	return out
}

// FormatRecord formats a single partitioned XML record
func FormatRecord(rcrd, prnt string, args FormatArgs) <-chan string {

	return xmlFormatter(rcrd, prnt, nil, 1, false, args)
}

// FormatTokens formats an XML token stream
func FormatTokens(inp <-chan XMLToken, args FormatArgs) <-chan string {

	return xmlFormatter("", "", inp, 0, true, args)
}
