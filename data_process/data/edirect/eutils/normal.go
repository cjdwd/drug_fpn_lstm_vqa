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
// File Name:  normal.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"fmt"
	"html"
	"os"
	"strings"
)

// NormalizeXML adjusts Entrez XML fields to conform to common conventions
func NormalizeXML(rdr <-chan XMLBlock, db string) <-chan string {

	if rdr == nil || db == "" {
		return nil
	}

	out := make(chan string, chanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "Unable to create normalize channel\n")
		os.Exit(1)
	}

	tknq := CreateTokenizer(rdr)

	if tknq == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create normalize tokenizer\n")
		os.Exit(1)
	}

	// force -strict cleanup flag for most databases (even after CreateReader and CreateTokenizer are called)
	switch db {
	case "bioc", "biocollections", "clinvar", "dbvar", "gap", "gapplus", "grasp", "pccompound", "pcsubstance":
		doMixed = true
		doStrict = false
	default:
		doStrict = true
		doMixed = false
	}
	allowEmbed = true
	contentMods = true

	normalizeXML := func(rdr <-chan XMLBlock, out chan<- string) {

		// close channel when all chunks have been sent
		defer close(out)

		var buffer strings.Builder

		count := 0

		uid := ""
		isDocsum := false
		prevName := ""

		// removes mixed content tags
		mfix := strings.NewReplacer(
			"<b>", "",
			"<i>", "",
			"<u>", "",
			"</b>", "",
			"</i>", "",
			"</u>", "",
			"<b/>", "",
			"<i/>", "",
			"<u/>", "",
			"<sup>", "",
			"<sub>", "",
			"</sup>", "",
			"</sub>", "",
			"<sup/>", "",
			"<sub/>", "",
		)

		// reencodes < and > to &lt and &gt
		rfix := strings.NewReplacer(
			"<", "&lt;",
			">", "&gt;",
		)

		for tkn := range tknq {

			tag := tkn.Tag
			name := tkn.Name
			attr := tkn.Attr

			switch tag {
			case STARTTAG:
				if name == "Id" && uid != "" {
					uid = ""
				}
				if uid != "" {
					// if object after DocumentSummary is not already Id, create Id from rescued attribute
					buffer.WriteString("<Id>\n")
					buffer.WriteString(uid)
					buffer.WriteString("\n</Id>\n")
					// clear until next docsum
					uid = ""
				}
				if name == "DocumentSummary" {
					isDocsum = true
					atts := ParseAttributes(attr)
					for i := 0; i < len(atts)-1; i += 2 {
						if atts[i] == "uid" {
							// store uid from DocumentSummary
							uid = atts[i+1]
							// if uid found, remove all attributes
							attr = ""
						}
					}
				}
				buffer.WriteString("<")
				buffer.WriteString(name)
				if attr != "" {
					attr = strings.TrimSpace(attr)
					attr = CompressRunsOfSpaces(attr)
					buffer.WriteString(" ")
					buffer.WriteString(attr)
				}
				buffer.WriteString(">\n")
				prevName = name
			case SELFTAG:
				buffer.WriteString("<")
				buffer.WriteString(name)
				if attr != "" {
					attr = strings.TrimSpace(attr)
					attr = CompressRunsOfSpaces(attr)
					buffer.WriteString(" ")
					buffer.WriteString(attr)
				}
				buffer.WriteString("/>\n")
			case STOPTAG:
				buffer.WriteString("</")
				buffer.WriteString(name)
				buffer.WriteString(">\n")
			case CONTENTTAG:
				if isDocsum {
					if db == "pubmed" && prevName == "Title" {
						if strings.Contains(name, "&") ||
							strings.Contains(name, "<") ||
							strings.Contains(name, ">") {
							ctype := tkn.Cont
							name = CleanupContents(name, (ctype&ASCII) != 0, (ctype&AMPER) != 0, (ctype&MIXED) != 0)
							if HasFlankingSpace(name) {
								name = strings.TrimSpace(name)
							}
							name = html.UnescapeString(name)
							// remove mixed content tags
							name = mfix.Replace(name)
							// reencode < and > to avoid breaking XML
							if strings.Contains(name, "<") || strings.Contains(name, ">") {
								name = rfix.Replace(name)
							}
						}
					} else if db == "gene" && prevName == "Summary" {
						if strings.Contains(name, "&amp;") {
							if HasFlankingSpace(name) {
								name = strings.TrimSpace(name)
							}
							name = html.UnescapeString(name)
							// reencode < and > to avoid breaking XML
							if strings.Contains(name, "<") || strings.Contains(name, ">") {
								name = rfix.Replace(name)
							}
						}
					} else if (db == "biosample" && prevName == "SampleData") ||
						(db == "medgen" && prevName == "ConceptMeta") ||
						(db == "sra" && prevName == "ExpXml") ||
						(db == "sra" && prevName == "Runs") {
						if strings.Contains(name, "&lt;") && strings.Contains(name, "&gt;") {
							if HasFlankingSpace(name) {
								name = strings.TrimSpace(name)
							}
							name = html.UnescapeString(name)
						}
					}
				} else {
					if db == "pubmed" {
						ctype := tkn.Cont
						name = CleanupContents(name, (ctype&ASCII) != 0, (ctype&AMPER) != 0, (ctype&MIXED) != 0)
						if HasFlankingSpace(name) {
							name = strings.TrimSpace(name)
						}
					} else if db == "bioc" {
						name = CleanupContents(name, true, true, true)
						if HasFlankingSpace(name) {
							name = strings.TrimSpace(name)
						}
					}
				}
				// content normally printed
				if HasFlankingSpace(name) {
					name = strings.TrimSpace(name)
				}
				buffer.WriteString(name)
				buffer.WriteString("\n")
			case CDATATAG:
				if isDocsum {
					if db == "assembly" && prevName == "Meta" {
						if strings.Contains(name, "<") && strings.Contains(name, ">") {
							// if CDATA contains embedded XML, simply remove CDATA wrapper
							if HasFlankingSpace(name) {
								name = strings.TrimSpace(name)
							}
							buffer.WriteString(name)
							buffer.WriteString("\n")
						}
					} else if db == "gtr" && prevName == "Extra" {
						// remove entire CDATA contents
					}
				}
			case COMMENTTAG:
				if !isDocsum {
					if db == "sra" {
						if strings.Contains(name, "<") && strings.Contains(name, ">") {
							// if comment contains embedded XML, remove comment wrapper and trim to leading < bracket
							pos := strings.Index(name, "<")
							if pos > 0 {
								name = name[pos:]
							}
							if HasFlankingSpace(name) {
								name = strings.TrimSpace(name)
							}
							buffer.WriteString(name)
							buffer.WriteString("\n")
						}
					}
				}
			case DOCTYPETAG:
				doctype := strings.TrimSpace(name)
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
			case NOTAG:
			case ISCLOSED:
				// now handled at end of calling function
			default:
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

		txt := buffer.String()
		if txt != "" {
			// send remaining result through output channel
			out <- txt
		}
	}

	// launch single normalize goroutine
	go normalizeXML(rdr, out)

	return out
}
