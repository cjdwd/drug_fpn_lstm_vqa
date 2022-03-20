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
// File Name:  parse.go
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

// XML token type
const (
	NOTAG = iota
	STARTTAG
	SELFTAG
	STOPTAG
	ATTRIBTAG
	CONTENTTAG
	CDATATAG
	COMMENTTAG
	DOCTYPETAG
	OBJECTTAG
	CONTAINERTAG
	ISCLOSED
	BADTAG
)

// content bit flags for performing special cleanup steps only when needed
const (
	NONE  = iota
	MIXED = 1 << iota
	AMPER
	ASCII
	LFTSPACE
	RGTSPACE
)

// internal tracking state for detecting unrecognized mixed content
const (
	_ = iota
	START
	STOP
	CHAR
	OTHER
)

// XMLNode is the node for an internal tree structure representing a single XML record
type XMLNode struct {
	Name       string
	Parent     string
	Contents   string
	Attributes string
	Attribs    []string
	Children   *XMLNode
	Next       *XMLNode
}

// XMLFind contains individual field values for finding a particular object
type XMLFind struct {
	Index  string
	Parent string
	Match  string
	Attrib string
	Versn  string
}

// XMLToken is the unit of XML parsing
type XMLToken struct {
	Tag   int
	Cont  int
	Name  string
	Attr  string
	Index int
	Line  int
}

// ParseAttributes produces tag/value pairs, only run on request.
func ParseAttributes(attrb string) []string {

	if attrb == "" {
		return nil
	}

	attlen := len(attrb)

	// count equal signs
	num := 0
	inQuote := false

	for i := 0; i < attlen; i++ {
		ch := attrb[i]
		if ch == '"' || ch == '\'' {
			// "
			inQuote = !inQuote
		}
		if ch == '=' && !inQuote {
			num += 2
		}
	}
	if num < 1 {
		return nil
	}

	// allocate array of proper size
	arry := make([]string, num)
	if arry == nil {
		return nil
	}

	start := 0
	idx := 0
	itm := 0
	inQuote = false

	// place tag and value in successive array slots
	for idx < attlen && itm < num {
		ch := attrb[idx]
		if ch == '"' || ch == '\'' {
			// "
			inQuote = !inQuote
		}
		if ch == '=' && !inQuote {
			inQuote = true
			// skip past possible leading blanks
			for start < attlen {
				ch = attrb[start]
				if inBlank[ch] {
					start++
				} else {
					break
				}
			}
			// =
			arry[itm] = strings.TrimSpace(attrb[start:idx])
			itm++
			// skip past equal sign
			idx++
			ch = attrb[idx]
			if ch != '"' && ch != '\'' {
				// "
				// skip past unexpected blanks
				for inBlank[ch] {
					idx++
					ch = attrb[idx]
				}
				if ch != '"' && ch != '\'' {
					// "
					fmt.Fprintf(os.Stderr, "\nAttribute in '%s' missing double quote\n", attrb)
				}
			}
			// skip past leading double quote
			idx++
			start = idx
		} else if ch == '"' || ch == '\'' {
			// "
			inQuote = false
			arry[itm] = strings.TrimSpace(attrb[start:idx])
			itm++
			// skip past trailing double quote and (possible) space
			idx += 2
			start = idx
		} else {
			idx++
		}
	}

	return arry
}

// parseXML calls XML parser on a partitioned string or on an XMLBlock channel of trimmed strings.
// It is optimized for maximum processing speed, sends tokens for CDATA and COMMENT sections (for
// unpacking by NormalizeXML), and optionally tracks line numbers (for ValidateXML).
func parseXML(record, parent string, inp <-chan XMLBlock, tokens func(XMLToken), find *XMLFind, ids func(string)) (*XMLNode, string) {

	if record == "" && (inp == nil || tokens == nil) {
		return nil, ""
	}

	// token parser variables
	recLen := len(record)
	Idx := 0

	// line tracking variables
	lineNum := 1
	lag := 0

	// variables to track COMMENT or CDATA sections that span reader blocks
	which := NOTAG
	skipTo := ""

	// updateLineCount is used to keep track of the correct line count for XML validation
	updateLineCount := func(max int) {
		// count lines
		for i := lag; i < max; i++ {
			if record[i] == '\n' {
				lineNum++
			}
		}
		lag = Idx
	}

	// currentLineCount calculates correct line for warning messages, does not update lineNum or lag variables
	currentLineCount := func(max int) int {
		line := lineNum
		for i := lag; i < max; i++ {
			if record[i] == '\n' {
				line++
			}
		}
		return line
	}

	// nextToken returns the type and content fields for the next XML token
	nextToken := func(idx int) (int, int, string, string, int) {

		if record == "" {
			// buffer is empty
			if inp != nil {
				// read next block if available
				record = string(<-inp)
				recLen = len(record)
				Idx = 0
				idx = 0
				lag = 0
			}

			if record == "" {
				// signal end of XML data
				return ISCLOSED, NONE, "", "", 0
			}

			if which != NOTAG && skipTo != "" {
				// previous block ended inside CDATA object or COMMENT
				text := record[:]
				txtlen := recLen
				whch := which
				start := idx
				found := strings.Index(text[idx:], skipTo)
				if found < 0 {
					// no stop signal found in next block
					str := text[start:]
					if HasFlankingSpace(str) {
						str = strings.TrimSpace(str)
					}

					if countLines {
						updateLineCount(txtlen)
					}

					// signal end of current block
					record = ""

					// leave which and skipTo values unchanged as another continuation signal
					// send CDATA or COMMENT contents
					return whch, NONE, str[:], "", 0
				}
				// otherwise adjust position past end of skipTo string and return to normal processing
				idx += found
				str := text[start:idx]
				if HasFlankingSpace(str) {
					str = strings.TrimSpace(str)
				}
				idx += len(skipTo)
				// clear tracking variables
				which = NOTAG
				skipTo = ""
				// send CDATA or COMMENT contents
				return whch, NONE, str[:], "", idx
			}
		}

		text := record[:]
		txtlen := recLen

		// XML string, and all blocks, end with > character, acts as sentinel to check if past end of text
		if idx >= txtlen {
			if inp != nil {

				if countLines {
					updateLineCount(txtlen)
				}

				// signal end of current block, will read next block on next call
				record = ""

				return NOTAG, NONE, "", "", 0
			}

			// signal end of XML string
			return ISCLOSED, NONE, "", "", 0
		}

		ctype := NONE

		// skip past leading blanks
		ch := text[idx]
		if inBlank[ch] {
			ctype |= LFTSPACE
			idx++
			ch = text[idx]
			for inBlank[ch] {
				idx++
				ch = text[idx]
			}
		}

		start := idx

		plainContent := true

		if doStrict && ch == '<' {
			// check to see if an HTML or MathML element is at the beginning of a content string
			if HTMLAhead(text, idx, txtlen) != 0 {
				plainContent = false
			}
		}

		if plainContent && ch == '<' {

			// at start of element
			idx++
			ch = text[idx]

			// check for legal first character of element
			if inFirst[ch] {

				// read element name
				start = idx
				idx++

				ch = text[idx]
				for inElement[ch] {
					idx++
					ch = text[idx]
				}

				str := text[start:idx]

				if ch == '>' {

					// end of element
					idx++

					return STARTTAG, NONE, str[:], "", idx

				} else if ch == '/' {

					// self-closing element without attributes
					idx++
					ch = text[idx]
					if ch != '>' {
						// skip past unexpected blanks
						for inBlank[ch] {
							idx++
							ch = text[idx]
						}
						if ch != '>' {
							fmt.Fprintf(os.Stderr, "\nSelf-closing element missing right angle bracket\n")
						}
					}
					idx++

					return SELFTAG, NONE, str[:], "", idx

				} else if inBlank[ch] {

					// attributes
					idx++
					ch = text[idx]
					// skip past unexpected blanks
					for inBlank[ch] {
						idx++
						ch = text[idx]
					}
					start = idx
					for ch != '<' && ch != '>' {
						idx++
						ch = text[idx]
					}
					if ch != '>' {
						fmt.Fprintf(os.Stderr, "\nAttributes not followed by right angle bracket\n")
					}
					// walk back past trailing blanks
					lst := idx - 1
					ch = text[lst]
					for inBlank[ch] && lst > start {
						lst--
						ch = text[lst]
					}
					if ch == '/' {
						// self-closing
						atr := text[start:lst]
						idx++

						return SELFTAG, NONE, str[:], atr[:], idx
					}
					atr := text[start:idx]
					idx++

					return STARTTAG, NONE, str[:], atr[:], idx

				} else {

					if countLines {
						fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element, line %d\n", ch, currentLineCount(idx))
					} else {
						fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element\n", ch)
					}

					return STARTTAG, NONE, str[:], "", idx
				}

				// other punctuation character immediately after first angle bracket

			} else if ch == '/' {

				// at start of end tag
				idx++
				start = idx
				ch = text[idx]
				// expect legal first character of element
				if inFirst[ch] {
					idx++
					ch = text[idx]
					for inElement[ch] {
						idx++
						ch = text[idx]
					}
					str := text[start:idx]
					if ch != '>' {
						// skip past unexpected blanks
						for inBlank[ch] {
							idx++
							ch = text[idx]
						}
						if ch != '>' {
							fmt.Fprintf(os.Stderr, "\nUnexpected characters after end element name\n")
						}
					}
					idx++

					return STOPTAG, NONE, str[:], "", idx
				}
				// legal character not found after slash
				if countLines {
					fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element, line %d\n", ch, currentLineCount(idx))
				} else {
					fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element\n", ch)
				}

			} else if ch == '!' {

				// skip !DOCTYPE, !COMMENT, and ![CDATA[
				idx++
				start = idx
				ch = text[idx]
				which = NOTAG
				skipTo = ""
				if ch == '[' && strings.HasPrefix(text[idx:], "[CDATA[") {
					which = CDATATAG
					skipTo = "]]>"
					start += 7
				} else if ch == '-' && strings.HasPrefix(text[idx:], "--") {
					which = COMMENTTAG
					skipTo = "-->"
					start += 2
				} else if ch == 'D' && strings.HasPrefix(text[idx:], "DOCTYPE") {
					which = DOCTYPETAG
					skipTo = ">"
				}
				if which != NOTAG && skipTo != "" {
					whch := which
					// CDATA or COMMENT block may contain internal angle brackets
					found := strings.Index(text[idx:], skipTo)
					if found < 0 {
						// string stops in middle of CDATA or COMMENT
						if inp != nil {
							str := text[start:]
							if HasFlankingSpace(str) {
								str = strings.TrimSpace(str)
							}

							if countLines {
								updateLineCount(txtlen)
							}

							// signal end of current block
							record = ""

							// leave which and skipTo values unchanged as another continuation signal
							// send CDATA or COMMENT contents
							return whch, NONE, str[:], "", 0
						}

						return ISCLOSED, NONE, "", "", idx
					}
					// adjust position past end of CDATA or COMMENT
					if inp != nil {
						idx += found
						str := text[start:idx]
						if HasFlankingSpace(str) {
							str = strings.TrimSpace(str)
						}
						idx += len(skipTo)
						// clear tracking variables
						which = NOTAG
						skipTo = ""
						// send CDATA or COMMENT contents
						return whch, NONE, str[:], "", idx
					}

					idx += found + len(skipTo)
					return NOTAG, NONE, "", "", idx
				}
				// otherwise just skip to next right angle bracket
				for ch != '>' {
					idx++
					ch = text[idx]
				}
				idx++
				return NOTAG, NONE, "", "", idx

			} else if ch == '?' {

				// skip ?xml and ?processing instructions
				idx++
				ch = text[idx]
				for ch != '>' {
					idx++
					ch = text[idx]
				}
				idx++
				return NOTAG, NONE, "", "", idx

			} else {

				if countLines {
					fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element, line %d\n", ch, currentLineCount(idx))
				} else {
					fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element\n", ch)
				}
			}

		} else if ch != '>' {

			// at start of contents
			start = idx

			hasMarkup := false
			hasNonASCII := false

			// find end of contents
			if allowEmbed {

				for {
					for inContent[ch] {
						idx++
						ch = text[idx]
					}
					// set flags to speed up conditional content processing
					if ch == '&' {
						idx++
						ch = text[idx]
						if ch == 'a' {
							if strings.HasPrefix(text[idx:], "amp;") {
								hasMarkup = true
							}
						} else if ch == 'g' {
							if strings.HasPrefix(text[idx:], "gt;") {
								hasMarkup = true
							}
						} else if ch == 'l' {
							if strings.HasPrefix(text[idx:], "lt;") {
								hasMarkup = true
							}
						}
						continue
					}
					if ch > 127 {
						hasNonASCII = true
						idx++
						ch = text[idx]
						continue
					}
					if ch == '<' && doStrict {
						// optionally allow HTML text formatting elements and super/subscripts
						advance := HTMLAhead(text, idx, txtlen)
						if advance > 0 {
							idx += advance
							if idx < txtlen {
								ch = text[idx]
							}
							plainContent = false
							continue
						}
					}
					break
				}

			} else {
				for ch != '<' && ch != '>' {
					idx++
					ch = text[idx]
				}
			}

			// trim back past trailing blanks
			lst := idx - 1
			ch = text[lst]
			if inBlank[ch] && lst > start {
				ctype |= RGTSPACE
				lst--
				ch = text[lst]
				for inBlank[ch] && lst > start {
					lst--
					ch = text[lst]
				}
			}

			str := text[start : lst+1]

			if allowEmbed {
				if !plainContent {
					ctype |= MIXED
				}
				if hasMarkup {
					ctype |= AMPER
				}
				if hasNonASCII {
					ctype |= ASCII
				}
			}

			return CONTENTTAG, ctype, str[:], "", idx
		}

		return BADTAG, NONE, "", "", idx
	}

	// node farm variables
	farmPos := 0
	farmMax := farmSize
	farmItems := make([]XMLNode, farmMax)

	// nextNode allocates multiple nodes in a large array for memory management efficiency
	nextNode := func(strt, attr, prnt string) *XMLNode {

		// if farm array slots used up, allocate new array
		if farmPos >= farmMax {
			farmItems = make([]XMLNode, farmMax)
			farmPos = 0
		}

		if farmItems == nil {
			return nil
		}

		// take node from next available slot in farm array
		node := &farmItems[farmPos]

		node.Name = strt[:]
		node.Attributes = attr[:]
		node.Parent = prnt[:]

		farmPos++

		return node
	}

	// Parse tokens into tree structure for exploration

	// parseSpecial recursive definition
	var parseSpecial func(string, string, string) (*XMLNode, bool)

	// parseSpecial parses XML tags into tree structure for searching, no contentMods flags set
	parseSpecial = func(strt, attr, prnt string) (*XMLNode, bool) {

		var obj *XMLNode
		ok := true

		// nextNode obtains next node from farm
		node := nextNode(strt, attr, prnt)
		if node == nil {
			return nil, false
		}

		var lastNode *XMLNode

		status := START
		for {
			tag, _, name, attr, idx := nextToken(Idx)
			Idx = idx

			if tag == BADTAG {
				if countLines {
					fmt.Fprintf(os.Stderr, "\n%sERROR: Unparsable XML element, line %d%s\n", RED, lineNum, INIT)
				} else {
					fmt.Fprintf(os.Stderr, "\n%sERROR: Unparsable XML element%s\n", RED, INIT)
				}
				break
			}
			if tag == ISCLOSED {
				break
			}

			switch tag {
			case STARTTAG:
				if status == CHAR {
					fmt.Fprintf(os.Stderr, "%s ERROR: %s UNEXPECTED MIXED CONTENT <%s> IN <%s>%s\n", INVT, LOUD, name, prnt, INIT)
				}
				// read sub tree
				obj, ok = parseSpecial(name, attr, node.Name)
				if !ok {
					break
				}

				// adding next child to end of linked list gives better performance than appending to slice of nodes
				if node.Children == nil {
					node.Children = obj
				}
				if lastNode != nil {
					lastNode.Next = obj
				}
				lastNode = obj
				status = STOP
			case STOPTAG:
				// pop out of recursive call
				return node, ok
			case CONTENTTAG:
				node.Contents = name
				status = CHAR
			case SELFTAG:
				if attr == "" && !doSelf {
					// ignore if self-closing tag has no attributes
					continue
				}

				// self-closing tag has no contents, just create child node
				obj = nextNode(name, attr, node.Name)

				if doSelf {
					// add default value for self-closing tag
					obj.Contents = "1"
				}

				if node.Children == nil {
					node.Children = obj
				}
				if lastNode != nil {
					lastNode.Next = obj
				}
				lastNode = obj
				status = OTHER
				// continue on same level
			default:
				status = OTHER
			}
		}

		return node, ok
	}

	// parseLevel recursive definition
	var parseLevel func(string, string, string) (*XMLNode, bool)

	// parseLevel parses XML tags into tree structure for searching, some contentMods flags set
	parseLevel = func(strt, attr, prnt string) (*XMLNode, bool) {

		var obj *XMLNode
		ok := true

		// obtain next node from farm
		node := nextNode(strt, attr, prnt)
		if node == nil {
			return nil, false
		}

		var lastNode *XMLNode

		status := START
		for {
			tag, ctype, name, attr, idx := nextToken(Idx)
			Idx = idx

			if countLines && Idx > 0 {
				updateLineCount(Idx)
			}

			if tag == BADTAG {
				if countLines {
					fmt.Fprintf(os.Stderr, "\n%sERROR: Unparsable XML element, line %d%s\n", RED, lineNum, INIT)
				} else {
					fmt.Fprintf(os.Stderr, "\n%sERROR: Unparsable XML element%s\n", RED, INIT)
				}
				break
			}
			if tag == ISCLOSED {
				break
			}

			switch tag {
			case STARTTAG:
				if status == CHAR {
					if doStrict {
						fmt.Fprintf(os.Stderr, "%s ERROR: %s UNRECOGNIZED MIXED CONTENT <%s> IN <%s>%s\n", INVT, LOUD, name, prnt, INIT)
					} else if !doMixed {
						fmt.Fprintf(os.Stderr, "%s ERROR: %s UNEXPECTED MIXED CONTENT <%s> IN <%s>%s\n", INVT, LOUD, name, prnt, INIT)
					}
				}
				// read sub tree
				obj, ok = parseLevel(name, attr, node.Name)
				if !ok {
					break
				}

				// adding next child to end of linked list gives better performance than appending to slice of nodes
				if node.Children == nil {
					node.Children = obj
				}
				if lastNode != nil {
					lastNode.Next = obj
				}
				lastNode = obj
				status = STOP
			case STOPTAG:
				// pop out of recursive call
				return node, ok
			case CONTENTTAG:
				if doMixed {
					// create unnamed child node for content string
					con := nextNode("", "", "")
					if con == nil {
						break
					}
					str := CleanupContents(name, (ctype&ASCII) != 0, (ctype&AMPER) != 0, (ctype&MIXED) != 0)
					if (ctype & LFTSPACE) != 0 {
						str = " " + str
					}
					if (ctype & RGTSPACE) != 0 {
						str += " "
					}
					con.Contents = str
					if node.Children == nil {
						node.Children = con
					}
					if lastNode != nil {
						lastNode.Next = con
					}
					lastNode = con
				} else {
					node.Contents = CleanupContents(name, (ctype&ASCII) != 0, (ctype&AMPER) != 0, (ctype&MIXED) != 0)
				}
				status = CHAR
			case SELFTAG:
				if attr == "" && !doSelf {
					// ignore if self-closing tag has no attributes
					continue
				}

				// self-closing tag has no contents, just create child node
				obj = nextNode(name, attr, node.Name)

				if doSelf {
					// add default value for self-closing tag
					obj.Contents = "1"
				}

				if node.Children == nil {
					node.Children = obj
				}
				if lastNode != nil {
					lastNode.Next = obj
				}
				lastNode = obj
				status = OTHER
				// continue on same level
			default:
				status = OTHER
			}
		}

		return node, ok
	}

	// parseIndex recursive definition
	var parseIndex func(string, string, string) string

	// parseIndex parses XML tags looking for trie index element
	parseIndex = func(strt, attr, prnt string) string {

		versn := ""

		// check for version attribute match
		if attr != "" && find.Versn != "" && strings.Contains(attr, find.Versn) {
			if strt == find.Match || find.Match == "" {
				if find.Parent == "" || prnt == find.Parent {
					attribs := ParseAttributes(attr)
					for i := 0; i < len(attribs)-1; i += 2 {
						if attribs[i] == find.Versn {
							versn = attribs[i+1]
						}
					}
				}
			}
		}

		// check for attribute index match
		if attr != "" && find.Attrib != "" && strings.Contains(attr, find.Attrib) {
			if strt == find.Match || find.Match == "" {
				if find.Parent == "" || prnt == find.Parent {
					attribs := ParseAttributes(attr)
					for i := 0; i < len(attribs)-1; i += 2 {
						if attribs[i] == find.Attrib {
							return attribs[i+1]
						}
					}
				}
			}
		}

		for {
			tag, _, name, attr, idx := nextToken(Idx)
			Idx = idx

			if tag == BADTAG {
				if countLines {
					fmt.Fprintf(os.Stderr, "\n%sERROR: Unparsable XML element, line %d%s\n", RED, lineNum, INIT)
				} else {
					fmt.Fprintf(os.Stderr, "\n%sERROR: Unparsable XML element%s\n", RED, INIT)
				}
				break
			}
			if tag == ISCLOSED {
				break
			}

			switch tag {
			case STARTTAG:
				id := parseIndex(name, attr, strt)
				if id != "" {
					return id
				}
			case SELFTAG:
			case STOPTAG:
				// break recursion
				return ""
			case CONTENTTAG:
				// check for content index match
				if strt == find.Match || find.Match == "" {
					if find.Parent == "" || prnt == find.Parent {
						// append version if specified as parent/element@attribute^version
						if versn != "" {
							name += "."
							name += versn
						}
						if ids != nil {
							ids(name)
						} else {
							return name
						}
					}
				}
			default:
			}
		}

		return ""
	}

	// main loops

	// stream all tokens through callback
	if tokens != nil {

		for {
			tag, ctype, name, attr, idx := nextToken(Idx)
			Idx = idx

			if countLines && Idx > 0 {
				updateLineCount(Idx)
			}

			if tag == BADTAG {
				if countLines {
					fmt.Fprintf(os.Stderr, "\n%sERROR: Unparsable XML element, line %d%s\n", RED, lineNum, INIT)
				} else {
					fmt.Fprintf(os.Stderr, "\n%sERROR: Unparsable XML element%s\n", RED, INIT)
				}
				break
			}

			tkn := XMLToken{tag, ctype, name, attr, idx, lineNum}

			tokens(tkn)

			if tag == ISCLOSED {
				break
			}
		}

		return nil, ""
	}

	// find value of index element
	if find != nil && find.Index != "" {

		// return indexed identifier

		tag, _, name, attr, idx := nextToken(Idx)

		// loop until start tag
		for {
			Idx = idx

			if tag == BADTAG {
				if countLines {
					fmt.Fprintf(os.Stderr, "\n%sERROR: Unparsable XML element, line %d%s\n", RED, lineNum, INIT)
				} else {
					fmt.Fprintf(os.Stderr, "\n%sERROR: Unparsable XML element%s\n", RED, INIT)
				}
				break
			}
			if tag == ISCLOSED {
				break
			}

			if tag == STARTTAG {
				break
			}

			tag, _, name, attr, idx = nextToken(Idx)
		}

		return nil, parseIndex(name, attr, parent)
	}

	// otherwise create node tree for general data extraction
	tag, _, name, attr, idx := nextToken(Idx)

	// loop until start tag
	for {
		Idx = idx

		if tag == BADTAG {
			if countLines {
				fmt.Fprintf(os.Stderr, "\n%sERROR: Unparsable XML element, line %d%s\n", RED, lineNum, INIT)
			} else {
				fmt.Fprintf(os.Stderr, "\n%sERROR: Unparsable XML element%s\n", RED, INIT)
			}
			break
		}
		if tag == ISCLOSED {
			break
		}

		if tag == STARTTAG {
			break
		}

		tag, _, name, attr, idx = nextToken(Idx)
	}

	if contentMods {
		// slower parser also handles mixed content
		top, ok := parseLevel(name, attr, parent)

		if !ok {
			return nil, ""
		}

		return top, ""
	}

	// fastest parsing with no contentMods flags
	top, ok := parseSpecial(name, attr, parent)

	if !ok {
		return nil, ""
	}

	return top, ""
}

// ParseRecord is the main public access to parseXML.
func ParseRecord(text, parent string) *XMLNode {

	pat, _ := parseXML(text, parent, nil, nil, nil, nil)

	return pat
}

// FindIdentifier returns a single identifier.
func FindIdentifier(text, parent string, find *XMLFind) string {

	_, id := parseXML(text, parent, nil, nil, find, nil)

	return id
}

// FindIdentifiers returns a set of identifiers through a callback.
func FindIdentifiers(text, parent string, find *XMLFind, ids func(string)) {

	parseXML(text, parent, nil, nil, find, ids)
}

// StreamTokens streams tokens from a reader through a callback.
func StreamTokens(inp <-chan XMLBlock, streamer func(tkn XMLToken)) {

	parseXML("", "", inp, streamer, nil, nil)
}

// StreamValues streams token values from a parsed record through a callback.
func StreamValues(text, parent string, stream func(string, string, string)) {

	elementName := ""
	attributeName := ""

	streamer := func(tkn XMLToken) {

		switch tkn.Tag {
		case STARTTAG:
			elementName = tkn.Name
			attributeName = tkn.Attr
		case CONTENTTAG:
			// send element name and content to callback
			stream(elementName, attributeName, tkn.Name)
		default:
		}
	}

	parseXML(text, parent, nil, streamer, nil, nil)
}

// CreateTokenizer streams tokens through a channel.
func CreateTokenizer(inp <-chan XMLBlock) <-chan XMLToken {

	if inp == nil {
		return nil
	}

	out := make(chan XMLToken, chanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\n%sERROR: Unable to create XML tokenizer channel%s\n", RED, INIT)
		os.Exit(1)
	}

	// xmlTokenizer sends XML tokens through channel
	xmlTokenizer := func(inp <-chan XMLBlock, out chan<- XMLToken) {

		// close channel when all records have been processed
		defer close(out)

		// parse XML and send tokens through channel
		parseXML("", "", inp, func(tkn XMLToken) { out <- tkn }, nil, nil)
	}

	// launch single tokenizer goroutine
	go xmlTokenizer(inp, out)

	return out
}

// ExploreElements returns matching element values to callback.
func ExploreElements(curr *XMLNode, mask, prnt, match, attrib string, wildcard, unescape bool, level int, proc func(string, int)) {

	if curr == nil || proc == nil {
		return
	}

	// **/Object performs deep exploration of recursive data (*/Object also supported)
	deep := false
	if prnt == "**" || prnt == "*" {
		prnt = ""
		deep = true
	}

	// exploreChildren recursive definition
	var exploreChildren func(curr *XMLNode, acc func(string))

	// exploreChildren handles mixed-content chains of embedded tags
	exploreChildren = func(curr *XMLNode, acc func(string)) {

		if curr.Contents != "" {
			acc(curr.Contents)
		}
		for chld := curr.Children; chld != nil; chld = chld.Next {
			if chld.Name != "" {
				acc("<" + chld.Name + ">")
			}
			exploreChildren(chld, acc)
			if chld.Name != "" {
				acc("</" + chld.Name + ">")
			}
		}
	}

	// exploreElements recursive definition
	var exploreElements func(curr *XMLNode, skip string, lev int)

	// exploreElements visits nodes looking for matches to requested object
	exploreElements = func(curr *XMLNode, skip string, lev int) {

		if !deep && curr.Name == skip {
			// do not explore within recursive object
			return
		}

		if curr.Name == match ||
			// parent/* matches any subfield
			(match == "*" && prnt != "") ||
			// wildcard (internal colon) matches any namespace prefix
			(wildcard && strings.HasPrefix(match, ":") && strings.HasSuffix(curr.Name, match)) ||
			(match == "" && attrib != "") {

			if prnt == "" ||
				curr.Parent == prnt ||
				(wildcard && strings.HasPrefix(prnt, ":") && strings.HasSuffix(curr.Parent, prnt)) {

				if attrib != "" {
					if curr.Attributes != "" && curr.Attribs == nil {
						// parse attributes on-the-fly if queried
						curr.Attribs = ParseAttributes(curr.Attributes)
					}
					for i := 0; i < len(curr.Attribs)-1; i += 2 {
						// attributes now parsed into array as [ tag, value, tag, value, tag, value, ... ]
						if curr.Attribs[i] == attrib ||
							(wildcard && strings.HasPrefix(attrib, ":") && strings.HasSuffix(curr.Attribs[i], attrib)) {
							proc(curr.Attribs[i+1], level)
							return
						}
					}

				} else if curr.Contents != "" {

					str := curr.Contents[:]

					if unescape && HasAmpOrNotASCII(str) {
						// processing of <, >, &, ", and ' characters is now delayed until element contents is requested
						str = html.UnescapeString(str)
					}

					proc(str, level)
					return

				} else if curr.Children != nil {

					if doMixed {
						// match with mixed contents - send all child strings
						var buffr strings.Builder
						exploreChildren(curr, func(str string) {
							if str != "" {
								buffr.WriteString(str)
							}
						})
						str := buffr.String()

						// clean up reconstructed mixed content
						str = DoTrimFlankingHTML(str)
						if HasBadSpace(str) {
							str = CleanupBadSpaces(str)
						}
						if HasAdjacentSpaces(str) {
							str = CompressRunsOfSpaces(str)
						}
						if NeedsTightening(str) {
							str = TightenParentheses(str)
						}
						if unescape && HasAmpOrNotASCII(str) {
							str = html.UnescapeString(str)
						}

						proc(str, level)
						return
					}

					// for XML container object, send empty string to callback to increment count
					proc("", level)
					// and continue exploring

				} else if curr.Attributes != "" {

					// for self-closing object, indicate presence by sending empty string to callback
					proc("", level)
					return
				}
			}
		}

		for chld := curr.Children; chld != nil; chld = chld.Next {
			// inner exploration is subject to recursive object exclusion
			exploreElements(chld, mask, lev+1)
		}
	}

	// start recursive exploration from current scope
	exploreElements(curr, "", level)
}

// ExploreNodes visits XML container nodes.
func ExploreNodes(curr *XMLNode, prnt, match string, index, level int, proc func(*XMLNode, int, int)) {

	if curr == nil || proc == nil {
		return
	}

	// leading colon indicates namespace prefix wildcard
	wildcard := false
	if strings.HasPrefix(prnt, ":") || strings.HasPrefix(match, ":") {
		wildcard = true
	}

	// **/Object performs deep exploration of recursive data
	deep := false
	if prnt == "**" {
		prnt = "*"
		deep = true
	}
	// Object/** performs exhaustive exploration of nodes
	tall := false
	if match == "**" {
		match = "*"
		tall = true
	}

	// exploreNodes recursive definition
	var exploreNodes func(*XMLNode, int, int, bool, func(*XMLNode, int, int)) int

	// exploreNodes visits all nodes that match the selection criteria
	exploreNodes = func(curr *XMLNode, indx, levl int, force bool, proc func(*XMLNode, int, int)) int {

		if curr == nil || proc == nil {
			return indx
		}

		// match is "*" for heterogeneous data constructs, e.g., -group PubmedArticleSet/*
		// wildcard matches any namespace prefix
		if curr.Name == match ||
			match == "*" ||
			(wildcard && strings.HasPrefix(match, ":") && strings.HasSuffix(curr.Name, match)) {

			if prnt == "" ||
				curr.Parent == prnt ||
				force ||
				(wildcard && strings.HasPrefix(prnt, ":") && strings.HasSuffix(curr.Parent, prnt)) {

				proc(curr, indx, levl)
				indx++

				if tall && prnt != "" {
					// exhaustive exploration of child nodes within region of parent match
					for chld := curr.Children; chld != nil; chld = chld.Next {
						indx = exploreNodes(chld, indx, levl+1, true, proc)
					}
				}

				if !deep {
					// do not explore within recursive object
					return indx
				}
			}
		}

		// clearing prnt "*" now allows nested exploration within recursive data, e.g., -pattern Taxon -block */Taxon
		if prnt == "*" {
			prnt = ""
		}

		// explore child nodes
		for chld := curr.Children; chld != nil; chld = chld.Next {
			indx = exploreNodes(chld, indx, levl+1, false, proc)
		}

		return indx
	}

	exploreNodes(curr, index, level, false, proc)
}
