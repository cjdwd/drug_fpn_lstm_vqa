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
// File Name:  transmute.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package main

import (
	"encoding/base64"
	"eutils"
	"fmt"
	"html"
	"io"
	"io/ioutil"
	"net/url"
	"os"
	"regexp"
	"runtime"
	"runtime/debug"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"sync"
	"unicode"
)

// XML FORMATTING FUNCTIONS

// createFormatters does concurrent reformatting, using flush-left to remove leading spaces
func createFormatters(parent string, format string, inp <-chan eutils.XMLRecord) <-chan eutils.XMLRecord {

	if inp == nil {
		return nil
	}

	out := make(chan eutils.XMLRecord, eutils.ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create formatter channel\n")
		os.Exit(1)
	}

	if format == "" {
		format = "flush"
	}

	// xmlFormatter reads partitioned XML from channel and formats on a per-record basis
	xmlFormatter := func(wg *sync.WaitGroup, parent string, inp <-chan eutils.XMLRecord, out chan<- eutils.XMLRecord) {

		// report when this formatter has no more records to process
		defer wg.Done()

		// read partitioned XML from producer channel
		for ext := range inp {

			idx := ext.Index
			text := ext.Text

			if text == "" {
				// should never see empty input data
				out <- eutils.XMLRecord{Index: idx, Text: text}
				continue
			}

			// str := doFormat(text[:], parent)

			frm := eutils.FormatRecord(text, parent, eutils.FormatArgs{Format: format})
			str := eutils.ChanToString(frm)

			// send even if empty to get all record counts for reordering
			out <- eutils.XMLRecord{Index: idx, Text: str}
		}
	}

	var wg sync.WaitGroup

	// launch multiple formatter goroutines
	for i := 0; i < eutils.NumServe(); i++ {
		wg.Add(1)
		go xmlFormatter(&wg, parent, inp, out)
	}

	// launch separate anonymous goroutine to wait until all formatters are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

// processFormat reformats XML for ease of reading
func processFormat(rdr <-chan eutils.XMLBlock, args []string) {

	if rdr == nil || args == nil {
		return
	}

	// skip past command name
	args = args[1:]

	format := ""
	xml := ""
	doctype := ""

	doCombine := false
	doSelf := false
	doComment := false
	doCdata := false

	if len(args) > 0 {
		// look for [compact|flush|indent|expand] specification
		format = args[0]
		if strings.HasPrefix(format, "-") {
			// ran into next argument, default to indent
			format = "indent"
		} else {
			// skip past first argument
			args = args[1:]
		}
	} else {
		format = "indent"
	}

	// look for remaining arguments
	for len(args) > 0 {

		switch args[0] {
		case "-xml":
			args = args[1:]
			// -xml argument must be followed by value to use in xml line
			if len(args) < 1 || strings.HasPrefix(args[0], "-") {
				fmt.Fprintf(os.Stderr, "\nERROR: -xml argument is missing\n")
				os.Exit(1)
			}
			xml = args[0]
			args = args[1:]
		case "-doctype":
			args = args[1:]
			if len(args) > 0 {
				// if -doctype argument followed by value, use instead of DOCTYPE line
				doctype = args[0]
				args = args[1:]
			}
		/*
			// allow setting of unicode, script, and mathml flags within -format
			case "-unicode":
				if len(args) < 2 {
					fmt.Fprintf(os.Stderr, "\nERROR: Unicode argument is missing\n")
					os.Exit(1)
				}
				// unicodePolicy = args[1]
				args = args[2:]
			case "-script":
				if len(args) < 2 {
					fmt.Fprintf(os.Stderr, "\nERROR: Script argument is missing\n")
					os.Exit(1)
				}
				// scriptPolicy = args[1]
				args = args[2:]
			case "-mathml":
				if len(args) < 2 {
					fmt.Fprintf(os.Stderr, "\nERROR: MathML argument is missing\n")
					os.Exit(1)
				}
				// mathmlPolicy = args[1]
				args = args[2:]
		*/

		// also allow setting additional processing flags within -format (undocumented)
		case "-combine", "-combined":
			// explicit flag to remove internal top-level tags, replaces -separate
			doCombine = true
			args = args[1:]
		case "-separate", "-separated":
			// deprecated, leaves default behavior of not combining internal top-level objects from multiple chunked reads
			args = args[1:]
		case "-self", "-self-closing":
			doSelf = true
			args = args[1:]
		case "-comment":
			doComment = true
			args = args[1:]
		case "-cdata":
			doCdata = true
			args = args[1:]
		default:
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized option after -format command\n")
			os.Exit(1)
		}
	}

	tknq := eutils.CreateTokenizer(rdr)

	frgs := eutils.FormatArgs{
		Format: format, XML: xml, Doctype: doctype,
		Combine: doCombine, Self: doSelf,
		Comment: doComment, Cdata: doCdata}

	frm := eutils.FormatTokens(tknq, frgs)

	eutils.ChanToStdout(frm)
}

// processTokens shows individual tokens in stream (undocumented)
func processTokens(rdr <-chan eutils.XMLBlock) {

	if rdr == nil {
		return
	}

	tknq := eutils.CreateTokenizer(rdr)

	if tknq == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create debug tokenizer\n")
		os.Exit(1)
	}

	var buffer strings.Builder

	count := 0
	indent := 0

	for tkn := range tknq {

		tag := tkn.Tag
		name := tkn.Name
		attr := tkn.Attr

		switch tag {
		case eutils.STARTTAG:
			buffer.WriteString("ST: ")
			for i := 0; i < indent; i++ {
				buffer.WriteString("  ")
			}
			buffer.WriteString(name)
			buffer.WriteString("\n")
			if attr != "" {
				buffer.WriteString("AT: ")
				for i := 0; i < indent; i++ {
					buffer.WriteString("  ")
				}
				buffer.WriteString(attr)
				buffer.WriteString("\n")
			}
			indent++
		case eutils.SELFTAG:
			buffer.WriteString("SL: ")
			for i := 0; i < indent; i++ {
				buffer.WriteString("  ")
			}
			buffer.WriteString(name)
			buffer.WriteString("/")
			buffer.WriteString("\n")
			if attr != "" {
				buffer.WriteString("AT: ")
				for i := 0; i < indent; i++ {
					buffer.WriteString("  ")
				}
				buffer.WriteString(attr)
				buffer.WriteString("\n")
			}
		case eutils.STOPTAG:
			indent--
			buffer.WriteString("SP: ")
			for i := 0; i < indent; i++ {
				buffer.WriteString("  ")
			}
			buffer.WriteString(name)
			buffer.WriteString("/")
			buffer.WriteString("\n")
		case eutils.CONTENTTAG:
			ctype := tkn.Cont
			if (ctype & eutils.LFTSPACE) != 0 {
				if (ctype & eutils.RGTSPACE) != 0 {
					buffer.WriteString("FL: ")
				} else {
					buffer.WriteString("LF: ")
				}
			} else if (ctype & eutils.RGTSPACE) != 0 {
				buffer.WriteString("RT: ")
			} else {
				buffer.WriteString("VL: ")
			}
			for i := 0; i < indent; i++ {
				buffer.WriteString("  ")
			}
			buffer.WriteString(name)
			buffer.WriteString("\n")
		case eutils.CDATATAG:
			buffer.WriteString("CD: ")
			for i := 0; i < indent; i++ {
				buffer.WriteString("  ")
			}
			buffer.WriteString(name)
			buffer.WriteString("\n")
		case eutils.COMMENTTAG:
			buffer.WriteString("CO: ")
			for i := 0; i < indent; i++ {
				buffer.WriteString("  ")
			}
			buffer.WriteString(name)
			buffer.WriteString("\n")
		case eutils.DOCTYPETAG:
			buffer.WriteString("DC: ")
			for i := 0; i < indent; i++ {
				buffer.WriteString("  ")
			}
			buffer.WriteString(name)
			buffer.WriteString("\n")
		case eutils.NOTAG:
			buffer.WriteString("NO:")
			if indent != 0 {
				buffer.WriteString(" (indent ")
				buffer.WriteString(strconv.Itoa(indent))
				buffer.WriteString(")")
			}
			buffer.WriteString("\n")
		case eutils.ISCLOSED:
			buffer.WriteString("CL:")
			if indent != 0 {
				buffer.WriteString(" (indent ")
				buffer.WriteString(strconv.Itoa(indent))
				buffer.WriteString(")")
			}
			buffer.WriteString("\n")
			txt := buffer.String()
			if txt != "" {
				// print final buffer
				fmt.Fprintf(os.Stdout, "%s", txt)
			}
			return
		default:
			buffer.WriteString("UNKONWN:")
			if indent != 0 {
				buffer.WriteString(" (indent ")
				buffer.WriteString(strconv.Itoa(indent))
				buffer.WriteString(")")
			}
			buffer.WriteString("\n")
		}

		count++
		if count > 1000 {
			count = 0
			txt := buffer.String()
			if txt != "" {
				// print current buffered output
				fmt.Fprintf(os.Stdout, "%s", txt)
			}
			buffer.Reset()
		}
	}
}

// processOutline displays outline of XML structure
func processOutline(rdr <-chan eutils.XMLBlock) {

	if rdr == nil {
		return
	}

	tknq := eutils.CreateTokenizer(rdr)

	if tknq == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create outline tokenizer\n")
		os.Exit(1)
	}

	var buffer strings.Builder

	count := 0
	indent := 0

	for tkn := range tknq {

		tag := tkn.Tag
		name := tkn.Name

		switch tag {
		case eutils.STARTTAG:
			if name == "eSummaryResult" ||
				name == "eLinkResult" ||
				name == "eInfoResult" ||
				name == "PubmedArticleSet" ||
				name == "DocumentSummarySet" ||
				name == "INSDSet" ||
				name == "Entrezgene-Set" ||
				name == "TaxaSet" {
				break
			}
			for i := 0; i < indent; i++ {
				buffer.WriteString("  ")
			}
			buffer.WriteString(name)
			buffer.WriteString("\n")
			indent++
		case eutils.SELFTAG:
			for i := 0; i < indent; i++ {
				buffer.WriteString("  ")
			}
			buffer.WriteString(name)
			buffer.WriteString("\n")
		case eutils.STOPTAG:
			indent--
		case eutils.DOCTYPETAG:
		case eutils.NOTAG:
		case eutils.ISCLOSED:
			txt := buffer.String()
			if txt != "" {
				// print final buffer
				fmt.Fprintf(os.Stdout, "%s", txt)
			}
			return
		default:
		}

		count++
		if count > 1000 {
			count = 0
			txt := buffer.String()
			if txt != "" {
				// print current buffered output
				fmt.Fprintf(os.Stdout, "%s", txt)
			}
			buffer.Reset()
		}
	}
}

// processSynopsis displays paths to XML elements
func processSynopsis(rdr <-chan eutils.XMLBlock, leaf bool, delim string) {

	if rdr == nil {
		return
	}

	tknq := eutils.CreateTokenizer(rdr)

	if tknq == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create synopsis tokenizer\n")
		os.Exit(1)
	}

	var buffer strings.Builder
	count := 0

	// synopsisLevel recursive definition
	var synopsisLevel func(string) bool

	synopsisLevel = func(parent string) bool {

		for tkn := range tknq {

			tag := tkn.Tag
			name := tkn.Name

			switch tag {
			case eutils.STARTTAG:
				if name == "eSummaryResult" ||
					name == "eLinkResult" ||
					name == "eInfoResult" ||
					name == "PubmedArticleSet" ||
					name == "DocumentSummarySet" ||
					name == "INSDSet" ||
					name == "Entrezgene-Set" ||
					name == "TaxaSet" {
					break
				}
				if leaf {
					if name == "root" ||
						name == "opt" ||
						name == "anon" {
						break
					}
				}
				if !leaf {
					// show all paths, including container objects
					if parent != "" {
						buffer.WriteString(parent)
						buffer.WriteString(delim)
					}
					buffer.WriteString(name)
					buffer.WriteString("\n")
				}
				path := parent
				if path != "" {
					path += delim
				}
				path += name
				if synopsisLevel(path) {
					return true
				}
			case eutils.SELFTAG:
				if parent != "" {
					buffer.WriteString(parent)
					buffer.WriteString(delim)
				}
				buffer.WriteString(name)
				buffer.WriteString("\n")
			case eutils.STOPTAG:
				// break recursion
				return false
			case eutils.CONTENTTAG:
				if leaf {
					// only show endpoint paths
					if parent != "" {
						buffer.WriteString(parent)
						buffer.WriteString("\n")
					}
				}
			case eutils.DOCTYPETAG:
			case eutils.NOTAG:
			case eutils.ISCLOSED:
				txt := buffer.String()
				if txt != "" {
					// print final buffer
					fmt.Fprintf(os.Stdout, "%s", txt)
				}
				return true
			default:
			}

			count++
			if count > 1000 {
				count = 0
				txt := buffer.String()
				if txt != "" {
					// print current buffered output
					fmt.Fprintf(os.Stdout, "%s", txt)
				}
				buffer.Reset()
			}
		}
		return true
	}

	for {
		// may have concatenated XMLs, loop through all
		if synopsisLevel("") {
			return
		}
	}
}

// processFilter modifies XML content, comments, or CDATA
func processFilter(rdr <-chan eutils.XMLBlock, args []string) {

	if rdr == nil || args == nil {
		return
	}

	tknq := eutils.CreateTokenizer(rdr)

	if tknq == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create filter tokenizer\n")
		os.Exit(1)
	}

	var buffer strings.Builder

	count := 0

	// skip past command name
	args = args[1:]

	max := len(args)
	if max < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to transmute -filter\n")
		os.Exit(1)
	}

	pttrn := args[0]

	args = args[1:]
	max--

	if max < 2 {
		fmt.Fprintf(os.Stderr, "\nERROR: No object name supplied to transmute -filter\n")
		os.Exit(1)
	}

	type ActionType int

	const (
		NOACTION ActionType = iota
		DORETAIN
		DOREMOVE
		DOENCODE
		DODECODE
		DOSHRINK
		DOEXPAND
		DOACCENT
	)

	action := args[0]

	what := NOACTION
	switch action {
	case "retain":
		what = DORETAIN
	case "remove":
		what = DOREMOVE
	case "encode":
		what = DOENCODE
	case "decode":
		what = DODECODE
	case "shrink":
		what = DOSHRINK
	case "expand":
		what = DOEXPAND
	case "accent":
		what = DOACCENT
	default:
		fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized action '%s' supplied to transmute -filter\n", action)
		os.Exit(1)
	}

	trget := args[1]

	which := eutils.NOTAG
	switch trget {
	case "attribute", "attributes":
		which = eutils.ATTRIBTAG
	case "content", "contents":
		which = eutils.CONTENTTAG
	case "cdata", "CDATA":
		which = eutils.CDATATAG
	case "comment", "comments":
		which = eutils.COMMENTTAG
	case "object":
		// object normally retained
		which = eutils.OBJECTTAG
	case "container":
		which = eutils.CONTAINERTAG
	default:
		fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized target '%s' supplied to transmute -filter\n", trget)
		os.Exit(1)
	}

	inPattern := false
	prevName := ""

	for tkn := range tknq {

		tag := tkn.Tag
		name := tkn.Name
		attr := tkn.Attr

		switch tag {
		case eutils.STARTTAG:
			prevName = name
			if name == pttrn {
				inPattern = true
				if which == eutils.CONTAINERTAG && what == DOREMOVE {
					continue
				}
			}
			if inPattern && which == eutils.OBJECTTAG && what == DOREMOVE {
				continue
			}
			buffer.WriteString("<")
			buffer.WriteString(name)
			if attr != "" {
				if which != eutils.ATTRIBTAG || what != DOREMOVE {
					attr = strings.TrimSpace(attr)
					attr = eutils.CompressRunsOfSpaces(attr)
					buffer.WriteString(" ")
					buffer.WriteString(attr)
				}
			}
			buffer.WriteString(">\n")
		case eutils.SELFTAG:
			if inPattern && which == eutils.OBJECTTAG && what == DOREMOVE {
				continue
			}
			buffer.WriteString("<")
			buffer.WriteString(name)
			if attr != "" {
				if which != eutils.ATTRIBTAG || what != DOREMOVE {
					attr = strings.TrimSpace(attr)
					attr = eutils.CompressRunsOfSpaces(attr)
					buffer.WriteString(" ")
					buffer.WriteString(attr)
				}
			}
			buffer.WriteString("/>\n")
		case eutils.STOPTAG:
			if name == pttrn {
				inPattern = false
				if which == eutils.OBJECTTAG && what == DOREMOVE {
					continue
				}
				if which == eutils.CONTAINERTAG && what == DOREMOVE {
					continue
				}
			}
			if inPattern && which == eutils.OBJECTTAG && what == DOREMOVE {
				continue
			}
			buffer.WriteString("</")
			buffer.WriteString(name)
			buffer.WriteString(">\n")
		case eutils.CONTENTTAG:
			if inPattern && which == eutils.OBJECTTAG && what == DOREMOVE {
				continue
			}
			if inPattern && which == eutils.CONTENTTAG && what == DOEXPAND {
				var words []string
				if strings.Contains(name, "|") {
					words = strings.FieldsFunc(name, func(c rune) bool {
						return c == '|'
					})
				} else if strings.Contains(name, ",") {
					words = strings.FieldsFunc(name, func(c rune) bool {
						return c == ','
					})
				} else {
					words = strings.Fields(name)
				}
				between := ""
				for _, item := range words {
					max := len(item)
					for max > 1 {
						ch := item[max-1]
						if ch != '.' && ch != ',' && ch != ':' && ch != ';' {
							break
						}
						// trim trailing punctuation
						item = item[:max-1]
						// continue checking for runs of punctuation at end
						max--
					}
					if eutils.HasFlankingSpace(item) {
						item = strings.TrimSpace(item)
					}
					if item != "" {
						if between != "" {
							buffer.WriteString(between)
						}
						buffer.WriteString(item)
						buffer.WriteString("\n")
						between = "</" + prevName + ">\n<" + prevName + ">\n"
					}
				}
				continue
			}
			if inPattern && which == tag {
				switch what {
				case DORETAIN:
					// default behavior for content - can use -filter X retain content as a no-op
				case DOREMOVE:
					continue
				case DOENCODE:
					name = html.EscapeString(name)
				case DODECODE:
					name = html.UnescapeString(name)
				case DOSHRINK:
					name = eutils.CompressRunsOfSpaces(name)
				case DOACCENT:
					if eutils.IsNotASCII(name) {
						name = eutils.DoAccentTransform(name)
					}
				default:
					continue
				}
			}
			// content normally printed
			if eutils.HasFlankingSpace(name) {
				name = strings.TrimSpace(name)
			}
			buffer.WriteString(name)
			buffer.WriteString("\n")
		case eutils.CDATATAG:
			if inPattern && which == eutils.OBJECTTAG && what == DOREMOVE {
				continue
			}
			if inPattern && which == tag {
				switch what {
				case DORETAIN:
					// cdata requires explicit retain command
				case DOREMOVE:
					continue
				case DOENCODE:
					name = html.EscapeString(name)
				case DODECODE:
					name = html.UnescapeString(name)
				case DOSHRINK:
					name = eutils.CompressRunsOfSpaces(name)
				case DOACCENT:
					if eutils.IsNotASCII(name) {
						name = eutils.DoAccentTransform(name)
					}
				default:
					continue
				}
				// cdata normally removed
				if eutils.HasFlankingSpace(name) {
					name = strings.TrimSpace(name)
				}
				buffer.WriteString(name)
				buffer.WriteString("\n")
			}
		case eutils.COMMENTTAG:
			if inPattern && which == eutils.OBJECTTAG && what == DOREMOVE {
				continue
			}
			if inPattern && which == tag {
				switch what {
				case DORETAIN:
					// comment requires explicit retain command
				case DOREMOVE:
					continue
				case DOENCODE:
					name = html.EscapeString(name)
				case DODECODE:
					name = html.UnescapeString(name)
				case DOSHRINK:
					name = eutils.CompressRunsOfSpaces(name)
				case DOACCENT:
					if eutils.IsNotASCII(name) {
						name = eutils.DoAccentTransform(name)
					}
				default:
					continue
				}
				// comment normally removed
				if eutils.HasFlankingSpace(name) {
					name = strings.TrimSpace(name)
				}
				buffer.WriteString(name)
				buffer.WriteString("\n")
			}
		case eutils.DOCTYPETAG:
		case eutils.NOTAG:
		case eutils.ISCLOSED:
			txt := buffer.String()
			if txt != "" {
				// print final buffer
				fmt.Fprintf(os.Stdout, "%s", txt)
			}
			return
		default:
		}

		count++
		if count > 1000 {
			count = 0
			txt := buffer.String()
			if txt != "" {
				// print current buffered output
				fmt.Fprintf(os.Stdout, "%s", txt)
			}
			buffer.Reset()
		}
	}
}

// STRING CONVERTERS

func encodeURL(inp io.Reader) {

	if inp == nil {
		return
	}

	data, _ := ioutil.ReadAll(inp)
	txt := string(data)
	txt = strings.TrimSuffix(txt, "\n")

	str := url.QueryEscape(txt)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

func decodeURL(inp io.Reader) {

	if inp == nil {
		return
	}

	byt, _ := ioutil.ReadAll(inp)
	txt := string(byt)
	txt = strings.TrimSuffix(txt, "\n")

	str, _ := url.QueryUnescape(txt)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

func encodeB64(inp io.Reader) {

	if inp == nil {
		return
	}

	data, _ := ioutil.ReadAll(inp)

	str := base64.StdEncoding.EncodeToString(data)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

func decodeB64(inp io.Reader) {

	if inp == nil {
		return
	}

	byt, _ := ioutil.ReadAll(inp)

	data, _ := base64.StdEncoding.DecodeString(string(byt))
	str := string(data)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

func decodeHGVS(inp io.Reader) {

	if inp == nil {
		return
	}

	byt, _ := ioutil.ReadAll(inp)
	txt := string(byt)

	os.Stdout.WriteString("<HGVS>\n")

	str := eutils.ParseHGVS(txt)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}

	os.Stdout.WriteString("</HGVS>\n")
}

// COLUMN ALIGNMENT FORMATTER

// processAlign aligns a tab-delimited table by individual column widths
func processAlign(inp io.Reader, args []string) {

	// tab-delimited-table to padded-by-spaces alignment inspired by
	// Steve Kinzler's align script - see http://kinzler.com/me/align/

	if inp == nil {
		return
	}

	mrg := 0
	pdg := 0
	mnw := 0
	aln := ""

	// skip past command name
	args = args[1:]

	for len(args) > 0 {

		switch args[0] {
		case "-g":
			pdg = eutils.GetNumericArg(args, "-g spacing between columns", 0, 1, 30)
			args = args[2:]
		case "-h":
			mrg = eutils.GetNumericArg(args, "-i indent before columns", 0, 1, 30)
			args = args[2:]
		case "-w":
			mnw = eutils.GetNumericArg(args, "-w minimum column width", 0, 1, 30)
			args = args[2:]
		case "-a":
			aln = eutils.GetStringArg(args, "-a column alignment code string")
			args = args[2:]
		default:
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized option after -align command\n")
			os.Exit(1)
		}
	}

	algn := eutils.AlignColumns(inp, mrg, pdg, mnw, aln)

	if algn == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create alignment function\n")
		os.Exit(1)
	}

	eutils.ChanToStdout(algn)

	return
}

// SEQUENCE EDITING

func readOneFastaSequence(inp io.Reader) string {

	fsta := eutils.FASTAConverter(inp, false)

	// return first FASTA sequence
	for fsa := range fsta {
		return fsa.Sequence
	}

	return ""
}

func sequenceRemove(inp io.Reader, args []string) {

	if inp == nil {
		return
	}

	first := ""
	last := ""

	// skip past command name
	args = args[1:]

	for len(args) > 0 {

		switch args[0] {
		case "-first":
			first = eutils.GetStringArg(args, "Bases to delete at beginning")
			first = strings.ToUpper(first)
			args = args[2:]
		case "-last":
			last = eutils.GetStringArg(args, "Bases to delete at end")
			last = strings.ToUpper(last)
			args = args[2:]
		default:
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized option after -remove command\n")
			os.Exit(1)
		}
	}

	str := readOneFastaSequence(inp)

	str = eutils.SequenceRemove(str, first, last)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

func sequenceRetain(inp io.Reader, args []string) {

	if inp == nil {
		return
	}

	lead := 0
	trail := 0

	// skip past command name
	args = args[1:]

	for len(args) > 0 {

		switch args[0] {
		case "-leading":
			lead = eutils.GetNumericArg(args, "Bases to keep at beginning", 0, -1, -1)
			args = args[2:]
		case "-trailing":
			trail = eutils.GetNumericArg(args, "Bases to keep at end", 0, -1, -1)
			args = args[2:]
		default:
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized option after -retain command\n")
			os.Exit(1)
		}
	}

	str := readOneFastaSequence(inp)

	str = eutils.SequenceRetain(str, lead, trail)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

func sequenceReplace(inp io.Reader, args []string) {

	if inp == nil {
		return
	}

	pos := 0
	del := ""
	ins := ""
	lower := false

	// skip past command name
	args = args[1:]

	for len(args) > 0 {

		switch args[0] {
		case "-offset":
			pos = eutils.GetNumericArg(args, "0-based position", 0, -1, -1)
			args = args[2:]
		case "-column":
			val := eutils.GetNumericArg(args, "1-based position", 1, -1, -1)
			pos = val - 1
			args = args[2:]
		case "-delete":
			del = eutils.GetStringArg(args, "Number to delete")
			del = strings.ToUpper(del)
			args = args[2:]
		case "-insert":
			ins = eutils.GetStringArg(args, "Bases to insert")
			ins = strings.ToUpper(ins)
			args = args[2:]
		case "-lower":
			lower = true
			args = args[1:]
		default:
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized option after -replace command\n")
			os.Exit(1)
		}
	}

	str := readOneFastaSequence(inp)

	if lower {
		str = strings.ToLower(str)
	}

	str = eutils.SequenceReplace(str, pos, del, ins)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

func sequenceExtract(inp io.Reader, args []string) {

	if inp == nil {
		return
	}

	featLoc := ""
	isOneBased := true
	lower := false

	// skip past command name
	args = args[1:]

	for len(args) > 0 {

		switch args[0] {
		case "-0-based":
			isOneBased = false
			args = args[1:]
		case "-1-based":
			isOneBased = true
			args = args[1:]
		case "-lower":
			lower = true
			args = args[1:]
		default:
			// read output of xtract -insd feat_location qualifier
			featLoc = args[0]
			args = args[1:]
		}
	}

	if featLoc == "" {
		fmt.Fprintf(os.Stderr, "\nERROR: Missing argument after -extract command\n")
		os.Exit(1)
	}

	str := readOneFastaSequence(inp)

	str = eutils.SequenceExtract(str, featLoc, isOneBased)

	if lower {
		str = strings.ToLower(str)
	}

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

func sequenceSearch(inp io.Reader, args []string) {

	if inp == nil {
		return
	}

	// skip past command name
	args = args[1:]

	protein := false
	circular := false
	topStrand := false

	for len(args) > 0 {
		if args[0] == "-protein" {
			protein = true
			args = args[1:]
		} else if args[0] == "-circular" {
			circular = true
			args = args[1:]
		} else if args[0] == "-top" {
			topStrand = true
			args = args[1:]
		} else {
			break
		}
	}

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Missing argument after -search command\n")
		os.Exit(1)
	}

	var arry []string

	// allow one or more patterns to be passed in each argument
	for len(args) > 0 {
		pat := args[0]
		args = args[1:]

		pat = strings.TrimSpace(pat)
		itms := strings.Split(pat, " ")
		for _, trm := range itms {
			arry = append(arry, trm)
		}
	}

	str := readOneFastaSequence(inp)

	srch := eutils.SequenceSearcher(arry, protein, circular, topStrand)

	res := srch.Search(str)

	txt := ""
	for _, hit := range res {
		txt = fmt.Sprintf("%d\t%s\n", hit.Point, hit.Match)
		os.Stdout.WriteString(txt)
	}
	if !strings.HasSuffix(txt, "\n") {
		os.Stdout.WriteString("\n")
	}
}

func readAllIntoString(inp io.Reader) string {

	if inp == nil {
		return ""
	}

	data, _ := ioutil.ReadAll(inp)
	txt := string(data)

	if txt == "" {
		return ""
	}

	// replace whitespace substrings with a single space
	re := regexp.MustCompile(`\s+`)
	txt = re.ReplaceAllString(txt, " ")

	txt = strings.TrimSpace(txt)

	return txt
}

func stringFind(inp io.Reader, args []string) {

	if inp == nil {
		return
	}

	// skip past command name
	args = args[1:]

	caseSensitive := false
	relaxed := false

	for len(args) > 0 {
		if args[0] == "-sensitive" {
			caseSensitive = true
			args = args[1:]
		} else if args[0] == "-relaxed" {
			relaxed = true
			args = args[1:]
		} else {
			break
		}
	}

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Missing argument after -find command\n")
		os.Exit(1)
	}

	str := readAllIntoString(inp)

	srch := eutils.StringSearcher(args, caseSensitive, relaxed)

	res := srch.Search(str)

	txt := ""
	for _, hit := range res {
		txt = fmt.Sprintf("%d\t%s\n", hit.Point, hit.Match)
		os.Stdout.WriteString(txt)
	}
	if !strings.HasSuffix(txt, "\n") {
		os.Stdout.WriteString("\n")
	}
}

func relaxString(inp io.Reader) {

	str := readAllIntoString(inp)

	str = eutils.RelaxString(str)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

func upperString(inp io.Reader) {

	str := readAllIntoString(inp)

	str = strings.ToUpper(str)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

func lowerString(inp io.Reader) {

	str := readAllIntoString(inp)

	str = strings.ToLower(str)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

// FASTA BASE COUNT

// baseCount prints a summary of base or residue counts
func baseCount(inp io.Reader) {

	if inp == nil {
		return
	}

	fsta := eutils.FASTAConverter(inp, false)

	countLetters := func(id, seq string) {

		counts := make(map[rune]int)

		for _, base := range seq {
			counts[base]++
		}

		var keys []rune
		for ky := range counts {
			keys = append(keys, ky)
		}
		sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })

		fmt.Fprintf(os.Stdout, "%s", id)
		for _, base := range keys {
			num := counts[base]
			fmt.Fprintf(os.Stdout, "\t%c %d", base, num)
		}
		fmt.Fprintf(os.Stdout, "\n")
	}

	for fsa := range fsta {
		countLetters(fsa.SeqID, fsa.Sequence)
	}
}

// REVERSE SEQUENCE

// seqFlip reverses without complementing - e.g., minus strand proteins translated in reverse order
func seqFlip(inp io.Reader) {

	if inp == nil {
		return
	}

	str := readOneFastaSequence(inp)

	str = eutils.SequenceReverse(str)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

// REVERSE COMPLEMENT

func nucRevComp(inp io.Reader) {

	if inp == nil {
		return
	}

	str := readOneFastaSequence(inp)

	str = eutils.ReverseComplement(str)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

/*
func fastaRevComp(inp io.Reader) {

	if inp == nil {
		return
	}

	fsta := eutils.FASTAConverter(inp, false)

	for fsa := range fsta {

		str := fsa.Sequence

		str = eutils.ReverseComplement(str)

		os.Stdout.WriteString(">")
		if fsa.SeqID != "" {
			os.Stdout.WriteString(fsa.SeqID)
			if fsa.Title != "" {
				os.Stdout.WriteString(" ")
				os.Stdout.WriteString(fsa.Title)
			}
		}
		os.Stdout.WriteString("\n")

		os.Stdout.WriteString(str)
		if !strings.HasSuffix(str, "\n") {
			os.Stdout.WriteString("\n")
		}
	}
}
*/

// FASTA DIFFERENCES

func printFastaPairs(frst, scnd string) {

	frst = strings.ToLower(frst)
	scnd = strings.ToLower(scnd)

	fst := frst[:]
	scd := scnd[:]

	// next functions return spaces after end of sequence
	nextF := func() rune {

		if len(fst) < 1 {
			return ' '
		}
		ch := fst[0]
		fst = fst[1:]

		return rune(ch)
	}

	nextS := func() rune {

		if len(scd) < 1 {
			return ' '
		}
		ch := scd[0]
		scd = scd[1:]

		return rune(ch)
	}

	var fs []rune
	var sc []rune
	mx := 0

	// populate output arrays
	for {

		f, s := nextF(), nextS()
		// if both spaces, end of both sequences
		if f == ' ' && s == ' ' {
			break
		}
		if f == s {
			fs = append(fs, f)
			sc = append(sc, ' ')
		} else {
			// show mismatches in upper case
			fs = append(fs, unicode.ToUpper(f))
			sc = append(sc, unicode.ToUpper(s))
		}
		mx++
	}

	// pad output to multiple of 50
	j := mx % 50
	if j > 0 {
		for j < 50 {
			fs = append(fs, ' ')
			sc = append(sc, ' ')
			j++
			mx++
		}
	}

	// print in blocks of 50 bases or residues
	for i := 0; i < mx; i += 50 {
		dl := 50
		if mx-i < 50 {
			dl = mx - i
		}
		lf := fs[:dl]
		rt := sc[:dl]
		fs = fs[dl:]
		sc = sc[dl:]
		tm := strings.TrimRight(string(lf), " ")
		fmt.Fprintf(os.Stdout, "%s %6d\n%s\n", string(lf), i+len(tm), string(rt))
	}
}

func fastaDiff(inp io.Reader, args []string) {

	if inp == nil {
		return
	}

	// skip past command name
	args = args[1:]

	if len(args) != 2 {
		fmt.Fprintf(os.Stderr, "\nERROR: Two files required by -diff command\n")
		os.Exit(1)
	}

	frst := args[0]
	scnd := args[1]

	readSeqFromFile := func(fname string) string {

		f, err := os.Open(fname)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Unable to open file %s - %s\n", fname, err.Error())
			os.Exit(1)
		}

		defer f.Close()

		seq := readOneFastaSequence(f)

		return seq
	}

	frstFasta := readSeqFromFile(frst)
	scndFasta := readSeqFromFile(scnd)

	if frstFasta == scndFasta {
		return
	}

	// sequences are assumed to be aligned, this code highlight mismatches
	printFastaPairs(frstFasta, scndFasta)
}

// PROTEIN WEIGHT

func protWeight(inp io.Reader, args []string) {

	if inp == nil {
		return
	}

	trimLeadingMet := true

	// skip past command name
	args = args[1:]

	for len(args) > 0 {

		switch args[0] {
		case "-met":
			trimLeadingMet = false
			args = args[1:]
		default:
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized option after -molwt command\n")
			os.Exit(1)
		}
	}

	str := readOneFastaSequence(inp)

	str = eutils.ProteinWeight(str, trimLeadingMet)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

// cdRegionToProtein reads all of stdin as sequence data
func cdRegionToProtein(inp io.Reader, args []string) {

	if inp == nil {
		return
	}

	genCode := 1
	frame := 0
	includeStop := false
	doEveryCodon := false
	removeTrailingX := false
	is5primeComplete := true
	is3primeComplete := true
	between := ""

	repeat := 1

	// skip past command name
	args = args[1:]

	for len(args) > 0 {

		switch args[0] {
		case "-code", "-gencode":
			genCode = eutils.GetNumericArg(args, "genetic code number", 0, 1, 30)
			args = args[2:]
		case "-frame":
			frame = eutils.GetNumericArg(args, "offset into coding sequence", 0, 1, 30)
			args = args[2:]
		case "-stop", "-stops":
			includeStop = true
			args = args[1:]
		case "-every", "-all":
			doEveryCodon = true
			args = args[1:]
		case "-trim", "-trailing":
			removeTrailingX = true
			args = args[1:]
		case "-part5", "-partial5", "-lt5":
			is5primeComplete = false
			args = args[1:]
		case "-part3", "-partial3", "-gt3":
			is3primeComplete = false
			args = args[1:]
		case "-between":
			between = eutils.GetStringArg(args, "separator between residues")
			args = args[2:]
		case "-repeat":
			repeat = eutils.GetNumericArg(args, "number of repetitions for testing", 1, 1, 100)
			args = args[2:]
		case "-":
			// lone dash is default for -every -trim
			doEveryCodon = true
			removeTrailingX = true
			args = args[1:]
		default:
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized option after -cds2prot command\n")
			os.Exit(1)
		}
	}

	txt := readOneFastaSequence(inp)

	for i := 0; i < repeat; i++ {

		// repeat multiple times for performance testing (undocumented)
		str := eutils.TranslateCdRegion(txt, genCode, frame, includeStop, doEveryCodon, removeTrailingX, is5primeComplete, is3primeComplete, between)

		os.Stdout.WriteString(str)
		if !strings.HasSuffix(str, "\n") {
			os.Stdout.WriteString("\n")
		}
	}
}

// nucProtCodonReport prints amino acid residues under nucleotide codons
func nucProtCodonReport(args []string) {

	nuc := ""
	prt := ""
	frame := 0
	threeLetter := false

	// skip past command name
	args = args[1:]

	for len(args) > 0 {

		switch args[0] {
		case "-nuc":
			nuc = eutils.GetStringArg(args, "separator between residues")
			args = args[2:]
		case "-prt":
			prt = eutils.GetStringArg(args, "separator between residues")
			args = args[2:]
		case "-frame":
			frame = eutils.GetNumericArg(args, "offset into coding sequence", 0, 1, 30)
			args = args[2:]
		case "-three", "-triple", "-triples", "-triplet", "-triplets":
			threeLetter = true
			args = args[1:]
		default:
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized option after -cds2prot command\n")
			os.Exit(1)
		}
	}

	str := eutils.NucProtCodonReport(nuc, prt, frame, threeLetter)

	os.Stdout.WriteString(str)
	if !strings.HasSuffix(str, "\n") {
		os.Stdout.WriteString("\n")
	}
}

// MAIN FUNCTION

func main() {

	// skip past executable name
	args := os.Args[1:]

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: No command-line arguments supplied to transmute\n")
		os.Exit(1)
	}

	// performance arguments
	chanDepth := 0
	farmSize := 0
	heapSize := 0
	numServe := 0
	goGc := 0

	// processing option arguments
	doCompress := false
	doCleanup := false
	doStrict := false
	doMixed := false
	doSelf := false
	deAccent := false
	doASCII := false

	/*
		doUnicode := false
		doScript := false
		doMathML := false
	*/

	// CONCURRENCY, CLEANUP, AND DEBUGGING FLAGS

	// do these first because -defcpu and -maxcpu can be sent from wrapper before other arguments

	ncpu := runtime.NumCPU()
	if ncpu < 1 {
		ncpu = 1
	}

	// wrapper can limit maximum number of processors to use (undocumented)
	maxProcs := ncpu
	defProcs := 0

	// concurrent performance tuning parameters, can be overridden by -proc and -cons
	numProcs := 0
	serverRatio := 4

	// -flag sets -strict or -mixed cleanup flags from argument
	flgs := ""

	/*
		unicodePolicy := ""
		scriptPolicy := ""
		mathmlPolicy := ""
	*/

	// read data from file instead of stdin
	fileName := ""

	// debugging
	stts := false
	timr := false

	// profiling
	prfl := false

	inSwitch := true

	// get concurrency, cleanup, and debugging flags in any order
	for {

		inSwitch = true

		switch args[0] {

		// concurrency override arguments can be passed in by local wrapper script (undocumented)
		case "-maxcpu":
			maxProcs = eutils.GetNumericArg(args, "Maximum number of processors", 1, 1, ncpu)
			args = args[1:]
		case "-defcpu":
			defProcs = eutils.GetNumericArg(args, "Default number of processors", ncpu, 1, ncpu)
			args = args[1:]
		// performance tuning flags
		case "-proc":
			numProcs = eutils.GetNumericArg(args, "Number of processors", ncpu, 1, ncpu)
			args = args[1:]
		case "-cons":
			serverRatio = eutils.GetNumericArg(args, "Parser to processor ratio", 4, 1, 32)
			args = args[1:]
		case "-serv":
			numServe = eutils.GetNumericArg(args, "Concurrent parser count", 0, 1, 128)
			args = args[1:]
		case "-chan":
			chanDepth = eutils.GetNumericArg(args, "Communication channel depth", 0, ncpu, 128)
			args = args[1:]
		case "-heap":
			heapSize = eutils.GetNumericArg(args, "Unshuffler heap size", 8, 8, 64)
			args = args[1:]
		case "-farm":
			farmSize = eutils.GetNumericArg(args, "Node buffer length", 4, 4, 2048)
			args = args[1:]
		case "-gogc":
			goGc = eutils.GetNumericArg(args, "Garbage collection percentage", 0, 50, 1000)
			args = args[1:]

		// read data from file
		case "-input":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Input file name is missing\n")
				os.Exit(1)
			}
			fileName = args[1]
			// skip past first of two arguments
			args = args[1:]

		// data cleanup flags
		case "-compress", "-compressed":
			doCompress = true
		case "-spaces", "-cleanup":
			doCleanup = true
		case "-strict":
			doStrict = true
		case "-mixed":
			doMixed = true
		case "-self":
			doSelf = true
		case "-accent":
			deAccent = true
		case "-ascii":
			doASCII = true

		// previously visible processing flags (undocumented)
		case "-stems", "-stem":
			// ignore
		case "-stops", "-stop":
			// ignore

		// allow setting of unicode, script, and mathml flags (undocumented)
		case "-unicode":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: -unicode argument is missing\n")
				os.Exit(1)
			}
			// unicodePolicy = eutils.GetStringArg(args, "Unicode argument")
			args = args[1:]
		case "-script":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: -script argument is missing\n")
				os.Exit(1)
			}
			// scriptPolicy = eutils.GetStringArg(args, "Script argument")
			args = args[1:]
		case "-mathml":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: -mathml argument is missing\n")
				os.Exit(1)
			}
			// mathmlPolicy = eutils.GetStringArg(args, "MathML argument")
			args = args[1:]

		case "-flag", "-flags":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: -flags argument is missing\n")
				os.Exit(1)
			}
			flgs = eutils.GetStringArg(args, "Flags argument")
			args = args[1:]

		// debugging flags
		case "-stats", "-stat":
			stts = true
		case "-timer":
			timr = true
		case "-profile":
			prfl = true

		default:
			// if not any of the controls, set flag to break out of for loop
			inSwitch = false
		}

		if !inSwitch {
			break
		}

		// skip past argument
		args = args[1:]

		if len(args) < 1 {
			break
		}
	}

	// -flag allows script to set -strict or -mixed (or -stems, or -stops) from argument
	switch flgs {
	case "strict":
		doStrict = true
	case "mixed":
		doMixed = true
	case "stems", "stem":
		// ignore
	case "stops", "stop":
		// ignore
	case "none", "default":
	default:
		if flgs != "" {
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized -flag value '%s'\n", flgs)
			os.Exit(1)
		}
	}

	/*
		UnicodeFix = ParseMarkup(unicodePolicy, "-unicode")
		ScriptFix = ParseMarkup(scriptPolicy, "-script")
		MathMLFix = ParseMarkup(mathmlPolicy, "-mathml")

		if UnicodeFix != NOMARKUP {
			doUnicode = true
		}

		if ScriptFix != NOMARKUP {
			doScript = true
		}

		if MathMLFix != NOMARKUP {
			doMathML = true
		}
	*/

	if numProcs == 0 {
		if defProcs > 0 {
			numProcs = defProcs
		} else if maxProcs > 0 {
			numProcs = maxProcs
		}
	}
	if numProcs > ncpu {
		numProcs = ncpu
	}
	if numProcs > maxProcs {
		numProcs = maxProcs
	}

	eutils.SetTunings(numProcs, numServe, serverRatio, chanDepth, farmSize, heapSize, goGc, false)

	eutils.SetOptions(doStrict, doMixed, doSelf, deAccent, doASCII, doCompress, doCleanup)

	// -stats prints number of CPUs and performance tuning values if no other arguments (undocumented)
	if stts && len(args) < 1 {

		eutils.PrintStats()

		return
	}

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to transmute\n")
		os.Exit(1)
	}

	// DOCUMENTATION COMMANDS

	inSwitch = true

	switch args[0] {
	case "-version":
		fmt.Printf("%s\n", eutils.EDirectVersion)
	case "-help", "help":
		eutils.PrintHelp("transmute", "transmute-help.txt")
	case "-extra", "-extras":
		eutils.PrintHelp("transmute", "transmute-extras.txt")
	case "-degenerate":
		// generate new genetic code data tables (undocumented)
		eutils.GenerateGeneticCodeMaps()
	default:
		// if not any of the documentation commands, keep going
		inSwitch = false
	}

	if inSwitch {
		return
	}

	// FILE NAME CAN BE SUPPLIED WITH -input COMMAND

	in := os.Stdin

	// check for data being piped into stdin
	isPipe := false
	fi, err := os.Stdin.Stat()
	if err == nil {
		isPipe = bool((fi.Mode() & os.ModeNamedPipe) != 0)
	}

	usingFile := false

	if fileName != "" {

		inFile, err := os.Open(fileName)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
			os.Exit(1)
		}

		defer inFile.Close()

		// use indicated file instead of stdin
		in = inFile
		usingFile = true

		if isPipe && runtime.GOOS != "windows" {
			mode := fi.Mode().String()
			fmt.Fprintf(os.Stderr, "\nERROR: Input data from both stdin and file '%s', mode is '%s'\n", fileName, mode)
			os.Exit(1)
		}
	}

	// check for -input command after extraction arguments
	for _, str := range args {
		if str == "-input" {
			fmt.Fprintf(os.Stderr, "\nERROR: Misplaced -input command\n")
			os.Exit(1)
		}
	}

	// START PROFILING IF REQUESTED

	if prfl {

		f, err := os.Create("cpu.pprof")
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create profile output file\n")
			os.Exit(1)
		}

		pprof.StartCPUProfile(f)

		defer pprof.StopCPUProfile()
	}

	// INITIALIZE RECORD COUNT

	recordCount := 0
	byteCount := 0

	// print processing rate and program duration
	printDuration := func(name string) {

		eutils.PrintDuration(name, recordCount, byteCount)
	}

	nextArg := func() (string, bool) {

		if len(args) < 1 {
			return "", false
		}

		// remove next token from slice
		nxt := args[0]
		args = args[1:]

		return nxt, true
	}

	// The several converter functions that follow must be called
	// before CreateXMLStreamer starts draining stdin

	// JSON TO XML CONVERTER

	if args[0] == "-j2x" || args[0] == "-json2xml" {

		// skip past command name
		args = args[1:]

		set := "root"
		rec := ""
		nest := "element"

		// look for optional arguments
		for {
			arg, ok := nextArg()
			if !ok {
				break
			}

			switch arg {
			case "-set":
				// override set wrapper
				set, ok = nextArg()
				if ok && set == "-" {
					set = ""
				}
			case "-rec":
				// override record wrapper
				rec, ok = nextArg()
				if ok && rec == "-" {
					rec = ""
				}
			case "-nest":
				// specify nested array naming policy
				nest, ok = nextArg()
				if !ok {
					fmt.Fprintf(os.Stderr, "Nested array naming policy is missing\n")
					os.Exit(1)
				}
				if ok && nest == "-" {
					nest = "flat"
				}
				switch nest {
				case "flat", "plural", "name", "recurse", "recursive", "same", "depth", "deep", "level", "element", "elem", "_E":
				default:
					fmt.Fprintf(os.Stderr, "Unrecognized nested array naming policy\n")
					os.Exit(1)
				}
			default:
				// alternative form uses positional arguments to override set and rec
				set = arg
				if set == "-" {
					set = ""
				}
				rec, ok = nextArg()
				if ok && rec == "-" {
					rec = ""
				}
			}
		}

		// use output channel of tokenizer as input channel of converter
		jcnv := eutils.JSONConverter(in, set, rec, nest)

		if jcnv == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create JSON to XML converter\n")
			os.Exit(1)
		}

		// drain output of channel
		for str := range jcnv {

			if str == "" {
				continue
			}

			// send result to output
			os.Stdout.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				os.Stdout.WriteString("\n")
			}

			recordCount++
			runtime.Gosched()
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("blocks")
		}

		return
	}

	// ASN.1 TO XML CONVERTER

	if args[0] == "-a2x" || args[0] == "-asn2xml" {

		// skip past command name
		args = args[1:]

		set := ""
		rec := ""

		// look for optional arguments
		for {
			arg, ok := nextArg()
			if !ok {
				break
			}

			switch arg {
			case "-set":
				// override set wrapper
				set, ok = nextArg()
				if ok && set == "-" {
					set = ""
				}
			case "-rec":
				// override record wrapper
				rec, ok = nextArg()
				if ok && rec == "-" {
					rec = ""
				}
			}
		}

		acnv := eutils.ASN1Converter(in, set, rec)

		if acnv == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create ASN.1 to XML converter\n")
			os.Exit(1)
		}

		// drain output of channel
		for str := range acnv {

			if str == "" {
				continue
			}

			// send result to output
			os.Stdout.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				os.Stdout.WriteString("\n")
			}

			recordCount++
			runtime.Gosched()
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("blocks")
		}

		return
	}

	// READ TAB-DELIMITED FILE AND WRAP IN XML FIELDS

	doTable := func(delim string) {

		// skip past command name
		args = args[1:]

		set := ""
		rec := ""

		skip := 0
		header := false
		lower := false
		upper := false
		indent := true

		var fields []string
		numFlds := 0

		for len(args) > 0 {
			str := args[0]
			switch str {
			case "-set":
				args = args[1:]
				if len(args) < 1 {
					fmt.Fprintf(os.Stderr, "\nERROR: No argument after -set\n")
					os.Exit(1)
				}
				set = args[0]
				args = args[1:]
			case "-rec":
				args = args[1:]
				if len(args) < 1 {
					fmt.Fprintf(os.Stderr, "\nERROR: No argument after -rec\n")
					os.Exit(1)
				}
				rec = args[0]
				args = args[1:]
			case "-skip":
				args = args[1:]
				if len(args) < 1 {
					fmt.Fprintf(os.Stderr, "\nERROR: No argument after -skip\n")
					os.Exit(1)
				}
				tmp := args[0]
				val, err := strconv.Atoi(tmp)
				if err != nil {
					fmt.Fprintf(os.Stderr, "\nERROR: -skip argument (%s) is not an integer\n", tmp)
					os.Exit(1)
				}
				skip = val
				args = args[1:]
			case "-header", "-headers", "-heading":
				header = true
				args = args[1:]
			case "-lower":
				lower = true
				args = args[1:]
			case "-upper":
				upper = true
				args = args[1:]
			case "-indent":
				indent = true
				args = args[1:]
			case "-flush":
				indent = false
				args = args[1:]
			default:
				// remaining arguments are names for columns
				if str != "" && str != "*" {
					fields = append(fields, str)
					numFlds++
				}
				args = args[1:]
			}
		}

		if numFlds < 1 && !header {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient arguments for table converter\n")
			os.Exit(1)
		}

		tble := eutils.TableConverter(in, delim, set, rec, skip, header, lower, upper, indent, fields)

		if tble == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create table to XML converter\n")
			os.Exit(1)
		}

		// drain output of channel
		for str := range tble {

			if str == "" {
				continue
			}

			// send result to output
			os.Stdout.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				os.Stdout.WriteString("\n")
			}

			recordCount++
			runtime.Gosched()
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("lines")
		}
	}

	if len(args) > 1 && args[0] == "-t2x" {

		doTable("\t")
		return
	}

	if len(args) > 1 && args[0] == "-c2x" {

		doTable(",")
		return
	}

	// READ GENBANK FLATFILE AND TRANSLATE TO INSDSEQ XML

	if len(args) > 0 && args[0] == "-g2x" {

		gbk := eutils.GenBankConverter(in)

		if gbk == nil {
			fmt.Fprintf(os.Stderr, "Unable to create GenBank to XML converter\n")
			os.Exit(1)
		}

		head := `<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE INSDSet PUBLIC "-//NCBI//INSD INSDSeq/EN" "https://www.ncbi.nlm.nih.gov/dtd/INSD_INSDSeq.dtd">
<INSDSet>
`
		tail := ""

		// drain output of last channel in service chain
		for str := range gbk {

			if str == "" {
				continue
			}

			if head != "" {
				os.Stdout.WriteString(head)
				head = ""
				tail = `</INSDSet>
`
			}

			// send result to stdout
			os.Stdout.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				os.Stdout.WriteString("\n")
			}

			recordCount++

			runtime.Gosched()
		}

		if tail != "" {
			os.Stdout.WriteString(tail)
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// STRING CONVERSION COMMANDS

	inSwitch = true

	switch args[0] {
	case "-encodeURL":
		encodeURL(in)
	case "-decodeURL":
		decodeURL(in)
	case "-encode64", "-encodeB64", "-encodeBase64":
		encodeB64(in)
	case "-decode64", "-decodeB64", "-decodeBase64":
		decodeB64(in)
	case "-hgvs":
		decodeHGVS(in)
	case "-align":
		processAlign(in, args)
	case "-remove":
		sequenceRemove(in, args)
	case "-retain":
		sequenceRetain(in, args)
	case "-replace":
		sequenceReplace(in, args)
	case "-extract":
		sequenceExtract(in, args)
	case "-search":
		sequenceSearch(in, args)
	case "-find":
		stringFind(in, args)
	case "-relax":
		relaxString(in)
	case "-upper":
		upperString(in)
	case "-lower":
		lowerString(in)
	case "-counts", "-basecount":
		baseCount(in)
	case "-revcomp":
		nucRevComp(in)
	case "-reverse":
		seqFlip(in)
	case "-molwt":
		protWeight(in, args)
	case "-cds2prot":
		cdRegionToProtein(in, args)
	case "-codons":
		nucProtCodonReport(args)
	case "-diff":
		fastaDiff(in, args)
	default:
		// if not any of the conversion commands, keep going
		inSwitch = false
	}

	if inSwitch {

		debug.FreeOSMemory()

		if timr {
			printDuration("bases")
		}

		return
	}

	// CREATE XML BLOCK READER FROM STDIN OR FILE

	rdr := eutils.CreateXMLStreamer(in)
	if rdr == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
		os.Exit(1)
	}

	// CONFIRM INPUT DATA AVAILABILITY AFTER RUNNING COMMAND GENERATORS

	if fileName == "" && runtime.GOOS != "windows" {

		fromStdin := bool((fi.Mode() & os.ModeCharDevice) == 0)
		if !isPipe || !fromStdin {
			mode := fi.Mode().String()
			fmt.Fprintf(os.Stderr, "\nERROR: No data supplied to transmute from stdin or file, mode is '%s'\n", mode)
			os.Exit(1)
		}
	}

	if !usingFile && !isPipe {

		fmt.Fprintf(os.Stderr, "\nERROR: No XML input data supplied to transmute\n")
		os.Exit(1)
	}

	// SPECIAL FORMATTING COMMANDS

	inSwitch = true
	leaf := false

	switch args[0] {
	case "-format":
		processFormat(rdr, args)
	case "-filter":
		processFilter(rdr, args)
	case "-normalize", "-normal":
		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: No database supplied to -normalize\n")
			os.Exit(1)
		}
		db := args[1]
		nrm := eutils.NormalizeXML(rdr, db)
		eutils.ChanToStdout(nrm)
	case "-outline":
		processOutline(rdr)
	case "-contour":
		leaf = true
		fallthrough
	case "-synopsis":
		args = args[1:]
		delim := "/"
		if len(args) > 0 {
			delim = args[0]
			if len(delim) > 3 {
				delim = "/"
			}
		}
		processSynopsis(rdr, leaf, delim)
	case "-tokens":
		processTokens(rdr)
	default:
		// if not any of the formatting commands, keep going
		inSwitch = false
	}

	if inSwitch {

		debug.FreeOSMemory()

		// suppress printing of lines if not properly counted
		if recordCount == 1 {
			recordCount = 0
		}

		if timr {
			printDuration("lines")
		}

		return
	}

	// SPECIFY STRINGS TO GO BEFORE AND AFTER ENTIRE OUTPUT OR EACH RECORD

	head := ""
	tail := ""

	hd := ""
	tl := ""

	for {

		inSwitch = true

		switch args[0] {
		case "-head":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -head command\n")
				os.Exit(1)
			}
			head = eutils.ConvertSlash(args[1])
		case "-tail":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -tail command\n")
				os.Exit(1)
			}
			tail = eutils.ConvertSlash(args[1])
		case "-hd":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -hd command\n")
				os.Exit(1)
			}
			hd = eutils.ConvertSlash(args[1])
		case "-tl":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -tl command\n")
				os.Exit(1)
			}
			tl = eutils.ConvertSlash(args[1])
		case "-wrp":
			// shortcut to wrap records in XML tags
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -wrp command\n")
				os.Exit(1)
			}
			tmp := eutils.ConvertSlash(args[1])
			lft, rgt := eutils.SplitInTwoLeft(tmp, ",")
			if lft != "" {
				head = "<" + lft + ">"
				tail = "</" + lft + ">"
			}
			if rgt != "" {
				hd = "<" + rgt + ">"
				tl = "</" + rgt + ">"
			}
		case "-set":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -set command\n")
				os.Exit(1)
			}
			tmp := eutils.ConvertSlash(args[1])
			if tmp != "" {
				head = "<" + tmp + ">"
				tail = "</" + tmp + ">"
			}
		case "-rec":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -rec command\n")
				os.Exit(1)
			}
			tmp := eutils.ConvertSlash(args[1])
			if tmp != "" {
				hd = "<" + tmp + ">"
				tl = "</" + tmp + ">"
			}
		default:
			// if not any of the controls, set flag to break out of for loop
			inSwitch = false
		}

		if !inSwitch {
			break
		}

		// skip past arguments
		args = args[2:]

		if len(args) < 1 {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to transmute\n")
			os.Exit(1)
		}
	}

	// ENSURE PRESENCE OF PATTERN ARGUMENT

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to transmute\n")
		os.Exit(1)
	}

	// allow -record as synonym of -pattern (undocumented)
	if args[0] == "-record" || args[0] == "-Record" {
		args[0] = "-pattern"
	}

	// make sure top-level -pattern command is next
	if args[0] != "-pattern" && args[0] != "-Pattern" {
		fmt.Fprintf(os.Stderr, "\nERROR: No -pattern in command-line arguments\n")
		os.Exit(1)
	}
	if len(args) < 2 {
		fmt.Fprintf(os.Stderr, "\nERROR: Item missing after -pattern command\n")
		os.Exit(1)
	}

	topPat := args[1]
	if topPat == "" {
		fmt.Fprintf(os.Stderr, "\nERROR: Item missing after -pattern command\n")
		os.Exit(1)
	}
	if strings.HasPrefix(topPat, "-") {
		fmt.Fprintf(os.Stderr, "\nERROR: Misplaced %s command\n", topPat)
		os.Exit(1)
	}

	// look for -pattern Parent/* construct for heterogeneous data, e.g., -pattern PubmedArticleSet/*
	topPattern, star := eutils.SplitInTwoLeft(topPat, "/")
	if topPattern == "" {
		return
	}

	// CONCURRENT REFORMATTING OF PARSED XML RECORDS

	// -pattern plus -format does concurrent flush-left reformatting
	if len(args) > 2 && args[2] == "-format" {

		format := "flush"
		if len(args) > 3 {
			format = args[3]
			if strings.HasPrefix(format, "-") {
				format = "flush"
			}
		}

		xmlq := eutils.CreateXMLProducer(topPattern, star, false, rdr)
		fchq := createFormatters(topPattern, format, xmlq)
		unsq := eutils.CreateXMLUnshuffler(fchq)

		if xmlq == nil || fchq == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create formatter\n")
			os.Exit(1)
		}

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		// drain output channel
		for curr := range unsq {

			str := curr.Text

			if str == "" {
				continue
			}

			if hd != "" {
				os.Stdout.WriteString(hd)
				os.Stdout.WriteString("\n")
			}

			// send result to output
			os.Stdout.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				os.Stdout.WriteString("\n")
			}

			if tl != "" {
				os.Stdout.WriteString(tl)
				os.Stdout.WriteString("\n")
			}

			recordCount++
			runtime.Gosched()
		}

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// REPORT UNRECOGNIZED COMMAND

	fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized transmute command\n")
	os.Exit(1)
}
