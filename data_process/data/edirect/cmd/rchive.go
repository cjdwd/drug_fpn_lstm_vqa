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
// File Name:  rchive.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"container/heap"
	"eutils"
	"fmt"
	"hash/crc32"
	"html"
	"io"
	"io/ioutil"
	"os"
	"os/user"
	"path/filepath"
	"runtime"
	"runtime/debug"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

func reportEncodedMarkup(typ, id, str string) {

	var buffer strings.Builder

	max := len(str)

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

	findContext := func(fr, to int) string {

		numSpaces := 0

		for fr > 0 {
			ch := str[fr]
			if ch == ' ' {
				numSpaces++
				if numSpaces > 1 {
					fr++
					break
				}
			} else if ch == '\n' || ch == '>' {
				fr++
				break
			}
			fr--
		}

		numSpaces = 0

		for to < max {
			ch := str[to]
			if ch == ' ' {
				numSpaces++
				if numSpaces > 1 {
					break
				}
			} else if ch == '\n' || ch == '<' {
				break
			}
			to++
		}

		return str[fr:to]
	}

	reportMarkup := func(lbl string, fr, to int, txt string) {

		if lbl == typ || typ == "ALL" {
			// extract XML of SELF, SINGLE, DOUBLE, or AMPER types, or ALL
			buffer.WriteString(str)
			buffer.WriteString("\n")
		} else if typ == "" {
			// print report
			buffer.WriteString(id)
			buffer.WriteString("\t")
			buffer.WriteString(lbl)
			buffer.WriteString("\t")
			buffer.WriteString(txt)
			buffer.WriteString("\t| ")
			ctx := findContext(fr, to)
			buffer.WriteString(ctx)
			if eutils.HasUnicodeMarkup(ctx) {
				ctx = eutils.RepairUnicodeMarkup(ctx, eutils.SPACE)
			}
			ctx = eutils.RepairEncodedMarkup(ctx)
			buffer.WriteString("\t| ")
			buffer.WriteString(ctx)
			if eutils.HasAmpOrNotASCII(ctx) {
				ctx = html.UnescapeString(ctx)
			}
			buffer.WriteString("\t| ")
			buffer.WriteString(ctx)
			buffer.WriteString("\n")
		}
	}

	/*
		badTags := [10]string{
			"<i/>",
			"<i />",
			"<b/>",
			"<b />",
			"<u/>",
			"<u />",
			"<sup/>",
			"<sup />",
			"<sub/>",
			"<sub />",
		}
	*/

	skip := 0

	/*
		var prev rune
	*/

	for i, ch := range str {
		if skip > 0 {
			skip--
			continue
		}
		/*
			if ch > 127 {
				if IsUnicodeSuper(ch) {
					if IsUnicodeSubsc(prev) {
						// reportMarkup("UNIUP", i, i+2, string(ch))
					}
				} else if IsUnicodeSubsc(ch) {
					if IsUnicodeSuper(prev) {
						// reportMarkup("UNIDN", i, i+2, string(ch))
					}
				} else if ch == '\u0038' || ch == '\u0039' {
					// reportMarkup("ANGLE", i, i+2, string(ch))
				}
				prev = ch
				continue
			} else {
				prev = ' '
			}
		*/
		if ch == '<' {
			/*
				j := i + 1
				if j < max {
					nxt := str[j]
					if nxt == 'i' || nxt == 'b' || nxt == 'u' || nxt == 's' {
						for _, tag := range badTags {
							if strings.HasPrefix(str, tag) {
								k := len(tag)
								reportMarkup("SELF", i, i+k, tag)
								break
							}
						}
					}
				}
				if strings.HasPrefix(str[i:], "</sup><sub>") {
					// reportMarkup("SUPSUB", i, i+11, "</sup><sub>")
				} else if strings.HasPrefix(str[i:], "</sub><sup>") {
					// reportMarkup("SUBSUP", i, i+11, "</sub><sup>")
				}
			*/
			continue
		} else if ch != '&' {
			continue
		} else if strings.HasPrefix(str[i:], "&lt;") {
			sub := lookAhead(str[i:], 14)
			_, ok := eutils.HTMLRepair(sub)
			if ok {
				skip = len(sub) - 1
				reportMarkup("SINGLE", i, i+skip+1, sub)
				continue
			}
		} else if strings.HasPrefix(str[i:], "&amp;lt;") {
			sub := lookAhead(str[i:], 22)
			_, ok := eutils.HTMLRepair(sub)
			if ok {
				skip = len(sub) - 1
				reportMarkup("DOUBLE", i, i+skip+1, sub)
				continue
			}
		} else if strings.HasPrefix(str[i:], "&amp;amp;") {
			reportMarkup("AMPER", i, i+9, "&amp;amp;")
			skip = 8
			continue
		}
	}

	res := buffer.String()

	os.Stdout.WriteString(res)
}

// CONCURRENT GOROUTINE SERVERS

// processes with single goroutine call defer close(out) so consumer(s) can range over channel
// processes with multiple instances call defer wg.Done(), separate goroutine uses wg.Wait() to delay close(out)

func createExternalIndexer(args []string, zipp bool, in io.Reader) int {

	recordCount := 0

	transform := make(map[string]string)

	readTransformTable := func(tf string) {

		inFile, err := os.Open(tf)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Unable to open transformation file %s - %s\n", tf, err.Error())
			os.Exit(1)
		}

		scant := bufio.NewScanner(inFile)

		// populate transformation map
		for scant.Scan() {

			line := scant.Text()
			frst, scnd := eutils.SplitInTwoLeft(line, "\t")

			transform[frst] = scnd
		}

		inFile.Close()
	}

	// BIOCONCEPTS INDEXER

	// create intermediate table for {chemical|disease|gene}2pubtatorcentral.gz (undocumented)
	if args[0] == "-bioconcepts" {

		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient arguments for -bioconcepts\n")
			os.Exit(1)
		}

		// read transformation file
		tf := args[1]
		readTransformTable(tf)

		var buffer strings.Builder
		count := 0
		okay := false

		wrtr := bufio.NewWriter(os.Stdout)

		scanr := bufio.NewScanner(in)

		currpmid := ""

		// read lines of PMIDs and extracted concepts
		for scanr.Scan() {

			line := scanr.Text()

			cols := strings.Split(line, "\t")
			if len(cols) != 5 {
				continue
			}

			pmid := cols[0]
			if currpmid != pmid {
				// end current block
				currpmid = pmid

				if pmid == "" {
					continue
				}

				recordCount++
				count++

				if count >= 1000 {
					count = 0
					txt := buffer.String()
					if txt != "" {
						// print current buffer
						wrtr.WriteString(txt[:])
					}
					buffer.Reset()
				}

				okay = true
			}

			typ := cols[1]
			val := cols[2]
			switch typ {
			case "Gene":
				genes := strings.Split(val, ";")
				for _, gene := range genes {
					if gene == "None" {
						continue
					}
					buffer.WriteString(pmid)
					buffer.WriteString("\t")
					buffer.WriteString("GENE")
					buffer.WriteString("\t")
					buffer.WriteString(gene)
					buffer.WriteString("\n")
					gn, ok := transform[gene]
					if !ok || gn == "" {
						continue
					}
					buffer.WriteString(pmid)
					buffer.WriteString("\t")
					buffer.WriteString("GENE")
					buffer.WriteString("\t")
					buffer.WriteString(gn)
					buffer.WriteString("\n")
				}
			case "Disease":
				if strings.HasPrefix(val, "MESH:") {
					diszs := strings.Split(val[5:], "|")
					for _, disz := range diszs {
						buffer.WriteString(pmid)
						buffer.WriteString("\t")
						buffer.WriteString("DISZ")
						buffer.WriteString("\t")
						buffer.WriteString(disz)
						buffer.WriteString("\n")
						dn, ok := transform[disz]
						if !ok || dn == "" {
							continue
						}
						buffer.WriteString(pmid)
						buffer.WriteString("\t")
						buffer.WriteString("DISZ")
						buffer.WriteString("\t")
						buffer.WriteString(dn)
						buffer.WriteString("\n")
					}
				} else if strings.HasPrefix(val, "OMIM:") {
					omims := strings.Split(val[5:], "|")
					for _, omim := range omims {
						// was OMIM, now fused with DISZ
						buffer.WriteString(pmid)
						buffer.WriteString("\t")
						buffer.WriteString("DISZ")
						buffer.WriteString("\t")
						// tag OMIM identifiers with M prefix
						buffer.WriteString("M" + omim)
						buffer.WriteString("\n")
					}
				}
			case "Chemical":
				if strings.HasPrefix(val, "MESH:") {
					chems := strings.Split(val[5:], "|")
					for _, chem := range chems {
						buffer.WriteString(pmid)
						buffer.WriteString("\t")
						buffer.WriteString("CHEM")
						buffer.WriteString("\t")
						buffer.WriteString(chem)
						buffer.WriteString("\n")
						ch, ok := transform[chem]
						if !ok || ch == "" {
							continue
						}
						buffer.WriteString(pmid)
						buffer.WriteString("\t")
						buffer.WriteString("CHEM")
						buffer.WriteString("\t")
						buffer.WriteString(ch)
						buffer.WriteString("\n")
					}
				} else if strings.HasPrefix(val, "CHEBI:") {
					chebs := strings.Split(val[6:], "|")
					for _, cheb := range chebs {
						// was CEBI, now fused with CHEM
						buffer.WriteString(pmid)
						buffer.WriteString("\t")
						buffer.WriteString("CHEM")
						buffer.WriteString("\t")
						// tag CHEBI identifiers with H prefix
						buffer.WriteString("H" + cheb)
						buffer.WriteString("\n")
					}
				}
			case "Species":
			case "Mutation":
			case "CellLine":
			default:
			}
		}

		if okay {
			txt := buffer.String()
			if txt != "" {
				// print current buffer
				wrtr.WriteString(txt[:])
			}
		}
		buffer.Reset()

		wrtr.Flush()

		return recordCount
	}

	// GENERIF INDEXER

	// create intermediate table for generifs_basic.gz (undocumented)
	if args[0] == "-generif" || args[0] == "-generifs" {

		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient arguments for -generif\n")
			os.Exit(1)
		}

		// read transformation file
		tf := args[1]
		readTransformTable(tf)

		var buffer strings.Builder
		count := 0
		okay := false

		wrtr := bufio.NewWriter(os.Stdout)

		scanr := bufio.NewScanner(in)

		currpmid := ""

		// skip first line with column heading names
		for scanr.Scan() {

			line := scanr.Text()
			cols := strings.Split(line, "\t")
			if len(cols) != 5 {
				fmt.Fprintf(os.Stderr, "Unexpected number of columns (%d) in generifs_basic.gz\n", len(cols))
				os.Exit(1)
			}
			if len(cols) != 5 || cols[0] != "#Tax ID" {
				fmt.Fprintf(os.Stderr, "Unrecognized contents in generifs_basic.gz\n")
				os.Exit(1)
			}
			break
		}

		// read lines of PMIDs and gene references
		for scanr.Scan() {

			line := scanr.Text()

			cols := strings.Split(line, "\t")
			if len(cols) != 5 {
				continue
			}

			val := cols[2]
			pmids := strings.Split(val, ",")
			for _, pmid := range pmids {
				if currpmid != pmid {
					// end current block
					currpmid = pmid

					if pmid == "" {
						continue
					}

					recordCount++
					count++

					if count >= 1000 {
						count = 0
						txt := buffer.String()
						if txt != "" {
							// print current buffer
							wrtr.WriteString(txt[:])
						}
						buffer.Reset()
					}

					okay = true
				}

				gene := cols[1]
				// was GRIF, now fused with GENE
				buffer.WriteString(pmid)
				buffer.WriteString("\t")
				buffer.WriteString("GENE")
				buffer.WriteString("\t")
				buffer.WriteString(gene)
				buffer.WriteString("\n")
				gn, ok := transform[gene]
				if !ok || gn == "" {
					continue
				}
				buffer.WriteString(pmid)
				buffer.WriteString("\t")
				buffer.WriteString("GENE")
				buffer.WriteString("\t")
				buffer.WriteString(gn)
				buffer.WriteString("\n")
			}
		}

		if okay {
			txt := buffer.String()
			if txt != "" {
				// print current buffer
				wrtr.WriteString(txt[:])
			}
		}
		buffer.Reset()

		wrtr.Flush()

		return recordCount
	}

	// THEME INDEXER

	/*
	   A+    Agonism, activation                      N     Inhibits
	   A-    Antagonism, blocking                     O     Transport, channels
	   B     Binding, ligand                          Pa    Alleviates, reduces
	   C     Inhibits cell growth                     Pr    Prevents, suppresses
	   D     Drug targets                             Q     Production by cell population
	   E     Affects expression/production            Rg    Regulation
	   E+    Increases expression/production          Sa    Side effect/adverse event
	   E-    Decreases expression/production          T     Treatment/therapy
	   G     Promotes progression                     Te    Possible therapeutic effect
	   H     Same protein or complex                  U     Causal mutations
	   I     Signaling pathway                        Ud    Mutations affecting disease course
	   J     Role in disease pathogenesis             V+    Activates, stimulates
	   K     Metabolism, pharmacokinetics             W     Enhances response
	   L     Improper regulation linked to disease    X     Overexpression in disease
	   Md    Biomarkers (diagnostic)                  Y     Polymorphisms alter risk
	   Mp    Biomarkers (progression)                 Z     Enzyme activity
	*/

	// create intermediate table for chemical-gene-disease themes (undocumented)
	if args[0] == "-theme" || args[0] == "-themes" {

		if len(args) < 4 {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient arguments for -theme\n")
			os.Exit(1)
		}

		one := args[1]
		two := args[2]
		tag := args[3]

		// for disambiguating B, E, E+, and J themes, in CHDI, CHGE, GEDI, and GEGE data sets
		sfx := ""
		if len(tag) > 0 {
			switch tag[0] {
			case 'C':
				sfx = "c"
			case 'G':
				sfx = "g"
			}
		}

		fl, err := os.Open(one)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", one)
			os.Exit(1)
		}

		scanr := bufio.NewScanner(fl)

		var columns []string
		numCols := 0

		// read first line with column heading names
		if scanr.Scan() {

			line := scanr.Text()
			line = strings.Replace(line, "+", "p", -1)
			line = strings.Replace(line, "-", "m", -1)
			columns = strings.Split(line, "\t")
			numCols = len(columns)

			if numCols < 3 {
				fmt.Fprintf(os.Stderr, "Unexpected number of columns (%d) in part-i file\n", numCols)
				os.Exit(1)
			}
			if columns[0] != "path" {
				fmt.Fprintf(os.Stderr, "Unrecognized contents in part-i file\n")
				os.Exit(1)
			}
		}

		var scores []int

		for i := 0; i < numCols; i++ {
			scores = append(scores, 0)
		}

		mapper := make(map[string]string)

		scorer := make(map[string]int)

		// read lines of dependency paths, scores for each theme
		for scanr.Scan() {

			line := scanr.Text()

			cols := strings.Split(line, "\t")
			if len(cols) != numCols {
				fmt.Fprintf(os.Stderr, "Mismatched columns in '%s'\n", line)
				continue
			}

			for i := 0; i < numCols; i++ {
				scores[i] = 0
			}

			sum := 0
			// increment by 2 to ignore flagship indicator fields
			for i := 1; i < numCols; i += 2 {
				str := cols[i]
				str, _ = eutils.SplitInTwoLeft(str, ".")
				val, err := strconv.Atoi(str)
				if err != nil {
					fmt.Fprintf(os.Stderr, "Unrecognized number '%s'\n", str)
					continue
				}
				scores[i] = val
				sum += val
			}
			if sum == 0 {
				continue
			}

			path := cols[0]
			path = strings.ToLower(path)
			themes := ""
			comma := ""
			for i := 1; i < numCols; i += 2 {
				// find scores over cutoff
				if scores[i]*3 > sum {
					theme := columns[i]
					themes += comma
					themes += theme
					comma = ","
					scorer[path+"_"+theme] = scores[i] * 100 / sum
				}
			}
			if themes == "" {
				continue
			}
			// populate theme lookup table
			mapper[path] = themes
		}

		fl.Close()

		fl, err = os.Open(two)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", two)
			os.Exit(1)
		}

		var buffer strings.Builder
		count := 0
		okay := false

		wrtr := bufio.NewWriter(os.Stdout)

		scanr = bufio.NewScanner(fl)

		// read lines of PMIDs and dependency paths
		for scanr.Scan() {

			line := scanr.Text()

			cols := strings.Split(line, "\t")
			if len(cols) != 14 {
				fmt.Fprintf(os.Stderr, "Mismatched columns in '%s'\n", line)
				continue
			}

			pmid := cols[0]
			path := cols[12]
			path = strings.ToLower(path)
			themes, ok := mapper[path]
			if !ok {
				continue
			}

			fwd := tag[0:2] + tag[2:4]
			rev := tag[2:4] + tag[0:2]

			cleanCol := func(str string) string {
				if str == "" || str == "null" {
					return ""
				}
				// remove database prefixes, tag OMIM and CHEBI identifiers with M and H prefixes
				if strings.HasPrefix(str, "MESH:") {
					str = strings.TrimPrefix(str, "MESH:")
				} else if strings.HasPrefix(str, "OMIM:") {
					str = "M" + strings.TrimPrefix(str, "OMIM:")
				} else if strings.HasPrefix(str, "CHEBI:") {
					str = "H" + strings.TrimPrefix(str, "CHEBI:")
				}
				idx := strings.Index(str, "(")
				if idx > 0 {
					// remove parenthetical Tax suffix
					str = str[:idx]
				}
				str = strings.ToLower(str)
				return str
			}

			splitCol := func(str string) []string {
				// multiple genes may be separated by semicolons
				if strings.Index(str, ";") >= 0 {
					return strings.Split(str, ";")
				}
				// mesh, omim, and chebi may be separated by vertical bars
				return strings.Split(str, "|")
			}

			frst := splitCol(cols[8])
			scnd := splitCol(cols[9])

			printItem := func(pmid, fld, item string) {
				if pmid == "" || fld == "" || item == "" {
					return
				}
				buffer.WriteString(pmid)
				buffer.WriteString("\t")
				buffer.WriteString(fld)
				buffer.WriteString("\t")
				buffer.WriteString(item)
				buffer.WriteString("\n")
			}

			thms := strings.Split(themes, ",")
			for _, theme := range thms {
				if theme == "" {
					continue
				}

				printItem(pmid, "THME", theme)

				// disambiguate B, E, E+, and J themes that appear in two data sets
				switch theme {
				case "B", "E", "J":
					printItem(pmid, "THME", theme+sfx)
				case "Ep":
					printItem(pmid, "THME", "E"+sfx+"p")
				case "Em":
					printItem(pmid, "THME", "E"+sfx+"m")
				}
			}

			for _, frs := range frst {
				fst := cleanCol(frs)
				if fst == "" {
					continue
				}

				for _, snd := range scnd {
					scd := cleanCol(snd)
					if scd == "" {
						continue
					}

					printItem(pmid, "CONV", fwd)

					printItem(pmid, "CONV", rev)

					printItem(pmid, "CONV", fst+" "+scd)

					printItem(pmid, "CONV", scd+" "+fst)

					printItem(pmid, "CONV", fwd+" "+fst+" "+scd)

					printItem(pmid, "CONV", rev+" "+scd+" "+fst)

					for _, theme := range thms {
						if theme == "" {
							continue
						}

						score := scorer[path+"_"+theme]
						pct := strconv.Itoa(score)

						printItem(pmid, "CONV", theme+" "+fwd)

						printItem(pmid, "CONV", theme+" "+rev)

						printItem(pmid, "CONV", theme+" "+fst+" "+scd+" "+pct)

						printItem(pmid, "CONV", theme+" "+scd+" "+fst+" "+pct)

						printItem(pmid, "CONV", theme+" "+fwd+" "+fst+" "+scd+" "+pct)

						printItem(pmid, "CONV", theme+" "+rev+" "+scd+" "+fst+" "+pct)
					}
				}
			}

			recordCount++
			count++

			if count >= 1000 {
				count = 0
				txt := buffer.String()
				if txt != "" {
					// print current buffer
					wrtr.WriteString(txt[:])
				}
				buffer.Reset()
			}

			okay = true
		}

		if okay {
			txt := buffer.String()
			if txt != "" {
				// print current buffer
				wrtr.WriteString(txt[:])
			}
		}
		buffer.Reset()

		wrtr.Flush()

		fl.Close()

		return recordCount
	}

	// DEPENDENCY PATH INDEXER

	// create intermediate table for chemical-gene-disease dependency paths (undocumented)
	if args[0] == "-dpath" || args[0] == "-dpaths" {

		var buffer strings.Builder
		count := 0
		okay := false

		replr := strings.NewReplacer(
			">", "_gtrthan_",
			"<", "_lssthan_",
			"/", "_slash_",
			"%", "_prcnt_",
			":", "_colln_",
			"+", "_pluss_",
			"!", "_exclam_",
			"?", "_qmark_",
			"'", "_squot_",
			"(", "_lparen_",
			")", "_rparen_",
		)
		if replr == nil {
			fmt.Fprintf(os.Stderr, "Unable to create replacer\n")
			os.Exit(1)
		}

		wrtr := bufio.NewWriter(os.Stdout)

		scanr := bufio.NewScanner(in)

		// read lines of PMIDs and dependency paths
		for scanr.Scan() {

			line := scanr.Text()

			cols := strings.Split(line, "\t")
			if len(cols) != 14 {
				fmt.Fprintf(os.Stderr, "Mismatched columns in '%s'\n", line)
				continue
			}

			pmid := cols[0]
			path := cols[12]
			path = strings.ToLower(path)

			// rescue known characters
			tmp := eutils.CompressRunsOfSpaces(path)
			tmp = strings.TrimSpace(tmp)

			tmp = " " + tmp + " "

			tmp = replr.Replace(tmp)

			tmp = eutils.CompressRunsOfSpaces(tmp)
			tmp = strings.TrimSpace(tmp)

			// final cleanup
			tmp = strings.Replace(tmp, "|", "_", -1)
			tmp = strings.Replace(tmp, "__", "_", -1)

			pths := strings.Split(tmp, " ")
			for _, pth := range pths {
				if pth == "" {
					continue
				}
				buffer.WriteString(pmid)
				buffer.WriteString("\t")
				buffer.WriteString("PATH")
				buffer.WriteString("\t")
				buffer.WriteString(pth)
				buffer.WriteString("\n")
			}

			recordCount++
			count++

			if count >= 1000 {
				count = 0
				txt := buffer.String()
				if txt != "" {
					// print current buffer
					wrtr.WriteString(txt[:])
				}
				buffer.Reset()
			}

			okay = true
		}

		if okay {
			txt := buffer.String()
			if txt != "" {
				// print current buffer
				wrtr.WriteString(txt[:])
			}
		}
		buffer.Reset()

		wrtr.Flush()

		return recordCount
	}

	// THESIS INDEXER

	// create -e2index file for bioconcepts, geneRIFs, and themes and their dependency paths (undocumented)
	if args[0] == "-thesis" {

		// e.g., -thesis 250000 "$target" "biocchem"
		if len(args) < 4 {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient arguments for -thesis\n")
			os.Exit(1)
		}

		chunk, err := strconv.Atoi(args[1])
		if err != nil {
			fmt.Fprintf(os.Stderr, "Unrecognized count - '%s'\n", err.Error())
			os.Exit(1)
		}
		target := strings.TrimSuffix(args[2], "/")
		prefix := args[3]

		suffix := "e2x"
		sfx := suffix
		if zipp {
			sfx += ".gz"
		}

		fnum := 0

		fld := ""

		scanr := bufio.NewScanner(os.Stdin)

		processChunk := func() bool {

			// map for combined index
			indexed := make(map[string][]string)

			writeChunk := func() {

				var (
					fl   *os.File
					wrtr *bufio.Writer
					zpr  *gzip.Writer
					err  error
				)

				fnum++
				fpath := fmt.Sprintf("%s/%s%03d.%s", target, prefix, fnum, sfx)
				fl, err = os.Create(fpath)
				if err != nil {
					fmt.Fprintf(os.Stderr, "%s\n", err.Error())
					return
				}
				defer fl.Close()

				pth := fmt.Sprintf("%s%03d.%s", prefix, fnum, suffix)
				os.Stderr.WriteString(pth + "\n")

				var out io.Writer

				out = fl

				if zipp {

					zpr, err = gzip.NewWriterLevel(fl, gzip.BestSpeed)
					if err != nil {
						fmt.Fprintf(os.Stderr, "%s\n", err.Error())
						return
					}

					out = zpr
				}

				wrtr = bufio.NewWriter(out)
				if wrtr == nil {
					fmt.Fprintf(os.Stderr, "Unable to create bufio.NewWriter\n")
					return
				}

				var buffer strings.Builder
				count := 0

				buffer.WriteString("<IdxDocumentSet>\n")

				// sort fields in alphabetical order
				var keys []string
				for ky := range indexed {
					keys = append(keys, ky)
				}

				if len(keys) > 1 {
					sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })
				}

				for _, idx := range keys {

					item, ok := indexed[idx]
					if !ok {
						continue
					}

					uid := item[0]
					data := item[1:]

					if uid == "" || len(data) < 1 {
						continue
					}

					if len(data) > 1 {
						sort.Slice(data, func(i, j int) bool { return data[i] < data[j] })
					}

					buffer.WriteString("  <IdxDocument>\n")
					buffer.WriteString("    <IdxUid>")
					buffer.WriteString(uid)
					buffer.WriteString("</IdxUid>\n")
					buffer.WriteString("    <IdxSearchFields>\n")

					prev := ""
					for len(data) > 0 {
						val := data[0]
						data = data[1:]

						if val == prev {
							continue
						}

						buffer.WriteString("      <")
						buffer.WriteString(fld)
						buffer.WriteString(">")
						buffer.WriteString(val)
						buffer.WriteString("</")
						buffer.WriteString(fld)
						buffer.WriteString(">\n")

						prev = val
					}

					buffer.WriteString("    </IdxSearchFields>\n")
					buffer.WriteString("  </IdxDocument>\n")

					recordCount++
					count++

					if count >= 1000 {
						count = 0
						txt := buffer.String()
						if txt != "" {
							// print current buffer
							wrtr.WriteString(txt[:])
						}
						buffer.Reset()
					}
				}

				buffer.WriteString("</IdxDocumentSet>\n")

				txt := buffer.String()
				if txt != "" {
					// print current buffer
					wrtr.WriteString(txt[:])
				}
				buffer.Reset()

				wrtr.Flush()

				if zpr != nil {
					zpr.Close()
				}
			}

			lineCount := 0
			okay := false

			// read lines of dependency paths, scores for each theme
			for scanr.Scan() {

				line := scanr.Text()

				cols := strings.Split(line, "\t")
				if len(cols) != 3 {
					fmt.Fprintf(os.Stderr, "Mismatched columns in '%s'\n", line)
					continue
				}

				uid := cols[0]
				fd := cols[1]
				val := cols[2]
				if uid == "" || fd == "" || val == "" {
					continue
				}
				if fld == "" {
					fld = fd
				}
				if fld != fd {
					fmt.Fprintf(os.Stderr, "Field '%s' expected, '%s' found\n", fld, fd)
					continue
				}

				val = strings.ToLower(val)
				// convert angle brackets in chemical names
				val = html.EscapeString(val)

				data, ok := indexed[uid]
				if !ok {
					data = make([]string, 0, 2)
					// first entry on new slice is uid
					data = append(data, uid)
				}
				data = append(data, val)
				// always need to update indexed, since data may be reallocated
				indexed[uid] = data

				okay = true

				lineCount++
				if lineCount > chunk {
					break
				}
			}

			if okay {
				writeChunk()
				return true
			}

			return false
		}

		for processChunk() {
			// loop until scanner runs out of lines
		}

		return recordCount
	}

	return 0
}

func createExternalArchive(stash string, args []string) <-chan string {

	makePresenters := func(args []string) []<-chan eutils.Plex {

		if args == nil {
			return nil
		}

		numFiles := len(args)
		if numFiles < 1 {
			fmt.Fprintf(os.Stderr, "\nERROR: Not enough indexed files to merge\n")
			os.Exit(1)
		}

		chns := make([]<-chan eutils.Plex, numFiles)
		if chns == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create presenter channel array\n")
			os.Exit(1)
		}

		// xmlPresenter sends partitioned XML strings through channel
		xmlPresenter := func(fileNum int, fileName string, out chan<- eutils.Plex) {

			// close channel when all records have been processed
			defer close(out)

			f, err := os.Open(fileName)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
				os.Exit(1)
			}

			// close input file when all records have been processed
			defer f.Close()

			var in io.Reader

			in = f

			// if suffix is ".gz", use decompressor
			iszip := false
			if strings.HasSuffix(fileName, ".gz") {
				iszip = true
			}

			if iszip {
				brd := bufio.NewReader(f)
				if brd == nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to create buffered reader on '%s'\n", fileName)
					os.Exit(1)
				}
				zpr, err := gzip.NewReader(brd)
				if err != nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to create decompressor on '%s'\n", fileName)
					os.Exit(1)
				}

				// close decompressor when all records have been processed
				defer zpr.Close()

				// use decompressor for reading file
				in = zpr
			}

			rdr := eutils.CreateXMLStreamer(in)

			if rdr == nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
				os.Exit(1)
			}

			find := eutils.ParseIndex("IdxUid")

			// partition all input by pattern and send XML substring through channel
			eutils.PartitionPattern("IdxDocument", "", false, rdr,
				func(str string) {
					id := eutils.FindIdentifier(str[:], "IdxDocument", find)

					out <- eutils.Plex{Which: fileNum, Ident: id, Text: str, Index: 0, Sibs: nil}
				})
		}

		// launch multiple presenter goroutines
		for i, str := range args {

			chn := make(chan eutils.Plex, eutils.ChanDepth())
			if chn == nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create presenter channel\n")
				os.Exit(1)
			}

			go xmlPresenter(i, str, chn)

			chns[i] = chn
		}

		// no need for separate anonymous goroutine to wait until all presenters are done

		return chns
	}

	makeManifold := func(inp []<-chan eutils.Plex) <-chan eutils.Plex {

		if inp == nil {
			return nil
		}

		out := make(chan eutils.Plex, eutils.ChanDepth())
		if out == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create manifold channel\n")
			os.Exit(1)
		}

		// xmlManifold restores alphabetical order of merged postings
		xmlManifold := func(inp []<-chan eutils.Plex, out chan<- eutils.Plex) {

			// close channel when all records have been processed
			defer close(out)

			// initialize empty heap
			hp := &eutils.PlexHeap{}
			heap.Init(hp)

			// read first object from all input channels in turn
			for _, chn := range inp {
				plx, ok := <-chn
				if ok {
					heap.Push(hp, plx)
				}
			}

			// array to collect strings with same identifier
			var arry []string

			prevIdent := ""
			rec := 0

			// reading from heap returns objects in alphabetical order
			for hp.Len() > 0 {

				// remove lowest item from heap, use interface type assertion
				curr := heap.Pop(hp).(eutils.Plex)

				// compare adjacent record identifiers
				if prevIdent == curr.Ident {

					// save next indexed object string in slice
					arry = append(arry, curr.Text)

				} else {

					if len(arry) > 0 {

						rec++
						// send set from previous identifier to output channel
						out <- eutils.Plex{Which: 0, Ident: prevIdent, Text: "", Index: rec, Sibs: arry}

						// empty the slice
						arry = nil
					}

					// remember new identifier
					prevIdent = curr.Ident

					// save first indexed object with this identifier
					arry = append(arry, curr.Text)
				}

				// read next object from channel that just supplied lowest item
				chn := inp[curr.Which]
				plx, ok := <-chn
				if ok {
					heap.Push(hp, plx)
				}
			}

			if len(arry) > 0 {

				rec++
				// send last record
				out <- eutils.Plex{Which: 0, Ident: prevIdent, Text: "", Index: rec, Sibs: arry}

				arry = nil
			}
		}

		// launch single manifold goroutine
		go xmlManifold(inp, out)

		return out
	}

	makeMergers := func(inp <-chan eutils.Plex) <-chan eutils.XMLRecord {

		if inp == nil {
			return nil
		}

		out := make(chan eutils.XMLRecord, eutils.ChanDepth())
		if out == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create merger channel\n")
			os.Exit(1)
		}

		// xmlMerger fuses adjacent IdxDocument records with the same identifier
		xmlMerger := func(wg *sync.WaitGroup, inp <-chan eutils.Plex, out chan<- eutils.XMLRecord) {

			defer wg.Done()

			var buffer strings.Builder

			fusePostings := func(key string, data []string) string {

				fields := make(map[string]map[string]string)

				addIdents := func(fld, pos, uid string) {

					// no need for mutex here
					positions, ok := fields[fld]
					if !ok {
						positions = make(map[string]string)
						fields[fld] = positions
					}

					positions[uid] = pos
				}

				addUID := func(tag, attr, content string) {

					if tag != "IdxUid" {

						addIdents(tag, attr, content)
					}
				}

				for _, str := range data {
					eutils.StreamValues(str[:], "IdxDocument", addUID)
				}

				buffer.Reset()

				buffer.WriteString("<IdxDocument>\n")
				buffer.WriteString("<IdxUid>")
				buffer.WriteString(key)
				buffer.WriteString("</IdxUid>\n")
				buffer.WriteString("<InvIDs>\n")

				// sort fields in alphabetical order
				var keys []string
				for ky := range fields {
					keys = append(keys, ky)
				}
				sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })

				for _, fld := range keys {

					positions := fields[fld]

					var arry []string

					for item := range positions {
						arry = append(arry, item)
					}

					if len(arry) > 1 {
						sort.Slice(arry, func(i, j int) bool {
							// numeric sort on strings checks lengths first
							lni := len(arry[i])
							lnj := len(arry[j])
							// shorter string is numerically less, assuming no leading zeros
							if lni < lnj {
								return true
							}
							if lni > lnj {
								return false
							}
							// same length, can now do string comparison on contents
							return arry[i] < arry[j]
						})
					}

					// print list of UIDs, skipping duplicates
					last := ""
					for _, uid := range arry {
						// detect duplicate UIDs, now in same list after conversion of one term entry from foreign alphabet
						if uid == last {
							continue
						}
						buffer.WriteString("<")
						buffer.WriteString(fld)
						atr := positions[uid]
						if atr != "" {
							buffer.WriteString(" ")
							buffer.WriteString(atr)
						}
						buffer.WriteString(">")
						buffer.WriteString(uid)
						buffer.WriteString("</")
						buffer.WriteString(fld)
						buffer.WriteString(">\n")

						last = uid
					}
				}

				buffer.WriteString("</InvIDs>\n")
				buffer.WriteString("</IdxDocument>\n")

				txt := buffer.String()

				return txt
			}

			for plx := range inp {

				rec := plx.Index
				key := plx.Ident
				data := plx.Sibs

				if len(data) < 1 {
					continue
				}

				str := fusePostings(key, data)

				out <- eutils.XMLRecord{Index: rec, Ident: key, Text: str}

				runtime.Gosched()
			}
		}

		var wg sync.WaitGroup

		// launch multiple merger goroutines
		for i := 0; i < eutils.NumServe(); i++ {
			wg.Add(1)
			go xmlMerger(&wg, inp, out)
		}

		// launch separate anonymous goroutine to wait until all mergers are done
		go func() {
			wg.Wait()
			close(out)
		}()

		return out
	}

	chns := makePresenters(args)
	mfld := makeManifold(chns)
	mrgr := makeMergers(mfld)
	stsq := eutils.CreateStashers(stash, "IdxDocument", "IdxDocument/IdxUid", ".e2x", false, true, 50000, mrgr)

	if chns == nil || mfld == nil || mrgr == nil || stsq == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create extra index stasher\n")
		os.Exit(1)
	}

	return stsq
}

// MAIN FUNCTION

func main() {

	// skip past executable name
	args := os.Args[1:]

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: No command-line arguments supplied to rchive\n")
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
	deStop := true

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

	// garbage collector control can be set by environment variable or default value with -gogc 0
	goGc = 200
	gcdefault := true

	// -flag sets -strict or -mixed cleanup flags from argument
	flgs := ""

	// read data from file instead of stdin
	fileName := ""

	// debugging
	stts := false
	timr := false

	// profiling
	prfl := false

	// element to use as local data index
	indx := ""

	// file of index values for removing duplicates
	unqe := ""

	// path for local data indexed as trie
	stsh := ""
	ftch := ""
	strm := ""

	// path for local extra link data
	smmn := ""

	// flag for inverted index
	nvrt := false

	// flag for combining sets of inverted files
	join := false

	// flag for combining sets of inverted files
	fuse := false

	// destination directory for merging and splitting inverted files
	merg := ""

	// base destination directory for promoting inverted index to retrieval indices
	prom := ""

	// fields for promoting inverted index files
	fild := ""

	// base for queries
	base := ""

	// query by phrase, normalized terms (with truncation wildcarding)
	phrs := ""
	rlxd := false
	xact := false
	titl := false
	mock := false
	btch := false

	// print term list with counts
	trms := ""
	plrl := false
	psns := false

	// use gzip compression on local data files
	zipp := false

	// print UIDs and hash values
	hshv := false

	// convert UIDs to archive trie
	trei := false

	// compare input record against stash
	cmpr := false
	cmprType := ""
	ignr := ""

	// flag missing identifiers
	msng := false

	// flag records with damaged embedded HTML tags
	dmgd := false
	dmgdType := ""

	// kludge to use non-threaded fetching for windows
	windows := false

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
			gcdefault = false

		// read data from file
		case "-input":
			fileName = eutils.GetStringArg(args, "Input file name")
			args = args[1:]

		// file with selected indexes for removing duplicates
		case "-unique":
			unqe = eutils.GetStringArg(args, "Unique identifier file")
			args = args[1:]

		// local directory path for indexing
		case "-archive", "-stash":
			stsh = eutils.GetStringArg(args, "Archive path")
			if stsh != "" && !strings.HasSuffix(stsh, "/") {
				stsh += "/"
			}
			args = args[1:]
		// local directory path for retrieval
		case "-fetch":
			ftch = eutils.GetStringArg(args, "Fetch path")
			if ftch != "" && !strings.HasSuffix(ftch, "/") {
				ftch += "/"
			}
			args = args[1:]
		// local directory path for retrieval of compressed XML
		case "-stream":
			strm = eutils.GetStringArg(args, "Stream path")
			if strm != "" && !strings.HasSuffix(strm, "/") {
				strm += "/"
			}
			args = args[1:]

		// local directory path for extra link retrieval
		case "-summon":
			smmn = eutils.GetStringArg(args, "Summon path")
			args = args[1:]

		// data element for indexing
		case "-index":
			indx = eutils.GetStringArg(args, "Index element")
			args = args[1:]

		// build inverted index
		case "-invert":
			nvrt = true

		// combine sets of inverted index files
		case "-join":
			join = true

		case "-fuse":
			fuse = true

		// merge inverted index files, distribute by prefix
		case "-merge":
			merg = eutils.GetStringArg(args, "Merge field")
			args = args[1:]

		// promote inverted index to term-specific postings files
		case "-promote":
			if len(args) < 3 {
				fmt.Fprintf(os.Stderr, "\nERROR: Promote path is missing\n")
				os.Exit(1)
			}
			prom = args[1]
			fild = args[2]
			// skip past first and second arguments
			args = args[2:]

		case "-path":
			base = eutils.GetStringArg(args, "Postings path")
			args = args[1:]

		case "-title":
			titl = true
			fallthrough
		case "-exact":
			xact = true
			fallthrough
		case "-search":
			rlxd = true
			fallthrough
		case "-query":
			if xact && rlxd {
				rlxd = false
			}
			phrs = eutils.GetStringArg(args, "Query argument")
			args = args[1:]

		case "-batch":
			btch = true

		case "-mockt":
			titl = true
			fallthrough
		case "-mockx":
			xact = true
			fallthrough
		case "-mocks":
			rlxd = true
			fallthrough
		case "-mock":
			if xact && rlxd {
				rlxd = false
			}
			phrs = eutils.GetStringArg(args, "Query argument")
			mock = true
			args = args[1:]

		// -countp tests the files containing positions of terms per UID (undocumented)
		case "-countp":
			psns = true
			fallthrough
		case "-counts":
			plrl = true
			fallthrough
		case "-countr":
			rlxd = true
			fallthrough
		case "-count":
			if plrl && rlxd {
				rlxd = false
			}
			trms = eutils.GetStringArg(args, "Count argument")
			args = args[1:]

		case "-gzip":
			zipp = true
		case "-hash":
			hshv = true
		case "-trie":
			trei = true
		// check for missing records
		case "-missing":
			msng = true

		// use non-threaded fetch function for windows (undocumented)
		case "-windows":
			windows = true

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
			// ignored, kept for backwards compatibility
		case "-stops", "-stop":
			deStop = false

		case "-unicode":
			// DoUnicode = true
		case "-script":
			// DoScript = true
		case "-mathml":
			// DoMathML = true

		case "-flag", "-flags":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: -flags argument is missing\n")
				os.Exit(1)
			}
			flgs = eutils.GetStringArg(args, "Flags argument")
			args = args[1:]

		// debugging flags
		case "-damaged", "-damage", "-broken":
			dmgd = true
			if len(args) > 1 {
				next := args[1]
				// if next argument is not another flag
				if next != "" && next[0] != '-' {
					// get optional extraction class (SELF, SINGLE, DOUBLE, AMPER, or ALL)
					dmgdType = next
					// skip past first of two arguments
					args = args[1:]
				}
			}
		case "-prepare":
			cmpr = true
			if len(args) > 1 {
				next := args[1]
				// if next argument is not another flag
				if next != "" && next[0] != '-' {
					// get optional data source specifier
					cmprType = next
					// skip past first of two arguments
					args = args[1:]
				}
			}
		case "-ignore":
			ignr = eutils.GetStringArg(args, "-ignore value")
			args = args[1:]

		// debugging flags
		case "-debug":
			// dbug = true
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

	// -flag allows script to set -strict or -mixed (or -stops) from argument
	switch flgs {
	case "strict":
		doStrict = true
	case "mixed":
		doMixed = true
	case "stems", "stem":
		// ignored, kept for backwards compatibility
	case "stops", "stop":
		deStop = false
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

	// if copying from local files accessed by identifier, add dummy argument to bypass length tests
	if stsh != "" && indx == "" {
		args = append(args, "-dummy")
	} else if ftch != "" || strm != "" || smmn != "" {
		args = append(args, "-dummy")
	} else if base != "" {
		args = append(args, "-dummy")
	} else if trei || dmgd || cmpr {
		args = append(args, "-dummy")
	}

	// expand -archive ~/ to home directory path
	if stsh != "" {

		if stsh[:2] == "~/" {
			cur, err := user.Current()
			if err == nil {
				hom := cur.HomeDir
				stsh = strings.Replace(stsh, "~/", hom+"/", 1)
			}
		}
	}

	// expand -fetch ~/ to home directory path
	if ftch != "" {

		if ftch[:2] == "~/" {
			cur, err := user.Current()
			if err == nil {
				hom := cur.HomeDir
				ftch = strings.Replace(ftch, "~/", hom+"/", 1)
			}
		}
	}

	// expand -stream ~/ to home directory path
	if strm != "" {

		if strm[:2] == "~/" {
			cur, err := user.Current()
			if err == nil {
				hom := cur.HomeDir
				strm = strings.Replace(strm, "~/", hom+"/", 1)
			}
		}
	}

	// expand -promote ~/ to home directory path
	if prom != "" {

		if prom[:2] == "~/" {
			cur, err := user.Current()
			if err == nil {
				hom := cur.HomeDir
				prom = strings.Replace(prom, "~/", hom+"/", 1)
			}
		}
	}

	// expand -summon ~/ to home directory path
	if smmn != "" {

		if smmn[:2] == "~/" {
			cur, err := user.Current()
			if err == nil {
				hom := cur.HomeDir
				smmn = strings.Replace(smmn, "~/", hom+"/", 1)
			}
		}
	}

	// DOCUMENTATION COMMANDS

	if len(args) > 0 {

		inSwitch = true

		switch args[0] {
		case "-version":
			fmt.Printf("%s\n", eutils.EDirectVersion)
		case "-help", "help":
			eutils.PrintHelp("rchive", "rchive-help.txt")
		case "-extras", "-extra", "-advanced":
			eutils.PrintHelp("rchive", "rchive-extras.txt")
		case "-internal", "-internals":
			eutils.PrintHelp("rchive", "rchive-internal.txt")
		default:
			// if not any of the documentation commands, keep going
			inSwitch = false
		}

		if inSwitch {
			return
		}
	}

	// FILE NAME CAN BE SUPPLIED WITH -input COMMAND

	in := os.Stdin

	// check for data being piped into stdin
	isPipe := false
	fi, staterr := os.Stdin.Stat()
	if staterr == nil {
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

	// EXTERNAL INDEXERS AND LINK ARCHIVER

	if len(args) > 0 {
		switch args[0] {
		case "-bioconcepts", "-generif", "-generifs":
			recordCount = createExternalIndexer(args, zipp, in)

			debug.FreeOSMemory()

			if timr {
				printDuration("records")
			}

			return
		case "-theme", "-themes", "-dpath", "-dpaths", "-thesis":
			recordCount = createExternalIndexer(args, zipp, in)

			debug.FreeOSMemory()

			if timr {
				printDuration("lines")
			}

			return
		default:
		}
	}

	if len(args) > 1 {
		switch args[0] {
		case "-distribute":
			args = args[1:]

			// first argument is archive path
			path := args[0]
			if path == "" {
				fmt.Fprintf(os.Stderr, "\nERROR: Need path in order to create extra index stasher\n")
				os.Exit(1)
			}
			if path[:2] == "~/" {
				cur, err := user.Current()
				if err == nil {
					hom := cur.HomeDir
					path = strings.Replace(path, "~/", hom+"/", 1)
				}
			}

			// remaining arguments are *.e2x files
			// e.g., rchive -timer -distribute archive_directory *.e2x
			args = args[1:]
			stsq := createExternalArchive(path, args)

			if stsq == nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create extra index stasher\n")
				os.Exit(1)
			}

			// drain output channel
			for str := range stsq {

				if hshv {
					// print table of UIDs and hash values
					os.Stdout.WriteString(str)
				}

				recordCount++
				runtime.Gosched()
			}

			debug.FreeOSMemory()

			if timr {
				printDuration("records")
			}

			return
		default:
		}
	}

	// JOIN SUBSETS OF INVERTED INDEX FILES

	// -join combines subsets of inverted files for subsequent -merge operation
	if join {

		// environment variable can override garbage collector (undocumented)
		gcEnv := os.Getenv("EDIRECT_JOIN_GOGC")
		if gcEnv != "" {
			val, err := strconv.Atoi(gcEnv)
			if err == nil {
				if val >= 50 && val <= 1000 {
					debug.SetGCPercent(val)
				} else {
					debug.SetGCPercent(100)
				}
			}
		} else if gcdefault {
			// default to 200 for join
			debug.SetGCPercent(200)
		}

		// environment variable can override number of servers (undocumented)
		svEnv := os.Getenv("EDIRECT_JOIN_SERV")
		if svEnv != "" {
			val, err := strconv.Atoi(svEnv)
			if err == nil {
				if val >= 1 && val <= 128 {
					numServe = val
				} else {
					numServe = 1
				}
				eutils.SetTunings(numProcs, numServe, serverRatio, chanDepth, farmSize, heapSize, goGc, false)
			}
		}

		chns := eutils.CreatePresenters(args)
		mfld := eutils.CreateManifold(chns)
		mrgr := eutils.CreateMergers(mfld)
		unsq := eutils.CreateXMLUnshuffler(mrgr)

		if chns == nil || mfld == nil || mrgr == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create inverted index joiner\n")
			os.Exit(1)
		}

		/*
			if dbug {

				// drain results, but suppress normal output
				for range unsq {
					recordCount++
					runtime.Gosched()
				}

				// force garbage collection, return memory to operating system
				debug.FreeOSMemory()

				// print processing parameters as XML object
				stopTime := time.Now()
				duration := stopTime.Sub(StartTime)
				seconds := float64(duration.Nanoseconds()) / 1e9

				// Threads is a more easily explained concept than GOMAXPROCS
				fmt.Printf("<Xtract>\n")
				fmt.Printf("  <Threads>%d</Threads>\n", numProcs)
				fmt.Printf("  <Parsers>%d</Parsers>\n", NumServe)
				fmt.Printf("  <Time>%.3f</Time>\n", seconds)
				if seconds >= 0.001 && recordCount > 0 {
					rate := int(float64(recordCount) / seconds)
					fmt.Printf("  <Rate>%d</Rate>\n", rate)
				}
				fmt.Printf("</Xtract>\n")

				return
			}
		*/

		var out io.Writer

		out = os.Stdout

		if zipp {

			zpr, err := gzip.NewWriterLevel(out, gzip.BestSpeed)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create compressor\n")
				os.Exit(1)
			}

			// close decompressor when all records have been processed
			defer zpr.Close()

			// use compressor for writing file
			out = zpr
		}

		// create buffered writer layer
		wrtr := bufio.NewWriter(out)

		wrtr.WriteString("<InvDocumentSet>\n")

		// drain channel of alphabetized results
		for curr := range unsq {

			str := curr.Text

			if str == "" {
				continue
			}

			// send result to output
			wrtr.WriteString(str)

			recordCount++
			runtime.Gosched()
		}

		wrtr.WriteString("</InvDocumentSet>\n\n")

		wrtr.Flush()

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	// MERGE INVERTED INDEX FILES AND GROUP BY TERM

	// -merge combines inverted files, distributes by prefix
	if merg != "" {

		// environment variable can override garbage collector (undocumented)
		gcEnv := os.Getenv("EDIRECT_MERGE_GOGC")
		if gcEnv != "" {
			val, err := strconv.Atoi(gcEnv)
			if err == nil {
				if val >= 50 && val <= 1000 {
					debug.SetGCPercent(val)
				} else {
					debug.SetGCPercent(100)
				}
			}
		} else if gcdefault {
			// default to 100 for merge
			debug.SetGCPercent(100)
		}

		// environment variable can override number of servers (undocumented)
		svEnv := os.Getenv("EDIRECT_MERGE_SERV")
		if svEnv != "" {
			val, err := strconv.Atoi(svEnv)
			if err == nil {
				if val >= 1 && val <= 128 {
					numServe = val
				} else {
					numServe = 1
				}
				eutils.SetTunings(numProcs, numServe, serverRatio, chanDepth, farmSize, heapSize, goGc, false)
			}
		}

		chns := eutils.CreatePresenters(args)
		mfld := eutils.CreateManifold(chns)
		mrgr := eutils.CreateMergers(mfld)
		unsq := eutils.CreateXMLUnshuffler(mrgr)
		sptr := eutils.CreateSplitter(merg, zipp, unsq)

		if chns == nil || mfld == nil || mrgr == nil || unsq == nil || sptr == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create inverted index merger\n")
			os.Exit(1)
		}

		/*
			if dbug {

				// drain results, but suppress normal output
				for range sptr {
					recordCount++
					runtime.Gosched()
				}

				// force garbage collection, return memory to operating system
				debug.FreeOSMemory()

				// print processing parameters as XML object
				stopTime := time.Now()
				duration := stopTime.Sub(StartTime)
				seconds := float64(duration.Nanoseconds()) / 1e9

				// Threads is a more easily explained concept than GOMAXPROCS
				fmt.Printf("<Xtract>\n")
				fmt.Printf("  <Threads>%d</Threads>\n", numProcs)
				fmt.Printf("  <Parsers>%d</Parsers>\n", NumServe)
				fmt.Printf("  <Time>%.3f</Time>\n", seconds)
				if seconds >= 0.001 && recordCount > 0 {
					rate := int(float64(recordCount) / seconds)
					fmt.Printf("  <Rate>%d</Rate>\n", rate)
				}
				fmt.Printf("</Xtract>\n")

				return
			}
		*/

		// drain channel, print two-to-four-character index name
		startTime := time.Now()
		first := true
		col := 0
		spaces := "       "

		for str := range sptr {

			stopTime := time.Now()
			duration := stopTime.Sub(startTime)
			seconds := float64(duration.Nanoseconds()) / 1e9

			if timr {
				if first {
					first = false
				} else {
					fmt.Fprintf(os.Stdout, "%.3f\n", seconds)
				}
				fmt.Fprintf(os.Stdout, "%s\t", str)
			} else {
				blank := 7 - len(str)
				if blank > 0 {
					fmt.Fprintf(os.Stdout, "%s", spaces[:blank])
				}
				fmt.Fprintf(os.Stdout, "%s", str)
				col++
				if col >= 10 {
					col = 0
					fmt.Fprintf(os.Stdout, "\n")
				}
			}

			recordCount++
			runtime.Gosched()

			startTime = time.Now()
		}

		stopTime := time.Now()
		duration := stopTime.Sub(startTime)
		seconds := float64(duration.Nanoseconds()) / 1e9

		if timr {
			fmt.Fprintf(os.Stdout, "%.3f\n", seconds)
		} else if col > 0 {
			fmt.Fprintf(os.Stdout, "\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("groups")
		}

		return
	}

	// PROMOTE MERGED INVERTED INDEX TO TERM LIST AND POSTINGS FILES

	if prom != "" && fild != "" {

		prmq := eutils.CreatePromoters(prom, fild, args)

		if prmq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create new postings file generator\n")
			os.Exit(1)
		}

		col := 0
		spaces := "       "

		// drain channel, print 2-4 character file prefix
		for str := range prmq {

			blank := 7 - len(str)
			if blank > 0 {
				fmt.Fprintf(os.Stdout, "%s", spaces[:blank])
			}
			fmt.Fprintf(os.Stdout, "%s", str)
			col++
			if col >= 10 {
				col = 0
				fmt.Fprintf(os.Stdout, "\n")
			}

			recordCount++
			runtime.Gosched()
		}

		if col > 0 {
			fmt.Fprintf(os.Stdout, "\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	// QUERY POSTINGS FILES

	if phrs != "" || trms != "" || btch {
		if base == "" {
			// obtain path from environment variable within rchive as a convenience
			base = os.Getenv("EDIRECT_PUBMED_MASTER")
			if base != "" {
				if !strings.HasSuffix(base, "/") {
					base += "/"
				}
				base += "Postings"
			}
		}
	}

	if base != "" && btch {

		// read query lines for exact match
		scanr := bufio.NewScanner(in)

		for scanr.Scan() {
			txt := scanr.Text()

			// deStop should match value used in building the indices
			recordCount += eutils.ProcessSearch(base, txt, true, false, false, deStop)
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	if base != "" && phrs != "" {

		// deStop should match value used in building the indices
		if mock {
			recordCount = eutils.ProcessMock(base, phrs, xact, titl, rlxd, deStop)
		} else {
			recordCount = eutils.ProcessSearch(base, phrs, xact, titl, rlxd, deStop)
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	if base != "" && trms != "" {

		// deStop should match value used in building the indices
		recordCount = eutils.ProcessCount(base, trms, plrl, psns, rlxd, deStop)

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	// CONFIRM INPUT DATA AVAILABILITY AFTER RUNNING COMMAND GENERATORS

	if fileName == "" && runtime.GOOS != "windows" {

		fromStdin := bool((fi.Mode() & os.ModeCharDevice) == 0)
		if !isPipe || !fromStdin {
			mode := fi.Mode().String()
			fmt.Fprintf(os.Stderr, "\nERROR: No data supplied to rchive from stdin or file, mode is '%s'\n", mode)
			os.Exit(1)
		}
	}

	if !usingFile && !isPipe {

		fmt.Fprintf(os.Stderr, "\nERROR: No XML input data supplied to rchive\n")
		os.Exit(1)
	}

	// SPECIFY STRINGS TO GO BEFORE AND AFTER ENTIRE OUTPUT OR EACH RECORD

	head := ""
	tail := ""

	hd := ""
	tl := ""

	for {

		if len(args) < 1 {
			break
		}

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
		default:
			// if not any of the controls, set flag to break out of for loop
			inSwitch = false
		}

		if !inSwitch {
			break
		}

		// skip past arguments
		args = args[2:]
	}

	// PRODUCE ARCHIVE SUBPATH FROM IDENTIFIER

	// -trie converts identifier to directory subpath plus file name (undocumented)
	if trei {

		scanr := bufio.NewScanner(in)

		sfx := ".xml"
		if zipp {
			sfx += ".gz"
		}

		// read lines of identifiers
		for scanr.Scan() {

			file := scanr.Text()

			var arry [132]rune
			trie := eutils.MakeArchiveTrie(file, arry)
			if trie == "" || file == "" {
				continue
			}

			fpath := filepath.Join(trie, file+sfx)
			if fpath == "" {
				continue
			}

			os.Stdout.WriteString(fpath)
			os.Stdout.WriteString("\n")
		}

		return
	}

	// CHECK FOR MISSING RECORDS IN LOCAL DIRECTORY INDEXED BY TRIE ON IDENTIFIER

	// -archive plus -missing checks for missing records
	if stsh != "" && msng {

		scanr := bufio.NewScanner(in)

		sfx := ".xml"
		if zipp {
			sfx += ".gz"
		}

		// read lines of identifiers
		for scanr.Scan() {

			file := scanr.Text()

			pos := strings.Index(file, ".")
			if pos >= 0 {
				// remove version suffix
				file = file[:pos]
			}

			var arry [132]rune
			trie := eutils.MakeArchiveTrie(file, arry)

			if file == "" || trie == "" {
				continue
			}

			fpath := filepath.Join(stsh, trie, file+sfx)
			if fpath == "" {
				continue
			}

			_, err := os.Stat(fpath)

			// if failed to find ".xml" file, try ".xml.gz" without requiring -gzip
			if err != nil && os.IsNotExist(err) && !zipp {
				fpath := filepath.Join(stsh, trie, file+".xml.gz")
				if fpath == "" {
					continue
				}
				_, err = os.Stat(fpath)
			}
			if err != nil && os.IsNotExist(err) {
				// record is missing from local file cache
				os.Stdout.WriteString(file)
				os.Stdout.WriteString("\n")
			}
		}

		return
	}

	// RETRIEVE XML COMPONENT RECORDS FROM LOCAL DIRECTORY INDEXED BY TRIE ON IDENTIFIER

	// alternative windows version limits memory by not using goroutines
	if ftch != "" && indx == "" && runtime.GOOS == "windows" && windows {

		scanr := bufio.NewScanner(in)
		if scanr == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create UID scanner\n")
			os.Exit(1)
		}

		sfx := ".xml"
		if zipp {
			sfx += ".gz"
		}

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		var buf bytes.Buffer

		for scanr.Scan() {

			// read next identifier
			file := scanr.Text()

			pos := strings.Index(file, ".")
			if pos >= 0 {
				// remove version suffix
				file = file[:pos]
			}

			var arry [132]rune
			trie := eutils.MakeArchiveTrie(file, arry)

			if file == "" || trie == "" {
				continue
			}

			fpath := filepath.Join(ftch, trie, file+sfx)
			if fpath == "" {
				continue
			}

			iszip := zipp

			inFile, err := os.Open(fpath)

			// if failed to find ".xml" file, try ".xml.gz" without requiring -gzip
			if err != nil && os.IsNotExist(err) && !zipp {
				iszip = true
				fpath := filepath.Join(ftch, trie, file+".xml.gz")
				if fpath == "" {
					continue
				}
				inFile, err = os.Open(fpath)
			}
			if err != nil {
				continue
			}

			buf.Reset()

			brd := bufio.NewReader(inFile)

			if iszip {

				zpr, err := gzip.NewReader(brd)

				if err == nil {
					// copy and decompress cached file contents
					buf.ReadFrom(zpr)
				}

				zpr.Close()

			} else {

				// copy cached file contents
				buf.ReadFrom(brd)
			}

			inFile.Close()

			str := buf.String()

			if str == "" {
				continue
			}

			recordCount++

			if hd != "" {
				os.Stdout.WriteString(hd)
				os.Stdout.WriteString("\n")
			}

			if hshv {
				// calculate hash code for verification table
				hsh := crc32.NewIEEE()
				hsh.Write([]byte(str))
				val := hsh.Sum32()
				res := strconv.FormatUint(uint64(val), 10)
				txt := file + "\t" + res + "\n"
				os.Stdout.WriteString(txt)
			} else {
				// send result to output
				os.Stdout.WriteString(str)
				if !strings.HasSuffix(str, "\n") {
					os.Stdout.WriteString("\n")
				}
			}

			if tl != "" {
				os.Stdout.WriteString(tl)
				os.Stdout.WriteString("\n")
			}

			debug.FreeOSMemory()
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

	// -fetch without -index retrieves XML files in trie-based directory structure
	if ftch != "" && indx == "" {

		uidq := eutils.CreateUIDReader(in)
		strq := eutils.CreateFetchers(ftch, ".xml", zipp, uidq)
		unsq := eutils.CreateXMLUnshuffler(strq)

		if uidq == nil || strq == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create archive reader\n")
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

			if hshv {
				// calculate hash code for verification table
				hsh := crc32.NewIEEE()
				hsh.Write([]byte(str))
				val := hsh.Sum32()
				res := strconv.FormatUint(uint64(val), 10)
				txt := curr.Ident + "\t" + res + "\n"
				os.Stdout.WriteString(txt)
			} else {
				// send result to output
				os.Stdout.WriteString(str)
				if !strings.HasSuffix(str, "\n") {
					os.Stdout.WriteString("\n")
				}
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

	// -stream without -index retrieves compressed XML files in trie-based directory structure
	if strm != "" && indx == "" {

		uidq := eutils.CreateUIDReader(in)
		strq := eutils.CreateCacheStreamers(strm, uidq)
		unsq := eutils.CreateXMLUnshuffler(strq)

		if uidq == nil || strq == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create archive reader\n")
			os.Exit(1)
		}

		// drain output channel
		for curr := range unsq {

			data := curr.Data

			if data == nil {
				continue
			}

			recordCount++
			runtime.Gosched()

			_, err := os.Stdout.Write(data)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			}
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// -summon retrieves link files in trie-based directory structure
	if smmn != "" && indx == "" {

		uidq := eutils.CreateUIDReader(in)
		strq := eutils.CreateFetchers(smmn, ".e2x", zipp, uidq)
		unsq := eutils.CreateXMLUnshuffler(strq)

		if uidq == nil || strq == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create link reader\n")
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

			if hshv {
				// calculate hash code for verification table
				hsh := crc32.NewIEEE()
				hsh.Write([]byte(str))
				val := hsh.Sum32()
				res := strconv.FormatUint(uint64(val), 10)
				txt := curr.Ident + "\t" + res + "\n"
				os.Stdout.WriteString(txt)
			} else {
				// send result to output
				os.Stdout.WriteString(str)
				if !strings.HasSuffix(str, "\n") {
					os.Stdout.WriteString("\n")
				}
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

	// CREATE XML BLOCK READER FROM STDIN OR FILE

	rdr := eutils.CreateXMLStreamer(in)
	if rdr == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
		os.Exit(1)
	}

	// ENTREZ INDEX INVERSION

	// -invert reads IdxDocumentSet XML and creates an inverted index
	if nvrt {

		// environment variable can override garbage collector (undocumented)
		gcEnv := os.Getenv("EDIRECT_INVERT_GOGC")
		if gcEnv != "" {
			val, err := strconv.Atoi(gcEnv)
			if err == nil {
				if val >= 50 && val <= 1000 {
					debug.SetGCPercent(val)
				} else {
					debug.SetGCPercent(100)
				}
			}
		}

		// environment variable can override number of servers (undocumented)
		svEnv := os.Getenv("EDIRECT_INVERT_SERV")
		if svEnv != "" {
			val, err := strconv.Atoi(svEnv)
			if err == nil {
				if val >= 1 && val <= 128 {
					numServe = val
				} else {
					numServe = 1
				}
				eutils.SetTunings(numProcs, numServe, serverRatio, chanDepth, farmSize, heapSize, goGc, false)
			}
		}

		colq := eutils.CreateXMLProducer("IdxDocument", "", false, rdr)
		dspq := eutils.CreateDispensers(colq)
		invq := eutils.CreateInverters(dspq)
		rslq := eutils.CreateResolver(invq)

		if colq == nil || dspq == nil || invq == nil || rslq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create inverter\n")
			os.Exit(1)
		}

		/*
			if dbug {

				// drain results, but suppress normal output
				for range rslq {
					recordCount++
					runtime.Gosched()
				}

				// force garbage collection, return memory to operating system
				debug.FreeOSMemory()

				// print processing parameters as XML object
				stopTime := time.Now()
				duration := stopTime.Sub(StartTime)
				seconds := float64(duration.Nanoseconds()) / 1e9

				// Threads is a more easily explained concept than GOMAXPROCS
				fmt.Printf("<Xtract>\n")
				fmt.Printf("  <Threads>%d</Threads>\n", numProcs)
				fmt.Printf("  <Parsers>%d</Parsers>\n", NumServe)
				fmt.Printf("  <Time>%.3f</Time>\n", seconds)
				if seconds >= 0.001 && recordCount > 0 {
					rate := int(float64(recordCount) / seconds)
					fmt.Printf("  <Rate>%d</Rate>\n", rate)
				}
				fmt.Printf("</Xtract>\n")

				return
			}
		*/

		var out io.Writer

		out = os.Stdout

		if zipp {

			zpr, err := gzip.NewWriterLevel(out, gzip.BestSpeed)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create compressor\n")
				os.Exit(1)
			}

			// close decompressor when all records have been processed
			defer zpr.Close()

			// use compressor for writing file
			out = zpr
		}

		// create buffered writer layer
		wrtr := bufio.NewWriter(out)

		wrtr.WriteString("<InvDocumentSet>\n")

		// drain channel of alphabetized results
		for str := range rslq {

			// send result to output
			wrtr.WriteString(str)

			recordCount++
			runtime.Gosched()
		}

		wrtr.WriteString("</InvDocumentSet>\n\n")

		wrtr.Flush()

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	// FUSE SUBSETS OF INVERTED INDEX FILES

	// -fuse combines subsets of inverted files for subsequent -merge operation
	if fuse {

		// environment variable can override garbage collector (undocumented)
		gcEnv := os.Getenv("EDIRECT_FUSE_GOGC")
		if gcEnv != "" {
			val, err := strconv.Atoi(gcEnv)
			if err == nil {
				if val >= 50 && val <= 1000 {
					debug.SetGCPercent(val)
				} else {
					debug.SetGCPercent(100)
				}
			}
		} else if gcdefault {
			// default to 100 for fuse and merge
			debug.SetGCPercent(100)
		}

		// environment variable can override number of servers (undocumented)
		svEnv := os.Getenv("EDIRECT_FUSE_SERV")
		if svEnv != "" {
			val, err := strconv.Atoi(svEnv)
			if err == nil {
				if val >= 1 && val <= 128 {
					numServe = val
				} else {
					numServe = 1
				}
				eutils.SetTunings(numProcs, numServe, serverRatio, chanDepth, farmSize, heapSize, goGc, false)
			}
		}

		chns := eutils.CreateXMLProducer("InvDocument", "", false, rdr)
		fusr := eutils.CreateFusers(chns)
		mrgr := eutils.CreateMergers(fusr)
		unsq := eutils.CreateXMLUnshuffler(mrgr)

		if chns == nil || fusr == nil || mrgr == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create inverted index fuser\n")
			os.Exit(1)
		}

		var out io.Writer

		out = os.Stdout

		if zipp {

			zpr, err := gzip.NewWriterLevel(out, gzip.BestSpeed)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create compressor\n")
				os.Exit(1)
			}

			// close decompressor when all records have been processed
			defer zpr.Close()

			// use compressor for writing file
			out = zpr
		}

		// create buffered writer layer
		wrtr := bufio.NewWriter(out)

		wrtr.WriteString("<InvDocumentSet>\n")

		// drain channel of alphabetized results
		for curr := range unsq {

			str := curr.Text

			if str == "" {
				continue
			}

			// send result to output
			wrtr.WriteString(str)

			recordCount++
			runtime.Gosched()
		}

		wrtr.WriteString("</InvDocumentSet>\n\n")

		wrtr.Flush()

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	// ENSURE PRESENCE OF PATTERN ARGUMENT

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to rchive\n")
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

	parent := ""
	if star == "*" {
		parent = topPattern
	} else if star != "" {
		fmt.Fprintf(os.Stderr, "\nERROR: -pattern Parent/Child construct is not supported\n")
		os.Exit(1)
	}

	// REPORT RECORDS THAT CONTAIN DAMAGED EMBEDDED HTML TAGS

	// -damaged plus -index plus -pattern reports records with multiply-encoded HTML tags
	if dmgd && indx != "" {

		find := eutils.ParseIndex(indx)

		eutils.PartitionPattern(topPattern, star, false, rdr,
			func(str string) {
				recordCount++

				id := eutils.FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}

				// remove default version suffix
				if strings.HasSuffix(id, ".1") {
					idlen := len(id)
					id = id[:idlen-2]
				}

				reportEncodedMarkup(dmgdType, id, str)
			})

		if timr {
			printDuration("records")
		}

		return
	}

	// COMPARE XML UPDATES TO LOCAL DIRECTORY, RETAIN NEW OR SUBSTANTIVELY CHANGED RECORDS

	// -prepare plus -archive plus -index plus -pattern compares XML files against stash
	if stsh != "" && indx != "" && cmpr {

		doReport := false
		if cmprType == "" || cmprType == "report" {
			doReport = true
		} else if cmprType != "release" {
			fmt.Fprintf(os.Stderr, "\nERROR: -prepare argument must be release or report\n")
			os.Exit(1)
		}

		find := eutils.ParseIndex(indx)

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		eutils.PartitionPattern(topPattern, star, false, rdr,
			func(str string) {
				recordCount++

				id := eutils.FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}

				pos := strings.Index(id, ".")
				if pos >= 0 {
					// remove version suffix
					id = id[:pos]
				}

				var arry [132]rune
				trie := eutils.MakeArchiveTrie(id, arry)

				if id == "" || trie == "" {
					return
				}

				fpath := filepath.Join(stsh, trie, id+".xml")
				if fpath == "" {
					return
				}

				// print new or updated XML record
				printRecord := func(stn string, isNew bool) {

					if stn == "" {
						return
					}

					if doReport {
						if isNew {
							os.Stdout.WriteString("NW ")
							os.Stdout.WriteString(id)
							os.Stdout.WriteString("\n")
						} else {
							os.Stdout.WriteString("UP ")
							os.Stdout.WriteString(id)
							os.Stdout.WriteString("\n")
						}
						return
					}

					if hd != "" {
						os.Stdout.WriteString(hd)
						os.Stdout.WriteString("\n")
					}

					os.Stdout.WriteString(stn)
					os.Stdout.WriteString("\n")

					if tl != "" {
						os.Stdout.WriteString(tl)
						os.Stdout.WriteString("\n")
					}
				}

				_, err := os.Stat(fpath)
				if err != nil && os.IsNotExist(err) {
					// new record
					printRecord(str, true)
					return
				}
				if err != nil {
					return
				}

				buf, err := ioutil.ReadFile(fpath)
				if err != nil {
					return
				}

				txt := string(buf[:])
				txt = strings.TrimSuffix(txt, "\n")

				// check for optional -ignore argument
				if ignr != "" {

					// ignore differences inside specified object
					ltag := "<" + ignr + ">"
					sleft, _ := eutils.SplitInTwoLeft(str, ltag)
					tleft, _ := eutils.SplitInTwoLeft(txt, ltag)

					rtag := "</" + ignr + ">"
					_, srght := eutils.SplitInTwoRight(str, rtag)
					_, trght := eutils.SplitInTwoRight(txt, rtag)

					if sleft == tleft && srght == trght {
						if doReport {
							os.Stdout.WriteString("NO ")
							os.Stdout.WriteString(id)
							os.Stdout.WriteString("\n")
						}
						return
					}

				} else {

					// compare entirety of objects
					if str == txt {
						if doReport {
							os.Stdout.WriteString("NO ")
							os.Stdout.WriteString(id)
							os.Stdout.WriteString("\n")
						}
						return
					}
				}

				// substantively modified record
				printRecord(str, false)
			})

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		if timr {
			printDuration("records")
		}

		return
	}

	// SAVE XML COMPONENT RECORDS TO LOCAL DIRECTORY INDEXED BY TRIE ON IDENTIFIER

	// -archive plus -index plus -pattern saves XML files in trie-based directory structure
	if stsh != "" && indx != "" {

		xmlq := eutils.CreateXMLProducer(topPattern, star, false, rdr)
		stsq := eutils.CreateStashers(stsh, parent, indx, ".xml", hshv, zipp, 1000, xmlq)

		if xmlq == nil || stsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create stash generator\n")
			os.Exit(1)
		}

		// drain output channel
		for str := range stsq {

			if hshv {
				// print table of UIDs and hash values
				os.Stdout.WriteString(str)
			}

			recordCount++
			runtime.Gosched()
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// READ FILE OF IDENTIFIERS AND EXTRACT SELECTED RECORDS FROM XML INPUT FILE

	// -index plus -unique [plus -head/-tail/-hd/-tl] plus -pattern with no other extraction arguments
	// takes an XML input file and a file of its UIDs and keeps only the last version of each record
	if indx != "" && unqe != "" && len(args) == 2 {

		// read file of identifiers to use for filtering
		fl, err := os.Open(unqe)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open identifier file '%s'\n", unqe)
			os.Exit(1)
		}

		// create map that counts instances of each UID
		order := make(map[string]int)

		scanr := bufio.NewScanner(fl)

		// read lines of identifiers
		for scanr.Scan() {

			id := scanr.Text()

			// map records count for given identifier
			val := order[id]
			val++
			order[id] = val
		}

		fl.Close()

		find := eutils.ParseIndex(indx)

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		eutils.PartitionPattern(topPattern, star, false, rdr,
			func(str string) {
				recordCount++

				id := eutils.FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}

				val, ok := order[id]
				if !ok {
					// not in identifier list, skip
					return
				}
				// decrement count in map
				val--
				order[id] = val
				if val > 0 {
					// only write last record with a given identifier
					return
				}

				if hd != "" {
					os.Stdout.WriteString(hd)
					os.Stdout.WriteString("\n")
				}

				// write selected record
				os.Stdout.WriteString(str[:])
				os.Stdout.WriteString("\n")

				if tl != "" {
					os.Stdout.WriteString(tl)
					os.Stdout.WriteString("\n")
				}
			})

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		if timr {
			printDuration("records")
		}

		return
	}

	// GENERATE RECORD INDEX ON XML INPUT FILE

	// -index plus -pattern prints record identifier and XML size
	if indx != "" {

		lbl := ""
		// check for optional filename label after -pattern argument (undocumented)
		if len(args) > 3 && args[2] == "-lbl" {
			lbl = args[3]

			lbl = strings.TrimSpace(lbl)
			if strings.HasPrefix(lbl, "pubmed") {
				lbl = lbl[7:]
			}
			if strings.HasSuffix(lbl, ".xml.gz") {
				xlen := len(lbl)
				lbl = lbl[:xlen-7]
			}
			lbl = strings.TrimSpace(lbl)
		}

		// legend := "ID\tREC\tSIZE"

		find := eutils.ParseIndex(indx)

		eutils.PartitionPattern(topPattern, star, false, rdr,
			func(str string) {
				recordCount++

				id := eutils.FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}
				if lbl != "" {
					fmt.Printf("%s\t%d\t%s\n", id, len(str), lbl)
				} else {
					fmt.Printf("%s\t%d\n", id, len(str))
				}
			})

		if timr {
			printDuration("records")
		}

		return
	}

	// SORT XML RECORDS BY IDENTIFIER

	// -pattern record_name -sort parent/element@attribute^version, strictly alphabetic sort order (undocumented)
	if len(args) == 4 && args[2] == "-sort" {

		indx := args[3]

		// create map that records each UID
		order := make(map[string][]string)

		find := eutils.ParseIndex(indx)

		eutils.PartitionPattern(topPattern, star, false, rdr,
			func(str string) {
				recordCount++

				id := eutils.FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}

				data, ok := order[id]
				if !ok {
					data = make([]string, 0, 1)
				}
				data = append(data, str)
				// always need to update order, since data may be reallocated
				order[id] = data
			})

		var keys []string
		for ky := range order {
			keys = append(keys, ky)
		}
		// sort fields in alphabetical order, unlike xtract version, which sorts numbers by numeric value
		sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		for _, id := range keys {

			strs := order[id]
			for _, str := range strs {
				os.Stdout.WriteString(str)
				os.Stdout.WriteString("\n")
			}
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

	fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized rchive command\n")
	os.Exit(1)
}
