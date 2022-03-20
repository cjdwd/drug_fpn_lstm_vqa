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
// File Name:  gbff.go
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
	"strings"
)

// GenBankConverter reads flatfiles and sends INSDSeq XML records down a channel
func GenBankConverter(inp io.Reader) <-chan string {

	if inp == nil {
		return nil
	}

	out := make(chan string, chanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "Unable to create GenBank converter channel\n")
		os.Exit(1)
	}

	const twelvespaces = "            "
	const twentyonespaces = "                     "

	var rec strings.Builder
	var alt strings.Builder
	var con strings.Builder
	var seq strings.Builder

	scanr := bufio.NewScanner(inp)

	convertGenBank := func(inp io.Reader, out chan<- string) {

		// close channel when all records have been sent
		defer close(out)

		row := 0

		nextLine := func() string {

			for scanr.Scan() {
				line := scanr.Text()
				if line == "" {
					continue
				}
				return line
			}
			return ""

		}

		for {

			rec.Reset()

			// read first line of next record
			line := nextLine()
			if line == "" {
				break
			}

			row++

			for {
				if line == "" {
					break
				}
				if !strings.HasPrefix(line, "LOCUS") {
					// skip release file header information
					line = nextLine()
					row++
					continue
				}
				break
			}

			readContinuationLines := func(str string) string {

				for {
					// read next line
					line = nextLine()
					row++
					if !strings.HasPrefix(line, twelvespaces) {
						// if not continuation line, break out of loop
						break
					}
					// append subsequent line and continue with loop
					txt := strings.TrimPrefix(line, twelvespaces)
					str += " " + txt
				}

				str = CompressRunsOfSpaces(str)
				str = strings.TrimSpace(str)

				return str
			}

			writeOneElement := func(spaces, tag, value string) {

				rec.WriteString(spaces)
				rec.WriteString("<")
				rec.WriteString(tag)
				rec.WriteString(">")
				value = html.EscapeString(value)
				rec.WriteString(value)
				rec.WriteString("</")
				rec.WriteString(tag)
				rec.WriteString(">\n")
			}

			// each section will exit with the next line ready to process

			if strings.HasPrefix(line, "LOCUS") {

				// start of record
				rec.WriteString("  <INSDSeq>\n")

				// do not break if given artificial multi-line LOCUS
				str := readContinuationLines(line)

				cols := strings.Fields(str)
				ln := len(cols)
				if ln == 8 {
					moleculetype := cols[4]
					strandedness := ""
					if strings.HasPrefix(moleculetype, "ds-") {
						moleculetype = strings.TrimPrefix(moleculetype, "ds-")
						strandedness = "double"
					} else if strings.HasPrefix(moleculetype, "ss-") {
						moleculetype = strings.TrimPrefix(moleculetype, "ss-")
						strandedness = "single"
					} else if strings.HasPrefix(moleculetype, "ms-") {
						moleculetype = strings.TrimPrefix(moleculetype, "ms-")
						strandedness = "mixed"
					} else if strings.HasSuffix(moleculetype, "DNA") {
						strandedness = "double"
					} else if strings.HasSuffix(moleculetype, "RNA") {
						strandedness = "single"
					}

					writeOneElement("    ", "INSDSeq_locus", cols[1])
					writeOneElement("    ", "INSDSeq_length", cols[2])

					if strandedness != "" {
						writeOneElement("    ", "INSDSeq_strandedness", strandedness)
					}

					writeOneElement("    ", "INSDSeq_moltype", moleculetype)
					writeOneElement("    ", "INSDSeq_topology", cols[5])
					writeOneElement("    ", "INSDSeq_division", cols[6])
					writeOneElement("    ", "INSDSeq_update-date", cols[7])

				} else if ln == 7 {

					writeOneElement("    ", "INSDSeq_locus", cols[1])
					writeOneElement("    ", "INSDSeq_length", cols[2])
					writeOneElement("    ", "INSDSeq_moltype", "AA")
					writeOneElement("    ", "INSDSeq_topology", cols[4])
					writeOneElement("    ", "INSDSeq_division", cols[5])
					writeOneElement("    ", "INSDSeq_update-date", cols[6])

				} else {
					fmt.Fprintf(os.Stderr, "ERROR: "+str+"\n")
				}

				// read next line and continue - handled by readContinuationLines above
				// line = nextLine()
				// row++
			}

			if strings.HasPrefix(line, "DEFINITION") {

				txt := strings.TrimPrefix(line, "DEFINITION")
				def := readContinuationLines(txt)
				def = strings.TrimSuffix(def, ".")

				writeOneElement("    ", "INSDSeq_definition", def)
			}

			var secondaries []string

			if strings.HasPrefix(line, "ACCESSION") {

				txt := strings.TrimPrefix(line, "ACCESSION")
				str := readContinuationLines(txt)
				accessions := strings.Fields(str)
				ln := len(accessions)
				if ln > 1 {

					writeOneElement("    ", "INSDSeq_primary-accession", accessions[0])

					// skip past primary accession, collect secondaries
					secondaries = accessions[1:]

				} else if ln == 1 {

					writeOneElement("    ", "INSDSeq_primary-accession", accessions[0])

				} else {
					fmt.Fprintf(os.Stderr, "ERROR: ACCESSION "+str+"\n")
				}
			}

			accnver := ""
			gi := ""

			if strings.HasPrefix(line, "VERSION") {

				cols := strings.Fields(line)
				if len(cols) == 2 {

					accnver = cols[1]
					writeOneElement("    ", "INSDSeq_accession-version", accnver)

				} else if len(cols) == 3 {

					accnver = cols[1]
					writeOneElement("    ", "INSDSeq_accession-version", accnver)

					// collect gi for other-seqids
					if strings.HasPrefix(cols[2], "GI:") {
						gi = strings.TrimPrefix(cols[2], "GI:")
					}

				} else {
					fmt.Fprintf(os.Stderr, "ERROR: "+line+"\n")
				}

				// read next line and continue
				line = nextLine()
				row++

			}

			if gi != "" {

				rec.WriteString("    <INSDSeq_other-seqids>\n")
				writeOneElement("      ", "INSDSeqid", "gi|"+gi)
				rec.WriteString("    </INSDSeq_other-seqids>\n")
			}

			if len(secondaries) > 0 {

				rec.WriteString("    <INSDSeq_secondary-accessions>\n")

				for _, secndry := range secondaries {

					if strings.HasPrefix(secndry, "REGION") {
						break
					}
					writeOneElement("      ", "INSDSecondary-accn", secndry)
				}

				rec.WriteString("    </INSDSeq_secondary-accessions>\n")

				// reset secondaries slice
				secondaries = nil
			}

			dblink := ""

			if strings.HasPrefix(line, "DBLINK") {

				txt := strings.TrimPrefix(line, "DBLINK")
				dblink = readContinuationLines(txt)
			}

			srcdb := ""

			if strings.HasPrefix(line, "DBSOURCE") {

				txt := strings.TrimPrefix(line, "DBSOURCE")
				srcdb = readContinuationLines(txt)
			}

			if strings.HasPrefix(line, "KEYWORDS") {

				txt := strings.TrimPrefix(line, "KEYWORDS")
				key := readContinuationLines(txt)
				key = strings.TrimSuffix(key, ".")

				if key != "" {
					rec.WriteString("    <INSDSeq_keywords>\n")
					kywds := strings.Split(key, ";")
					for _, kw := range kywds {
						kw = strings.TrimSpace(kw)
						if kw == "" || kw == "." {
							continue
						}

						writeOneElement("      ", "INSDKeyword", kw)
					}
					rec.WriteString("    </INSDSeq_keywords>\n")
				}
			}

			if strings.HasPrefix(line, "SOURCE") {

				txt := strings.TrimPrefix(line, "SOURCE")
				src := readContinuationLines(txt)

				writeOneElement("    ", "INSDSeq_source", src)
			}

			if strings.HasPrefix(line, "  ORGANISM") {

				org := strings.TrimPrefix(line, "  ORGANISM")
				org = CompressRunsOfSpaces(org)
				org = strings.TrimSpace(org)

				writeOneElement("    ", "INSDSeq_organism", org)

				line = nextLine()
				row++
				if strings.HasPrefix(line, twelvespaces) {
					txt := strings.TrimPrefix(line, twelvespaces)
					tax := readContinuationLines(txt)
					tax = strings.TrimSuffix(tax, ".")

					writeOneElement("    ", "INSDSeq_taxonomy", tax)
				}
			}

			rec.WriteString("    <INSDSeq_references>\n")
			for {
				if !strings.HasPrefix(line, "REFERENCE") {
					// exit out of reference section
					break
				}

				ref := "0"

				rec.WriteString("      <INSDReference>\n")

				txt := strings.TrimPrefix(line, "REFERENCE")
				str := readContinuationLines(txt)
				str = CompressRunsOfSpaces(str)
				str = strings.TrimSpace(str)
				idx := strings.Index(str, "(")
				if idx > 0 {
					ref = strings.TrimSpace(str[:idx])

					writeOneElement("        ", "INSDReference_reference", ref)

					posn := str[idx+1:]
					posn = strings.TrimSuffix(posn, ")")
					posn = strings.TrimSpace(posn)
					if posn == "sites" {

						writeOneElement("        ", "INSDReference_position", posn)

					} else {
						var arry []string
						cls := strings.Split(posn, ";")
						for _, item := range cls {
							item = strings.TrimPrefix(item, "bases ")
							item = strings.TrimPrefix(item, "residues ")
							item = strings.TrimSpace(item)
							cols := strings.Fields(item)
							if len(cols) == 3 && cols[1] == "to" {
								arry = append(arry, cols[0]+".."+cols[2])
							}
						}
						if len(arry) > 0 {
							posit := strings.Join(arry, ",")
							writeOneElement("        ", "INSDReference_position", posit)
						} else {
							fmt.Fprintf(os.Stderr, "ERROR: "+posn+"\n")
						}
					}
				} else {
					ref = strings.TrimSpace(str)

					writeOneElement("        ", "INSDReference_reference", ref)
				}
				row++

				if strings.HasPrefix(line, "  AUTHORS") {

					txt := strings.TrimPrefix(line, "  AUTHORS")
					auths := readContinuationLines(txt)

					rec.WriteString("        <INSDReference_authors>\n")
					authors := strings.Split(auths, ", ")
					for _, auth := range authors {
						auth = strings.TrimSpace(auth)
						if auth == "" {
							continue
						}
						pair := strings.Split(auth, " and ")
						for _, name := range pair {

							writeOneElement("          ", "INSDAuthor", name)
						}
					}
					rec.WriteString("        </INSDReference_authors>\n")
				}

				if strings.HasPrefix(line, "  CONSRTM") {

					txt := strings.TrimPrefix(line, "  CONSRTM")
					cons := readContinuationLines(txt)

					writeOneElement("        ", "INSDReference_consortium", cons)
				}

				if strings.HasPrefix(line, "  TITLE") {

					txt := strings.TrimPrefix(line, "  TITLE")
					titl := readContinuationLines(txt)

					writeOneElement("        ", "INSDReference_title", titl)
				}

				if strings.HasPrefix(line, "  JOURNAL") {

					txt := strings.TrimPrefix(line, "  JOURNAL")
					jour := readContinuationLines(txt)

					writeOneElement("        ", "INSDReference_journal", jour)
				}

				if strings.HasPrefix(line, "   PUBMED") {

					txt := strings.TrimPrefix(line, "   PUBMED")
					pmid := readContinuationLines(txt)

					writeOneElement("        ", "INSDReference_pubmed", pmid)
				}

				if strings.HasPrefix(line, "  MEDLINE") {

					txt := strings.TrimPrefix(line, "  MEDLINE")
					// old MEDLINE uid not supported
					readContinuationLines(txt)
				}

				if strings.HasPrefix(line, "  REMARK") {

					txt := strings.TrimPrefix(line, "  REMARK")
					rem := readContinuationLines(txt)

					writeOneElement("        ", "INSDReference_remark", rem)
				}

				// end of this reference
				rec.WriteString("      </INSDReference>\n")
				// continue to next reference
			}
			rec.WriteString("    </INSDSeq_references>\n")

			if strings.HasPrefix(line, "COMMENT") {

				txt := strings.TrimPrefix(line, "COMMENT")
				com := readContinuationLines(txt)

				writeOneElement("    ", "INSDSeq_comment", com)
			}

			if strings.HasPrefix(line, "PRIMARY") {

				txt := strings.TrimPrefix(line, "PRIMARY")
				pmy := readContinuationLines(txt)

				writeOneElement("    ", "INSDSeq_primary", pmy)
			}

			if srcdb != "" {
				writeOneElement("    ", "INSDSeq_source-db", srcdb)
			}

			rec.WriteString("    <INSDSeq_feature-table>\n")
			if strings.HasPrefix(line, "FEATURES") {

				line = nextLine()
				row++

				for {
					if !strings.HasPrefix(line, "     ") {
						// exit out of features section
						break
					}
					if len(line) < 22 {
						fmt.Fprintf(os.Stderr, "ERROR: "+line+"\n")
						line = nextLine()
						row++
						continue
					}

					rec.WriteString("      <INSDFeature>\n")

					// read feature key and start of location
					fkey := line[5:21]
					fkey = strings.TrimSpace(fkey)

					writeOneElement("        ", "INSDFeature_key", fkey)

					loc := line[21:]
					loc = strings.TrimSpace(loc)
					for {
						line = nextLine()
						row++
						if !strings.HasPrefix(line, twentyonespaces) {
							break
						}
						txt := strings.TrimPrefix(line, twentyonespaces)
						if strings.HasPrefix(txt, "/") {
							// if not continuation of location, break out of loop
							break
						}
						// append subsequent line and continue with loop
						loc += strings.TrimSpace(txt)
					}

					writeOneElement("        ", "INSDFeature_location", loc)

					locationOperator := ""
					isComp := false
					prime5 := false
					prime3 := false

					// parseloc recursive definition
					var parseloc func(string) []string

					parseloc = func(str string) []string {

						var acc []string

						if strings.HasPrefix(str, "join(") && strings.HasSuffix(str, ")") {

							locationOperator = "join"

							str = strings.TrimPrefix(str, "join(")
							str = strings.TrimSuffix(str, ")")
							items := strings.Split(str, ",")

							for _, thisloc := range items {
								inner := parseloc(thisloc)
								for _, sub := range inner {
									acc = append(acc, sub)
								}
							}

						} else if strings.HasPrefix(str, "order(") && strings.HasSuffix(str, ")") {

							locationOperator = "order"

							str = strings.TrimPrefix(str, "order(")
							str = strings.TrimSuffix(str, ")")
							items := strings.Split(str, ",")

							for _, thisloc := range items {
								inner := parseloc(thisloc)
								for _, sub := range inner {
									acc = append(acc, sub)
								}
							}

						} else if strings.HasPrefix(str, "complement(") && strings.HasSuffix(str, ")") {

							isComp = true

							str = strings.TrimPrefix(str, "complement(")
							str = strings.TrimSuffix(str, ")")
							items := parseloc(str)

							// reverse items
							for i, j := 0, len(items)-1; i < j; i, j = i+1, j-1 {
								items[i], items[j] = items[j], items[i]
							}

							// reverse from and to positions, flip direction of angle brackets (partial flags)
							for _, thisloc := range items {
								pts := strings.Split(thisloc, "..")
								ln := len(pts)
								if ln == 2 {
									fst := pts[0]
									scd := pts[1]
									lf := ""
									rt := ""
									if strings.HasPrefix(fst, "<") {
										fst = strings.TrimPrefix(fst, "<")
										rt = ">"
									}
									if strings.HasPrefix(scd, ">") {
										scd = strings.TrimPrefix(scd, ">")
										lf = "<"
									}
									acc = append(acc, lf+scd+".."+rt+fst)
								} else if ln > 0 {
									acc = append(acc, pts[0])
								}
							}

						} else {

							// save individual interval or point if no leading accession
							if strings.Index(str, ":") < 0 {
								acc = append(acc, str)
							}
						}

						return acc
					}

					items := parseloc(loc)

					rec.WriteString("        <INSDFeature_intervals>\n")

					numIvals := 0

					// report individual intervals
					for _, thisloc := range items {
						if thisloc == "" {
							continue
						}

						numIvals++

						rec.WriteString("          <INSDInterval>\n")
						pts := strings.Split(thisloc, "..")
						if len(pts) == 2 {

							// fr..to
							fr := pts[0]
							to := pts[1]
							if strings.HasPrefix(fr, "<") {
								fr = strings.TrimPrefix(fr, "<")
								prime5 = true
							}
							if strings.HasPrefix(to, ">") {
								to = strings.TrimPrefix(to, ">")
								prime3 = true
							}
							writeOneElement("            ", "INSDInterval_from", fr)
							writeOneElement("            ", "INSDInterval_to", to)
							if isComp {
								rec.WriteString("            <INSDInterval_iscomp value=\"true\"/>\n")
							}
							writeOneElement("            ", "INSDInterval_accession", accnver)

						} else {

							crt := strings.Split(thisloc, "^")
							if len(crt) == 2 {

								// fr^to
								fr := crt[0]
								to := crt[1]
								writeOneElement("            ", "INSDInterval_from", fr)
								writeOneElement("            ", "INSDInterval_to", to)
								if isComp {
									rec.WriteString("            <INSDInterval_iscomp value=\"true\"/>\n")
								}
								rec.WriteString("            <INSDInterval_interbp value=\"true\"/>\n")
								writeOneElement("            ", "INSDInterval_accession", accnver)

							} else {

								// pt
								pt := pts[0]
								if strings.HasPrefix(pt, "<") {
									pt = strings.TrimPrefix(pt, "<")
									prime5 = true
								}
								if strings.HasPrefix(pt, ">") {
									pt = strings.TrimPrefix(pt, ">")
									prime3 = true
								}
								writeOneElement("            ", "INSDInterval_point", pt)
								writeOneElement("            ", "INSDInterval_accession", accnver)
							}
						}
						rec.WriteString("          </INSDInterval>\n")
					}

					rec.WriteString("        </INSDFeature_intervals>\n")

					if numIvals > 1 {
						writeOneElement("        ", "INSDFeature_operator", locationOperator)
					}
					if prime5 {
						rec.WriteString("        <INSDFeature_partial5 value=\"true\"/>\n")
					}
					if prime3 {
						rec.WriteString("        <INSDFeature_partial3 value=\"true\"/>\n")
					}

					hasQual := false
					for {
						if !strings.HasPrefix(line, twentyonespaces) {
							// if not qualifier line, break out of loop
							break
						}
						txt := strings.TrimPrefix(line, twentyonespaces)
						qual := ""
						val := ""
						if strings.HasPrefix(txt, "/") {
							if !hasQual {
								hasQual = true
								rec.WriteString("        <INSDFeature_quals>\n")
							}
							// read new qualifier and start of value
							qual = strings.TrimPrefix(txt, "/")
							qual = strings.TrimSpace(qual)
							idx := strings.Index(qual, "=")
							if idx > 0 {
								val = qual[idx+1:]
								qual = qual[:idx]
							}

							for {
								line = nextLine()
								row++
								if !strings.HasPrefix(line, twentyonespaces) {
									break
								}
								txt := strings.TrimPrefix(line, twentyonespaces)
								if strings.HasPrefix(txt, "/") {
									// if not continuation of qualifier, break out of loop
									break
								}
								// append subsequent line to value and continue with loop
								if qual == "transcription" || qual == "translation" || qual == "peptide" || qual == "anticodon" {
									val += strings.TrimSpace(txt)
								} else {
									val += " " + strings.TrimSpace(txt)
								}
							}

							rec.WriteString("          <INSDQualifier>\n")

							writeOneElement("            ", "INSDQualifier_name", qual)

							val = strings.TrimPrefix(val, "\"")
							val = strings.TrimSuffix(val, "\"")
							val = strings.TrimSpace(val)
							if val != "" {

								writeOneElement("            ", "INSDQualifier_value", val)
							}

							rec.WriteString("          </INSDQualifier>\n")
						}
					}
					if hasQual {
						rec.WriteString("        </INSDFeature_quals>\n")
					}

					// end of this feature
					rec.WriteString("      </INSDFeature>\n")
					// continue to next feature
				}
			}
			rec.WriteString("    </INSDSeq_feature-table>\n")

			// TSA, TLS, WGS, or CONTIG lines may be next

			altName := ""

			if strings.HasPrefix(line, "TSA") ||
				strings.HasPrefix(line, "TLS") ||
				strings.HasPrefix(line, "WGS") {

				alt.Reset()

				altName = line[:3]
				line = line[3:]
			}

			if strings.HasPrefix(line, "WGS_CONTIG") ||
				strings.HasPrefix(line, "WGS_SCAFLD") {

				alt.Reset()

				altName = line[:3]
				line = line[10:]
			}

			if altName != "" {

				altName = strings.ToLower(altName)
				txt := strings.TrimSpace(line)
				alt.WriteString(txt)
				for {
					// read next line
					line = nextLine()
					row++
					if !strings.HasPrefix(line, twelvespaces) {
						// if not continuation of contig, break out of loop
						break
					}
					// append subsequent line and continue with loop
					txt = strings.TrimPrefix(line, twelvespaces)
					txt = strings.TrimSpace(txt)
					alt.WriteString(txt)
				}
			}

			if strings.HasPrefix(line, "CONTIG") {

				// pathological records can have over 90,000 components, use strings.Builder
				con.Reset()

				txt := strings.TrimPrefix(line, "CONTIG")
				txt = strings.TrimSpace(txt)
				con.WriteString(txt)
				for {
					// read next line
					line = nextLine()
					row++
					if !strings.HasPrefix(line, twelvespaces) {
						// if not continuation of contig, break out of loop
						break
					}
					// append subsequent line and continue with loop
					txt = strings.TrimPrefix(line, twelvespaces)
					txt = strings.TrimSpace(txt)
					con.WriteString(txt)
				}
			}

			if strings.HasPrefix(line, "BASE COUNT") {

				txt := strings.TrimPrefix(line, "BASE COUNT")
				readContinuationLines(txt)
				// not supported
			}

			if strings.HasPrefix(line, "ORIGIN") {

				line = nextLine()
				row++
			}

			// remainder should be sequence

			// sequence can be millions of bases, use strings.Builder
			seq.Reset()

			for line != "" {

				if strings.HasPrefix(line, "//") {

					// end of record, print collected sequence
					str := seq.String()
					if str != "" {

						writeOneElement("    ", "INSDSeq_sequence", str)
					}
					seq.Reset()

					// print contig section
					str = con.String()
					str = strings.TrimSpace(str)
					if str != "" {
						writeOneElement("    ", "INSDSeq_contig", str)
					}
					con.Reset()

					if altName != "" {
						rec.WriteString("    <INSDSeq_alt-seq>\n")
						rec.WriteString("      <INSDAltSeqData>\n")
						str = alt.String()
						str = strings.TrimSpace(str)
						if str != "" {
							writeOneElement("        ", "INSDAltSeqData_name", altName)
							rec.WriteString("        <INSDAltSeqData_items>\n")
							fr, to := SplitInTwoLeft(str, "-")
							if fr != "" && to != "" {
								rec.WriteString("          <INSDAltSeqItem>\n")
								writeOneElement("            ", "INSDAltSeqItem_first-accn", fr)
								writeOneElement("            ", "INSDAltSeqItem_last-accn", to)
								rec.WriteString("          </INSDAltSeqItem>\n")
							} else {
								writeOneElement("          ", "INSDAltSeqItem_value", str)
							}
							rec.WriteString("        </INSDAltSeqData_items>\n")
						}
						alt.Reset()
						rec.WriteString("      </INSDAltSeqData>\n")
						rec.WriteString("    </INSDSeq_alt-seq>\n")
					}

					if dblink != "" {
						rec.WriteString("    <INSDSeq_xrefs>\n")
						// collect for database-reference
						flds := strings.Fields(dblink)
						for len(flds) > 1 {
							tag := flds[0]
							val := flds[1]
							flds = flds[2:]
							tag = strings.TrimSuffix(tag, ":")
							rec.WriteString("      <INSDXref>\n")
							writeOneElement("        ", "INSDXref_dbname", tag)
							writeOneElement("        ", "INSDXref_id", val)
							rec.WriteString("      </INSDXref>\n")
						}
						rec.WriteString("    </INSDSeq_xrefs>\n")
					}

					// end of record
					rec.WriteString("  </INSDSeq>\n")

					// send formatted record down channel
					txt := rec.String()
					out <- txt
					rec.Reset()
					// go to top of loop for next record
					break
				}

				// read next sequence line

				cols := strings.Fields(line)
				if len(cols) > 0 && !IsAllDigits(cols[0]) {
					fmt.Fprintf(os.Stderr, "ERROR: Unrecognized section "+cols[0]+"\n")
				}

				for _, str := range cols {

					if IsAllDigits(str) {
						continue
					}

					// append letters to sequence
					seq.WriteString(str)
				}

				// read next line and continue
				line = nextLine()
				row++

			}

			// continue to next record
		}
	}

	// launch single converter goroutine
	go convertGenBank(inp, out)

	return out
}
