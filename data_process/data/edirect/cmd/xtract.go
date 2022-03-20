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
// File Name:  xtract.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package main

import (
	"bufio"
	"encoding/base64"
	"eutils"
	"fmt"
	"github.com/fatih/color"
	"github.com/surgebase/porter2"
	"html"
	"io"
	"math"
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
	"time"
	"unicode"
)

// GLOBAL VARIABLES

var (
	doStem bool
	deStop bool
)

// TYPED CONSTANTS

// LevelType is the integer type for exploration arguments
type LevelType int

// LevelType keys for exploration arguments
const (
	_ LevelType = iota
	UNIT
	SUBSET
	SECTION
	BLOCK
	BRANCH
	GROUP
	DIVISION
	PATH
	PATTERN
)

// IndentType is the integer type for XML formatting
type IndentType int

// IndentType keys for XML formatting
const (
	SINGULARITY IndentType = iota
	COMPACT
	FLUSH
	INDENT
	SUBTREE
	WRAPPED
)

// OpType is the integer type for operations
type OpType int

// OpType keys for operations
const (
	UNSET OpType = iota
	ELEMENT
	FIRST
	LAST
	BACKWARD
	ENCODE
	DECODE
	PLAIN
	UPPER
	LOWER
	CHAIN
	TITLE
	AUTHOR
	ORDER
	YEAR
	DOI
	TRANSLATE
	REPLACE
	TERMS
	WORDS
	PAIRS
	REVERSE
	LETTERS
	CLAUSES
	INDICES
	ARTICLE
	MESHCODE
	MATRIX
	HISTOGRAM
	ACCENTED
	PFX
	SFX
	SEP
	TAB
	RET
	LBL
	CLR
	PFC
	DEQ
	PLG
	ELG
	FWD
	AWD
	WRP
	ENC
	PKG
	RST
	DEF
	REG
	EXP
	COLOR
	POSITION
	SELECT
	IF
	UNLESS
	MATCH
	AVOID
	AND
	OR
	EQUALS
	CONTAINS
	INCLUDES
	ISWITHIN
	STARTSWITH
	ENDSWITH
	ISNOT
	ISBEFORE
	ISAFTER
	MATCHES
	RESEMBLES
	ISEQUALTO
	DIFFERSFROM
	GT
	GE
	LT
	LE
	EQ
	NE
	NUM
	LEN
	SUM
	MIN
	MAX
	INC
	DEC
	SUB
	AVG
	DEV
	MED
	MUL
	DIV
	MOD
	BIN
	BIT
	ZEROBASED
	ONEBASED
	UCSCBASED
	REVCOMP
	NUCLEIC
	FASTA
	NCBI2NA
	NCBI4NA
	MOLWT
	HGVS
	ELSE
	VARIABLE
	ACCUMULATOR
	VALUE
	QUESTION
	TILDE
	STAR
	DOLLAR
	ATSIGN
	COUNT
	LENGTH
	DEPTH
	INDEX
	UNRECOGNIZED
)

// ArgumentType is the integer type for argument classification
type ArgumentType int

// ArgumentType keys for argument classification
const (
	_ ArgumentType = iota
	EXPLORATION
	CONDITIONAL
	EXTRACTION
	CUSTOMIZATION
)

// RangeType is the integer type for element range choices
type RangeType int

// RangeType keys for element range choices
const (
	NORANGE RangeType = iota
	STRINGRANGE
	VARIABLERANGE
	INTEGERRANGE
)

// SeqEndType is used for -ucsc-based decisions
type SeqEndType int

// SeqEndType keys for -ucsc-based decisions
const (
	_ SeqEndType = iota
	ISSTART
	ISSTOP
	ISPOS
)

// SequenceType is used to record XML tag and position for -ucsc-based
type SequenceType struct {
	Based int
	Which SeqEndType
}

// MUTEXES

var hlock sync.Mutex

var slock sync.RWMutex

// ARGUMENT MAPS

var argTypeIs = map[string]ArgumentType{
	"-unit":         EXPLORATION,
	"-Unit":         EXPLORATION,
	"-subset":       EXPLORATION,
	"-Subset":       EXPLORATION,
	"-section":      EXPLORATION,
	"-Section":      EXPLORATION,
	"-block":        EXPLORATION,
	"-Block":        EXPLORATION,
	"-branch":       EXPLORATION,
	"-Branch":       EXPLORATION,
	"-group":        EXPLORATION,
	"-Group":        EXPLORATION,
	"-division":     EXPLORATION,
	"-Division":     EXPLORATION,
	"-path":         EXPLORATION,
	"-Path":         EXPLORATION,
	"-pattern":      EXPLORATION,
	"-Pattern":      EXPLORATION,
	"-position":     CONDITIONAL,
	"-select":       CONDITIONAL,
	"-if":           CONDITIONAL,
	"-unless":       CONDITIONAL,
	"-match":        CONDITIONAL,
	"-avoid":        CONDITIONAL,
	"-and":          CONDITIONAL,
	"-or":           CONDITIONAL,
	"-equals":       CONDITIONAL,
	"-contains":     CONDITIONAL,
	"-includes":     CONDITIONAL,
	"-is-within":    CONDITIONAL,
	"-starts-with":  CONDITIONAL,
	"-ends-with":    CONDITIONAL,
	"-is-not":       CONDITIONAL,
	"-is-before":    CONDITIONAL,
	"-is-after":     CONDITIONAL,
	"-matches":      CONDITIONAL,
	"-resembles":    CONDITIONAL,
	"-is-equal-to":  CONDITIONAL,
	"-differs-from": CONDITIONAL,
	"-gt":           CONDITIONAL,
	"-ge":           CONDITIONAL,
	"-lt":           CONDITIONAL,
	"-le":           CONDITIONAL,
	"-eq":           CONDITIONAL,
	"-ne":           CONDITIONAL,
	"-element":      EXTRACTION,
	"-first":        EXTRACTION,
	"-last":         EXTRACTION,
	"-backward":     EXTRACTION,
	"-encode":       EXTRACTION,
	"-decode":       EXTRACTION,
	"-decode64":     EXTRACTION,
	"-plain":        EXTRACTION,
	"-upper":        EXTRACTION,
	"-lower":        EXTRACTION,
	"-chain":        EXTRACTION,
	"-title":        EXTRACTION,
	"-author":       EXTRACTION,
	"-order":        EXTRACTION,
	"-year":         EXTRACTION,
	"-doi":          EXTRACTION,
	"-translate":    EXTRACTION,
	"-replace":      EXTRACTION,
	"-terms":        EXTRACTION,
	"-words":        EXTRACTION,
	"-pairs":        EXTRACTION,
	"-reverse":      EXTRACTION,
	"-letters":      EXTRACTION,
	"-clauses":      EXTRACTION,
	"-indices":      EXTRACTION,
	"-article":      EXTRACTION,
	"-meshcode":     EXTRACTION,
	"-matrix":       EXTRACTION,
	"-histogram":    EXTRACTION,
	"-accented":     EXTRACTION,
	"-num":          EXTRACTION,
	"-len":          EXTRACTION,
	"-sum":          EXTRACTION,
	"-min":          EXTRACTION,
	"-max":          EXTRACTION,
	"-inc":          EXTRACTION,
	"-dec":          EXTRACTION,
	"-sub":          EXTRACTION,
	"-avg":          EXTRACTION,
	"-dev":          EXTRACTION,
	"-med":          EXTRACTION,
	"-mul":          EXTRACTION,
	"-div":          EXTRACTION,
	"-mod":          EXTRACTION,
	"-bin":          EXTRACTION,
	"-bit":          EXTRACTION,
	"-0-based":      EXTRACTION,
	"-zero-based":   EXTRACTION,
	"-1-based":      EXTRACTION,
	"-one-based":    EXTRACTION,
	"-ucsc":         EXTRACTION,
	"-ucsc-based":   EXTRACTION,
	"-ucsc-coords":  EXTRACTION,
	"-bed-based":    EXTRACTION,
	"-bed-coords":   EXTRACTION,
	"-revcomp":      EXTRACTION,
	"-nucleic":      EXTRACTION,
	"-fasta":        EXTRACTION,
	"-ncbi2na":      EXTRACTION,
	"-ncbi4na":      EXTRACTION,
	"-molwt":        EXTRACTION,
	"-hgvs":         EXTRACTION,
	"-else":         EXTRACTION,
	"-pfx":          CUSTOMIZATION,
	"-sfx":          CUSTOMIZATION,
	"-sep":          CUSTOMIZATION,
	"-tab":          CUSTOMIZATION,
	"-ret":          CUSTOMIZATION,
	"-lbl":          CUSTOMIZATION,
	"-clr":          CUSTOMIZATION,
	"-pfc":          CUSTOMIZATION,
	"-deq":          CUSTOMIZATION,
	"-plg":          CUSTOMIZATION,
	"-elg":          CUSTOMIZATION,
	"-fwd":          CUSTOMIZATION,
	"-awd":          CUSTOMIZATION,
	"-wrp":          CUSTOMIZATION,
	"-enc":          CUSTOMIZATION,
	"-pkg":          CUSTOMIZATION,
	"-rst":          CUSTOMIZATION,
	"-def":          CUSTOMIZATION,
	"-reg":          CUSTOMIZATION,
	"-exp":          CUSTOMIZATION,
	"-color":        CUSTOMIZATION,
}

var opTypeIs = map[string]OpType{
	"-element":      ELEMENT,
	"-first":        FIRST,
	"-last":         LAST,
	"-backward":     BACKWARD,
	"-encode":       ENCODE,
	"-decode":       DECODE,
	"-decode64":     DECODE,
	"-plain":        PLAIN,
	"-upper":        UPPER,
	"-lower":        LOWER,
	"-chain":        CHAIN,
	"-title":        TITLE,
	"-author":       AUTHOR,
	"-order":        ORDER,
	"-year":         YEAR,
	"-doi":          DOI,
	"-translate":    TRANSLATE,
	"-replace":      REPLACE,
	"-terms":        TERMS,
	"-words":        WORDS,
	"-pairs":        PAIRS,
	"-reverse":      REVERSE,
	"-letters":      LETTERS,
	"-clauses":      CLAUSES,
	"-indices":      INDICES,
	"-article":      ARTICLE,
	"-meshcode":     MESHCODE,
	"-matrix":       MATRIX,
	"-histogram":    HISTOGRAM,
	"-accented":     ACCENTED,
	"-pfx":          PFX,
	"-sfx":          SFX,
	"-sep":          SEP,
	"-tab":          TAB,
	"-ret":          RET,
	"-lbl":          LBL,
	"-clr":          CLR,
	"-pfc":          PFC,
	"-deq":          DEQ,
	"-plg":          PLG,
	"-elg":          ELG,
	"-fwd":          FWD,
	"-awd":          AWD,
	"-wrp":          WRP,
	"-enc":          ENC,
	"-pkg":          PKG,
	"-rst":          RST,
	"-def":          DEF,
	"-reg":          REG,
	"-exp":          EXP,
	"-color":        COLOR,
	"-position":     POSITION,
	"-select":       SELECT,
	"-if":           IF,
	"-unless":       UNLESS,
	"-match":        MATCH,
	"-avoid":        AVOID,
	"-and":          AND,
	"-or":           OR,
	"-equals":       EQUALS,
	"-contains":     CONTAINS,
	"-includes":     INCLUDES,
	"-is-within":    ISWITHIN,
	"-starts-with":  STARTSWITH,
	"-ends-with":    ENDSWITH,
	"-is-not":       ISNOT,
	"-is-before":    ISBEFORE,
	"-is-after":     ISAFTER,
	"-matches":      MATCHES,
	"-resembles":    RESEMBLES,
	"-is-equal-to":  ISEQUALTO,
	"-differs-from": DIFFERSFROM,
	"-gt":           GT,
	"-ge":           GE,
	"-lt":           LT,
	"-le":           LE,
	"-eq":           EQ,
	"-ne":           NE,
	"-num":          NUM,
	"-len":          LEN,
	"-sum":          SUM,
	"-min":          MIN,
	"-max":          MAX,
	"-inc":          INC,
	"-dec":          DEC,
	"-sub":          SUB,
	"-avg":          AVG,
	"-dev":          DEV,
	"-med":          MED,
	"-mul":          MUL,
	"-div":          DIV,
	"-mod":          MOD,
	"-bin":          BIN,
	"-bit":          BIT,
	"-0-based":      ZEROBASED,
	"-zero-based":   ZEROBASED,
	"-1-based":      ONEBASED,
	"-one-based":    ONEBASED,
	"-ucsc":         UCSCBASED,
	"-ucsc-based":   UCSCBASED,
	"-ucsc-coords":  UCSCBASED,
	"-bed-based":    UCSCBASED,
	"-bed-coords":   UCSCBASED,
	"-revcomp":      REVCOMP,
	"-nucleic":      NUCLEIC,
	"-fasta":        FASTA,
	"-ncbi2na":      NCBI2NA,
	"-ncbi4na":      NCBI4NA,
	"-molwt":        MOLWT,
	"-hgvs":         HGVS,
	"-else":         ELSE,
}

var sequenceTypeIs = map[string]SequenceType{
	"INSDSeq:INSDInterval_from":       {1, ISSTART},
	"INSDSeq:INSDInterval_to":         {1, ISSTOP},
	"DocumentSummary:ChrStart":        {0, ISSTART},
	"DocumentSummary:ChrStop":         {0, ISSTOP},
	"DocumentSummary:Chr_start":       {1, ISSTART},
	"DocumentSummary:Chr_end":         {1, ISSTOP},
	"DocumentSummary:Chr_inner_start": {1, ISSTART},
	"DocumentSummary:Chr_inner_end":   {1, ISSTOP},
	"DocumentSummary:Chr_outer_start": {1, ISSTART},
	"DocumentSummary:Chr_outer_end":   {1, ISSTOP},
	"DocumentSummary:start":           {1, ISSTART},
	"DocumentSummary:stop":            {1, ISSTOP},
	"DocumentSummary:display_start":   {1, ISSTART},
	"DocumentSummary:display_stop":    {1, ISSTOP},
	"Entrezgene:Seq-interval_from":    {0, ISSTART},
	"Entrezgene:Seq-interval_to":      {0, ISSTOP},
	"GenomicInfoType:ChrStart":        {0, ISSTART},
	"GenomicInfoType:ChrStop":         {0, ISSTOP},
	"RS:position":                     {0, ISPOS},
	"RS:@asnFrom":                     {0, ISSTART},
	"RS:@asnTo":                       {0, ISSTOP},
	"RS:@end":                         {0, ISSTOP},
	"RS:@leftContigNeighborPos":       {0, ISSTART},
	"RS:@physMapInt":                  {0, ISPOS},
	"RS:@protLoc":                     {0, ISPOS},
	"RS:@rightContigNeighborPos":      {0, ISSTOP},
	"RS:@start":                       {0, ISSTART},
	"RS:@structLoc":                   {0, ISPOS},
}

// DATA OBJECTS

// Step contains parameters for executing a single command step
type Step struct {
	Type   OpType
	Value  string
	Parent string
	Match  string
	Attrib string
	TypL   RangeType
	StrL   string
	IntL   int
	TypR   RangeType
	StrR   string
	IntR   int
	Norm   bool
	Wild   bool
}

// Operation breaks commands into sequential steps
type Operation struct {
	Type   OpType
	Value  string
	Stages []*Step
}

// Block contains nested instructions for executing commands
type Block struct {
	Visit      string
	Parent     string
	Match      string
	Path       []string
	Working    []string
	Parsed     []string
	Position   string
	Foreword   string
	Afterword  string
	Conditions []*Operation
	Commands   []*Operation
	Failure    []*Operation
	Subtasks   []*Block
}

// Limiter is used for collecting specific nodes (e.g., first and last)
type Limiter struct {
	Obj *eutils.XMLNode
	Idx int
	Lvl int
}

// UTILITIES

func hasSpaceOrHyphen(str string) bool {

	for _, ch := range str {
		if ch == ' ' || ch == '-' {
			return true
		}
	}

	return false
}

func isAllCapsOrDigits(str string) bool {

	for _, ch := range str {
		if !unicode.IsUpper(ch) && !unicode.IsDigit(ch) {
			return false
		}
	}

	return true
}

// hasCommaOrSemicolon reports on comma, semicolon, or hyphen
func hasCommaOrSemicolon(str string) bool {

	for _, ch := range str {
		if ch == ',' || ch == ';' || ch == '-' {
			return true
		}
	}

	return false
}

// removeCommaOrSemicolon replaces comma or semicolon with space
func removeCommaOrSemicolon(str string) string {

	str = strings.ToLower(str)

	if hasCommaOrSemicolon(str) {
		str = strings.Replace(str, ",", " ", -1)
		str = strings.Replace(str, ";", " ", -1)
		str = eutils.CompressRunsOfSpaces(str)
	}
	str = strings.TrimSpace(str)
	str = strings.TrimRight(str, ".?:")

	return str
}

// sortStringByWords sorts the individual words in a string
func sortStringByWords(str string) string {

	str = removeCommaOrSemicolon(str)

	// check for multiple words
	if hasSpaceOrHyphen(str) {
		flds := strings.Fields(str)
		sort.Slice(flds, func(i, j int) bool { return flds[i] < flds[j] })
		str = strings.Join(flds, " ")
		str = strings.Replace(str, "-", " ", -1)
		str = eutils.CompressRunsOfSpaces(str)
		str = strings.TrimRight(str, ".?:")
	}

	return str
}

func parseFlag(str string) (OpType, bool) {

	op, ok := opTypeIs[str]
	if ok {
		if argTypeIs[str] == EXTRACTION {
			return op, true
		}
		return op, false
	}

	if len(str) > 1 && str[0] == '-' && isAllCapsOrDigits(str[1:]) {
		return VARIABLE, true
	}

	if len(str) > 2 && strings.HasPrefix(str, "--") && isAllCapsOrDigits(str[2:]) {
		return ACCUMULATOR, true
	}

	if len(str) > 0 && str[0] == '-' {
		return UNRECOGNIZED, false
	}

	return UNSET, false
}

func parseMarkup(str, cmd string) int {

	switch str {
	case "fuse", "fused":
		return eutils.FUSE
	case "space", "spaces":
		return eutils.SPACE
	case "period", "periods":
		return eutils.PERIOD
	case "bracket", "brackets":
		return eutils.BRACKETS
	case "markdown":
		return eutils.MARKDOWN
	case "slash":
		return eutils.SLASH
	case "tag", "tags":
		return eutils.TAGS
	case "terse":
		return eutils.TERSE
	default:
		if str != "" {
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized %s value '%s'\n", cmd, str)
			os.Exit(1)
		}
	}
	return eutils.NOMARKUP
}

// DebugBlock examines structure of parsed arguments (undocumented)
/*
func DebugBlock(blk *Block, depth int) {

	doIndent := func(indt int) {
		for i := 1; i < indt; i++ {
			fmt.Fprintf(os.Stderr, "  ")
		}
	}

	doIndent(depth)

	if blk.Visit != "" {
		doIndent(depth + 1)
		fmt.Fprintf(os.Stderr, "<Visit> %s </Visit>\n", blk.Visit)
	}
	if len(blk.Parsed) > 0 {
		doIndent(depth + 1)
		fmt.Fprintf(os.Stderr, "<Parsed>")
		for _, str := range blk.Parsed {
			fmt.Fprintf(os.Stderr, " %s", str)
		}
		fmt.Fprintf(os.Stderr, " </Parsed>\n")
	}

	if len(blk.Subtasks) > 0 {
		for _, sub := range blk.Subtasks {
			DebugBlock(sub, depth+1)
		}
	}
}
*/

// PARSE COMMAND-LINE ARGUMENTS

// parseArguments parses nested exploration instruction from command-line arguments
func parseArguments(cmdargs []string, pttrn string) *Block {

	// different names of exploration control arguments allow multiple levels of nested "for" loops in a linear command line
	// (capitalized versions for backward-compatibility with original Perl implementation handling of recursive definitions)
	var (
		lcname = []string{
			"",
			"-unit",
			"-subset",
			"-section",
			"-block",
			"-branch",
			"-group",
			"-division",
			"-path",
			"-pattern",
		}

		ucname = []string{
			"",
			"-Unit",
			"-Subset",
			"-Section",
			"-Block",
			"-Branch",
			"-Group",
			"-Division",
			"-Path",
			"-Pattern",
		}
	)

	// parseCommands recursive definition
	var parseCommands func(parent *Block, startLevel LevelType)

	// parseCommands does initial parsing of exploration command structure
	parseCommands = func(parent *Block, startLevel LevelType) {

		// find next highest level exploration argument
		findNextLevel := func(args []string, level LevelType) (LevelType, string, string) {

			if len(args) > 1 {

				for {

					if level < UNIT {
						break
					}

					lctag := lcname[level]
					uctag := ucname[level]

					for _, txt := range args {
						if txt == lctag || txt == uctag {
							return level, lctag, uctag
						}
					}

					level--
				}
			}

			return 0, "", ""
		}

		arguments := parent.Working

		level, lctag, uctag := findNextLevel(arguments, startLevel)

		if level < UNIT {

			// break recursion
			return
		}

		// group arguments at a given exploration level
		subsetCommands := func(args []string) *Block {

			max := len(args)

			visit := ""

			// extract name of object to visit
			if max > 1 {
				visit = args[1]
				args = args[2:]
				max -= 2
			}

			partition := 0
			for cur, str := range args {

				// record point of next exploration command
				partition = cur + 1

				// skip if not a command
				if len(str) < 1 || str[0] != '-' {
					continue
				}

				if argTypeIs[str] == EXPLORATION {
					partition = cur
					break
				}
			}

			// convert slashes (e.g., parent/child construct) to periods (e.g., dotted exploration path)
			if strings.Contains(visit, "/") {
				if !strings.Contains(visit, ".") {
					visit = strings.Replace(visit, "/", ".", -1)
				}
			}

			// parse parent.child or dotted path construct
			// colon indicates a namespace prefix in any or all of the components
			prnt, rmdr := eutils.SplitInTwoRight(visit, ".")
			match, rest := eutils.SplitInTwoLeft(rmdr, ".")

			if rest != "" {

				// exploration match on first component, then search remainder one level at a time with subsequent components
				dirs := strings.Split(rmdr, ".")

				// signal with "path" position
				return &Block{Visit: visit, Parent: "", Match: prnt, Path: dirs, Position: "path", Parsed: args[0:partition], Working: args[partition:]}
			}

			// promote arguments parsed at this level
			return &Block{Visit: visit, Parent: prnt, Match: match, Parsed: args[0:partition], Working: args[partition:]}
		}

		cur := 0

		// search for positions of current exploration command

		for idx, txt := range arguments {
			if txt == lctag || txt == uctag {
				if idx == 0 {
					continue
				}

				blk := subsetCommands(arguments[cur:idx])
				parseCommands(blk, level-1)
				parent.Subtasks = append(parent.Subtasks, blk)

				cur = idx
			}
		}

		if cur < len(arguments) {
			blk := subsetCommands(arguments[cur:])
			parseCommands(blk, level-1)
			parent.Subtasks = append(parent.Subtasks, blk)
		}

		// clear execution arguments from parent after subsetting
		parent.Working = nil
	}

	// parse optional [min:max], [&VAR:&VAR], or [after|before] range specification
	parseRange := func(item, rnge string) (typL RangeType, strL string, intL int, typR RangeType, strR string, intR int) {

		typL = NORANGE
		typR = NORANGE
		strL = ""
		strR = ""
		intL = 0
		intR = 0

		if rnge == "" {
			// no range specification, return default values
			return
		}

		// check if last character is right square bracket
		if !strings.HasSuffix(rnge, "]") {
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized range %s\n", rnge)
			os.Exit(1)
		}

		rnge = strings.TrimSuffix(rnge, "]")

		if rnge == "" {
			fmt.Fprintf(os.Stderr, "\nERROR: Empty range %s[]\n", item)
			os.Exit(1)
		}

		// check for [after|before] variant
		if strings.Contains(rnge, "|") {

			strL, strR = eutils.SplitInTwoLeft(rnge, "|")
			// spacing matters, so do not call TrimSpace

			if strL == "" && strR == "" {
				fmt.Fprintf(os.Stderr, "\nERROR: Empty range %s[|]\n", item)
				os.Exit(1)
			}

			typL = STRINGRANGE
			typR = STRINGRANGE

			// return statement returns named variables
			return
		}

		// otherwise must have colon within brackets
		if !strings.Contains(rnge, ":") {
			fmt.Fprintf(os.Stderr, "\nERROR: Colon missing in range %s[%s]\n", item, rnge)
			os.Exit(1)
		}

		// split at colon
		lft, rgt := eutils.SplitInTwoLeft(rnge, ":")

		lft = strings.TrimSpace(lft)
		rgt = strings.TrimSpace(rgt)

		if lft == "" && rgt == "" {
			fmt.Fprintf(os.Stderr, "\nERROR: Empty range %s[:]\n", item)
			os.Exit(1)
		}

		// for variable, parse optional +/- offset suffix
		parseOffset := func(str string) (string, int) {

			if str == "" || str[0] == ' ' {
				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized variable '&%s'\n", str)
				os.Exit(1)
			}

			pls := ""
			mns := ""

			ofs := 0

			// check for &VAR+1 or &VAR-1 integer adjustment
			str, pls = eutils.SplitInTwoLeft(str, "+")
			str, mns = eutils.SplitInTwoLeft(str, "-")

			if pls != "" {
				val, err := strconv.Atoi(pls)
				if err != nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized range adjustment &%s+%s\n", str, pls)
					os.Exit(1)
				}
				ofs = val
			} else if mns != "" {
				val, err := strconv.Atoi(mns)
				if err != nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized range adjustment &%s-%s\n", str, mns)
					os.Exit(1)
				}
				ofs = -val
			}

			return str, ofs
		}

		// parse integer position, 1-based coordinate must be greater than 0
		parseInteger := func(str string, mustBePositive bool) int {
			if str == "" {
				return 0
			}

			val, err := strconv.Atoi(str)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized range component %s[%s:]\n", item, str)
				os.Exit(1)
			}
			if mustBePositive {
				if val < 1 {
					fmt.Fprintf(os.Stderr, "\nERROR: Range component %s[%s:] must be positive\n", item, str)
					os.Exit(1)
				}
			} else {
				if val == 0 {
					fmt.Fprintf(os.Stderr, "\nERROR: Range component %s[%s:] must not be zero\n", item, str)
					os.Exit(1)
				}
			}

			return val
		}

		if lft != "" {
			if lft[0] == '&' {
				lft = lft[1:]
				strL, intL = parseOffset(lft)
				typL = VARIABLERANGE
			} else {
				intL = parseInteger(lft, true)
				typL = INTEGERRANGE
			}
		}

		if rgt != "" {
			if rgt[0] == '&' {
				rgt = rgt[1:]
				strR, intR = parseOffset(rgt)
				typR = VARIABLERANGE
			} else {
				intR = parseInteger(rgt, false)
				typR = INTEGERRANGE
			}
		}

		// return statement required to return named variables
		return
	}

	parseConditionals := func(cmds *Block, arguments []string) []*Operation {

		max := len(arguments)
		if max < 1 {
			return nil
		}

		// check for missing condition command
		txt := arguments[0]
		if txt != "-if" && txt != "-unless" && txt != "-select" && txt != "-match" && txt != "-avoid" && txt != "-position" {
			fmt.Fprintf(os.Stderr, "\nERROR: Missing -if command before '%s'\n", txt)
			os.Exit(1)
		}
		if txt == "-position" && max > 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: Cannot combine -position with -if or -unless commands\n")
			os.Exit(1)
		}
		// check for missing argument after last condition
		txt = arguments[max-1]
		if len(txt) > 0 && txt[0] == '-' {
			fmt.Fprintf(os.Stderr, "\nERROR: Item missing after %s command\n", txt)
			os.Exit(1)
		}

		cond := make([]*Operation, 0, max)

		// parse conditional clause into execution step
		parseStep := func(op *Operation, elementColonValue bool) {

			if op == nil {
				return
			}

			str := op.Value

			status := ELEMENT

			// isolate and parse optional [min:max], [&VAR:&VAR], or [after|before] range specification
			str, rnge := eutils.SplitInTwoLeft(str, "[")

			str = strings.TrimSpace(str)
			rnge = strings.TrimSpace(rnge)

			if str == "" && rnge != "" {
				fmt.Fprintf(os.Stderr, "\nERROR: Variable missing in range specification [%s\n", rnge)
				os.Exit(1)
			}

			typL, strL, intL, typR, strR, intR := parseRange(str, rnge)

			// check for pound, percent, or caret character at beginning of name
			if len(str) > 1 {
				switch str[0] {
				case '&':
					if isAllCapsOrDigits(str[1:]) {
						status = VARIABLE
						str = str[1:]
					} else if strings.Contains(str, ":") {
						fmt.Fprintf(os.Stderr, "\nERROR: Unsupported construct '%s', use -if &VARIABLE -equals VALUE instead\n", str)
						os.Exit(1)
					} else {
						fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized variable '%s'\n", str)
						os.Exit(1)
					}
				case '#':
					status = COUNT
					str = str[1:]
				case '%':
					status = LENGTH
					str = str[1:]
				case '^':
					status = DEPTH
					str = str[1:]
				default:
				}
			} else if str == "+" {
				status = INDEX
			}

			// parse parent/element@attribute construct
			// colon indicates a namespace prefix in any or all of the components
			prnt, match := eutils.SplitInTwoRight(str, "/")
			match, attrib := eutils.SplitInTwoLeft(match, "@")
			val := ""

			// leading colon indicates namespace prefix wildcard
			wildcard := false
			if strings.HasPrefix(prnt, ":") || strings.HasPrefix(match, ":") || strings.HasPrefix(attrib, ":") {
				wildcard = true
			}

			if elementColonValue {

				// allow parent/element@attribute:value construct for deprecated -match and -avoid, and for subsequent -and and -or commands
				match, val = eutils.SplitInTwoLeft(str, ":")
				prnt, match = eutils.SplitInTwoRight(match, "/")
				match, attrib = eutils.SplitInTwoLeft(match, "@")
			}

			norm := true
			if rnge != "" {
				if typL != NORANGE || typR != NORANGE || strL != "" || strR != "" || intL != 0 || intR != 0 {
					norm = false
				}
			}

			tsk := &Step{Type: status, Value: str, Parent: prnt, Match: match, Attrib: attrib,
				TypL: typL, StrL: strL, IntL: intL, TypR: typR, StrR: strR, IntR: intR,
				Norm: norm, Wild: wildcard}

			op.Stages = append(op.Stages, tsk)

			// transform old -match "element:value" to -match element -equals value
			if val != "" {
				tsk := &Step{Type: EQUALS, Value: val}
				op.Stages = append(op.Stages, tsk)
			}
		}

		idx := 0

		// conditionals should alternate between command and object/value
		expectDash := true
		last := ""

		var op *Operation

		// flag to allow element-colon-value for deprecated -match and -avoid commands, otherwise colon is for namespace prefixes
		elementColonValue := false

		status := UNSET

		// parse command strings into operation structure
		for idx < max {
			str := arguments[idx]
			idx++

			// conditionals should alternate between command and object/value
			if expectDash {
				if len(str) < 1 || str[0] != '-' {
					fmt.Fprintf(os.Stderr, "\nERROR: Unexpected '%s' argument after '%s'\n", str, last)
					os.Exit(1)
				}
				expectDash = false
			} else {
				if len(str) > 0 && str[0] == '-' {
					fmt.Fprintf(os.Stderr, "\nERROR: Unexpected '%s' command after '%s'\n", str, last)
					os.Exit(1)
				}
				expectDash = true
			}
			last = str

			switch status {
			case UNSET:
				status, _ = parseFlag(str)
			case POSITION:
				if cmds.Position != "" {
					fmt.Fprintf(os.Stderr, "\nERROR: -position '%s' conflicts with existing '%s'\n", str, cmds.Position)
					os.Exit(1)
				}
				cmds.Position = str
				status = UNSET
			case MATCH, AVOID:
				elementColonValue = true
				fallthrough
			case SELECT, IF, UNLESS, AND, OR:
				op = &Operation{Type: status, Value: str}
				cond = append(cond, op)
				parseStep(op, elementColonValue)
				status = UNSET
			case EQUALS, CONTAINS, INCLUDES, ISWITHIN, STARTSWITH, ENDSWITH, ISNOT, ISBEFORE, ISAFTER:
				if op != nil {
					if len(str) > 1 && str[0] == '\\' {
						// first character may be backslash protecting dash (undocumented)
						str = str[1:]
					}
					tsk := &Step{Type: status, Value: str}
					op.Stages = append(op.Stages, tsk)
					op = nil
				} else {
					fmt.Fprintf(os.Stderr, "\nERROR: Unexpected adjacent string match constraints\n")
					os.Exit(1)
				}
				status = UNSET
			case MATCHES:
				if op != nil {
					if len(str) > 1 && str[0] == '\\' {
						// first character may be backslash protecting dash (undocumented)
						str = str[1:]
					}
					str = removeCommaOrSemicolon(str)
					tsk := &Step{Type: status, Value: str}
					op.Stages = append(op.Stages, tsk)
					op = nil
				} else {
					fmt.Fprintf(os.Stderr, "\nERROR: Unexpected adjacent string match constraints\n")
					os.Exit(1)
				}
				status = UNSET
			case RESEMBLES:
				if op != nil {
					if len(str) > 1 && str[0] == '\\' {
						// first character may be backslash protecting dash (undocumented)
						str = str[1:]
					}
					str = sortStringByWords(str)
					tsk := &Step{Type: status, Value: str}
					op.Stages = append(op.Stages, tsk)
					op = nil
				} else {
					fmt.Fprintf(os.Stderr, "\nERROR: Unexpected adjacent string match constraints\n")
					os.Exit(1)
				}
				status = UNSET
			case ISEQUALTO, DIFFERSFROM:
				if op != nil {
					if len(str) < 1 {
						fmt.Fprintf(os.Stderr, "\nERROR: Empty conditional argument\n")
						os.Exit(1)
					}
					ch := str[0]
					// uses element as second argument
					orig := str
					if ch == '#' || ch == '%' || ch == '^' {
						// check for pound, percent, or caret character at beginning of element (undocumented)
						str = str[1:]
						if len(str) < 1 {
							fmt.Fprintf(os.Stderr, "\nERROR: Unexpected conditional constraints\n")
							os.Exit(1)
						}
						ch = str[0]
					}
					if (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') {
						prnt, match := eutils.SplitInTwoRight(str, "/")
						match, attrib := eutils.SplitInTwoLeft(match, "@")
						wildcard := false
						if strings.HasPrefix(prnt, ":") || strings.HasPrefix(match, ":") || strings.HasPrefix(attrib, ":") {
							wildcard = true
						}
						tsk := &Step{Type: status, Value: orig, Parent: prnt, Match: match, Attrib: attrib, Wild: wildcard}
						op.Stages = append(op.Stages, tsk)
					} else {
						fmt.Fprintf(os.Stderr, "\nERROR: Unexpected conditional constraints\n")
						os.Exit(1)
					}
					op = nil
				}
				status = UNSET
			case GT, GE, LT, LE, EQ, NE:
				if op != nil {
					if len(str) > 1 && str[0] == '\\' {
						// first character may be backslash protecting minus sign (undocumented)
						str = str[1:]
					}
					if len(str) < 1 {
						fmt.Fprintf(os.Stderr, "\nERROR: Empty numeric match constraints\n")
						os.Exit(1)
					}
					ch := str[0]
					if (ch >= '0' && ch <= '9') || ch == '-' || ch == '+' {
						// literal numeric constant
						tsk := &Step{Type: status, Value: str}
						op.Stages = append(op.Stages, tsk)
					} else {
						// numeric test allows element as second argument
						orig := str
						if ch == '#' || ch == '%' || ch == '^' {
							// check for pound, percent, or caret character at beginning of element (undocumented)
							str = str[1:]
							if len(str) < 1 {
								fmt.Fprintf(os.Stderr, "\nERROR: Unexpected numeric match constraints\n")
								os.Exit(1)
							}
							ch = str[0]
						}
						if (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') {
							prnt, match := eutils.SplitInTwoRight(str, "/")
							match, attrib := eutils.SplitInTwoLeft(match, "@")
							wildcard := false
							if strings.HasPrefix(prnt, ":") || strings.HasPrefix(match, ":") || strings.HasPrefix(attrib, ":") {
								wildcard = true
							}
							tsk := &Step{Type: status, Value: orig, Parent: prnt, Match: match, Attrib: attrib, Wild: wildcard}
							op.Stages = append(op.Stages, tsk)
						} else {
							fmt.Fprintf(os.Stderr, "\nERROR: Unexpected numeric match constraints\n")
							os.Exit(1)
						}
					}
					op = nil
				} else {
					fmt.Fprintf(os.Stderr, "\nERROR: Unexpected adjacent numeric match constraints\n")
					os.Exit(1)
				}
				status = UNSET
			case UNRECOGNIZED:
				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized argument '%s'\n", str)
				os.Exit(1)
			default:
				fmt.Fprintf(os.Stderr, "\nERROR: Unexpected argument '%s'\n", str)
				os.Exit(1)
			}
		}

		return cond
	}

	parseExtractions := func(cmds *Block, arguments []string) []*Operation {

		max := len(arguments)
		if max < 1 {
			return nil
		}

		// check for missing -element (or -first, etc.) command
		txt := arguments[0]
		if len(txt) < 1 || txt[0] != '-' {
			fmt.Fprintf(os.Stderr, "\nERROR: Missing -element command before '%s'\n", txt)
			os.Exit(1)
		}
		// check for missing argument after last -element (or -first, etc.) command
		txt = arguments[max-1]
		if len(txt) > 0 && txt[0] == '-' {
			if txt == "-rst" {
				fmt.Fprintf(os.Stderr, "\nERROR: Unexpected position for %s command\n", txt)
				os.Exit(1)
			} else if txt == "-clr" {
				// main loop runs out after trailing -clr, add another so this one will be executed
				arguments = append(arguments, "-clr")
				max++
			} else if max < 2 || arguments[max-2] != "-lbl" {
				fmt.Fprintf(os.Stderr, "\nERROR: Item missing after %s command\n", txt)
				os.Exit(1)
			}
		}

		comm := make([]*Operation, 0, max)

		// parse next argument
		nextStatus := func(str string) (OpType, bool) {

			status, isExtraction := parseFlag(str)

			switch status {
			case VARIABLE:
				op := &Operation{Type: status, Value: str[1:]}
				comm = append(comm, op)
				status = VALUE
			case ACCUMULATOR:
				op := &Operation{Type: status, Value: str[2:]}
				comm = append(comm, op)
				status = VALUE
			case CLR, RST:
				op := &Operation{Type: status, Value: ""}
				comm = append(comm, op)
				status = UNSET
			case ELEMENT:
			case TAB, RET, PFX, SFX, SEP, LBL, PFC, DEQ, PLG, ELG, WRP, ENC, DEF, REG, EXP, COLOR:
			case FWD, AWD, PKG:
			case UNSET:
				fmt.Fprintf(os.Stderr, "\nERROR: No -element before '%s'\n", str)
				os.Exit(1)
			case UNRECOGNIZED:
				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized argument '%s'\n", str)
				os.Exit(1)
			default:
				if !isExtraction {
					// not ELEMENT through HGVS
					fmt.Fprintf(os.Stderr, "\nERROR: Misplaced %s command\n", str)
					os.Exit(1)
				}
			}

			return status, isExtraction
		}

		// parse extraction clause into individual steps
		parseSteps := func(op *Operation, pttrn string) {

			if op == nil {
				return
			}

			stat := op.Type
			str := op.Value

			// element names combined with commas are treated as a prefix-separator-suffix group
			comma := strings.Split(str, ",")

			rnge := ""
			for _, item := range comma {
				status := stat

				// isolate and parse optional [min:max], [&VAR:&VAR], or [after|before] range specification
				item, rnge = eutils.SplitInTwoLeft(item, "[")

				item = strings.TrimSpace(item)
				rnge = strings.TrimSpace(rnge)

				if item == "" && rnge != "" {
					fmt.Fprintf(os.Stderr, "\nERROR: Variable missing in range specification [%s\n", rnge)
					os.Exit(1)
				}

				typL, strL, intL, typR, strR, intR := parseRange(item, rnge)

				// check for special character at beginning of name
				if len(item) > 1 {
					switch item[0] {
					case '&':
						if isAllCapsOrDigits(item[1:]) {
							status = VARIABLE
							item = item[1:]
						} else {
							fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized variable '%s'\n", item)
							os.Exit(1)
						}
					case '#':
						status = COUNT
						item = item[1:]
					case '%':
						status = LENGTH
						item = item[1:]
					case '^':
						status = DEPTH
						item = item[1:]
					case '*':
						for _, ch := range item {
							if ch != '*' {
								break
							}
						}
						status = STAR
					default:
					}
				} else {
					switch item {
					case "?":
						status = QUESTION
					case "~":
						status = TILDE
					case "*":
						status = STAR
					case "$":
						status = DOLLAR
					case "@":
						status = ATSIGN
					case "+":
						status = INDEX
					default:
					}
				}

				// parse parent/element@attribute construct
				// colon indicates a namespace prefix in any or all of the components
				prnt, match := eutils.SplitInTwoRight(item, "/")
				match, attrib := eutils.SplitInTwoLeft(match, "@")

				// leading colon indicates namespace prefix wildcard
				wildcard := false
				if strings.HasPrefix(prnt, ":") || strings.HasPrefix(match, ":") || strings.HasPrefix(attrib, ":") {
					wildcard = true
				}

				// sequence coordinate adjustments
				switch status {
				case ZEROBASED, ONEBASED, UCSCBASED:
					seq := pttrn + ":"
					if attrib != "" {
						seq += "@"
						seq += attrib
					} else if match != "" {
						seq += match
					}
					// confirm -0-based or -1-based arguments are known sequence position elements or attributes
					slock.RLock()
					seqtype, ok := sequenceTypeIs[seq]
					slock.RUnlock()
					if !ok {
						fmt.Fprintf(os.Stderr, "\nERROR: Element '%s' is not suitable for sequence coordinate conversion\n", item)
						os.Exit(1)
					}
					switch status {
					case ZEROBASED:
						status = ELEMENT
						// if 1-based coordinates, decrement to get 0-based value
						if seqtype.Based == 1 {
							status = DEC
						}
					case ONEBASED:
						status = ELEMENT
						// if 0-based coordinates, increment to get 1-based value
						if seqtype.Based == 0 {
							status = INC
						}
					case UCSCBASED:
						status = ELEMENT
						// half-open intervals, start is 0-based, stop is 1-based
						if seqtype.Based == 0 && seqtype.Which == ISSTOP {
							status = INC
						} else if seqtype.Based == 1 && seqtype.Which == ISSTART {
							status = DEC
						}
					default:
						status = ELEMENT
					}
				default:
				}

				norm := true
				if rnge != "" {
					if typL != NORANGE || typR != NORANGE || strL != "" || strR != "" || intL != 0 || intR != 0 {
						norm = false
					}
				}

				tsk := &Step{Type: status, Value: item, Parent: prnt, Match: match, Attrib: attrib,
					TypL: typL, StrL: strL, IntL: intL, TypR: typR, StrR: strR, IntR: intR,
					Norm: norm, Wild: wildcard}

				op.Stages = append(op.Stages, tsk)
			}
		}

		idx := 0

		status := UNSET
		isExtraction := false

		// parse command strings into operation structure
		for idx < max {
			str := arguments[idx]
			idx++

			if argTypeIs[str] == CONDITIONAL {
				fmt.Fprintf(os.Stderr, "\nERROR: Misplaced %s command\n", str)
				os.Exit(1)
			}

			switch status {
			case UNSET:
				status, isExtraction = nextStatus(str)
			case TAB, RET, PFX, SFX, SEP, LBL, PFC, DEQ, PLG, ELG, WRP, ENC, DEF, REG, EXP, COLOR:
				op := &Operation{Type: status, Value: eutils.ConvertSlash(str)}
				comm = append(comm, op)
				status = UNSET
			case FWD:
				cmds.Foreword = eutils.ConvertSlash(str)
				status = UNSET
			case AWD:
				cmds.Afterword = eutils.ConvertSlash(str)
				status = UNSET
			case PKG:
				pkg := eutils.ConvertSlash(str)
				cmds.Foreword = ""
				cmds.Afterword = ""
				if pkg != "" && pkg != "-" {
					items := strings.Split(pkg, "/")
					for i := 0; i < len(items); i++ {
						cmds.Foreword += "<" + items[i] + ">"
					}
					for i := len(items) - 1; i >= 0; i-- {
						cmds.Afterword += "</" + items[i] + ">"
					}
				}
				status = UNSET
			case VARIABLE:
				op := &Operation{Type: status, Value: str[1:]}
				comm = append(comm, op)
				status = VALUE
			case ACCUMULATOR:
				op := &Operation{Type: status, Value: str[2:]}
				comm = append(comm, op)
				status = VALUE
			case VALUE:
				op := &Operation{Type: status, Value: str}
				comm = append(comm, op)
				parseSteps(op, pttrn)
				status = UNSET
			case UNRECOGNIZED:
				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized argument '%s'\n", str)
				os.Exit(1)
			default:
				if isExtraction {
					// ELEMENT through HGVS
					for !strings.HasPrefix(str, "-") {
						// create one operation per argument, even if under a single -element statement
						op := &Operation{Type: status, Value: str}
						comm = append(comm, op)
						parseSteps(op, pttrn)
						if idx >= max {
							break
						}
						str = arguments[idx]
						idx++
					}
					status = UNSET
					if idx < max {
						status, isExtraction = nextStatus(str)
					}
				}
			}
		}

		return comm
	}

	// parseOperations recursive definition
	var parseOperations func(parent *Block)

	// parseOperations converts parsed arguments to operations lists
	parseOperations = func(parent *Block) {

		args := parent.Parsed

		partition := 0
		for cur, str := range args {

			// record junction between conditional and extraction commands
			partition = cur + 1

			// skip if not a command
			if len(str) < 1 || str[0] != '-' {
				continue
			}

			if argTypeIs[str] != CONDITIONAL {
				partition = cur
				break
			}
		}

		// split arguments into conditional tests and extraction or customization commands
		conditionals := args[0:partition]
		args = args[partition:]

		partition = 0
		foundElse := false
		for cur, str := range args {

			// record junction at -else command
			partition = cur + 1

			// skip if not a command
			if len(str) < 1 || str[0] != '-' {
				continue
			}

			if str == "-else" {
				partition = cur
				foundElse = true
				break
			}
		}

		extractions := args[0:partition]
		alternative := args[partition:]

		if len(alternative) > 0 && alternative[0] == "-else" {
			alternative = alternative[1:]
		}

		// validate argument structure and convert to operations lists
		parent.Conditions = parseConditionals(parent, conditionals)
		parent.Commands = parseExtractions(parent, extractions)
		parent.Failure = parseExtractions(parent, alternative)

		// reality checks on placement of -else command
		if foundElse {
			if len(conditionals) < 1 {
				fmt.Fprintf(os.Stderr, "\nERROR: Misplaced -else command\n")
				os.Exit(1)
			}
			if len(alternative) < 1 {
				fmt.Fprintf(os.Stderr, "\nERROR: Misplaced -else command\n")
				os.Exit(1)
			}
			if len(parent.Subtasks) > 0 {
				fmt.Fprintf(os.Stderr, "\nERROR: Misplaced -else command\n")
				os.Exit(1)
			}
		}

		for _, sub := range parent.Subtasks {
			parseOperations(sub)
		}
	}

	// parseArguments

	head := &Block{}

	for _, txt := range cmdargs {
		head.Working = append(head.Working, txt)
	}

	// initial parsing of exploration command structure
	parseCommands(head, PATTERN)

	if len(head.Subtasks) != 1 {
		return nil
	}

	// skip past empty placeholder
	head = head.Subtasks[0]

	// convert command strings to array of operations for faster processing
	parseOperations(head)

	// check for no -element or multiple -pattern commands
	noElement := true
	numPatterns := 0
	for _, txt := range cmdargs {
		if argTypeIs[txt] == EXTRACTION {
			noElement = false
		}
		if txt == "-pattern" || txt == "-Pattern" {
			numPatterns++
		} else if txt == "-select" {
			noElement = false
			head.Position = "select"
		}
	}

	if numPatterns < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: No -pattern in command-line arguments\n")
		os.Exit(1)
	}

	if numPatterns > 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Only one -pattern command is permitted\n")
		os.Exit(1)
	}

	if noElement {
		fmt.Fprintf(os.Stderr, "\nERROR: No -element statement in argument list\n")
		os.Exit(1)
	}

	return head
}

// printSubtree supports compression styles selected by -element "*" through "****"
func printSubtree(node *eutils.XMLNode, style IndentType, printAttrs bool, proc func(string)) {

	if node == nil || proc == nil {
		return
	}

	// WRAPPED is SUBTREE plus each attribute on its own line
	wrapped := false
	if style == WRAPPED {
		style = SUBTREE
		wrapped = true
	}

	// INDENT is offset by two spaces to allow for parent tag, SUBTREE is not offset
	initial := 1
	if style == SUBTREE {
		style = INDENT
		initial = 0
	}

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

	// indent a specified number of spaces
	doIndent := func(indt int) {
		i := indt
		for i > 9 {
			proc("                    ")
			i -= 10
		}
		if i < 0 {
			return
		}
		proc(indentSpaces[i])
	}

	// doSubtree recursive definition
	var doSubtree func(*eutils.XMLNode, int)

	doSubtree = func(curr *eutils.XMLNode, depth int) {

		// suppress if it would be an empty self-closing tag
		if !eutils.IsNotJustWhitespace(curr.Attributes) && curr.Contents == "" && curr.Children == nil {
			return
		}

		if style == INDENT {
			doIndent(depth)
		}

		if curr.Name != "" {
			proc("<")
			proc(curr.Name)

			if printAttrs {

				attr := strings.TrimSpace(curr.Attributes)
				attr = eutils.CompressRunsOfSpaces(attr)

				if attr != "" {

					if wrapped {

						start := 0
						idx := 0

						attlen := len(attr)

						for idx < attlen {
							ch := attr[idx]
							if ch == '=' {
								str := attr[start:idx]
								proc("\n")
								doIndent(depth)
								proc(" ")
								proc(str)
								// skip past equal sign and leading double quote
								idx += 2
								start = idx
							} else if ch == '"' || ch == '\'' {
								str := attr[start:idx]
								proc("=\"")
								proc(str)
								proc("\"")
								// skip past trailing double quote and (possible) space
								idx += 2
								start = idx
							} else {
								idx++
							}
						}

						proc("\n")
						doIndent(depth)

					} else {

						proc(" ")
						proc(attr)
					}
				}
			}

			// see if suitable for for self-closing tag
			if curr.Contents == "" && curr.Children == nil {
				proc("/>")
				if style != COMPACT {
					proc("\n")
				}
				return
			}

			proc(">")
		}

		if curr.Contents != "" {

			proc(curr.Contents[:])

		} else {

			if style != COMPACT {
				proc("\n")
			}

			for chld := curr.Children; chld != nil; chld = chld.Next {
				doSubtree(chld, depth+1)
			}

			if style == INDENT {
				i := depth
				for i > 9 {
					proc("                    ")
					i -= 10
				}
				proc(indentSpaces[i])
			}
		}

		if curr.Name != "" {
			proc("<")
			proc("/")
			proc(curr.Name)
			proc(">")
		}

		if style != COMPACT {
			proc("\n")
		}
	}

	doSubtree(node, initial)
}

var (
	xlock sync.Mutex
	replx map[string]*regexp.Regexp
)

// processClause handles comma-separated -element arguments
func processClause(curr *eutils.XMLNode, stages []*Step, mask, prev, pfx, sfx, plg, sep, def, reg, exp string, wrp bool, status OpType, index, level int, variables map[string]string, transform map[string]string, histogram map[string]int) (string, bool) {

	if curr == nil || stages == nil {
		return "", false
	}

	if replx == nil {
		xlock.Lock()
		if replx == nil {
			replx = make(map[string]*regexp.Regexp)
		}
		xlock.Unlock()
	}

	// processElement handles individual -element constructs
	processElement := func(acc func(string)) {

		if acc == nil {
			return
		}

		// element names combined with commas are treated as a prefix-separator-suffix group
		for _, stage := range stages {

			stat := stage.Type
			item := stage.Value
			prnt := stage.Parent
			match := stage.Match
			attrib := stage.Attrib
			typL := stage.TypL
			strL := stage.StrL
			intL := stage.IntL
			typR := stage.TypR
			strR := stage.StrR
			intR := stage.IntR
			norm := stage.Norm
			wildcard := stage.Wild
			unescape := (stat != INDICES && stat != ARTICLE)

			// exploreElements is a wrapper for eutils.ExploreElements, obtaining most arguments as closures
			exploreElements := func(proc func(string, int)) {
				eutils.ExploreElements(curr, mask, prnt, match, attrib, wildcard, unescape, level, proc)
			}

			// sendSlice applies optional [min:max] range restriction and sends result to accumulator
			sendSlice := func(str string) {

				// handle usual situation with no range first
				if norm {
					if wrp {
						str = html.EscapeString(str)
					}
					acc(str)
					return
				}

				// check for [after|before] variant
				if typL == STRINGRANGE || typR == STRINGRANGE {
					if strL != "" {
						// use case-insensitive test
						strL = strings.ToUpper(strL)
						idx := strings.Index(strings.ToUpper(str), strL)
						if idx < 0 {
							// specified substring must be present in original string
							return
						}
						ln := len(strL)
						// remove leading text
						str = str[idx+ln:]
					}
					if strR != "" {
						strR = strings.ToUpper(strR)
						idx := strings.Index(strings.ToUpper(str), strR)
						if idx < 0 {
							// specified substring must be present in remaining string
							return
						}
						// remove trailing text
						str = str[:idx]
					}
					if str != "" {
						if wrp {
							str = html.EscapeString(str)
						}
						acc(str)
					}
					return
				}

				min := 0
				max := 0

				// slice arguments use variable value +- adjustment or integer constant
				if typL == VARIABLERANGE {
					if strL == "" {
						return
					}
					lft, ok := variables[strL]
					if !ok {
						return
					}
					val, err := strconv.Atoi(lft)
					if err != nil {
						return
					}
					// range argument values are inclusive and 1-based, decrement variable start +- offset to use in slice
					min = val + intL - 1
				} else if typL == INTEGERRANGE {
					// range argument values are inclusive and 1-based, decrement literal start to use in slice
					min = intL - 1
				}
				if typR == VARIABLERANGE {
					if strR == "" {
						return
					}
					rgt, ok := variables[strR]
					if !ok {
						return
					}
					val, err := strconv.Atoi(rgt)
					if err != nil {
						return
					}
					if val+intR < 0 {
						// negative value is 1-based inset from end of string (undocumented)
						max = len(str) + val + intR + 1
					} else {
						max = val + intR
					}
				} else if typR == INTEGERRANGE {
					if intR < 0 {
						// negative max is inset from end of string (undocumented)
						max = len(str) + intR + 1
					} else {
						max = intR
					}
				}

				doRevComp := false
				doUpCase := false
				if status == NUCLEIC {
					// -nucleic uses direction of range to decide between forward strand or reverse complement
					if min+1 > max {
						min, max = max-1, min+1
						doRevComp = true
					}
					doUpCase = true
				}

				// numeric range now calculated, apply slice to string
				if min == 0 && max == 0 {
					if doRevComp {
						str = eutils.ReverseComplement(str)
					}
					if doUpCase {
						str = strings.ToUpper(str)
					}
					if wrp {
						str = html.EscapeString(str)
					}
					acc(str)
				} else if max == 0 {
					if min > 0 && min < len(str) {
						str = str[min:]
						if str != "" {
							if doRevComp {
								str = eutils.ReverseComplement(str)
							}
							if doUpCase {
								str = strings.ToUpper(str)
							}
							if wrp {
								str = html.EscapeString(str)
							}
							acc(str)
						}
					}
				} else if min == 0 {
					if max > 0 && max <= len(str) {
						str = str[:max]
						if str != "" {
							if doRevComp {
								str = eutils.ReverseComplement(str)
							}
							if doUpCase {
								str = strings.ToUpper(str)
							}
							if wrp {
								str = html.EscapeString(str)
							}
							acc(str)
						}
					}
				} else {
					if min < max && min > 0 && max <= len(str) {
						str = str[min:max]
						if str != "" {
							if doRevComp {
								str = eutils.ReverseComplement(str)
							}
							if doUpCase {
								str = strings.ToUpper(str)
							}
							if wrp {
								str = html.EscapeString(str)
							}
							acc(str)
						}
					}
				}
			}

			switch stat {
			case ELEMENT:
				exploreElements(func(str string, lvl int) {
					if str != "" {
						sendSlice(str)
					}
				})
			case VARIABLE, ACCUMULATOR:
				// use value of stored variable
				val, ok := variables[match]
				if ok {
					sendSlice(val)
				}
			case NUM, COUNT:
				count := 0

				exploreElements(func(str string, lvl int) {
					count++
				})

				// number of element objects
				val := strconv.Itoa(count)
				acc(val)
			case LENGTH:
				length := 0

				exploreElements(func(str string, lvl int) {
					length += len(str)
				})

				// length of element strings
				val := strconv.Itoa(length)
				acc(val)
			case DEPTH:
				exploreElements(func(str string, lvl int) {
					// depth of each element in scope
					val := strconv.Itoa(lvl)
					acc(val)
				})
			case INDEX:
				// -element "+" prints index of current XML object
				val := strconv.Itoa(index)
				acc(val)
			case INC:
				// -inc, or component of -0-based, -1-based, or -ucsc-based
				exploreElements(func(str string, lvl int) {
					if str != "" {
						num, err := strconv.Atoi(str)
						if err == nil {
							// increment value
							num++
							val := strconv.Itoa(num)
							acc(val)
						}
					}
				})
			case DEC:
				// -dec, or component of -0-based, -1-based, or -ucsc-based
				exploreElements(func(str string, lvl int) {
					if str != "" {
						num, err := strconv.Atoi(str)
						if err == nil {
							// decrement value
							num--
							val := strconv.Itoa(num)
							acc(val)
						}
					}
				})
			case QUESTION:
				acc(curr.Name)
			case TILDE:
				acc(curr.Contents)
			case STAR:
				// -element "*" prints current XML subtree on a single line
				style := SINGULARITY
				printAttrs := true

				for _, ch := range item {
					if ch == '*' {
						style++
					} else if ch == '@' {
						printAttrs = false
					}
				}
				if style > WRAPPED {
					style = WRAPPED
				}
				if style < COMPACT {
					style = COMPACT
				}

				var buffer strings.Builder

				printSubtree(curr, style, printAttrs,
					func(str string) {
						if str != "" {
							buffer.WriteString(str)
						}
					})

				txt := buffer.String()
				if txt != "" {
					acc(txt)
				}
			case DOLLAR:
				for chld := curr.Children; chld != nil; chld = chld.Next {
					acc(chld.Name)
				}
			case ATSIGN:
				if curr.Attributes != "" && curr.Attribs == nil {
					curr.Attribs = eutils.ParseAttributes(curr.Attributes)
				}
				for i := 0; i < len(curr.Attribs)-1; i += 2 {
					acc(curr.Attribs[i])
				}
			default:
				exploreElements(func(str string, lvl int) {
					if str != "" {
						sendSlice(str)
					}
				})
			}
		}
	}

	ok := false

	// format results in buffer
	var buffer strings.Builder

	buffer.WriteString(prev)
	buffer.WriteString(plg)
	buffer.WriteString(pfx)
	between := ""

	switch status {
	case ELEMENT:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case FIRST:
		single := ""

		processElement(func(str string) {
			ok = true
			if single == "" {
				single = str
			}
		})

		if single != "" {
			buffer.WriteString(between)
			buffer.WriteString(single)
			between = sep
		}

	case LAST:
		single := ""

		processElement(func(str string) {
			ok = true
			single = str
		})

		if single != "" {
			buffer.WriteString(between)
			buffer.WriteString(single)
			between = sep
		}

	case BACKWARD:
		var arry []string

		processElement(func(str string) {
			if str != "" {
				ok = true
				arry = append(arry, str)
			}
		})

		if ok {
			for i := len(arry) - 1; i >= 0; i-- {
				buffer.WriteString(between)
				buffer.WriteString(arry[i])
				between = sep
			}
		}

	case ENCODE:
		processElement(func(str string) {
			if str != "" {
				ok = true
				if !wrp {
					str = html.EscapeString(str)
				}
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case DECODE:
		// superseded by transmute -decode64 (undocumented)
		processElement(func(str string) {
			if str != "" {
				txt, err := base64.StdEncoding.DecodeString(str)
				if err == nil {
					ok = true
					buffer.WriteString(between)
					buffer.WriteString(string(txt))
					between = sep
				}
			}
		})

	case PLAIN:
		processElement(func(str string) {
			if str != "" {
				ok = true
				if eutils.IsNotASCII(str) {
					str = eutils.DoAccentTransform(str)
					if eutils.HasUnicodeMarkup(str) {
						str = eutils.RepairUnicodeMarkup(str, eutils.SPACE)
					}
				}
				if eutils.HasBadSpace(str) {
					str = eutils.CleanupBadSpaces(str)
				}
				if eutils.HasAngleBracket(str) {
					str = eutils.RepairTableMarkup(str, eutils.SPACE)
					str = eutils.RemoveEmbeddedMarkup(str)
					str = eutils.CompressRunsOfSpaces(str)
				}
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case UPPER:
		processElement(func(str string) {
			if str != "" {
				ok = true
				str = strings.ToUpper(str)
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case LOWER:
		processElement(func(str string) {
			if str != "" {
				ok = true
				str = strings.ToLower(str)
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case CHAIN:
		processElement(func(str string) {
			if str != "" {
				ok = true
				str = strings.Replace(str, " ", "_", -1)
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case TITLE:
		processElement(func(str string) {
			if str != "" {
				ok = true
				str = strings.ToLower(str)
				str = strings.Title(str)
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case AUTHOR:
		processElement(func(str string) {
			if str != "" {
				ok = true
				str = strings.Replace(str, ",", " ", -1)
				str = strings.Replace(str, ".", "", -1)
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case ORDER:
		processElement(func(str string) {
			if str != "" {
				ok = true
				str = sortStringByWords(str)
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case YEAR:
		processElement(func(str string) {
			if str != "" {
				words := strings.FieldsFunc(str, func(c rune) bool {
					return !unicode.IsDigit(c)
				})
				for _, item := range words {
					if len(item) == 4 {
						str = item
						ok = true
						// only print first year, e.g., PubDate/MedlineDate "2008 Dec-2009 Jan"
						break
					}
				}
				if ok {
					buffer.WriteString(between)
					buffer.WriteString(str)
					between = sep
				}
			}
		})

	case DOI:
		processElement(func(str string) {
			if str != "" {
				ok = true
				str = strings.TrimPrefix(str, "doi:")
				str = strings.TrimSpace(str)
				str = strings.TrimPrefix(str, "/")
				str = strings.TrimPrefix(str, "https://doi.org/")
				str = strings.TrimPrefix(str, "http://dx.doi.org/")
				str = url.QueryEscape(str)
				str = "https://doi.org/" + str
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case TRANSLATE:
		processElement(func(str string) {
			if str != "" {
				txt, found := transform[str]
				if found {
					// require successful mapping
					ok = true
					buffer.WriteString(between)
					buffer.WriteString(txt)
					between = sep
				}
			}
		})

	case REPLACE:
		processElement(func(str string) {
			if str != "" {
				re, found := replx[str]
				if !found {
					xlock.Lock()
					re, found = replx[str]
					if !found {
						nw, err := regexp.Compile(reg)
						if err == nil {
							replx[str] = nw
							re = nw
						}
					}
					xlock.Unlock()
				}
				if re != nil {
					txt := re.ReplaceAllString(str, exp)
					if txt != "" {
						ok = true
						buffer.WriteString(between)
						buffer.WriteString(txt)
						between = sep
					}
				}
			}
		})

	case VALUE:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case NUM:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case INC:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case DEC:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case ZEROBASED:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case ONEBASED:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case UCSCBASED:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case NUCLEIC:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})

	case LEN:
		length := 0

		processElement(func(str string) {
			length += len(str)
			ok = true
		})

		if ok {
			// length of element strings
			val := strconv.Itoa(length)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}

	case SUM:
		sum := 0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				sum += value
				ok = true
			}
		})

		if ok {
			// sum of element values
			val := strconv.Itoa(sum)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}

	case MIN:
		min := 0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				if !ok || value < min {
					min = value
				}
				ok = true
			}
		})

		if ok {
			// minimum of element values
			val := strconv.Itoa(min)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}

	case MAX:
		max := 0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				if !ok || value > max {
					max = value
				}
				ok = true
			}
		})

		if ok {
			// maximum of element values
			val := strconv.Itoa(max)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}

	case SUB:
		first := 0
		second := 0
		count := 0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				count++
				if count == 1 {
					first = value
				} else if count == 2 {
					second = value
				}
			}
		})

		if count == 2 {
			// must have exactly 2 elements
			ok = true
			// difference of element values
			val := strconv.Itoa(first - second)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}

	case AVG:
		sum := 0
		count := 0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				sum += value
				count++
				ok = true
			}
		})

		if ok {
			// average of element values
			avg := int(float64(sum) / float64(count))
			val := strconv.Itoa(avg)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}

	case DEV:
		count := 0
		mean := 0.0
		m2 := 0.0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				// Welford algorithm for one-pass standard deviation
				count++
				x := float64(value)
				delta := x - mean
				mean += delta / float64(count)
				m2 += delta * (x - mean)
			}
		})

		if count > 1 {
			// must have at least 2 elements
			ok = true
			// standard deviation of element values
			vrc := m2 / float64(count-1)
			dev := int(math.Sqrt(vrc))
			val := strconv.Itoa(dev)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}

	case MED:
		var arry []int
		count := 0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				arry = append(arry, value)
				count++
				ok = true
			}
		})

		if ok {
			// median of element values
			sort.Slice(arry, func(i, j int) bool { return arry[i] < arry[j] })
			med := arry[count/2]
			val := strconv.Itoa(med)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}

	case MUL:
		first := 0
		second := 0
		count := 0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				count++
				if count == 1 {
					first = value
				} else if count == 2 {
					second = value
				}
			}
		})

		if count == 2 {
			// must have exactly 2 elements
			ok = true
			// product of element values
			val := strconv.Itoa(first * second)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}

	case DIV:
		first := 0
		second := 0
		count := 0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				count++
				if count == 1 {
					first = value
				} else if count == 2 {
					second = value
				}
			}
		})

		if count == 2 {
			// must have exactly 2 elements
			ok = true
			// quotient of element values
			val := strconv.Itoa(first / second)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}

	case MOD:
		first := 0
		second := 0
		count := 0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				count++
				if count == 1 {
					first = value
				} else if count == 2 {
					second = value
				}
			}
		})

		if count == 2 {
			// must have exactly 2 elements
			ok = true
			// modulus of element values
			val := strconv.Itoa(first % second)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}

	case BIN:
		processElement(func(str string) {
			num, err := strconv.Atoi(str)
			if err == nil {
				// convert to binary representation
				val := strconv.FormatInt(int64(num), 2)
				buffer.WriteString(between)
				buffer.WriteString(val)
				between = sep
				ok = true
			}
		})

	case BIT:
		processElement(func(str string) {
			num, err := strconv.Atoi(str)
			if err == nil {
				// Kernighan algorithm for counting set bits
				count := 0
				for num != 0 {
					num &= num - 1
					count++
				}
				val := strconv.Itoa(count)
				buffer.WriteString(between)
				buffer.WriteString(val)
				between = sep
				ok = true
			}
		})

	case REVCOMP:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				str = eutils.ReverseComplement(str)
				buffer.WriteString(str)
				between = sep
			}
		})

	case FASTA:
		processElement(func(str string) {
			for str != "" {
				mx := len(str)
				if mx > 70 {
					mx = 70
				}
				item := str[:mx]
				str = str[mx:]
				ok = true
				item = strings.ToUpper(item)
				buffer.WriteString(between)
				buffer.WriteString(item)
				between = sep
			}
		})

	case NCBI2NA:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				str = eutils.Ncbi2naToIupac(str)
				buffer.WriteString(str)
				between = sep
			}
		})

	case NCBI4NA:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				str = eutils.Ncbi4naToIupac(str)
				buffer.WriteString(str)
				between = sep
			}
		})

	case MOLWT:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				str = eutils.ProteinWeight(str, true)
				buffer.WriteString(str)
				between = sep
			}
		})

	case HGVS:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				str = eutils.ParseHGVS(str)
				buffer.WriteString(str)
				between = sep
			}
		})

	case INDICES, ARTICLE:
		tiab := make(map[string][]string)
		stem := make(map[string][]string)
		titl := make(map[string][]string)

		cumulative := 0

		// mutex for inverted index
		var ilock sync.Mutex

		addItem := func(field map[string][]string, term string, position int) {

			// protect with mutex
			ilock.Lock()

			arry, found := field[term]
			if !found {
				arry = make([]string, 0, 1)
			}
			arry = append(arry, strconv.Itoa(position))
			field[term] = arry

			ilock.Unlock()
		}

		processElement(func(str string) {

			if str == "" {
				return
			}

			if str == "[Not Available]." {
				return
			}

			if eutils.IsNotASCII(str) {
				str = eutils.DoAccentTransform(str)
				if eutils.HasUnicodeMarkup(str) {
					str = eutils.RepairUnicodeMarkup(str, eutils.SPACE)
				}
			}

			str = strings.ToLower(str)

			if eutils.HasBadSpace(str) {
				str = eutils.CleanupBadSpaces(str)
			}
			if eutils.HasAngleBracket(str) {
				str = eutils.RepairEncodedMarkup(str)
				str = eutils.RepairTableMarkup(str, eutils.SPACE)
				str = eutils.RepairScriptMarkup(str, eutils.SPACE)
				str = eutils.RepairMathMLMarkup(str, eutils.SPACE)
				// RemoveEmbeddedMarkup must be called before UnescapeString, which was suppressed in eutils.ExploreElements
				str = eutils.RemoveEmbeddedMarkup(str)
			}

			if eutils.HasAmpOrNotASCII(str) {
				str = html.UnescapeString(str)
				str = strings.ToLower(str)
			}

			if eutils.IsNotASCII(str) {
				if eutils.HasGreek(str) {
					str = eutils.SpellGreek(str)
					str = eutils.CompressRunsOfSpaces(str)
				}
			}

			str = strings.Replace(str, "(", " ", -1)
			str = strings.Replace(str, ")", " ", -1)

			str = strings.Replace(str, "_", " ", -1)

			if eutils.HasHyphenOrApostrophe(str) {
				str = eutils.FixSpecialCases(str)
			}

			str = strings.Replace(str, "-", " ", -1)

			// remove trailing punctuation from each word
			var arry []string

			terms := strings.Fields(str)
			for _, item := range terms {
				max := len(item)
				for max > 1 {
					ch := item[max-1]
					if ch != '.' && ch != ',' && ch != ':' && ch != ';' {
						break
					}
					// trim trailing period, comma, colon, and semicolon
					item = item[:max-1]
					// continue checking for runs of punctuation at end
					max--
				}
				if item == "" {
					continue
				}
				arry = append(arry, item)
			}

			// rejoin into string
			cleaned := strings.Join(arry, " ")

			// break clauses at punctuation other than space or underscore, and at non-ASCII characters
			clauses := strings.FieldsFunc(cleaned, func(c rune) bool {
				return (!unicode.IsLetter(c) && !unicode.IsDigit(c)) && c != ' ' && c != '_' || c > 127
			})

			// space replaces plus sign to separate runs of unpunctuated words
			phrases := strings.Join(clauses, " ")

			// break phrases into individual words
			words := strings.Fields(phrases)

			for _, item := range words {

				cumulative++

				// skip at site of punctuation break
				if item == "+" {
					continue
				}

				// skip terms that are all digits
				if eutils.IsAllDigitsOrPeriod(item) {
					continue
				}

				// optional stop word removal
				if deStop && eutils.IsStopWord(item) {
					continue
				}

				if status == ARTICLE {
					addItem(titl, item, cumulative)
					ok = true
					continue
				}

				// index single normalized term
				addItem(tiab, item, cumulative)
				ok = true

				// apply stemming algorithm
				item = porter2.Stem(item)
				item = strings.TrimSpace(item)
				addItem(stem, item, cumulative)
			}

			// pad to avoid false positive proximity match of words in adjacent paragraphs
			rounded := ((cumulative + 99) / 100) * 100
			if rounded-cumulative < 20 {
				rounded += 100
			}
			cumulative = rounded
		})

		prepareIndices := func(field map[string][]string, label string) {

			if len(field) < 1 {
				return
			}

			var arry []string

			for item := range field {
				arry = append(arry, item)
			}

			sort.Slice(arry, func(i, j int) bool { return arry[i] < arry[j] })

			last := ""
			for _, item := range arry {
				item = strings.TrimSpace(item)
				if item == "" {
					continue
				}
				if item == last {
					// skip duplicate entry
					continue
				}
				buffer.WriteString("      <")
				buffer.WriteString(label)
				if len(field[item]) > 0 {
					buffer.WriteString(" pos=\"")
					attr := strings.Join(field[item], ",")
					buffer.WriteString(attr)
					buffer.WriteString("\"")
				}
				buffer.WriteString(">")
				buffer.WriteString(item)
				buffer.WriteString("</")
				buffer.WriteString(label)
				buffer.WriteString(">\n")
				last = item
			}
		}

		if ok {
			if status == ARTICLE {
				prepareIndices(titl, "TITL")
			} else {
				prepareIndices(tiab, "TIAB")
				prepareIndices(stem, "STEM")
			}
		}

	case TERMS:
		processElement(func(str string) {
			if str != "" {

				terms := strings.Fields(str)
				for _, item := range terms {
					max := len(item)
					for max > 1 {
						ch := item[max-1]
						if ch != '.' && ch != ',' && ch != ':' && ch != ';' {
							break
						}
						// trim trailing period, comma, colon, and semicolon
						item = item[:max-1]
						// continue checking for runs of punctuation at end
						max--
					}
					if item == "" {
						continue
					}
					ok = true
					buffer.WriteString(between)
					buffer.WriteString(item)
					between = sep
				}
			}
		})

	case WORDS:
		processElement(func(str string) {
			if str != "" {

				words := strings.FieldsFunc(str, func(c rune) bool {
					return !unicode.IsLetter(c) && !unicode.IsDigit(c)
				})
				for _, item := range words {
					item = strings.ToLower(item)
					if deStop {
						if eutils.IsStopWord(item) {
							continue
						}
					}
					if doStem {
						item = porter2.Stem(item)
						item = strings.TrimSpace(item)
					}
					if item == "" {
						continue
					}
					ok = true
					buffer.WriteString(between)
					buffer.WriteString(item)
					between = sep
				}
			}
		})

	case PAIRS:
		processElement(func(str string) {
			if str != "" {

				// break clauses at punctuation other than space or underscore, and at non-ASCII characters
				clauses := strings.FieldsFunc(str, func(c rune) bool {
					return (!unicode.IsLetter(c) && !unicode.IsDigit(c)) && c != ' ' || c > 127
				})

				// plus sign separates runs of unpunctuated words
				phrases := strings.Join(clauses, " + ")

				// break phrases into individual words
				words := strings.FieldsFunc(phrases, func(c rune) bool {
					return !unicode.IsLetter(c) && !unicode.IsDigit(c)
				})

				if len(words) > 1 {
					past := ""
					for _, item := range words {
						if item == "+" {
							past = ""
							continue
						}
						item = strings.ToLower(item)
						if deStop {
							if eutils.IsStopWord(item) {
								past = ""
								continue
							}
						}
						if doStem {
							item = porter2.Stem(item)
							item = strings.TrimSpace(item)
						}
						if item == "" {
							past = ""
							continue
						}
						if past != "" {
							ok = true
							buffer.WriteString(between)
							buffer.WriteString(past + " " + item)
							between = sep
						}
						past = item
					}
				}
			}
		})

	case REVERSE:
		processElement(func(str string) {
			if str != "" {

				words := strings.FieldsFunc(str, func(c rune) bool {
					return !unicode.IsLetter(c) && !unicode.IsDigit(c)
				})
				for lf, rt := 0, len(words)-1; lf < rt; lf, rt = lf+1, rt-1 {
					words[lf], words[rt] = words[rt], words[lf]
				}
				for _, item := range words {
					item = strings.ToLower(item)
					if deStop {
						if eutils.IsStopWord(item) {
							continue
						}
					}
					if doStem {
						item = porter2.Stem(item)
						item = strings.TrimSpace(item)
					}
					if item == "" {
						continue
					}
					ok = true
					buffer.WriteString(between)
					buffer.WriteString(item)
					between = sep
				}
			}
		})

	case LETTERS:
		processElement(func(str string) {
			if str != "" {
				for _, ch := range str {
					ok = true
					buffer.WriteString(between)
					buffer.WriteRune(ch)
					between = sep
				}
			}
		})

	case CLAUSES:
		processElement(func(str string) {
			if str != "" {

				clauses := strings.FieldsFunc(str, func(c rune) bool {
					return c == '.' || c == ',' || c == ';' || c == ':'
				})
				for _, item := range clauses {
					item = strings.ToLower(item)
					item = strings.TrimSpace(item)
					if item == "" {
						continue
					}
					ok = true
					buffer.WriteString(between)
					buffer.WriteString(item)
					between = sep
				}
			}
		})

	case MESHCODE:
		var code []string
		var tree []string

		processElement(func(str string) {
			if str != "" {
				txt, found := transform[str]
				str = strings.ToLower(str)
				code = append(code, str)
				ok = true

				if !found {
					return
				}
				txt = strings.ToLower(txt)
				txt = strings.Replace(txt, ".", "_", -1)
				codes := strings.FieldsFunc(txt, func(c rune) bool {
					return c == ','
				})
				for _, item := range codes {
					ch := item[0]
					if item == "" {
						continue
					}
					switch ch {
					case 'a', 'c', 'd', 'e', 'f', 'g', 'z':
						tree = append(tree, item)
					default:
					}
				}
			}
		})

		if len(code) > 1 {
			sort.Slice(code, func(i, j int) bool { return code[i] < code[j] })
		}
		if len(tree) > 1 {
			sort.Slice(tree, func(i, j int) bool { return tree[i] < tree[j] })
		}

		last := ""
		for _, item := range code {
			if item == last {
				// skip duplicate entry
				continue
			}
			buffer.WriteString("      <CODE>")
			buffer.WriteString(item)
			buffer.WriteString("</CODE>\n")
			last = item
		}

		last = ""
		for _, item := range tree {
			if item == last {
				// skip duplicate entry
				continue
			}
			buffer.WriteString("      <TREE>")
			buffer.WriteString(item)
			buffer.WriteString("</TREE>\n")
			last = item
		}

	case MATRIX:
		var arry []string

		processElement(func(str string) {
			if str != "" {
				txt, found := transform[str]
				if found {
					str = txt
				}
				arry = append(arry, str)
				ok = true
			}
		})

		if len(arry) > 1 {
			sort.Slice(arry, func(i, j int) bool { return arry[i] < arry[j] })

			for i, frst := range arry {
				for j, scnd := range arry {
					if i == j {
						continue
					}
					buffer.WriteString(between)
					buffer.WriteString(frst)
					buffer.WriteString("\t")
					buffer.WriteString(scnd)
					between = "\n"
				}
			}
		}

	case HISTOGRAM:
		processElement(func(str string) {
			if str != "" {
				ok = true

				hlock.Lock()

				val := histogram[str]
				val++
				histogram[str] = val

				hlock.Unlock()
			}
		})

	case ACCENTED:
		processElement(func(str string) {
			if str != "" {
				found := false
				for _, ch := range str {
					if ch > 127 {
						found = true
						break
					}
				}
				if found {
					ok = true
					buffer.WriteString(between)
					buffer.WriteString(str)
					between = sep
				}
			}
		})

	default:
	}

	// use default value if nothing written
	if !ok && def != "" {
		ok = true
		buffer.WriteString(def)
	}

	buffer.WriteString(sfx)

	if !ok {
		return "", false
	}

	txt := buffer.String()

	return txt, true
}

// processInstructions performs extraction commands on a subset of XML
func processInstructions(commands []*Operation, curr *eutils.XMLNode, mask, tab, ret string, index, level int, variables map[string]string, transform map[string]string, histogram map[string]int, accum func(string)) (string, string) {

	if accum == nil {
		return tab, ret
	}

	sep := "\t"
	pfx := ""
	sfx := ""
	plg := ""
	elg := ""
	lst := ""

	def := ""

	reg := ""
	exp := ""

	col := "\t"
	lin := "\n"

	varname := ""
	isAccum := false

	wrp := false

	plain := true
	var currColor *color.Color

	// handles color, e.g., -color "red,bold", reset to plain by -color "-" (undocumented)
	printInColor := func(str string) {
		if plain || currColor == nil {
			accum(str)
		} else {
			tx := currColor.SprintFunc()
			tmp := fmt.Sprintf("%s", tx(str))
			accum(tmp)
		}
	}

	// process commands
	for _, op := range commands {

		str := op.Value

		switch op.Type {
		case ELEMENT:
			txt, ok := processClause(curr, op.Stages, mask, tab, pfx, sfx, plg, sep, def, reg, exp, wrp, op.Type, index, level, variables, transform, histogram)
			if ok {
				plg = ""
				lst = elg
				tab = col
				ret = lin
				if plain {
					accum(txt)
				} else {
					printInColor(txt)
				}
			}
		case HISTOGRAM:
			txt, ok := processClause(curr, op.Stages, mask, "", "", "", "", "", "", "", "", wrp, op.Type, index, level, variables, transform, histogram)
			if ok {
				accum(txt)
			}
		case TAB:
			col = str
		case RET:
			lin = str
		case PFX:
			pfx = str
		case SFX:
			sfx = str
		case SEP:
			sep = str
		case LBL:
			lbl := str
			accum(tab)
			accum(plg)
			accum(pfx)
			if plain {
				accum(lbl)
			} else {
				printInColor(lbl)
			}
			accum(sfx)
			plg = ""
			lst = elg
			tab = col
			ret = lin
		case PFC:
			// preface clears previous tab and sets prefix in one command
			pfx = str
			fallthrough
		case CLR:
			// clear previous tab after the fact
			tab = ""
		case DEQ:
			// set queued tab after the fact
			tab = str
		case PLG:
			plg = str
		case ELG:
			elg = str
		case WRP:
			// shortcut to wrap elements in XML tags
			if str == "" || str == "-" {
				sep = "\t"
				pfx = ""
				sfx = ""
				plg = ""
				elg = ""
				wrp = false
				break
			}
			// -wrp with comma-separated arguments is deprecated, but supported for backward compatibility
			lft, rgt := eutils.SplitInTwoRight(str, ",")
			if lft != "" {
				plg = "<" + lft + ">"
				elg = "</" + lft + ">"
			}
			if rgt != "" && rgt != "-" {
				pfx = "<" + rgt + ">"
				sfx = "</" + rgt + ">"
				sep = "</" + rgt + "><" + rgt + ">"
			}
			wrp = true
		case ENC:
			// shortcut to mark unexpanded instances with XML tags
			plg = ""
			elg = ""
			if str != "" && str != "-" {
				items := strings.Split(str, "/")
				for i := 0; i < len(items); i++ {
					plg += "<" + items[i] + ">"
				}
				for i := len(items) - 1; i >= 0; i-- {
					elg += "</" + items[i] + ">"
				}
			}
		case RST:
			pfx = ""
			sfx = ""
			plg = ""
			elg = ""
			sep = "\t"
			def = ""
			wrp = false
		case DEF:
			def = str
		case REG:
			reg = str
		case EXP:
			exp = str
		case COLOR:
			currColor = color.New()
			if str == "-" || str == "reset" || str == "clear" {
				plain = true
				break
			}
			plain = false
			items := strings.Split(str, ",")
			for _, itm := range items {
				switch itm {
				case "red":
					currColor.Add(color.FgRed)
				case "grn", "green":
					currColor.Add(color.FgGreen)
				case "blu", "blue":
					currColor.Add(color.FgBlue)
				case "blk", "black":
					currColor.Add(color.FgBlack)
				case "bld", "bold":
					currColor.Add(color.Bold)
				case "ital", "italic", "italics":
					currColor.Add(color.Italic)
				case "blink", "flash":
					currColor.Add(color.BlinkSlow)
				default:
					fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized color argument '%s'\n", itm)
					os.Exit(1)
				}
			}
		case ACCUMULATOR:
			isAccum = true
			varname = str
		case VARIABLE:
			isAccum = false
			varname = str
		case VALUE:
			length := len(str)
			if length > 1 && str[0] == '(' && str[length-1] == ')' {
				// set variable from literal text inside parentheses, e.g., -COM "(, )"
				variables[varname] = str[1 : length-1]
				// -if "&VARIABLE" will succeed if set to blank with empty parentheses "()"
			} else if str == "" {
				// -if "&VARIABLE" will fail if initialized with empty string ""
				delete(variables, varname)
			} else {
				txt, ok := processClause(curr, op.Stages, mask, "", pfx, sfx, plg, sep, def, reg, exp, wrp, op.Type, index, level, variables, transform, histogram)
				if ok {
					plg = ""
					lst = elg
					if isAccum {
						if variables[varname] == "" {
							variables[varname] = txt
						} else {
							variables[varname] += sep + txt
						}
					} else {
						variables[varname] = txt
					}
				}
			}
			varname = ""
			isAccum = false
		default:
			txt, ok := processClause(curr, op.Stages, mask, tab, pfx, sfx, plg, sep, def, reg, exp, wrp, op.Type, index, level, variables, transform, histogram)
			if ok {
				plg = ""
				lst = elg
				tab = col
				ret = lin
				if plain {
					accum(txt)
				} else {
					printInColor(txt)
				}
			}
		}
	}

	if plain {
		accum(lst)
	} else {
		printInColor(lst)
	}

	return tab, ret
}

// CONDITIONAL EXECUTION USES -if AND -unless STATEMENT, WITH SUPPORT FOR DEPRECATED -match AND -avoid STATEMENTS

// conditionsAreSatisfied tests a set of conditions to determine if extraction should proceed
func conditionsAreSatisfied(conditions []*Operation, curr *eutils.XMLNode, mask string, index, level int, variables map[string]string) bool {

	if curr == nil {
		return false
	}

	required := 0
	observed := 0
	forbidden := 0
	isMatch := false
	isAvoid := false

	// matchFound tests individual conditions
	matchFound := func(stages []*Step) bool {

		if stages == nil || len(stages) < 1 {
			return false
		}

		stage := stages[0]

		var constraint *Step

		if len(stages) > 1 {
			constraint = stages[1]
		}

		status := stage.Type
		prnt := stage.Parent
		match := stage.Match
		attrib := stage.Attrib
		typL := stage.TypL
		strL := stage.StrL
		intL := stage.IntL
		typR := stage.TypR
		strR := stage.StrR
		intR := stage.IntR
		norm := stage.Norm
		wildcard := stage.Wild
		unescape := true

		found := false
		number := ""

		// exploreElements is a wrapper for eutils.ExploreElements, obtaining most arguments as closures
		exploreElements := func(proc func(string, int)) {
			eutils.ExploreElements(curr, mask, prnt, match, attrib, wildcard, unescape, level, proc)
		}

		// test string or numeric constraints
		testConstraint := func(str string) bool {

			if str == "" || constraint == nil {
				return false
			}

			val := constraint.Value
			stat := constraint.Type

			switch stat {
			case EQUALS, CONTAINS, INCLUDES, ISWITHIN, STARTSWITH, ENDSWITH, ISNOT, ISBEFORE, ISAFTER, MATCHES, RESEMBLES:
				// substring test on element values
				str = strings.ToUpper(str)
				val = strings.ToUpper(val)

				switch stat {
				case EQUALS:
					if str == val {
						return true
					}
				case CONTAINS:
					if strings.Contains(str, val) {
						return true
					}
				case INCLUDES:
					str = strings.TrimSpace(str)
					val = strings.TrimSpace(val)
					if strings.Contains(" "+str+" ", " "+val+" ") {
						return true
					}
				case ISWITHIN:
					if strings.Contains(val, str) {
						return true
					}
				case STARTSWITH:
					if strings.HasPrefix(str, val) {
						return true
					}
				case ENDSWITH:
					if strings.HasSuffix(str, val) {
						return true
					}
				case ISNOT:
					if str != val {
						return true
					}
				case ISBEFORE:
					if str < val {
						return true
					}
				case ISAFTER:
					if str > val {
						return true
					}
				case MATCHES:
					if removeCommaOrSemicolon(str) == strings.ToLower(val) {
						return true
					}
				case RESEMBLES:
					if sortStringByWords(str) == strings.ToLower(val) {
						return true
					}
				default:
				}
			case ISEQUALTO, DIFFERSFROM:
				// conditional argument is element specifier
				if constraint.Parent != "" || constraint.Match != "" || constraint.Attrib != "" {
					ch := val[0]
					// pound, percent, and caret prefixes supported (undocumented)
					switch ch {
					case '#':
						count := 0
						eutils.ExploreElements(curr, mask, constraint.Parent, constraint.Match, constraint.Attrib, constraint.Wild, true, level, func(stn string, lvl int) {
							count++
						})
						val = strconv.Itoa(count)
					case '%':
						length := 0
						eutils.ExploreElements(curr, mask, constraint.Parent, constraint.Match, constraint.Attrib, constraint.Wild, true, level, func(stn string, lvl int) {
							if stn != "" {
								length += len(stn)
							}
						})
						val = strconv.Itoa(length)
					case '^':
						depth := 0
						eutils.ExploreElements(curr, mask, constraint.Parent, constraint.Match, constraint.Attrib, constraint.Wild, true, level, func(stn string, lvl int) {
							depth = lvl
						})
						val = strconv.Itoa(depth)
					default:
						eutils.ExploreElements(curr, mask, constraint.Parent, constraint.Match, constraint.Attrib, constraint.Wild, true, level, func(stn string, lvl int) {
							if stn != "" {
								val = stn
							}
						})
					}
				}
				str = strings.ToUpper(str)
				val = strings.ToUpper(val)

				switch stat {
				case ISEQUALTO:
					if str == val {
						return true
					}
				case DIFFERSFROM:
					if str != val {
						return true
					}
				default:
				}
			case GT, GE, LT, LE, EQ, NE:
				// second argument of numeric test can be element specifier
				if constraint.Parent != "" || constraint.Match != "" || constraint.Attrib != "" {
					ch := val[0]
					// pound, percent, and caret prefixes supported as potentially useful for data QA (undocumented)
					switch ch {
					case '#':
						count := 0
						eutils.ExploreElements(curr, mask, constraint.Parent, constraint.Match, constraint.Attrib, constraint.Wild, true, level, func(stn string, lvl int) {
							count++
						})
						val = strconv.Itoa(count)
					case '%':
						length := 0
						eutils.ExploreElements(curr, mask, constraint.Parent, constraint.Match, constraint.Attrib, constraint.Wild, true, level, func(stn string, lvl int) {
							if stn != "" {
								length += len(stn)
							}
						})
						val = strconv.Itoa(length)
					case '^':
						depth := 0
						eutils.ExploreElements(curr, mask, constraint.Parent, constraint.Match, constraint.Attrib, constraint.Wild, true, level, func(stn string, lvl int) {
							depth = lvl
						})
						val = strconv.Itoa(depth)
					default:
						eutils.ExploreElements(curr, mask, constraint.Parent, constraint.Match, constraint.Attrib, constraint.Wild, true, level, func(stn string, lvl int) {
							if stn != "" {
								_, errz := strconv.Atoi(stn)
								if errz == nil {
									val = stn
								}
							}
						})
					}
				}

				// numeric tests on element values
				x, errx := strconv.Atoi(str)
				y, erry := strconv.Atoi(val)

				// both arguments must resolve to integers
				if errx != nil || erry != nil {
					return false
				}

				switch stat {
				case GT:
					if x > y {
						return true
					}
				case GE:
					if x >= y {
						return true
					}
				case LT:
					if x < y {
						return true
					}
				case LE:
					if x <= y {
						return true
					}
				case EQ:
					if x == y {
						return true
					}
				case NE:
					if x != y {
						return true
					}
				default:
				}
			default:
			}

			return false
		}

		// checkConstraint applies optional [min:max] range restriction and sends result to testConstraint
		checkConstraint := func(str string) bool {

			// handle usual situation with no range first
			if norm {
				return testConstraint(str)
			}

			// check for [after|before] variant
			if typL == STRINGRANGE || typR == STRINGRANGE {
				if strL != "" {
					// use case-insensitive test
					strL = strings.ToUpper(strL)
					idx := strings.Index(strings.ToUpper(str), strL)
					if idx < 0 {
						// specified substring must be present in original string
						return false
					}
					ln := len(strL)
					// remove leading text
					str = str[idx+ln:]
				}
				if strR != "" {
					strR = strings.ToUpper(strR)
					idx := strings.Index(strings.ToUpper(str), strR)
					if idx < 0 {
						// specified substring must be present in remaining string
						return false
					}
					// remove trailing text
					str = str[:idx]
				}
				if str != "" {
					return testConstraint(str)
				}
				return false
			}

			min := 0
			max := 0

			// slice arguments use variable value +- adjustment or integer constant
			if typL == VARIABLERANGE {
				if strL == "" {
					return false
				}
				lft, ok := variables[strL]
				if !ok {
					return false
				}
				val, err := strconv.Atoi(lft)
				if err != nil {
					return false
				}
				// range argument values are inclusive and 1-based, decrement variable start +- offset to use in slice
				min = val + intL - 1
			} else if typL == INTEGERRANGE {
				// range argument values are inclusive and 1-based, decrement literal start to use in slice
				min = intL - 1
			}
			if typR == VARIABLERANGE {
				if strR == "" {
					return false
				}
				rgt, ok := variables[strR]
				if !ok {
					return false
				}
				val, err := strconv.Atoi(rgt)
				if err != nil {
					return false
				}
				if val+intR < 0 {
					// negative value is 1-based inset from end of string (undocumented)
					max = len(str) + val + intR + 1
				} else {
					max = val + intR
				}
			} else if typR == INTEGERRANGE {
				if intR < 0 {
					// negative max is inset from end of string (undocumented)
					max = len(str) + intR + 1
				} else {
					max = intR
				}
			}

			// numeric range now calculated, apply slice to string
			if min == 0 && max == 0 {
				return testConstraint(str)
			} else if max == 0 {
				if min > 0 && min < len(str) {
					str = str[min:]
					if str != "" {
						return testConstraint(str)
					}
				}
			} else if min == 0 {
				if max > 0 && max <= len(str) {
					str = str[:max]
					if str != "" {
						return testConstraint(str)
					}
				}
			} else {
				if min < max && min > 0 && max <= len(str) {
					str = str[min:max]
					if str != "" {
						return testConstraint(str)
					}
				}
			}

			return false
		}

		switch status {
		case ELEMENT:
			exploreElements(func(str string, lvl int) {
				// match to XML container object sends empty string, so do not check for str != "" here
				// test every selected element individually if value is specified
				if constraint == nil || checkConstraint(str) {
					found = true
				}
			})
		case VARIABLE:
			// use value of stored variable
			str, ok := variables[match]
			if ok {
				//  -if &VARIABLE -equals VALUE is the supported construct
				if constraint == nil || checkConstraint(str) {
					found = true
				}
			}
		case COUNT:
			count := 0

			exploreElements(func(str string, lvl int) {
				count++
				found = true
			})

			// number of element objects
			number = strconv.Itoa(count)
		case LENGTH:
			length := 0

			exploreElements(func(str string, lvl int) {
				length += len(str)
				found = true
			})

			// length of element strings
			number = strconv.Itoa(length)
		case DEPTH:
			depth := 0

			exploreElements(func(str string, lvl int) {
				depth = lvl
				found = true
			})

			// depth of last element in scope
			number = strconv.Itoa(depth)
		case INDEX:
			// index of explored parent object
			number = strconv.Itoa(index)
			found = true
		default:
		}

		if number == "" {
			return found
		}

		if constraint == nil || checkConstraint(number) {
			return true
		}

		return false
	}

	// test conditional arguments
	for _, op := range conditions {

		switch op.Type {
		// -if tests for presence of element (deprecated -match can test element:value)
		case SELECT, IF, MATCH:
			// checking for failure here allows for multiple -if [ -and / -or ] clauses
			if isMatch && observed < required {
				return false
			}
			if isAvoid && forbidden > 0 {
				return false
			}
			required = 0
			observed = 0
			forbidden = 0
			isMatch = true
			isAvoid = false
			// continue on to next two cases
			fallthrough
		case AND:
			required++
			// continue on to next case
			fallthrough
		case OR:
			if matchFound(op.Stages) {
				observed++
				// record presence of forbidden element if in -unless clause
				forbidden++
			}
		// -unless tests for absence of element, or presence but with failure of subsequent value test (deprecated -avoid can test element:value)
		case UNLESS, AVOID:
			if isMatch && observed < required {
				return false
			}
			if isAvoid && forbidden > 0 {
				return false
			}
			required = 0
			observed = 0
			forbidden = 0
			isMatch = false
			isAvoid = true
			if matchFound(op.Stages) {
				forbidden++
			}
		default:
		}
	}

	if isMatch && observed < required {
		return false
	}
	if isAvoid && forbidden > 0 {
		return false
	}

	return true
}

// RECURSIVELY PROCESS EXPLORATION COMMANDS AND XML DATA STRUCTURE

// processCommands visits XML nodes, performs conditional tests, and executes data extraction instructions
func processCommands(cmds *Block, curr *eutils.XMLNode, tab, ret string, index, level int, variables map[string]string, transform map[string]string, histogram map[string]int, accum func(string)) (string, string) {

	if accum == nil {
		return tab, ret
	}

	prnt := cmds.Parent
	match := cmds.Match

	// closure passes local variables to callback, which can modify caller tab and ret values
	processNode := func(node *eutils.XMLNode, idx, lvl int) {

		// apply -if or -unless tests
		if conditionsAreSatisfied(cmds.Conditions, node, match, idx, lvl, variables) {

			// execute data extraction commands
			if len(cmds.Commands) > 0 {
				tab, ret = processInstructions(cmds.Commands, node, match, tab, ret, idx, lvl, variables, transform, histogram, accum)
			}

			// process sub commands on child node
			for _, sub := range cmds.Subtasks {
				tab, ret = processCommands(sub, node, tab, ret, 1, lvl, variables, transform, histogram, accum)
			}

		} else {

			// execute commands after -else statement
			if len(cmds.Failure) > 0 {
				tab, ret = processInstructions(cmds.Failure, node, match, tab, ret, idx, lvl, variables, transform, histogram, accum)
			}
		}
	}

	// explorePath recursive definition
	var explorePath func(*eutils.XMLNode, []string, int, int, func(*eutils.XMLNode, int, int)) int

	// explorePath visits child nodes and matches against next entry in path
	explorePath = func(curr *eutils.XMLNode, path []string, indx, levl int, proc func(*eutils.XMLNode, int, int)) int {

		if curr == nil || proc == nil {
			return indx
		}

		if len(path) < 1 {
			proc(curr, indx, levl)
			indx++
			return indx
		}

		name := path[0]
		rest := path[1:]

		// explore next level of child nodes
		for chld := curr.Children; chld != nil; chld = chld.Next {
			if chld.Name == name {
				// recurse only if child matches next component in path
				indx = explorePath(chld, rest, indx, levl+1, proc)
			}
		}

		return indx
	}

	if cmds.Foreword != "" {
		accum(cmds.Foreword)
	}

	// apply -position test

	if cmds.Position == "" || cmds.Position == "all" {

		eutils.ExploreNodes(curr, prnt, match, index, level, processNode)

	} else if cmds.Position == "path" {

		eutils.ExploreNodes(curr, prnt, match, index, level,
			func(node *eutils.XMLNode, idx, lvl int) {
				// exploreNodes callback has matched first path component, now explore remainder one level and component at a time
				explorePath(node, cmds.Path, idx, lvl, processNode)
			})

	} else {

		var single *eutils.XMLNode
		lev := 0
		ind := 0

		if cmds.Position == "first" {

			eutils.ExploreNodes(curr, prnt, match, index, level,
				func(node *eutils.XMLNode, idx, lvl int) {
					if single == nil {
						single = node
						ind = idx
						lev = lvl
					}
				})

		} else if cmds.Position == "last" {

			eutils.ExploreNodes(curr, prnt, match, index, level,
				func(node *eutils.XMLNode, idx, lvl int) {
					single = node
					ind = idx
					lev = lvl
				})

		} else if cmds.Position == "outer" {

			// print only first and last nodes
			var beg *Limiter
			var end *Limiter

			eutils.ExploreNodes(curr, prnt, match, index, level,
				func(node *eutils.XMLNode, idx, lvl int) {
					if beg == nil {
						beg = &Limiter{node, idx, lvl}
					} else {
						end = &Limiter{node, idx, lvl}
					}
				})

			if beg != nil {
				processNode(beg.Obj, beg.Idx, beg.Lvl)
			}
			if end != nil {
				processNode(end.Obj, end.Idx, end.Lvl)
			}

		} else if cmds.Position == "inner" {

			// print all but first and last nodes
			var prev *Limiter
			var next *Limiter
			first := true

			eutils.ExploreNodes(curr, prnt, match, index, level,
				func(node *eutils.XMLNode, idx, lvl int) {
					if first {
						first = false
						return
					}

					prev = next
					next = &Limiter{node, idx, lvl}

					if prev != nil {
						processNode(prev.Obj, prev.Idx, prev.Lvl)
					}
				})

		} else if cmds.Position == "even" {

			okay := false

			eutils.ExploreNodes(curr, prnt, match, index, level,
				func(node *eutils.XMLNode, idx, lvl int) {
					if okay {
						processNode(node, idx, lvl)
					}
					okay = !okay
				})

		} else if cmds.Position == "odd" {

			okay := true

			eutils.ExploreNodes(curr, prnt, match, index, level,
				func(node *eutils.XMLNode, idx, lvl int) {
					if okay {
						processNode(node, idx, lvl)
					}
					okay = !okay
				})

		} else {

			// use numeric position
			number, err := strconv.Atoi(cmds.Position)
			if err == nil {

				pos := 0

				eutils.ExploreNodes(curr, prnt, match, index, level,
					func(node *eutils.XMLNode, idx, lvl int) {
						pos++
						if pos == number {
							single = node
							ind = idx
							lev = lvl
						}
					})

			} else {

				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized position '%s'\n", cmds.Position)
				os.Exit(1)
			}
		}

		if single != nil {
			processNode(single, ind, lev)
		}
	}

	if cmds.Afterword != "" {
		accum(cmds.Afterword)
	}

	return tab, ret
}

// PROCESS ONE XML COMPONENT RECORD

// processQuery perform data extraction driven by command-line arguments
func processQuery(text, parent string, index int, hd, tl string, transform map[string]string, histogram map[string]int, cmds *Block) string {

	if text == "" || cmds == nil {
		return ""
	}

	// exit from function will collect garbage of node structure for current XML object
	pat := eutils.ParseRecord(text, parent)

	if pat == nil {
		return ""
	}

	// exit from function will also free map of recorded variables for current -pattern
	variables := make(map[string]string)

	var buffer strings.Builder

	ok := false

	if hd != "" {
		buffer.WriteString(hd)
	}

	ret := ""

	if cmds.Position == "select" {

		if conditionsAreSatisfied(cmds.Conditions, pat, cmds.Match, index, 1, variables) {
			ok = true
			buffer.WriteString(text)
			ret = "\n"
		}

	} else {

		// start processing at top of command tree and top of XML subregion selected by -pattern
		_, ret = processCommands(cmds, pat, "", "", index, 1, variables, transform, histogram,
			func(str string) {
				if str != "" {
					ok = true
					buffer.WriteString(str)
				}
			})
	}

	if tl != "" {
		buffer.WriteString(tl)
	}

	if ret != "" {
		ok = true
		buffer.WriteString(ret)
	}

	txt := buffer.String()

	// remove leading newline (-insd -pfx artifact)
	if txt != "" && txt[0] == '\n' {
		txt = txt[1:]
	}

	if !ok {
		return ""
	}

	// return consolidated result string
	return txt
}

// INSDSEQ EXTRACTION COMMAND GENERATOR

// e.g., xtract -insd complete mat_peptide "%peptide" product peptide

// processINSD generates extraction commands for GenBank/RefSeq records in INSDSet format
func processINSD(args []string, isPipe, addDash, doIndex bool) []string {

	// legal GenBank / GenPept / RefSeq features

	features := []string{
		"-10_signal",
		"-35_signal",
		"3'clip",
		"3'UTR",
		"5'clip",
		"5'UTR",
		"allele",
		"assembly_gap",
		"attenuator",
		"Bond",
		"C_region",
		"CAAT_signal",
		"CDS",
		"centromere",
		"conflict",
		"D_segment",
		"D-loop",
		"enhancer",
		"exon",
		"gap",
		"GC_signal",
		"gene",
		"iDNA",
		"intron",
		"J_segment",
		"LTR",
		"mat_peptide",
		"misc_binding",
		"misc_difference",
		"misc_feature",
		"misc_recomb",
		"misc_RNA",
		"misc_signal",
		"misc_structure",
		"mobile_element",
		"modified_base",
		"mRNA",
		"mutation",
		"N_region",
		"ncRNA",
		"old_sequence",
		"operon",
		"oriT",
		"polyA_signal",
		"polyA_site",
		"precursor_RNA",
		"prim_transcript",
		"primer_bind",
		"promoter",
		"propeptide",
		"protein_bind",
		"Protein",
		"RBS",
		"Region",
		"regulatory",
		"rep_origin",
		"repeat_region",
		"repeat_unit",
		"rRNA",
		"S_region",
		"satellite",
		"scRNA",
		"sig_peptide",
		"Site",
		"snoRNA",
		"snRNA",
		"source",
		"stem_loop",
		"STS",
		"TATA_signal",
		"telomere",
		"terminator",
		"tmRNA",
		"transit_peptide",
		"tRNA",
		"unsure",
		"V_region",
		"V_segment",
		"variation",
	}

	// legal GenBank / GenPept / RefSeq qualifiers

	qualifiers := []string{
		"allele",
		"altitude",
		"anticodon",
		"artificial_location",
		"bio_material",
		"bond_type",
		"bound_moiety",
		"breed",
		"calculated_mol_wt",
		"cell_line",
		"cell_type",
		"chloroplast",
		"chromoplast",
		"chromosome",
		"circular_RNA",
		"citation",
		"clone_lib",
		"clone",
		"coded_by",
		"codon_start",
		"codon",
		"collected_by",
		"collection_date",
		"compare",
		"cons_splice",
		"country",
		"cultivar",
		"culture_collection",
		"cyanelle",
		"db_xref",
		"derived_from",
		"dev_stage",
		"direction",
		"EC_number",
		"ecotype",
		"encodes",
		"endogenous_virus",
		"environmental_sample",
		"estimated_length",
		"evidence",
		"exception",
		"experiment",
		"focus",
		"frequency",
		"function",
		"gap_type",
		"gdb_xref",
		"gene_synonym",
		"gene",
		"germline",
		"haplogroup",
		"haplotype",
		"host",
		"identified_by",
		"inference",
		"insertion_seq",
		"isolate",
		"isolation_source",
		"kinetoplast",
		"lab_host",
		"label",
		"lat_lon",
		"linkage_evidence",
		"locus_tag",
		"macronuclear",
		"map",
		"mating_type",
		"metagenome_source",
		"metagenomic",
		"mitochondrion",
		"mobile_element_type",
		"mobile_element",
		"mod_base",
		"mol_type",
		"name",
		"nat_host",
		"ncRNA_class",
		"non_functional",
		"note",
		"number",
		"old_locus_tag",
		"operon",
		"organelle",
		"organism",
		"partial",
		"PCR_conditions",
		"PCR_primers",
		"peptide",
		"phenotype",
		"plasmid",
		"pop_variant",
		"product",
		"protein_id",
		"proviral",
		"pseudo",
		"pseudogene",
		"rearranged",
		"recombination_class",
		"region_name",
		"regulatory_class",
		"replace",
		"ribosomal_slippage",
		"rpt_family",
		"rpt_type",
		"rpt_unit_range",
		"rpt_unit_seq",
		"rpt_unit",
		"satellite",
		"segment",
		"sequenced_mol",
		"serotype",
		"serovar",
		"sex",
		"site_type",
		"specific_host",
		"specimen_voucher",
		"standard_name",
		"strain",
		"structural_class",
		"sub_clone",
		"sub_species",
		"sub_strain",
		"submitter_seqid",
		"tag_peptide",
		"tissue_lib",
		"tissue_type",
		"trans_splicing",
		"transcript_id",
		"transcription",
		"transgenic",
		"transl_except",
		"transl_table",
		"translation",
		"transposon",
		"type_material",
		"UniProtKB_evidence",
		"usedin",
		"variety",
		"virion",
	}

	// legal INSDSeq XML fields

	insdtags := []string{
		"INSDAltSeqData_items",
		"INSDAltSeqData",
		"INSDAltSeqItem_first-accn",
		"INSDAltSeqItem_gap-comment",
		"INSDAltSeqItem_gap-length",
		"INSDAltSeqItem_gap-linkage",
		"INSDAltSeqItem_gap-type",
		"INSDAltSeqItem_interval",
		"INSDAltSeqItem_isgap",
		"INSDAltSeqItem_isgap@value",
		"INSDAltSeqItem_last-accn",
		"INSDAltSeqItem_value",
		"INSDAltSeqItem",
		"INSDAuthor",
		"INSDComment_paragraphs",
		"INSDComment_type",
		"INSDComment",
		"INSDCommentParagraph",
		"INSDFeature_intervals",
		"INSDFeature_key",
		"INSDFeature_location",
		"INSDFeature_operator",
		"INSDFeature_partial3",
		"INSDFeature_partial3@value",
		"INSDFeature_partial5",
		"INSDFeature_partial5@value",
		"INSDFeature_quals",
		"INSDFeature_xrefs",
		"INSDFeature",
		"INSDFeatureSet_annot-source",
		"INSDFeatureSet_features",
		"INSDFeatureSet",
		"INSDInterval_accession",
		"INSDInterval_from",
		"INSDInterval_interbp",
		"INSDInterval_interbp@value",
		"INSDInterval_iscomp",
		"INSDInterval_iscomp@value",
		"INSDInterval_point",
		"INSDInterval_to",
		"INSDInterval",
		"INSDKeyword",
		"INSDQualifier_name",
		"INSDQualifier_value",
		"INSDQualifier",
		"INSDReference_authors",
		"INSDReference_consortium",
		"INSDReference_journal",
		"INSDReference_position",
		"INSDReference_pubmed",
		"INSDReference_reference",
		"INSDReference_remark",
		"INSDReference_title",
		"INSDReference_xref",
		"INSDReference",
		"INSDSecondary-accn",
		"INSDSeq_accession-version",
		"INSDSeq_alt-seq",
		"INSDSeq_comment-set",
		"INSDSeq_comment",
		"INSDSeq_contig",
		"INSDSeq_create-date",
		"INSDSeq_create-release",
		"INSDSeq_database-reference",
		"INSDSeq_definition",
		"INSDSeq_division",
		"INSDSeq_entry-version",
		"INSDSeq_feature-set",
		"INSDSeq_feature-table",
		"INSDSeq_keywords",
		"INSDSeq_length",
		"INSDSeq_locus",
		"INSDSeq_moltype",
		"INSDSeq_organism",
		"INSDSeq_other-seqids",
		"INSDSeq_primary-accession",
		"INSDSeq_primary",
		"INSDSeq_project",
		"INSDSeq_references",
		"INSDSeq_secondary-accessions",
		"INSDSeq_segment",
		"INSDSeq_sequence",
		"INSDSeq_source-db",
		"INSDSeq_source",
		"INSDSeq_strandedness",
		"INSDSeq_struc-comments",
		"INSDSeq_taxonomy",
		"INSDSeq_topology",
		"INSDSeq_update-date",
		"INSDSeq_update-release",
		"INSDSeq_xrefs",
		"INSDSeq",
		"INSDSeqid",
		"INSDSet",
		"INSDStrucComment_items",
		"INSDStrucComment_name",
		"INSDStrucComment",
		"INSDStrucCommentItem_tag",
		"INSDStrucCommentItem_url",
		"INSDStrucCommentItem_value",
		"INSDStrucCommentItem",
		"INSDXref_dbname",
		"INSDXref_id",
		"INSDXref",
	}

	checkAgainstVocabulary := func(str, objtype string, arry []string) {

		if str == "" || arry == nil {
			return
		}

		// skip past pound, percent, or caret character at beginning of string
		if len(str) > 1 {
			switch str[0] {
			case '#', '%', '^':
				str = str[1:]
			default:
			}
		}

		for _, txt := range arry {
			if str == txt {
				return
			}
			if strings.ToUpper(str) == strings.ToUpper(txt) {
				fmt.Fprintf(os.Stderr, "\nERROR: Incorrect capitalization of '%s' %s, change to '%s'\n", str, objtype, txt)
				os.Exit(1)
			}
		}

		fmt.Fprintf(os.Stderr, "\nERROR: Item '%s' is not a legal -insd %s\n", str, objtype)
		os.Exit(1)
	}

	var acc []string

	max := len(args)
	if max < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract -insd\n")
		os.Exit(1)
	}

	if doIndex {
		if isPipe {
			acc = append(acc, "-head", "<IdxDocumentSet>", "-tail", "</IdxDocumentSet>")
			acc = append(acc, "-hd", "  <IdxDocument>\n", "-tl", "  </IdxDocument>")
			acc = append(acc, "-pattern", "INSDSeq", "-pfx", "    <IdxUid>", "-sfx", "</IdxUid>\n")
			acc = append(acc, "-element", "INSDSeq_accession-version", "-clr", "-rst", "-tab", "\n")
		} else {
			acc = append(acc, "-head", "\"<IdxDocumentSet>\"", "-tail", "\"</IdxDocumentSet>\"")
			acc = append(acc, "-hd", "\"  <IdxDocument>\\n\"", "-tl", "\"  </IdxDocument>\"")
			acc = append(acc, "-pattern", "INSDSeq", "-pfx", "\"    <IdxUid>\"", "-sfx", "\"</IdxUid>\\n\"")
			acc = append(acc, "-element", "INSDSeq_accession-version", "-clr", "-rst", "-tab", "\\n")
		}
	} else {
		acc = append(acc, "-pattern", "INSDSeq", "-ACCN", "INSDSeq_accession-version")
		acc = append(acc, "-LCUS", "INSDSeq_locus", "-SEQ", "INSDSeq_sequence")
	}

	if doIndex {
		if isPipe {
			acc = append(acc, "-group", "INSDSeq", "-lbl", "    <IdxSearchFields>\n")
		} else {
			acc = append(acc, "-group", "INSDSeq", "-lbl", "\"    <IdxSearchFields>\\n\"")
		}
	}

	printAccn := true

	// collect descriptors

	if strings.HasPrefix(args[0], "INSD") {

		if doIndex {
			acc = append(acc, "-clr", "-indices")
		} else {
			if isPipe {
				acc = append(acc, "-clr", "-pfx", "\\n", "-element", "&ACCN")
				acc = append(acc, "-group", "INSDSeq", "-sep", "|", "-element")
			} else {
				acc = append(acc, "-clr", "-pfx", "\"\\n\"", "-element", "\"&ACCN\"")
				acc = append(acc, "-group", "INSDSeq", "-sep", "\"|\"", "-element")
			}
			printAccn = false
		}

		for {
			if len(args) < 1 {
				return acc
			}
			str := args[0]
			if !strings.HasPrefix(args[0], "INSD") {
				break
			}
			checkAgainstVocabulary(str, "element", insdtags)
			acc = append(acc, str)
			args = args[1:]
		}

	} else if strings.HasPrefix(strings.ToUpper(args[0]), "INSD") {

		// report capitalization or vocabulary failure
		checkAgainstVocabulary(args[0], "element", insdtags)

		// program should not get to this point, but warn and exit anyway
		fmt.Fprintf(os.Stderr, "\nERROR: Item '%s' is not a legal -insd %s\n", args[0], "element")
		os.Exit(1)
	}

	// collect qualifiers

	partial := false
	complete := false

	if args[0] == "+" || args[0] == "complete" {
		complete = true
		args = args[1:]
		max--
	} else if args[0] == "-" || args[0] == "partial" {
		partial = true
		args = args[1:]
		max--
	}

	if max < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: No feature key supplied to xtract -insd\n")
		os.Exit(1)
	}

	acc = append(acc, "-group", "INSDFeature")

	// limit to designated features

	feature := args[0]

	fcmd := "-if"

	// can specify multiple features separated by plus sign (e.g., CDS+mRNA) or comma (e.g., CDS,mRNA)
	plus := strings.Split(feature, "+")
	for _, pls := range plus {
		comma := strings.Split(pls, ",")
		for _, cma := range comma {

			checkAgainstVocabulary(cma, "feature", features)
			acc = append(acc, fcmd, "INSDFeature_key", "-equals", cma)

			fcmd = "-or"
		}
	}

	if max < 2 {
		// still need at least one qualifier even on legal feature
		fmt.Fprintf(os.Stderr, "\nERROR: Feature '%s' must be followed by at least one qualifier\n", feature)
		os.Exit(1)
	}

	args = args[1:]

	if complete {
		acc = append(acc, "-unless", "INSDFeature_partial5", "-or", "INSDFeature_partial3")
	} else if partial {
		acc = append(acc, "-if", "INSDFeature_partial5", "-or", "INSDFeature_partial3")
	}

	if printAccn {
		if doIndex {
		} else {
			if isPipe {
				acc = append(acc, "-clr", "-pfx", "\\n", "-first", "&ACCN,&LCUS")
			} else {
				acc = append(acc, "-clr", "-pfx", "\"\\n\"", "-first", "\"&ACCN,&LCUS\"")
			}
		}
	}

	for _, str := range args {

		if str == "mol_wt" {
			str = "calculated_mol_wt"
		}

		if strings.HasPrefix(str, "INSD") {

			checkAgainstVocabulary(str, "element", insdtags)
			if doIndex {
				acc = append(acc, "-block", "INSDFeature", "-clr", "-indices")
			} else {
				if isPipe {
					acc = append(acc, "-block", "INSDFeature", "-sep", "|", "-element")
				} else {
					acc = append(acc, "-block", "INSDFeature", "-sep", "\"|\"", "-element")
				}
			}
			acc = append(acc, str)
			if addDash {
				acc = append(acc, "-block", "INSDFeature", "-unless", str)
				if strings.HasSuffix(str, "@value") {
					if isPipe {
						acc = append(acc, "-lbl", "false")
					} else {
						acc = append(acc, "-lbl", "\"false\"")
					}
				} else {
					if isPipe {
						acc = append(acc, "-lbl", "\\-")
					} else {
						acc = append(acc, "-lbl", "\"\\-\"")
					}
				}
			}

		} else if strings.HasPrefix(str, "#INSD") {

			checkAgainstVocabulary(str, "element", insdtags)
			if doIndex {
				acc = append(acc, "-block", "INSDFeature", "-clr", "-indices")
			} else {
				if isPipe {
					acc = append(acc, "-block", "INSDFeature", "-sep", "|", "-element")
					acc = append(acc, str)
				} else {
					acc = append(acc, "-block", "INSDFeature", "-sep", "\"|\"", "-element")
					ql := fmt.Sprintf("\"%s\"", str)
					acc = append(acc, ql)
				}
			}

		} else if strings.HasPrefix(strings.ToUpper(str), "#INSD") {

			// report capitalization or vocabulary failure
			checkAgainstVocabulary(str, "element", insdtags)

		} else if str == "sub_sequence" {

			// special sub_sequence qualifier shows sequence under feature intervals
			acc = append(acc, "-block", "INSDFeature_intervals")

			acc = append(acc, "-subset", "INSDInterval", "-FR", "INSDInterval_from", "-TO", "INSDInterval_to")
			if isPipe {
				acc = append(acc, "-pfx", "", "-tab", "", "-nucleic", "&SEQ[&FR:&TO]")
			} else {
				acc = append(acc, "-pfx", "\"\"", "-tab", "\"\"", "-nucleic", "\"&SEQ[&FR:&TO]\"")
			}

			acc = append(acc, "-subset", "INSDFeature_intervals")
			if isPipe {
				acc = append(acc, "-deq", "\\t")
			} else {
				acc = append(acc, "-deq", "\"\\t\"")
			}

		} else if str == "feat_location" {

			// special feat_location qualifier shows feature intervals, in 1-based GenBank convention
			acc = append(acc, "-block", "INSDFeature_intervals")

			acc = append(acc, "-subset", "INSDInterval", "-FR", "INSDInterval_from", "-TO", "INSDInterval_to")
			if isPipe {
				acc = append(acc, "-pfx", "", "-tab", "..", "-element", "&FR")
				acc = append(acc, "-pfx", "", "-tab", ",", "-element", "&TO")
			} else {
				acc = append(acc, "-pfx", "\"\"", "-tab", "\"..\"", "-element", "\"&FR\"")
				acc = append(acc, "-pfx", "\"\"", "-tab", "\",\"", "-element", "\"&TO\"")
			}

			acc = append(acc, "-subset", "INSDFeature_intervals")
			if isPipe {
				acc = append(acc, "-deq", "\\t")
			} else {
				acc = append(acc, "-deq", "\"\\t\"")
			}

		} else if str == "feat_intervals" {

			// special feat_intervals qualifier shows feature intervals, decremented to 0-based
			acc = append(acc, "-block", "INSDFeature_intervals")

			acc = append(acc, "-subset", "INSDInterval")
			if isPipe {
				acc = append(acc, "-pfx", "", "-tab", "..", "-dec", "INSDInterval_from")
				acc = append(acc, "-pfx", "", "-tab", ",", "-dec", "INSDInterval_to")
			} else {
				acc = append(acc, "-pfx", "\"\"", "-tab", "\"..\"", "-dec", "\"INSDInterval_from\"")
				acc = append(acc, "-pfx", "\"\"", "-tab", "\",\"", "-dec", "\"INSDInterval_to\"")
			}

			acc = append(acc, "-subset", "INSDFeature_intervals")
			if isPipe {
				acc = append(acc, "-deq", "\\t")
			} else {
				acc = append(acc, "-deq", "\"\\t\"")
			}

		} else if str == "chloroplast" ||
			str == "chromoplast" ||
			str == "cyanelle" ||
			str == "environmental_sample" ||
			str == "focus" ||
			str == "germline" ||
			str == "kinetoplast" ||
			str == "macronuclear" ||
			str == "metagenomic" ||
			str == "mitochondrion" ||
			str == "partial" ||
			str == "proviral" ||
			str == "pseudo" ||
			str == "rearranged" ||
			str == "ribosomal_slippage" ||
			str == "trans_splicing" ||
			str == "transgenic" ||
			str == "virion" {

			acc = append(acc, "-block", "INSDQualifier")

			checkAgainstVocabulary(str, "qualifier", qualifiers)
			if doIndex {
				acc = append(acc, "-if", "INSDQualifier_name", "-equals", str)
				acc = append(acc, "-clr", "-indices", "INSDQualifier_name")
			} else {
				acc = append(acc, "-if", "INSDQualifier_name", "-equals", str)
				acc = append(acc, "-lbl", str)
			}
			if addDash {
				acc = append(acc, "-block", "INSDFeature", "-unless", "INSDQualifier_name", "-equals", str)
				if isPipe {
					acc = append(acc, "-lbl", "\\-")
				} else {
					acc = append(acc, "-lbl", "\"\\-\"")
				}
			}

		} else {

			acc = append(acc, "-block", "INSDQualifier")

			checkAgainstVocabulary(str, "qualifier", qualifiers)
			if len(str) > 2 && str[0] == '%' {
				acc = append(acc, "-if", "INSDQualifier_name", "-equals", str[1:])
				if doIndex {
					if isPipe {
						acc = append(acc, "-clr", "-indices", "%INSDQualifier_value")
					} else {
						acc = append(acc, "-clr", "-indices", "\"%INSDQualifier_value\"")
					}
				} else {
					if isPipe {
						acc = append(acc, "-element", "%INSDQualifier_value")
					} else {
						acc = append(acc, "-element", "\"%INSDQualifier_value\"")
					}
				}
				if addDash {
					acc = append(acc, "-block", "INSDFeature", "-unless", "INSDQualifier_name", "-equals", str[1:])
					if isPipe {
						acc = append(acc, "-lbl", "\\-")
					} else {
						acc = append(acc, "-lbl", "\"\\-\"")
					}
				}
			} else {
				if doIndex {
					acc = append(acc, "-if", "INSDQualifier_name", "-equals", str)
					acc = append(acc, "-clr", "-indices", "INSDQualifier_value")
				} else {
					acc = append(acc, "-if", "INSDQualifier_name", "-equals", str)
					acc = append(acc, "-element", "INSDQualifier_value")
				}
				if addDash {
					acc = append(acc, "-block", "INSDFeature", "-unless", "INSDQualifier_name", "-equals", str)
					if isPipe {
						acc = append(acc, "-lbl", "\\-")
					} else {
						acc = append(acc, "-lbl", "\"\\-\"")
					}
				}
			}
		}
	}

	if doIndex {
		if isPipe {
			acc = append(acc, "-group", "INSDSeq", "-clr", "-lbl", "    </IdxSearchFields>\n")
		} else {
			acc = append(acc, "-group", "INSDSeq", "-clr", "-lbl", "\"    </IdxSearchFields>\\n\"")
		}
	}

	return acc
}

// BIOTHINGS EXTRACTION COMMAND GENERATOR

// processBiopath generates extraction commands for BioThings resources (undocumented)
func processBiopath(args []string, isPipe bool) []string {

	// nquire -get "http://myvariant.info/v1/variant/chr6:g.26093141G>A" \
	//   -fields clinvar.rcv.conditions.identifiers \
	//   -always_list clinvar.rcv.conditions.identifiers |
	// transmute -j2x |
	// xtract -biopath opt clinvar.rcv.conditions.identifiers.omim

	var acc []string

	max := len(args)
	if max < 2 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract -biopath\n")
		os.Exit(1)
	}

	obj := args[0]
	args = args[1:]

	acc = append(acc, "-pattern", obj)

	paths := args[0]

	items := strings.Split(paths, ",")

	for _, path := range items {

		dirs := strings.Split(path, ".")
		max = len(dirs)
		if max < 1 {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient path arguments supplied to xtract -biopath\n")
			os.Exit(1)
		}
		if max > 7 {
			fmt.Fprintf(os.Stderr, "\nERROR: Too many nodes in argument supplied to xtract -biopath\n")
			os.Exit(1)
		}

		str := dirs[max-1]

		acc = append(acc, "-path")
		if isPipe {
			acc = append(acc, path)
			acc = append(acc, "-tab", "\\n")
			acc = append(acc, "-element", str)
		} else {
			acc = append(acc, "\""+path+"\"")
			acc = append(acc, "-tab", "\"\\n\"")
			acc = append(acc, "-element", "\""+str+"\"")
		}
	}

	return acc
}

// HYDRA CITATION MATCHER COMMAND GENERATOR

// processHydra generates extraction commands for NCBI's in-house citation matcher (undocumented)
func processHydra(isPipe bool) []string {

	var acc []string

	// acceptable scores are 0.8 or higher, exact match on "1" rejects low value in scientific notation with minus sign present

	acc = append(acc, "-pattern", "Id")
	acc = append(acc, "-if", "@score", "-equals", "1")
	acc = append(acc, "-or", "@score", "-starts-with", "0.9")
	acc = append(acc, "-or", "@score", "-starts-with", "0.8")
	acc = append(acc, "-element", "Id")

	return acc
}

// ENTREZ2INDEX COMMAND GENERATOR

// processE2Index generates extraction commands to create input for Entrez2Index
func processE2Index(args []string, tform string, isPipe bool) []string {

	var acc []string

	max := len(args)
	if max < 3 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to -e2index\n")
		os.Exit(1)
	}

	year := ""
	patrn := args[0]
	args = args[1:]

	isAllNumeric := func(str string) bool {

		for _, ch := range str {
			if !unicode.IsDigit(ch) &&
				ch != '.' &&
				ch != '+' &&
				ch != '-' &&
				ch != '*' &&
				ch != '/' &&
				ch != ',' &&
				ch != '$' &&
				ch != '#' &&
				ch != '%' &&
				ch != '(' &&
				ch != ')' {
				return false
			}
		}

		return true
	}

	if isAllNumeric(patrn) {
		year = patrn
		patrn = args[0]
		args = args[1:]
	}

	ident := args[0]
	args = args[1:]

	if !isPipe {
		if !deStop {
			acc = append(acc, "-stops")
		}
		if doStem {
			acc = append(acc, "-stems")
		}
	}

	if isPipe {
		acc = append(acc, "-head", "<IdxDocumentSet>", "-tail", "</IdxDocumentSet>")
		acc = append(acc, "-hd", "  <IdxDocument>\\n", "-tl", "  </IdxDocument>")
		acc = append(acc, "-pattern")
		acc = append(acc, patrn)
		if year != "" {
			acc = append(acc, "-if", "PubDate/Year", "-ge", year)
			acc = append(acc, "-or", "PubDate/MedlineDate[1:4]", "-ge", year)
		}
		acc = append(acc, "-pfx", "    <IdxUid>", "-sfx", "</IdxUid>\\n")
		acc = append(acc, "-element")
		acc = append(acc, ident)
		acc = append(acc, "-clr", "-rst", "-tab", "")
		acc = append(acc, "-lbl", "    <IdxSearchFields>\\n")
		acc = append(acc, "-pfx", "      <YEAR>", "-sfx", "</YEAR>\\n")
		acc = append(acc, "-year", "PubDate/*")
		acc = append(acc, "-clr", "-rst", "-tab", "")
		for _, str := range args {
			if str == "ArticleTitle" {
				acc = append(acc, "-article")
			} else {
				acc = append(acc, "-indices")
			}
			acc = append(acc, str)
		}
		if tform != "" {
			acc = append(acc, "-clr", "-rst", "-tab", "\"\"")
			acc = append(acc, "-sep", ",", "-meshcode")
			acc = append(acc, "MeshHeading/DescriptorName@UI,Chemical/NameOfSubstance@UI,SupplMeshName@UI")
		}
		acc = append(acc, "-clr", "-lbl", "    </IdxSearchFields>\\n")
	} else {
		acc = append(acc, "-head", "\"<IdxDocumentSet>\"", "-tail", "\"</IdxDocumentSet>\"")
		acc = append(acc, "-hd", "\"  <IdxDocument>\\n\"", "-tl", "\"  </IdxDocument>\"")
		acc = append(acc, "-pattern")
		ql := fmt.Sprintf("\"%s\"", patrn)
		acc = append(acc, ql)
		if year != "" {
			acc = append(acc, "-if", "PubDate/Year", "-ge", year)
			acc = append(acc, "-or", "PubDate/MedlineDate[1:4]", "-ge", year)
		}
		acc = append(acc, "-pfx", "\"    <IdxUid>\"", "-sfx", "\"</IdxUid>\\n\"")
		acc = append(acc, "-element")
		ql = fmt.Sprintf("\"%s\"", ident)
		acc = append(acc, ql)
		acc = append(acc, "-clr", "-rst", "-tab", "\"\"")
		acc = append(acc, "-lbl", "\"    <IdxSearchFields>\\n\"")
		acc = append(acc, "-pfx", "\"      <YEAR>\"", "-sfx", "\"</YEAR>\\n\"")
		acc = append(acc, "-year", "\"PubDate/*\"")
		acc = append(acc, "-clr", "-rst", "-tab", "\"\"")
		for _, str := range args {
			if str == "ArticleTitle" {
				acc = append(acc, "-article")
			} else {
				acc = append(acc, "-indices")
			}
			ql = fmt.Sprintf("\"%s\"", str)
			acc = append(acc, ql)
		}
		if tform != "" {
			acc = append(acc, "-clr", "-rst", "-tab", "\"\"")
			acc = append(acc, "-sep", "\",\"", "-meshcode")
			acc = append(acc, "\"MeshHeading/DescriptorName@UI,Chemical/NameOfSubstance@UI,SupplMeshName@UI\"")
		}
		acc = append(acc, "-clr", "-lbl", "\"    </IdxSearchFields>\\n\"")
	}

	return acc
}

// CONCURRENT CONSUMER GOROUTINES PARSE AND PROCESS PARTITIONED XML OBJECTS

// StreamBlocks -> SplitPattern => XmlParse => StreamTokens => ProcessQuery -> MergeResults

// processes with single goroutine call defer close(out) so consumer(s) can range over channel
// processes with multiple instances call defer wg.Done(), separate goroutine uses wg.Wait() to delay close(out)

func createConsumers(cmds *Block, parent, hd, tl string, transform map[string]string, histogram map[string]int, inp <-chan eutils.XMLRecord) <-chan eutils.XMLRecord {

	if inp == nil {
		return nil
	}

	out := make(chan eutils.XMLRecord, eutils.ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create consumer channel\n")
		os.Exit(1)
	}

	// xmlConsumer reads partitioned XML from channel and calls parser for processing
	xmlConsumer := func(cmds *Block, parent string, wg *sync.WaitGroup, inp <-chan eutils.XMLRecord, out chan<- eutils.XMLRecord) {

		// report when this consumer has no more records to process
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

			str := processQuery(text[:], parent, idx, hd, tl, transform, histogram, cmds)

			// send even if empty to get all record counts for reordering
			out <- eutils.XMLRecord{Index: idx, Text: str}
		}
	}

	var wg sync.WaitGroup

	// launch multiple consumer goroutines
	for i := 0; i < eutils.NumServe(); i++ {
		wg.Add(1)
		go xmlConsumer(cmds, parent, &wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all consumers are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

func createSelectors(parent, indx string, order map[string]bool, inp <-chan eutils.XMLRecord) <-chan eutils.XMLRecord {

	if parent == "" || indx == "" || order == nil || inp == nil {
		return nil
	}

	find := eutils.ParseIndex(indx)

	out := make(chan eutils.XMLRecord, eutils.ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create selector channel\n")
		os.Exit(1)
	}

	// xmlSelector reads partitioned XML from channel and matches identifiers of records to keep
	xmlSelector := func(wg *sync.WaitGroup, inp <-chan eutils.XMLRecord, out chan<- eutils.XMLRecord) {

		// report when this selector has no more records to process
		defer wg.Done()

		// read partitioned XML from producer channel
		for ext := range inp {

			text := ext.Text

			found := false

			eutils.FindIdentifiers(text[:], parent, find,
				func(id string) {
					id = sortStringByWords(id)
					_, ok := order[id]
					if ok {
						found = true
					}
				})

			if !found {
				// identifier field not found or not in identifier list, send empty placeholder for unshuffler
				out <- eutils.XMLRecord{Index: ext.Index}
				continue
			}

			// send selected record
			out <- eutils.XMLRecord{Index: ext.Index, Text: text}
		}
	}

	var wg sync.WaitGroup

	// launch multiple selector goroutines
	for i := 0; i < eutils.NumServe(); i++ {
		wg.Add(1)
		go xmlSelector(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all selectors are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

// MAIN FUNCTION

// e.g., xtract -pattern PubmedArticle -element MedlineCitation/PMID -block Author -sep " " -element Initials,LastName

func main() {

	// skip past executable name
	args := os.Args[1:]

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: No command-line arguments supplied to xtract\n")
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
	doStem = false
	deStop = true

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

	// flag for indexed input file
	turbo := false

	// debugging
	mpty := false
	idnt := false
	stts := false
	timr := false

	// profiling
	prfl := false

	// repeat the specified extraction 5 times for each -proc from 1 to nCPU
	trial := false

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
			fileName = eutils.GetStringArg(args, "Input file name")
			args = args[1:]

		// input is indexed with <NEXT_RECORD_SIZE> objects
		case "-turbo":
			turbo = true

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
			doStem = true
		case "-stops", "-stop":
			deStop = false

		// allow setting of unicode, script, and mathml flags (undocumented)
		case "-unicode":
			// unicodePolicy = GetStringArg(args, "Unicode argument")
			args = args[1:]
		case "-script":
			// scriptPolicy = GetStringArg(args, "Script argument")
			args = args[1:]
		case "-mathml":
			// mathmlPolicy = GetStringArg(args, "MathML argument")
			args = args[1:]

		case "-flag", "-flags":
			flgs = eutils.GetStringArg(args, "Flags argument")
			args = args[1:]

		// debugging flags
		case "-debug":
			// dbug = true
		case "-empty":
			mpty = true
		case "-ident":
			idnt = true
		case "-stats", "-stat":
			stts = true
		case "-timer":
			timr = true
		case "-profile":
			prfl = true
		case "-trial", "-trials":
			trial = true

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
		doStem = true
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
		UnicodeFix = parseMarkup(unicodePolicy, "-unicode")
		ScriptFix = parseMarkup(scriptPolicy, "-script")
		MathMLFix = parseMarkup(mathmlPolicy, "-mathml")

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

	eutils.SetTunings(numProcs, numServe, serverRatio, chanDepth, farmSize, heapSize, goGc, turbo)

	eutils.SetOptions(doStrict, doMixed, doSelf, deAccent, doASCII, doCompress, doCleanup)

	// -stats prints number of CPUs and performance tuning values if no other arguments (undocumented)
	if stts && len(args) < 1 {

		eutils.PrintStats()

		return
	}

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract\n")
		os.Exit(1)
	}

	// DOCUMENTATION COMMANDS

	inSwitch = true

	switch args[0] {
	case "-version":
		fmt.Printf("%s\n", eutils.EDirectVersion)
	case "-help", "help":
		eutils.PrintHelp("xtract", "xtract-help.txt")
	case "-examples", "-example":
		eutils.PrintHelp("xtract", "xtract-examples.txt")
	case "-extras", "-extra", "-advanced":
		fmt.Printf("Please run rchive -help for local record indexing information\n")
	case "-internal", "-internals":
		eutils.PrintHelp("xtract", "xtract-internal.txt")
	case "-keys":
		eutils.PrintHelp("xtract", "xtract-keys.txt")
	case "-unix":
		eutils.PrintHelp("xtract", "xtract-unix.txt")
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

	// NAME OF OUTPUT STRING TRANSFORMATION FILE

	tform := ""
	transform := make(map[string]string)

	populateTx := func(tf string) {

		inFile, err := os.Open(tf)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Unable to open transformation file %s\n", err.Error())
			os.Exit(1)
		}
		defer inFile.Close()

		scanr := bufio.NewScanner(inFile)

		// populate transformation map for -translate (and -matrix) output
		for scanr.Scan() {

			line := scanr.Text()
			frst, scnd := eutils.SplitInTwoLeft(line, "\t")

			transform[frst] = scnd
		}
	}

	if len(args) > 2 && args[0] == "-transform" {
		tform = args[1]
		args = args[2:]
		if tform != "" {
			populateTx(tform)
		}
	}

	// CREATE XML BLOCK READER FROM STDIN OR FILE

	const FIRST_BUFF_SIZE = 4096

	getFirstBlock := func() string {

		buffer := make([]byte, FIRST_BUFF_SIZE)
		n, err := in.Read(buffer)
		if err != nil && err != io.EOF {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to read first block: %s\n", err.Error())
			// os.Exit(1)
		}
		bufr := buffer[:n]
		return string(bufr)
	}

	first := getFirstBlock()

	mlt := io.MultiReader(strings.NewReader(first), in)

	isJsn := false
	isAsn := false
	isGbf := false
	matched := 0

	// auto-detect XML, JSON, or ASN.1 format
	if first != "" {
		posJ1 := strings.Index(first, "{")
		posJ2 := strings.Index(first, "\":")
		if posJ1 >= 0 && posJ2 >= 0 && posJ1 < posJ2 {
			isJsn = true
			matched++
		} else {
			posJ1 = FIRST_BUFF_SIZE
			posJ2 = FIRST_BUFF_SIZE
		}
		posA1 := strings.Index(first, "::=")
		posA2 := strings.Index(first, "{")
		if posA1 >= 0 && posA2 >= 0 && posA1 < posA2 {
			isAsn = true
			matched++
		} else {
			posA1 = FIRST_BUFF_SIZE
			posA2 = FIRST_BUFF_SIZE
		}
		posG1 := strings.Index(first, "LOCUS")
		posG2 := strings.Index(first, "DEFINITION")
		if posG1 >= 0 && posG2 >= 0 && posG1 < posG2 {
			isGbf = true
			matched++
		} else {
			posG1 = FIRST_BUFF_SIZE
			posG2 = FIRST_BUFF_SIZE
		}
		posX1 := strings.Index(first, "<")
		posX2 := strings.Index(first, ">")
		if posX1 >= 0 && posX2 >= 0 && posX1 < posX2 {
			matched++
		} else {
			posX1 = FIRST_BUFF_SIZE
			posX2 = FIRST_BUFF_SIZE
		}
		if matched > 1 {
			if posX1 < posJ1 && posX1 < posA1 && posX1 < posG1 {
				isJsn = false
				isAsn = false
				isGbf = false
			} else if posJ1 < posA1 && posJ1 < posG1 {
				isAsn = false
				isGbf = false
			} else if posA1 < posJ1 && posA1 < posG1 {
				isJsn = false
				isGbf = false
			} else if posG1 < posJ1 && posG1 < posA1 {
				isJsn = false
				isAsn = false
			}
		}
	}

	if isJsn {
		jrdr := eutils.JSONConverter(mlt, "root", "opt", "element")
		mlt = eutils.ChanToReader(jrdr)
	} else if isAsn {
		ardr := eutils.ASN1Converter(mlt, "", "")
		mlt = eutils.ChanToReader(ardr)
	} else if isGbf {
		grdr := eutils.GenBankConverter(mlt)
		mlt = eutils.ChanToReader(grdr)
	}

	rdr := eutils.CreateXMLStreamer(mlt)
	if rdr == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
		os.Exit(1)
	}

	// SEQUENCE RECORD EXTRACTION COMMAND GENERATOR

	// -insd simplifies extraction of INSDSeq qualifiers
	if args[0] == "-insd" || args[0] == "-insd-" || args[0] == "-insd-idx" {

		addDash := true
		doIndex := false
		// -insd- variant suppresses use of dash as placeholder for missing qualifiers (undocumented)
		if args[0] == "-insd-" {
			addDash = false
		}
		// -insd-idx variant creates word index using -indices command (undocumented)
		if args[0] == "-insd-idx" {
			doIndex = true
			addDash = false
		}

		args = args[1:]

		insd := processINSD(args, isPipe || usingFile, addDash, doIndex)

		if !isPipe && !usingFile {
			// no piped input, so write output instructions
			fmt.Printf("xtract")
			for _, str := range insd {
				fmt.Printf(" %s", str)
			}
			fmt.Printf("\n")
			return
		}

		// data in pipe, so replace arguments, execute dynamically
		args = insd
	}

	// CITATION MATCHER EXTRACTION COMMAND GENERATOR

	// -hydra filters HydraResponse output by relevance score (undocumented)
	if args[0] == "-hydra" {

		hydra := processHydra(isPipe || usingFile)

		if !isPipe && !usingFile {
			// no piped input, so write output instructions
			fmt.Printf("xtract")
			for _, str := range hydra {
				fmt.Printf(" %s", str)
			}
			fmt.Printf("\n")
			return
		}

		// data in pipe, so replace arguments, execute dynamically
		args = hydra
	}

	// BIOTHINGS EXTRACTION COMMAND GENERATOR

	// -biopath takes a parent object and a dotted exploration path for BioThings resources (undocumented)
	if args[0] == "-biopath" {

		args = args[1:]

		biopath := processBiopath(args, isPipe || usingFile)

		if !isPipe && !usingFile {
			// no piped input, so write output instructions
			fmt.Printf("xtract")
			for _, str := range biopath {
				fmt.Printf(" %s", str)
			}
			fmt.Printf("\n")
			return
		}

		// data in pipe, so replace arguments, execute dynamically
		args = biopath
	}

	// ENTREZ2INDEX COMMAND GENERATOR

	// -e2index shortcut for experimental indexing code (documented in rchive.go)
	if args[0] == "-e2index" {

		// e.g., xtract -transform "$EDIRECT_MESH_TREE" -e2index

		args = args[1:]

		if len(args) == 0 {
			// if no arguments, use default values
			args = []string{"PubmedArticle", "MedlineCitation/PMID", "ArticleTitle", "ArticleTitle,Abstract/AbstractText"}
		}

		// environment variable can override garbage collector (undocumented)
		gcEnv := os.Getenv("EDIRECT_INDEX_GOGC")
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
		svEnv := os.Getenv("EDIRECT_INDEX_SERV")
		if svEnv != "" {
			val, err := strconv.Atoi(svEnv)
			if err == nil {
				if val >= 1 && val <= 128 {
					numServe = val
				} else {
					numServe = 1
				}
				eutils.SetTunings(numProcs, numServe, serverRatio, chanDepth, farmSize, heapSize, goGc, turbo)
			}
		}

		res := processE2Index(args, tform, isPipe || usingFile)

		if !isPipe && !usingFile {
			// no piped input, so write output instructions
			fmt.Printf("xtract")
			if tform != "" {
				fmt.Printf(" -transform %s", tform)
			}
			for _, str := range res {
				fmt.Printf(" %s", str)
			}
			fmt.Printf("\n")
			return
		}

		// data in pipe, so replace arguments, execute dynamically
		args = res
	}

	// CONFIRM INPUT DATA AVAILABILITY AFTER RUNNING COMMAND GENERATORS

	if fileName == "" && runtime.GOOS != "windows" {

		fromStdin := bool((fi.Mode() & os.ModeCharDevice) == 0)
		if !isPipe || !fromStdin {
			mode := fi.Mode().String()
			fmt.Fprintf(os.Stderr, "\nERROR: No data supplied to xtract from stdin or file, mode is '%s'\n", mode)
			os.Exit(1)
		}
	}

	if !usingFile && !isPipe {

		fmt.Fprintf(os.Stderr, "\nERROR: No XML input data supplied to xtract\n")
		os.Exit(1)
	}

	// XML VALIDATION

	nextArg := func() (string, bool) {

		if len(args) < 1 {
			return "", false
		}

		// remove next token from slice
		nxt := args[0]
		args = args[1:]

		return nxt, true
	}

	if args[0] == "-verify" || args[0] == "-validate" {

		// skip past command name
		args = args[1:]

		find := ""
		html := false

		// look for optional arguments
		for {
			arg, ok := nextArg()
			if !ok {
				break
			}

			switch arg {
			case "-find":
				// override set wrapper
				find, ok = nextArg()
			case "-html":
				html = true
			}
		}

		recordCount = eutils.ValidateXML(rdr, find, html)

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

	// MISCELLANEOUS TIMING COMMANDS

	if args[0] == "-chunk" {

		for str := range rdr {
			recordCount++
			byteCount += len(str)
		}

		printDuration("blocks")

		return
	}

	if args[0] == "-split" {

		if len(args) > 1 {
			if args[1] == "-pattern" {
				// skip past -split if followed by -pattern
				args = args[1:]
			}
		}
		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -split command\n")
			os.Exit(1)
		}
		pat := args[1]

		eutils.PartitionPattern(pat, "", turbo, rdr,
			func(str string) {
				recordCount++
				byteCount += len(str)
			})

		printDuration("patterns")

		return
	}

	if args[0] == "-token" {

		eutils.StreamTokens(rdr,
			func(tkn eutils.XMLToken) {
				recordCount++
				byteCount += len(tkn.Name) + len(tkn.Attr)
			})

		printDuration("tokens")

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
			// allow splitting of -head argument, keep appending until next command (undocumented)
			ofs, nxt := 0, args[2:]
			for {
				if len(nxt) < 1 {
					break
				}
				tmp := nxt[0]
				if strings.HasPrefix(tmp, "-") {
					break
				}
				ofs++
				txt := eutils.ConvertSlash(tmp)
				if head != "" && !strings.HasSuffix(head, "\t") {
					head += "\t"
				}
				head += txt
				nxt = nxt[1:]
			}
			if ofs > 0 {
				args = args[ofs:]
			}
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
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract\n")
			os.Exit(1)
		}
	}

	// INDEXED XML FILE PREPARATION

	// cat carotene.xml | xtract -timer -index -pattern PubmedArticle > carindex.txt
	// xtract -timer -turbo -input carindex.txt -pattern PubmedArticle -element LastName
	if args[0] == "-index" {

		if len(args) > 1 {
			if args[1] == "-pattern" {
				// skip past -index if followed by -pattern
				args = args[1:]
			}
		}
		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -index command\n")
			os.Exit(1)
		}
		pat := args[1]

		retlength := len("\n")

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		eutils.PartitionPattern(pat, "", false, rdr,
			func(str string) {
				recordCount++
				nxt := len(str)
				byteCount += nxt
				newln := false

				if !strings.HasSuffix(str, "\n") {
					nxt += retlength
					newln = true
				}

				os.Stdout.WriteString("<NEXT_RECORD_SIZE>")
				val := strconv.Itoa(nxt)
				os.Stdout.WriteString(val)
				os.Stdout.WriteString("</NEXT_RECORD_SIZE>\n")

				os.Stdout.WriteString(str)
				if newln {
					os.Stdout.WriteString("\n")
				}
			})

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

	// ENSURE PRESENCE OF PATTERN ARGUMENT

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract\n")
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

	// READ FILE OF IDENTIFIERS AND CONCURRENTLY EXTRACT SELECTED RECORDS

	// -pattern record_name -select parent/element@attribute^version -in file_of_identifiers
	if len(args) == 6 && args[2] == "-select" && (args[4] == "-in" || args[4] == "-retaining") {

		indx := args[3]
		unqe := args[5]

		// read file of identifiers to use for filtering
		fl, err := os.Open(unqe)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open identifier file '%s'\n", unqe)
			os.Exit(1)
		}

		// create map that records each UID
		order := make(map[string]bool)

		scanr := bufio.NewScanner(fl)

		// read lines of identifiers
		for scanr.Scan() {

			line := scanr.Text()
			id, _ := eutils.SplitInTwoLeft(line, "\t")

			id = sortStringByWords(id)

			// add identifier to map
			order[id] = true
		}

		fl.Close()

		xmlq := eutils.CreateXMLProducer(topPattern, star, false, rdr)
		fchq := createSelectors(topPattern, indx, order, xmlq)
		unsq := eutils.CreateXMLUnshuffler(fchq)

		if xmlq == nil || fchq == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create selector\n")
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

	// READ FILE OF IDENTIFIERS AND EXCLUDE SELECTED RECORDS

	// -pattern record_name -exclude element -excluding file_of_identifiers (undocumented)
	if len(args) == 6 && args[2] == "-select" && args[4] == "-excluding" {

		indx := args[3]
		unqe := args[5]

		// read file of identifiers to use for filtering
		fl, err := os.Open(unqe)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open identifier file '%s'\n", unqe)
			os.Exit(1)
		}

		// create map that records each UID
		order := make(map[string]bool)

		scanr := bufio.NewScanner(fl)

		// read lines of identifiers
		for scanr.Scan() {

			line := scanr.Text()
			id, _ := eutils.SplitInTwoLeft(line, "\t")
			id = strings.ToLower(id)

			// add identifier to map
			order[id] = true
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
				if id != "" {
					id = strings.ToLower(id)
					_, ok := order[id]
					if ok {
						// in exclusion list, skip
						return
					}
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

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// READ ORDERED FILE OF IDENTIFIERS AND XML STRINGS, APPEND XML JUST INSIDE CLOSING TAG OF APPROPRIATE RECORD

	// -pattern record_name -select element -appending file_of_identifiers_and_metadata (undocumented)
	if len(args) == 6 && args[2] == "-select" && args[4] == "-appending" {

		indx := args[3]
		apnd := args[5]

		fl, err := os.Open(apnd)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open transformation file '%s'\n", apnd)
			os.Exit(1)
		}

		scanr := bufio.NewScanner(fl)

		find := eutils.ParseIndex(indx)

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		rgt := "</" + topPattern + ">"

		eutils.PartitionPattern(topPattern, star, false, rdr,
			func(str string) {
				recordCount++

				id := eutils.FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}
				id = strings.ToLower(id)

				for scanr.Scan() {

					line := scanr.Text()
					frst, scnd := eutils.SplitInTwoLeft(line, "\t")
					frst = strings.ToLower(frst)

					if id != frst {
						return
					}
					if !strings.HasSuffix(str, rgt) {
						return
					}

					lft := strings.TrimSuffix(str, rgt)
					str = lft + "  " + scnd + "\n" + rgt

					if hd != "" {
						os.Stdout.WriteString(hd)
						os.Stdout.WriteString("\n")
					}

					os.Stdout.WriteString(str[:])
					os.Stdout.WriteString("\n")

					if tl != "" {
						os.Stdout.WriteString(tl)
						os.Stdout.WriteString("\n")
					}

					break
				}
			})

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		fl.Close()

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// SORT XML RECORDS BY IDENTIFIER

	// -pattern record_name -sort parent/element@attribute^version
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
		// sort fields in alphabetical or numeric order
		sort.Slice(keys, func(i, j int) bool {
			// numeric sort on strings checks lengths first
			if eutils.IsAllDigits(keys[i]) && eutils.IsAllDigits(keys[j]) {
				lni := len(keys[i])
				lnj := len(keys[j])
				// shorter string is numerically less, assuming no leading zeros
				if lni < lnj {
					return true
				}
				if lni > lnj {
					return false
				}
			}
			// same length or non-numeric, can now do string comparison on contents
			return keys[i] < keys[j]
		})

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

	// SPLIT FILE BY BY RECORD COUNT

	// split XML record into subfiles by count
	if len(args) == 8 && args[2] == "-split" && args[4] == "-prefix" && args[6] == "-suffix" {

		// e.g., -head "<IdxDocumentSet>" -tail "</IdxDocumentSet>" -pattern IdxDocument -split 250000 -prefix "biocon" -suffix "e2x"
		count := 0
		fnum := 0
		var (
			fl  *os.File
			err error
		)
		chunk, err := strconv.Atoi(args[3])
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return
		}
		prefix := args[5]
		suffix := args[7]

		eutils.PartitionPattern(topPattern, star, false, rdr,
			func(str string) {
				recordCount++

				if count >= chunk {
					if tail != "" {
						fl.WriteString(tail)
						fl.WriteString("\n")
					}
					fl.Close()
					count = 0
				}
				if count == 0 {
					fpath := fmt.Sprintf("%s%03d.%s", prefix, fnum, suffix)
					fl, err = os.Create(fpath)
					if err != nil {
						fmt.Fprintf(os.Stderr, "%s\n", err.Error())
						return
					}
					os.Stderr.WriteString(fpath + "\n")
					fnum++
					if head != "" {
						fl.WriteString(head)
						fl.WriteString("\n")
					}
				}
				count++

				fl.WriteString(str[:])
				fl.WriteString("\n")
			})

		if count >= chunk {
			if tail != "" {
				fl.WriteString(tail)
				fl.WriteString("\n")
			}
			fl.Close()
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// PARSE AND VALIDATE EXTRACTION ARGUMENTS

	// parse nested exploration instruction from command-line arguments
	cmds := parseArguments(args, topPattern)
	if cmds == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Problem parsing command-line arguments\n")
		os.Exit(1)
	}

	// GLOBAL MAP FOR SORT-UNIQ-COUNT HISTOGRAM ARGUMENT

	histogram := make(map[string]int)

	// PERFORMANCE TIMING COMMAND

	// -stats with an extraction command prints XML size and processing time for each record
	if stts {

		legend := "REC\tOFST\tSIZE\tTIME"

		rec := 0

		eutils.PartitionPattern(topPattern, star, turbo, rdr,
			func(str string) {
				rec++
				beginTime := time.Now()
				processQuery(str[:], parent, rec, hd, tl, transform, histogram, cmds)
				endTime := time.Now()
				duration := endTime.Sub(beginTime)
				micro := int(float64(duration.Nanoseconds()) / 1e3)
				if legend != "" {
					fmt.Printf("%s\n", legend)
					legend = ""
				}
				fmt.Printf("%d\t%d\t%d\n", rec, len(str), micro)
			})

		return
	}

	// PERFORMANCE OPTIMIZATION FUNCTION

	// -trial -input fileName runs the specified extraction for each -proc from 1 to nCPU
	if trial && fileName != "" {

		legend := "CPU\tRATE\tDEV"

		for numServ := 1; numServ <= ncpu; numServ++ {

			numServe = numServ

			eutils.SetTunings(numProcs, numServe, serverRatio, chanDepth, farmSize, heapSize, goGc, turbo)

			runtime.GOMAXPROCS(numServ)

			sum := 0
			count := 0
			mean := 0.0
			m2 := 0.0

			// calculate mean and standard deviation of processing rate
			for trials := 0; trials < 5; trials++ {

				inFile, err := os.Open(fileName)
				if err != nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
					os.Exit(1)
				}

				trdr := eutils.CreateXMLStreamer(inFile)
				if trdr == nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to read input file\n")
					os.Exit(1)
				}

				xmlq := eutils.CreateXMLProducer(topPattern, star, turbo, trdr)
				tblq := createConsumers(cmds, parent, hd, tl, transform, histogram, xmlq)

				if xmlq == nil || tblq == nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to create servers\n")
					os.Exit(1)
				}

				begTime := time.Now()
				recordCount = 0

				for range tblq {
					recordCount++
					runtime.Gosched()
				}

				inFile.Close()

				debug.FreeOSMemory()

				endTime := time.Now()
				expended := endTime.Sub(begTime)
				secs := float64(expended.Nanoseconds()) / 1e9

				if secs >= 0.000001 && recordCount > 0 {
					speed := int(float64(recordCount) / secs)
					sum += speed
					count++
					x := float64(speed)
					delta := x - mean
					mean += delta / float64(count)
					m2 += delta * (x - mean)
				}
			}

			if legend != "" {
				fmt.Printf("%s\n", legend)
				legend = ""
			}
			if count > 1 {
				vrc := m2 / float64(count-1)
				dev := int(math.Sqrt(vrc))
				fmt.Printf("%d\t%d\t%d\n", numServ, sum/count, dev)
			}
		}

		return
	}

	// PROCESS SINGLE SELECTED RECORD IF -pattern ARGUMENT IS IMMEDIATELY FOLLOWED BY -position COMMAND

	posn := ""
	if cmds.Visit == topPat {
		if cmds.Position == "outer" ||
			cmds.Position == "inner" ||
			cmds.Position == "even" ||
			cmds.Position == "odd" ||
			cmds.Position == "all" {
			// filter by record position when draining unshuffler channel
			posn = cmds.Position
			cmds.Position = ""
		}
	}

	if cmds.Visit == topPat && cmds.Position != "" && cmds.Position != "select" {

		qry := ""
		idx := 0
		rec := 0

		if cmds.Position == "first" {

			eutils.PartitionPattern(topPattern, star, turbo, rdr,
				func(str string) {
					rec++
					if rec == 1 {
						qry = str
						idx = rec
					}
				})

		} else if cmds.Position == "last" {

			eutils.PartitionPattern(topPattern, star, turbo, rdr,
				func(str string) {
					qry = str
					idx = rec
				})

		} else {

			// use numeric position
			number, err := strconv.Atoi(cmds.Position)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized position '%s'\n", cmds.Position)
				os.Exit(1)
			}

			eutils.PartitionPattern(topPattern, star, turbo, rdr,
				func(str string) {
					rec++
					if rec == number {
						qry = str
						idx = rec
					}
				})
		}

		if qry == "" {
			return
		}

		// clear position on top node to prevent condition test failure
		cmds.Position = ""

		// process single selected record
		res := processQuery(qry[:], parent, idx, hd, tl, transform, histogram, cmds)

		if res != "" {
			fmt.Printf("%s", res)
		}

		return
	}

	// LAUNCH PRODUCER, CONSUMER, AND UNSHUFFLER GOROUTINES

	// launch producer goroutine to partition XML by pattern
	xmlq := eutils.CreateXMLProducer(topPattern, star, turbo, rdr)

	// launch consumer goroutines to parse and explore partitioned XML objects
	tblq := createConsumers(cmds, parent, hd, tl, transform, histogram, xmlq)

	// launch unshuffler goroutine to restore order of results
	unsq := eutils.CreateXMLUnshuffler(tblq)

	if xmlq == nil || tblq == nil || unsq == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create servers\n")
		os.Exit(1)
	}

	// PERFORMANCE SUMMARY

	/*
		if dbug {

			// drain results, but suppress extraction output
			for ext := range unsq {
				byteCount += len(ext.Text)
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

	// DRAIN OUTPUT CHANNEL TO EXECUTE EXTRACTION COMMANDS, RESTORE OUTPUT ORDER WITH HEAP

	var buffer strings.Builder
	count := 0
	okay := false
	lastTime := time.Now()

	wrtr := bufio.NewWriter(os.Stdout)

	// printResult prints output for current pattern, handles -empty and -ident flags, and periodically flushes buffer
	printResult := func(curr eutils.XMLRecord) {

		str := curr.Text

		if mpty {

			if str == "" {

				okay = true

				idx := curr.Index
				val := strconv.Itoa(idx)
				buffer.WriteString(val[:])
				buffer.WriteString("\n")

				count++
			}

		} else if str != "" {

			okay = true

			if idnt {
				idx := curr.Index
				val := strconv.Itoa(idx)
				buffer.WriteString(val[:])
				buffer.WriteString("\t")
			}

			// save output to byte buffer
			buffer.WriteString(str[:])

			count++
		}

		thisTime := time.Now()
		duration := thisTime.Sub(lastTime)
		milliSeconds := duration.Milliseconds()

		if count > 1000 || milliSeconds > 4999 {
			count = 0
			lastTime = thisTime
			txt := buffer.String()
			if txt != "" {
				// print current buffer
				wrtr.WriteString(txt[:])
			}
			buffer.Reset()
		}
	}

	if head != "" {
		buffer.WriteString(head[:])
		buffer.WriteString("\n")
	}

	// drain unshuffler channel

	if posn == "outer" {

		// print only first and last records
		var beg *eutils.XMLRecord
		var end *eutils.XMLRecord

		for curr := range unsq {

			if beg == nil {
				beg = &eutils.XMLRecord{Index: curr.Index, Ident: curr.Ident, Text: curr.Text}
			} else {
				end = &eutils.XMLRecord{Index: curr.Index, Ident: curr.Ident, Text: curr.Text}
			}

			recordCount++
		}

		if beg != nil {
			printResult(*beg)
		}
		if end != nil {
			printResult(*end)
		}

	} else if posn == "inner" {

		// print all but first and last records
		var prev *eutils.XMLRecord
		var next *eutils.XMLRecord
		first := true

		for curr := range unsq {

			if first {
				first = false
			} else {
				prev = next
				next = &eutils.XMLRecord{Index: curr.Index, Ident: curr.Ident, Text: curr.Text}
			}

			if prev != nil {
				printResult(*prev)
			}

			recordCount++
		}

	} else if posn == "even" {

		even := false

		for curr := range unsq {

			if even {
				printResult(curr)
			}
			even = !even

			recordCount++
		}

	} else if posn == "odd" {

		odd := true

		for curr := range unsq {

			if odd {
				printResult(curr)
			}
			odd = !odd

			recordCount++
		}

	} else {

		// default or -position all
		for curr := range unsq {

			// send result to output
			printResult(curr)

			recordCount++
			runtime.Gosched()
		}
	}

	if tail != "" {
		buffer.WriteString(tail[:])
		buffer.WriteString("\n")
	}

	// do not print head or tail if no extraction output
	if okay {
		txt := buffer.String()
		if txt != "" {
			// print final buffer
			wrtr.WriteString(txt[:])
		}
	}
	buffer.Reset()

	wrtr.Flush()

	// print -histogram results, if populated
	var keys []string
	for ky := range histogram {
		keys = append(keys, ky)
	}
	if len(keys) > 0 {
		// sort fields in alphabetical or numeric order
		sort.Slice(keys, func(i, j int) bool {
			// numeric sort on strings checks lengths first
			if eutils.IsAllDigits(keys[i]) && eutils.IsAllDigits(keys[j]) {
				lni := len(keys[i])
				lnj := len(keys[j])
				// shorter string is numerically less, assuming no leading zeros
				if lni < lnj {
					return true
				}
				if lni > lnj {
					return false
				}
			}
			// same length or non-numeric, can now do string comparison on contents
			return keys[i] < keys[j]
		})

		for _, str := range keys {

			count := histogram[str]
			val := strconv.Itoa(count)
			os.Stdout.WriteString(val)
			os.Stdout.WriteString("\t")
			os.Stdout.WriteString(str)
			os.Stdout.WriteString("\n")
		}
	}

	// force garbage collection and return memory before calculating processing rate
	debug.FreeOSMemory()

	if timr {
		printDuration("records")
	}
}
