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
// File Name:  utils.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"bufio"
	"fmt"
	"github.com/klauspost/cpuid"
	"github.com/pbnjay/memory"
	"os"
	"path/filepath"
	"runtime"
	"runtime/debug"
	"strconv"
	"time"
)

// EDirectVersion is the current EDirect release number
const EDirectVersion = "16.2"

// ASCII terminal color constants
const (
	RED  = "\033[31m"
	BLUE = "\033[34m"
	BOLD = "\033[1m"
	FLIP = "\033[7m"
	INIT = "\033[0m"
	LOUD = INIT + RED + BOLD
	INVT = LOUD + FLIP
)

// parser character type lookup tables
var (
	inBlank     [256]bool
	inFirst     [256]bool
	inElement   [256]bool
	inLower     [256]bool
	inContent   [256]bool
	inAsnTag    [256]bool
	inAsnString [256]bool
	inAsnBits   [256]bool
)

// program execution timer
var (
	startTime time.Time
)

// performance tuning variables
var (
	chanDepth   int
	farmSize    int
	heapSize    int
	numServe    int
	goGc        int
	nCPU        int
	numProcs    int
	serverRatio int
)

// reading and cleaning options
var (
	doStrict   bool
	doMixed    bool
	doSelf     bool
	deAccent   bool
	doASCII    bool
	doCompress bool
	doCleanup  bool
)

// additional options
var (
	doUnicode bool
	doScript  bool
	doMathML  bool
)

// derived policy variables
var (
	allowEmbed  bool
	contentMods bool
	countLines  bool
)

var (
	unicodeFix = NOMARKUP
	scriptFix  = NOMARKUP
	mathMLFix  = NOMARKUP
)

// SetTunings sets performance parameters
func SetTunings(nmProcs, nmServe, svRatio, chnDepth, frmSize, hepSize, gogc int, turbo bool) {

	// set default values
	if frmSize < 1 {
		frmSize = 64
	}

	if hepSize < 8 || hepSize > 64 {
		hepSize = 16
	}

	if gogc < 50 || gogc > 1000 {
		gogc = 600
	}

	farmSize = frmSize
	heapSize = hepSize
	goGc = gogc

	// calculate number of simultaneous threads for multiplexed goroutines

	nCPU = runtime.NumCPU()
	if nCPU < 1 {
		nCPU = 1
	}

	// reality checks on number of processors to use
	// performance degrades if capacity is above maximum number of partitions per second (context switching?)
	// best performance varies slightly among PubmedArticle, gene DocumentSummary, and INSDSeq sequence records
	if nmProcs < 1 {
		if turbo {
			// best performance measurement with size-indexed data is obtained when 12 hyperthreads are assigned
			nmProcs = 12
			threads := nCPU
			if threads > 4 && threads < 12 {
				nmProcs = threads
			}
		} else {
			// best performance measurement with Boyer-Moore-Horspool when 6 to 8 processors are assigned
			nmProcs = 8
			if cpuid.CPU.ThreadsPerCore > 1 {
				cores := nCPU / cpuid.CPU.ThreadsPerCore
				if cores > 4 && cores < 8 {
					nmProcs = cores
				}
			}
		}
	}

	if nmProcs > nCPU {
		nmProcs = nCPU
	}

	numProcs = nmProcs

	// allow simultaneous threads for multiplexed goroutines
	runtime.GOMAXPROCS(numProcs)

	// adjust garbage collection target percentage
	debug.SetGCPercent(goGc)

	if svRatio < 1 || svRatio > 32 {
		svRatio = 4
	}

	serverRatio = svRatio

	if nmServe > 0 {
		serverRatio = nmServe / numProcs
	} else {
		nmServe = numProcs * serverRatio
	}

	if nmServe > 128 {
		nmServe = 128
	} else if nmServe < 1 {
		nmServe = numProcs
	}

	numServe = nmServe

	// number of channels usually equals number of servers
	if chnDepth < nCPU || chnDepth > 128 {
		chnDepth = numServe
	}

	chanDepth = chnDepth
}

// SetOptions sets processing options
func SetOptions(strict, mixed, self, accent, ascii, compress, cleanup bool) {

	doStrict = strict
	doMixed = mixed

	doSelf = self

	deAccent = accent
	doASCII = ascii

	doCompress = compress
	doCleanup = cleanup

	countLines = false

	// set dependent flags
	countLines = doMixed
	allowEmbed = doStrict || doMixed
	contentMods = allowEmbed || doCompress || doUnicode || doScript || doMathML || deAccent || doASCII
}

// ChanDepth returns the communication channel depth
func ChanDepth() int {

	return chanDepth
}

// NumServe returns the number of concurrent servers
func NumServe() int {

	return numServe
}

// GetTunings returns performance parameter values
func GetTunings() (nmProcs, nmServe, svRatio, chnDepth, frmSize, hepSize, gogc int) {

	return numProcs, numServe, serverRatio, chanDepth, farmSize, heapSize, goGc
}

// GetOptions returns processing option values
func GetOptions() (strict, mixed, self, accent, ascii, compress, cleanup bool) {

	return doStrict, doMixed, doSelf, deAccent, doASCII, doCompress, doCleanup
}

// GetNumericArg returns an integer argument, reporting an error if no remaining arguments
func GetNumericArg(args []string, name string, zer, min, max int) int {

	if len(args) < 2 {
		fmt.Fprintf(os.Stderr, "\nERROR: %s is missing\n", name)
		os.Exit(1)
	}
	value, err := strconv.Atoi(args[1])
	if err != nil {
		fmt.Fprintf(os.Stderr, "\nERROR: %s (%s) is not an integer\n", name, args[1])
		os.Exit(1)
	}

	// special case for argument value of 0
	if value < 1 {
		return zer
	}
	// limit value to between specified minimum and maximum
	if value < min && min > 0 {
		return min
	}
	if value > max && max > 0 {
		return max
	}
	return value
}

// GetStringArg returns a string argument, reporting an error if no remaining arguments
func GetStringArg(args []string, name string) string {

	if len(args) < 2 {
		fmt.Fprintf(os.Stderr, "\nERROR: %s is missing\n", name)
		os.Exit(1)
	}
	return args[1]
}

// PrintDuration prints processing rate and program duration
func PrintDuration(name string, recordCount, byteCount int) {

	stopTime := time.Now()
	duration := stopTime.Sub(startTime)
	seconds := float64(duration.Nanoseconds()) / 1e9

	prec := 3
	if seconds >= 100 {
		prec = 1
	} else if seconds >= 10 {
		prec = 2
	}

	if recordCount >= 1000000 {
		throughput := float64(recordCount/100000) / 10.0
		fmt.Fprintf(os.Stderr, "\nProcessed %.1f million %s in %.*f seconds", throughput, name, prec, seconds)
	} else if recordCount > 0 {
		fmt.Fprintf(os.Stderr, "\nProcessed %d %s in %.*f seconds", recordCount, name, prec, seconds)
	} else {
		fmt.Fprintf(os.Stderr, "\nProcessing completed in %.*f seconds", prec, seconds)
	}

	if seconds >= 0.001 && recordCount > 0 {
		rate := int(float64(recordCount) / seconds)
		if rate >= 1000000 {
			fmt.Fprintf(os.Stderr, " (%d million %s/second", rate/1000000, name)
		} else {
			fmt.Fprintf(os.Stderr, " (%d %s/second", rate, name)
		}
		if byteCount > 0 {
			rate := int(float64(byteCount) / seconds)
			if rate >= 1000000 {
				fmt.Fprintf(os.Stderr, ", %d megabytes/second", rate/1000000)
			} else if rate >= 1000 {
				fmt.Fprintf(os.Stderr, ", %d kilobytes/second", rate/1000)
			} else {
				fmt.Fprintf(os.Stderr, ", %d bytes/second", rate)
			}
		}
		fmt.Fprintf(os.Stderr, ")")
	}

	fmt.Fprintf(os.Stderr, "\n\n")
}

// PrintMemory is adapted from PrintMemUsage in: https://golangcode.com/print-the-current-memory-usage/
func PrintMemory() {

	var m runtime.MemStats
	runtime.ReadMemStats(&m)

	bToMb := func(b uint64) uint64 {
		return b / 1024 / 1024
	}

	// For info on each, see: https://golang.org/pkg/runtime/#MemStats
	fmt.Fprintf(os.Stderr, "Alloc = %v MiB", bToMb(m.Alloc))
	fmt.Fprintf(os.Stderr, "\tTotalAlloc = %v MiB", bToMb(m.TotalAlloc))
	fmt.Fprintf(os.Stderr, "\tSys = %v MiB", bToMb(m.Sys))
	fmt.Fprintf(os.Stderr, "\tNumGC = %v\n", m.NumGC)
}

// PrintStats prints performance tuning parameters
func PrintStats() {

	fmt.Fprintf(os.Stderr, "Thrd %d\n", nCPU)
	if cpuid.CPU.ThreadsPerCore > 0 {
		fmt.Fprintf(os.Stderr, "Core %d\n", nCPU/cpuid.CPU.ThreadsPerCore)
	}
	if cpuid.CPU.LogicalCores > 0 {
		fmt.Fprintf(os.Stderr, "Sock %d\n", nCPU/cpuid.CPU.LogicalCores)
	}
	fmt.Fprintf(os.Stderr, "Mmry %d\n", memory.TotalMemory()/(1024*1024*1024))

	fmt.Fprintf(os.Stderr, "Proc %d\n", numProcs)
	fmt.Fprintf(os.Stderr, "Serv %d\n", numServe)
	fmt.Fprintf(os.Stderr, "Chan %d\n", chanDepth)
	fmt.Fprintf(os.Stderr, "Heap %d\n", heapSize)
	fmt.Fprintf(os.Stderr, "Farm %d\n", farmSize)
	fmt.Fprintf(os.Stderr, "Gogc %d\n", goGc)

	fi, err := os.Stdin.Stat()
	if err == nil {
		mode := fi.Mode().String()
		fmt.Fprintf(os.Stderr, "Mode %s\n", mode)
	}

	fmt.Fprintf(os.Stderr, "\n")
}

// PrintHelp finds and prints indicated help file
func PrintHelp(programName, helpFile string) {

	ex, eerr := os.Executable()
	if eerr == nil {
		fmt.Printf("%s %s\n\n", programName, EDirectVersion)
		exPath := filepath.Dir(ex)
		fpath := filepath.Join(exPath, "help", helpFile)
		file, ferr := os.Open(fpath)
		if file != nil && ferr == nil {
			scanner := bufio.NewScanner(file)
			for scanner.Scan() {
				fmt.Println(scanner.Text())
			}
		}
		file.Close()
		fmt.Printf("\n")
	}
}

// creates lookup tables to simplify the tokenizer
func init() {

	startTime = time.Now()

	for i := range inBlank {
		inBlank[i] = false
	}
	inBlank[' '] = true
	inBlank['\t'] = true
	inBlank['\n'] = true
	inBlank['\r'] = true
	inBlank['\f'] = true

	// first character of element cannot be a digit, dash, or period
	for i := range inFirst {
		inFirst[i] = false
	}
	for ch := 'A'; ch <= 'Z'; ch++ {
		inFirst[ch] = true
	}
	for ch := 'a'; ch <= 'z'; ch++ {
		inFirst[ch] = true
	}
	inFirst['_'] = true
	// extend legal initial letter with at sign and digits to handle biological data converted from JSON
	inFirst['@'] = true
	for ch := '0'; ch <= '9'; ch++ {
		inFirst[ch] = true
	}

	// remaining characters also includes colon for namespace
	for i := range inElement {
		inElement[i] = false
	}
	for ch := 'A'; ch <= 'Z'; ch++ {
		inElement[ch] = true
	}
	for ch := 'a'; ch <= 'z'; ch++ {
		inElement[ch] = true
	}
	for ch := '0'; ch <= '9'; ch++ {
		inElement[ch] = true
	}
	inElement['_'] = true
	inElement['-'] = true
	inElement['.'] = true
	inElement[':'] = true

	// embedded markup and math tags are lower case
	for i := range inLower {
		inLower[i] = false
	}
	for ch := 'a'; ch <= 'z'; ch++ {
		inLower[ch] = true
	}
	for ch := '0'; ch <= '9'; ch++ {
		inLower[ch] = true
	}
	inLower['_'] = true
	inLower['-'] = true
	inLower['.'] = true
	inLower[':'] = true

	// shortcut to find <, >, or &, or non-ASCII
	for i := range inContent {
		inContent[i] = false
	}
	for i := 0; i <= 127; i++ {
		inContent[i] = true
	}
	inContent['<'] = false
	inContent['>'] = false
	inContent['&'] = false

	// ASN.1 tag can have letters, digits, and dashes
	for i := range inAsnTag {
		inAsnTag[i] = false
	}
	for ch := 'A'; ch <= 'Z'; ch++ {
		inAsnTag[ch] = true
	}
	for ch := 'a'; ch <= 'z'; ch++ {
		inAsnTag[ch] = true
	}
	for ch := '0'; ch <= '9'; ch++ {
		inAsnTag[ch] = true
	}
	inAsnTag['-'] = true

	// ASN.1 string starts and ends with double quote
	for i := range inAsnString {
		inAsnString[i] = false
	}
	for ch := ' '; ch <= '~'; ch++ {
		inAsnString[ch] = true
	}
	inAsnString['"'] = false
	// "

	// ASN.1 bit string flanked by apostrophe, followed by H or B
	for i := range inAsnBits {
		inAsnBits[i] = false
	}
	for ch := ' '; ch <= '~'; ch++ {
		inAsnBits[ch] = true
	}
	inAsnBits['\''] = false

	// initialize reading and cleaning options with default values
	SetOptions(false, false, false, false, false, false, false)

	// initialize performance tuning variables with default values
	SetTunings(0, 0, 0, 0, 0, 0, 0, false)
}
