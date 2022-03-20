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
// File Name:  cache.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"hash/crc32"
	"io"
	"os"
	"path"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"
)

// CreateUIDReader sends PMIDs and their numeric orders down a channel.
// This allows detection of updated records that appear shortly after
// earlier versions, preventing the wrong version from being saved.
func CreateUIDReader(in io.Reader) <-chan XMLRecord {

	if in == nil {
		return nil
	}

	out := make(chan XMLRecord, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create uid reader channel\n")
		os.Exit(1)
	}

	// uidReader reads uids from input stream and sends through channel
	uidReader := func(in io.Reader, out chan<- XMLRecord) {

		// close channel when all records have been processed
		defer close(out)

		scanr := bufio.NewScanner(in)

		idx := 0
		for scanr.Scan() {

			// read lines of identifiers
			file := scanr.Text()
			idx++

			pos := strings.Index(file, ".")
			if pos >= 0 {
				// remove version suffix
				file = file[:pos]
			}

			out <- XMLRecord{Index: idx, Text: file}
		}
	}

	// launch single uid reader goroutine
	go uidReader(in, out)

	return out
}

// CreateStashers saves records to archive, multithreaded for performance, use of UID
// position index allows it to prevent earlier version from overwriting later version
func CreateStashers(stash, parent, indx, sfx string, hash, zipp bool, report int, inp <-chan XMLRecord) <-chan string {

	if inp == nil {
		return nil
	}

	out := make(chan string, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create stasher channel\n")
		os.Exit(1)
	}

	find := ParseIndex(indx)

	if zipp {
		sfx += ".gz"
	}

	type StasherType int

	const (
		OKAY StasherType = iota
		WAIT
		BAIL
	)

	// mutex to protect access to inUse map
	var wlock sync.Mutex

	// map to track files currently being written
	inUse := make(map[string]int)

	// lockFile function prevents colliding writes
	lockFile := func(id string, index int) StasherType {

		// map is non-reentrant, protect with mutex
		wlock.Lock()

		// multiple return paths, schedule the unlock command up front
		defer wlock.Unlock()

		idx, ok := inUse[id]

		if ok {
			if idx < index {
				// later version is being written by another goroutine, skip this
				return BAIL
			}
			// earlier version is being written by another goroutine, wait
			return WAIT
		}

		// okay to write file, mark in use to prevent collision
		inUse[id] = index
		return OKAY
	}

	// freeFile function removes entry from inUse map
	freeFile := func(id string) {

		wlock.Lock()

		// free entry in map, later versions of same record can now be written
		delete(inUse, id)

		wlock.Unlock()
	}

	// mutex to protect access to rollingCount variable
	var tlock sync.Mutex

	rollingCount := 0

	countSuccess := func() {

		tlock.Lock()

		rollingCount++
		if rollingCount >= report {
			rollingCount = 0
			// print dot (progress monitor)
			fmt.Fprintf(os.Stderr, ".")
		}

		tlock.Unlock()
	}

	// stashRecord saves individual XML record to archive file accessed by trie
	stashRecord := func(str, id string, index int) string {

		pos := strings.Index(id, ".")
		if pos >= 0 {
			// remove version from UID
			id = id[:pos]
		}

		var arry [132]rune
		trie := MakeArchiveTrie(id, arry)
		if trie == "" {
			return ""
		}

		attempts := 5
		keepChecking := true

		for keepChecking {
			// check if file is not being written by another goroutine
			switch lockFile(id, index) {
			case OKAY:
				// okay to save this record now
				keepChecking = false
			case WAIT:
				// earlier version is being saved, wait one second and try again
				time.Sleep(time.Second)
				attempts--
				if attempts < 1 {
					// could not get lock after several attempts
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to save '%s'\n", id)
					return ""
				}
			case BAIL:
				// later version is being saved, skip this one
				return ""
			default:
			}
		}

		// delete lock after writing file
		defer freeFile(id)

		dpath := path.Join(stash, trie)
		if dpath == "" {
			return ""
		}
		err := os.MkdirAll(dpath, os.ModePerm)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return ""
		}
		fpath := path.Join(dpath, id+sfx)
		if fpath == "" {
			return ""
		}

		// overwrites and truncates existing file
		fl, err := os.Create(fpath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return ""
		}

		res := ""

		if hash {
			// calculate hash code for verification table
			hsh := crc32.NewIEEE()
			hsh.Write([]byte(str))
			val := hsh.Sum32()
			res = strconv.FormatUint(uint64(val), 10)
		}

		if zipp {

			zpr, err := gzip.NewWriterLevel(fl, gzip.DefaultCompression)

			if err == nil {

				wrtr := bufio.NewWriter(zpr)

				// compress and copy record to file
				wrtr.WriteString(str)
				if !strings.HasSuffix(str, "\n") {
					wrtr.WriteString("\n")
				}

				err = wrtr.Flush()
				if err != nil {
					fmt.Fprintf(os.Stderr, "%s\n", err.Error())
					return ""
				}

				err = zpr.Close()
				if err != nil {
					fmt.Fprintf(os.Stderr, "%s\n", err.Error())
					return ""
				}

			} else {

				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			}

		} else {

			// copy uncompressed record to file
			fl.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				fl.WriteString("\n")
			}
		}

		// fl.Sync()

		err = fl.Close()
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return ""
		}

		// progress monitor prints dot every 1000 (.xml) or 50000 (.e2x) records
		countSuccess()

		return res
	}

	// xmlStasher reads from channel and calls stashRecord
	xmlStasher := func(wg *sync.WaitGroup, inp <-chan XMLRecord, out chan<- string) {

		defer wg.Done()

		for ext := range inp {

			ext.Ident = FindIdentifier(ext.Text, parent, find)

			hsh := stashRecord(ext.Text, ext.Ident, ext.Index)

			res := ext.Ident
			if hash {
				res += "\t" + hsh
			}
			res += "\n"

			runtime.Gosched()

			out <- res
		}
	}

	var wg sync.WaitGroup

	// launch multiple stasher goroutines
	for i := 0; i < NumServe(); i++ {
		wg.Add(1)
		go xmlStasher(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all stashers are done
	go func() {
		wg.Wait()
		close(out)
		// print newline after rows of dots (progress monitor)
		fmt.Fprintf(os.Stderr, "\n")
	}()

	return out
}

// CreateFetchers returns uncompressed records from archive, multithreaded for speed
func CreateFetchers(stash, sfx string, zipp bool, inp <-chan XMLRecord) <-chan XMLRecord {

	if inp == nil {
		return nil
	}

	out := make(chan XMLRecord, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create fetcher channel\n")
		os.Exit(1)
	}

	if zipp {
		sfx += ".gz"
	}

	fetchRecord := func(file string, buf bytes.Buffer) string {

		var arry [132]rune
		trie := MakeArchiveTrie(file, arry)

		if file == "" || trie == "" {
			return ""
		}

		fpath := path.Join(stash, trie, file+sfx)
		if fpath == "" {
			return ""
		}

		iszip := zipp

		inFile, err := os.Open(fpath)

		// if failed to find ".xml" or ".e2x" file, try ".xml.gz" or ".e2x.gz" without requiring -gzip
		if err != nil && os.IsNotExist(err) && !zipp {
			iszip = true
			fpath := path.Join(stash, trie, file+sfx+".gz")
			if fpath == "" {
				return ""
			}
			inFile, err = os.Open(fpath)
		}
		if err != nil {
			msg := err.Error()
			if !strings.HasSuffix(msg, "no such file or directory") && !strings.HasSuffix(msg, "cannot find the path specified.") {
				fmt.Fprintf(os.Stderr, "%s\n", msg)
			}
			return ""
		}

		defer inFile.Close()

		brd := bufio.NewReader(inFile)

		if iszip {

			zpr, err := gzip.NewReader(brd)

			defer zpr.Close()

			if err == nil {
				// copy and decompress cached file contents
				buf.ReadFrom(zpr)
			}

		} else {

			// copy cached file contents
			buf.ReadFrom(brd)
		}

		str := buf.String()

		return str
	}

	// xmlFetcher reads XML from file
	xmlFetcher := func(wg *sync.WaitGroup, inp <-chan XMLRecord, out chan<- XMLRecord) {

		// report when more records to process
		defer wg.Done()

		var buf bytes.Buffer

		for ext := range inp {

			buf.Reset()

			str := fetchRecord(ext.Text, buf)

			runtime.Gosched()

			out <- XMLRecord{Index: ext.Index, Text: str}
		}
	}

	var wg sync.WaitGroup

	// launch multiple fetcher goroutines
	for i := 0; i < NumServe(); i++ {
		wg.Add(1)
		go xmlFetcher(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all fetchers are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

// CreateCacheStreamers returns compressed records from archive, multithreaded for speed,
// could be used for sending records over network to be decompressed later by client
func CreateCacheStreamers(stash string, inp <-chan XMLRecord) <-chan XMLRecord {

	if inp == nil {
		return nil
	}

	out := make(chan XMLRecord, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create streamer channel\n")
		os.Exit(1)
	}

	sfx := ".xml.gz"

	getRecord := func(file string, buf bytes.Buffer) []byte {

		var arry [132]rune
		trie := MakeArchiveTrie(file, arry)

		if file == "" || trie == "" {
			return nil
		}

		fpath := path.Join(stash, trie, file+sfx)
		if fpath == "" {
			return nil
		}

		inFile, err := os.Open(fpath)

		if err != nil {
			msg := err.Error()
			if !strings.HasSuffix(msg, "no such file or directory") && !strings.HasSuffix(msg, "cannot find the path specified.") {
				fmt.Fprintf(os.Stderr, "%s\n", msg)
			}
			return nil
		}

		brd := bufio.NewReader(inFile)

		// copy cached file contents
		buf.ReadFrom(brd)

		data := buf.Bytes()

		inFile.Close()

		return data
	}

	// xmlStreamer reads compressed XML from file
	xmlStreamer := func(wg *sync.WaitGroup, inp <-chan XMLRecord, out chan<- XMLRecord) {

		// report when more records to process
		defer wg.Done()

		var buf bytes.Buffer

		for ext := range inp {

			buf.Reset()

			data := getRecord(ext.Text, buf)

			runtime.Gosched()

			out <- XMLRecord{Index: ext.Index, Data: data}
		}
	}

	var wg sync.WaitGroup

	// launch multiple streamer goroutines
	for i := 0; i < NumServe(); i++ {
		wg.Add(1)
		go xmlStreamer(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all streamers are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}
