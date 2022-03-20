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
// File Name:  chan.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"io"
	"os"
	"strings"
)

// stringChanReader connect a string output channel to an io.Reader interface
type stringChanReader struct {
	c <-chan string
	s string
}

func (r *stringChanReader) Read(b []byte) (n int, err error) {

	if r.s != "" {
		n = copy(b, []byte(r.s))
		r.s = r.s[n:]
		return
	}

	for str := range r.c {
		r.s = str
		n = copy(b, []byte(r.s))
		r.s = r.s[n:]
		return
	}

	return 0, io.EOF
}

// ChanToReader converts a string channel to an ioReader
func ChanToReader(inp <-chan string) io.Reader {

	if inp == nil {
		return nil
	}

	return &stringChanReader{c: inp}
}

// ChanToStdout sends a string channel to stdout
func ChanToStdout(inp <-chan string) {

	if inp == nil {
		return
	}

	last := ""

	for str := range inp {
		last = str
		os.Stdout.WriteString(str)
	}

	if !strings.HasSuffix(last, "\n") {
		os.Stdout.WriteString("\n")
	}
}

// ChanToString converts a string channel to a string
func ChanToString(inp <-chan string) string {

	if inp == nil {
		return ""
	}

	var buffer strings.Builder

	last := ""

	for str := range inp {
		last = str
		buffer.WriteString(str)
	}

	if !strings.HasSuffix(last, "\n") {
		buffer.WriteString("\n")
	}

	return buffer.String()
}
