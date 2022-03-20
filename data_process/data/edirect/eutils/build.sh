#!/bin/sh

# Public domain notice for all NCBI EDirect scripts is located at:
# https://www.ncbi.nlm.nih.gov/books/NBK179288/#chapter6.Public_Domain_Notice

if [ ! -f "go.mod" ]
then
  go mod init eutils
fi
if [ ! -f "go.sum" ]
then
  go mod tidy
fi

while [ "$#" -ne 0 ]
do
  case "$1" in
    -vendor | vendor )
      shift
      if [ ! -d "vendor" ]
      then
        go mod vendor -e
      fi
      ;;
    * )
      break
      ;;
  esac
done

go build
