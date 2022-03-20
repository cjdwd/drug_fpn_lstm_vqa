#!/bin/sh

# Public domain notice for all NCBI EDirect scripts is located at:
# https://www.ncbi.nlm.nih.gov/books/NBK179288/#chapter6.Public_Domain_Notice

# determine current platform
platform=""
osname=`uname -s`
cputype=`uname -m`
case "$osname-$cputype" in
  Linux-x86_64 )           platform=Linux ;;
  Darwin-x86_64 )          platform=Darwin ;;
  Darwin-*arm* )           platform=Silicon ;;
  CYGWIN_NT-* | MINGW*-* ) platform=CYGWIN_NT ;;
  Linux-*arm* )            platform=ARM ;;
  * )                      platform=UNSUPPORTED ;;
esac

canBuildSilicon=false

gv=$( go version )
case "$gv" in
  # only build for Apple Silicon if using Go compiler 1.16 or later
  *go1.16* | *go1.17* | *go1.18* )
    canBuildSilicon=true
    ;;
  * )
    ;;
esac

crossCompileAll=false
install=false
cleanup=true
vendor=false
target=""

# process optional command-line arguments
while [ "$#" -ne 0 ]
do
  case "$1" in
    -install | install )
      # install native executables on development machine
      install=true
      # default executable path
      target="$HOME/Misc/scripts/"
      shift
      ;;
    -desktop | desktop )
      # place native executables on desktop
      install=true
      target="$HOME/Desktop/"
      shift
      ;;
    -silicon | silicon )
      # coerce platform to create Silicon executables
      platform=Silicon
      install=true
      target="$HOME/Desktop/"
      # but do not remove existing native binaries on deskop
      cleanup=false
      shift
      ;;
    -distrib | -distribute | distrib | distribute )
      # cross-compile all versions for ftp distribution
      crossCompileAll=true
      install=true
      # default distribution path
      target="$HOME/goxtract/"
      shift
      ;;
    -vendor | vendor )
      vendor=true
      shift
      ;;
    * )
      if [ -n "$1" ]
      then
        # allow override of default target path
        install=true
        target="$1"
      fi
      # break out of loop
      break
      ;;
  esac
done

# create module files
if [ ! -f "go.mod" ]
then
  go mod init edirect
  # add explicit location to find local helper package
  echo "replace eutils => ../eutils" >> go.mod
fi
if [ ! -f "go.sum" ]
then
  go mod tidy
fi

# cache external dependencies
if [ "$vendor" = true ] && [ ! -d "vendor" ]
then
  go mod vendor -e
fi

# erase any existing executables in current directory
for plt in Darwin Silicon Linux CYGWIN_NT ARM
do
  rm -f *.$plt
done

# platform-specific compiler environment variable values
mods="darwin amd64 Darwin \
      darwin arm64 Silicon \
      linux amd64 Linux \
      windows 386 CYGWIN_NT \
      linux arm ARM"

# build all executables for each selected platform
echo "$mods" |
xargs -n 3 sh -c 'echo "$0 $1 $2"' |
while read os ar pl
do
  if [ "$pl" != "$platform" ] && [ "$crossCompileAll" = false ]
  then
    continue
  fi
  if [ "$pl" = "Silicon" ] && [ "$canBuildSilicon" = false ]
  then
    continue
  fi
  for exc in xtract rchive transmute
  do
    env GOOS="$os" GOARCH="$ar" go build -o "$exc.$pl" "$exc.go"
  done
done

if [ "$install" = true ] && [ -n "$target" ]
then
  if [ "$cleanup" = true ]
  then
    # remove old executables from target
    for plt in Darwin Silicon Linux CYGWIN_NT ARM
    do
      rm -f $target/*.$plt
    done
  fi
  # copy new executables to target
  for plt in Darwin Silicon Linux CYGWIN_NT ARM
  do
    for exc in xtract rchive transmute
    do
      if [ -f "$exc.$plt" ]
      then
        mv -f "$exc.$plt" "$target"
      fi
    done
  done
fi

# erase any remaining executables after compiling
for plt in Darwin Silicon Linux CYGWIN_NT ARM
do
  rm -f *.$plt
done
