#!/bin/bash

## Source required scripts
source $(dirname "$0")/script/colors.sh

base_dir=""
if [ -n "$1" ]
then
  if [ -d "$1" ]; then base_dir="$1"
  else
    (>&2 echo "[ERROR] Given directory ($1) doesn't exist.")
    exit -1
  fi
else
  (>&2 echo "[ERROR] No directory given.")
  exit -1
fi

PROGRAM_NAME="gsl"
VERSION="-2.4"
ARCH_EXT=".tar.gz"
SRC_URL_PREFIX="http://ftp-stud.hs-esslingen.de/pub/Mirrors/ftp.gnu.org/gsl"

cd $base_dir
base_dir_abs=$(pwd)
INSTALL_DIR="$base_dir_abs/$PROGRAM_NAME"
if [ -d "$INSTALL_DIR" ]; then rm -rf $INSTALL_DIR/*; fi
BUILD_DIR="$INSTALL_DIR/build"
SRC_DOWN_DIR="$INSTALL_DIR/src"

for dir in $INSTALL_DIR $BUILD_DIR $SRC_DOWN_DIR
do
  if [ ! -d "$dir" ]; then mkdir $dir; fi
done

SRC_NAME="$PROGRAM_NAME$VERSION"
SRC_TARBALL="$SRC_NAME$ARCH_EXT"
SRC_URL="$SRC_URL_PREFIX/$SRC_TARBALL"
SRC_DIR="$SRC_DOWN_DIR/$SRC_NAME"

echo "SRC_NAME: $SRC_NAME"
echo "SRC_DIR: $SRC_DIR"
echo "BUILD_DIR: $BUILD_DIR"

cd $SRC_DOWN_DIR
wget $SRC_URL
if [ "$?" -ne "0" ]
then
  (>&2 echo "${ERROR} Failed to wget during downloading '$SRC_NAME' from '$SRC_URL'")
  echo "${LOG} Trying download using curl . . ."
  curl -LO --fail "$SRC_URL"
  if [ "$?" -ne 0 ]
    then (>&2 echo "${ERROR} Failed to run 'curl' for downloading '$SRC_NAME' from '$SRC_URL'")
    exit -1
  fi
fi

tar xzvf $SRC_TARBALL
cd $BUILD_DIR
$SRC_DIR/configure --prefix=$INSTALL_DIR
make -j4
make check
make install
make installcheck

return_code="$?"
if [ "$return_code" -ne "0" ]
then
  (>&2 echo "[ERROR] Failed to install GSL. Aborting . . .")
  exit -1
fi

rm -rf $SRC_DOWN_DIR $BUILD_DIR

exit 0

