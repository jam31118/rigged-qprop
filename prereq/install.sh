#!/bin/bash

if [ -n "$1" ]; then base_dir="$1"; else base_dir=$HOME; fi

#base_dir=$HOME
#if [ -n "$1" ]; then base_dir="$1"; fi

echo "[ LOG ] Entered base directory: $base_dir"
original_dir=$(pwd)
cd $base_dir
base_dir_abs=$(pwd)
cd $original_dir
echo "[ LOG ] Absolute base directory: $base_dir_abs"

#printf "[  Q  ] Is this directory correct? [y/n] "
#
#read reply
#echo "[ LOG ] Accepted answer: $reply"
#if [ "$reply" == "y" ]; then echo "[ LOG ] Installation is starting . . .";
#else
#  if [ "$reply" == "n" ]; then echo "[ LOG ] Please select right install directory. Aborting . . ."; exit 1;
#  else echo "[ERROR] Select either 'y' or 'n'"; exit 1; fi
#fi

#echo base_dir is $base_dir

home_dir="$base_dir_abs/tool"
if [ ! -d "$home_dir" ]; then mkdir -p $home_dir; fi
srcs_dir="$home_dir/src"
if [ ! -d "$srcs_dir" ]; then mkdir -p $srcs_dir; fi
installs_dir="$home_dir/install"
if [ ! -d "$installs_dir" ]; then mkdir -p $installs_dir; fi
builds_dir="$home_dir/build"
if [ ! -d "$builds_dir" ]; then mkdir -p $builds_dir; fi

echo "home_dir: $home_dir"
echo "srcs_dir: $srcs_dir"
echo "installs_dir: $installs_dir"
echo "builds_dir: $builds_dir"

## Install GSL 2.4 ##
SRC_NAME="gsl-2.4"
SRC_DIRNAME="$SRC_NAME"

SRC_DIR="$srcs_dir/$SRC_NAME"
BUILD_DIR="$builds_dir/$SRC_NAME"
if [ ! -d "$BUILD_DIR" ]; then mkdir -p $BUILD_DIR; fi
INSTALL_DIR="$installs_dir/$SRC_NAME"
if [ ! -d "$INSTALL_DIR" ]; then mkdir -p $INSTALL_DIR; fi
echo "Installation started for $SRC_NAME at $INSTALL_DIR"

SRC_TARBALL="$SRC_NAME.tar.gz"
SRC_URL="http://ftp-stud.hs-esslingen.de/pub/Mirrors/ftp.gnu.org/gsl/$SRC_TARBALL"

cd $srcs_dir
wget $SRC_URL
tar xzvf $SRC_TARBALL
cd $BUILD_DIR
$SRC_DIR/configure --prefix=$INSTALL_DIR
make -j4
make check
make install
make installcheck


## Install OpenMPI 1.10.7
SRC_NAME="openmpi-1.10.7"
SRC_DIRNAME="$SRC_NAME"

SRC_DIR="$srcs_dir/$SRC_NAME"
BUILD_DIR="$builds_dir/$SRC_NAME"
if [ ! -d "$BUILD_DIR" ]; then mkdir -p $BUILD_DIR; fi
INSTALL_DIR="$installs_dir/$SRC_NAME"
if [ ! -d "$INSTALL_DIR" ]; then mkdir -p $INSTALL_DIR; fi
echo "Installation started for $SRC_NAME at $INSTALL_DIR"

SRC_TARBALL="$SRC_NAME.tar.gz"
SRC_URL="http://www.open-mpi.de/software/ompi/v1.10/downloads/$SRC_TARBALL"

cd $srcs_dir
wget $SRC_URL
tar xzvf $SRC_TARBALL
cd $BUILD_DIR
$SRC_DIR/configure --prefix=$INSTALL_DIR
make all -j4
make install


## Install Boost 1.65.1 ##

SRC_NAME="boost_1_65_1"
SRC_DIRNAME="$SRC_NAME"

SRC_DIR="$srcs_dir/$SRC_NAME"
#BUILD_DIR="$builds_dir/$SRC_NAME"
#if [ ! -d "$BUILD_DIR" ]; then mkdir -p $BUILD_DIR; fi
INSTALL_DIR="$installs_dir/$SRC_NAME"
if [ ! -d "$INSTALL_DIR" ]; then mkdir -p $INSTALL_DIR; fi
echo "Installation started for $SRC_NAME at $INSTALL_DIR"

SRC_TARBALL="$SRC_NAME.tar.gz"
SRC_URL="https://dl.bintray.com/boostorg/release/1.65.1/source/$SRC_TARBALL"

cd $srcs_dir
wget $SRC_URL
tar xzvf $SRC_TARBALL
#cd $BUILD_DIR
cd $SRC_DIR
./bootstrap.sh --prefix=$INSTALL_DIR --with-libraries=timer
./b2 install

