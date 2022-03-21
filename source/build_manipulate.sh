#!/bin/bash -e
#Set default variables:
BUILD_TYPE=Release
CLEAN_INST="false"

#Parse input options
while [ "$1" != "" ]; do
  case $1 in
    -d | --debug)
      BUILD_TYPE=Debug
      shift;;
    -C | --Clean-install)
      CLEAN_INST="true"
      shift;;
    *)
      echo "Unrecognised build option: $1"
      return 1;;
  esac
done

#Set bin directory based on operating system
SRC_DIR=$(pwd)
if [ "$OSTYPE" == "msys" ]; then
  BIN_DIR="${SRC_DIR}/../bin/win"
else
  BIN_DIR="${SRC_DIR}/../bin/unix"
fi

#If desired wipe build directory
if [ $CLEAN_INST == "true" ]; then
  rm -rf ./build
fi

#Check if build directory is present and create if not
if [ ! -d ./build ]; then
  mkdir build
  #Configure
  CMAKE_OPTS="
      -D CMAKE_BUILD_TYPE=${BUILD_TYPE}
      -D CMAKE_INSTALL_DIR=${BIN_DIR}
      -D CMAKE_Fortran_COMPILER="gfortran"
      -D CMAKE_C_COMPILER="gcc"
      -D CMAKE_CXX_COMPILER="gcc"
  "
  cd ./build

  #Call cmake to create build files
  cmake $CMAKE_OPTS ../ -G "Unix Makefiles"
else
  cd ./build
fi

#Execute Build
if [ "$OSTYPE" == "msys" ]; then
  make -j install SHELL=cmd
else
  make -j install
fi

#Move back out of build directory
cd ../
