#!/bin/bash

# cmake build script for IMS processing code 
# uses GDASApp environment for consistency with workflow(s) 
# only actually needs intel, impi, and netcdf
# eg: intel/2020.2 impi/2018.0.4 netcdf/4.7.0

if [ $# == 1 ]; then 
        echo "setting env from $1"
        env_file=$1
else 
        # assume installed a subdir of DA_update
        env_file="../env_GDASApp" 
fi 
source $env_file

if [[ -d exec ]]; then 
  rm -rf exec 
fi
mkdir exec

if [[ -d build ]]; then 
  rm -rf build 
fi
mkdir build
cd build 

# configure 
cmake .. -DCMAKE_INSTALL_PREFIX=../exec

# build 
cmake --build  .

# install 
cmake --build . --target install

cd .. 

exit 0

