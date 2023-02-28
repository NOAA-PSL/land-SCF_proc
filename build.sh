#! /usr/bin/env bash
set -eux
machine='aws'

# check if part of workflow. If so, use those modules.
if [ $machine == "hera" ]; then
   source ./hera_modules
elif [ $machine == 'aws' ]; then
  source ./aws_modules
fi 

export FCMP=ifort

# Check final exec folder exists
if [ ! -d "./exec" ]; then
  mkdir ./exec
fi

cd ./sorc/
./Makefile.sh
