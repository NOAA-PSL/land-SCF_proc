#! /usr/bin/env bash
set -eux

# check if part of workflow. If so, use those modules.
if [ -f ../workflow_mods_bash ]; then 
  echo 'using workflow modules'
  source ../workflow_mods_bash
else
  echo 'using own modules'
  source hera_modules
fi 

export FCMP=mpiifort

# Check final exec folder exists
if [ ! -d "./exec" ]; then
  mkdir ./exec
fi

cd ./sorc/
./makefile.sh
