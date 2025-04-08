#!/bin/bash

ulimit -S -s unlimited


YYYY=2016
MM=05
DD=01
HH=12

RES=384

EXEC_DIR=../exec/bin/
RESTART_DIR=/scratch2/NCEPDEV/stmp1/Yuan.Xue/pull_request/restart_files/

############# NO CHANGES NEEDED BELOW THIS LINE
source ../env_GDASApp 

export SNOW_OBS_DIR=/scratch2/NCEPDEV/land/data/DA/

if [ -e fscf.nml ]; then 
rm fscf.nml
fi
if [ -e ${YYYY}${MM}${DD}.${HH}0000.sfc_data.tile1.nc ]; then 
rm ${YYYY}${MM}${DD}.${HH}0000.sfc_data.tile*
fi 
if [ -e calcfSCF.exe ]; then 
rm calcfSCF.exe
fi

DOY=$(date -d "${YYYY}-${MM}-${DD}" +%j)
JDATE=$YYYY$DOY

TSTUB=C${RES}.mx025_oro_data

# Users can specify viirs_threshold below when constructing the namelist file
# Tested thresholds are: 0.1, 0.3(default), 0.5
cat >> fscf.nml << EOF
 &fSCF_nml
  source=2,
  idim=$RES, jdim=$RES,
  otype=${TSTUB},
  jdate=${YYYY}${DOY},
  yyyymmddhh=${YYYY}${MM}${DD}.${HH},
  viirsversion="002",
  viirs_threshold=0.3,
  VIIRS_OBS_PATH="${SNOW_OBS_DIR}/snow_ice_cover/VIIRS/${YYYY}/",
  VIIRS_IND_PATH="${SNOW_OBS_DIR}/snow_ice_cover/VIIRS/index_files/"
  /
EOF


# stage restarts
for tt in 1 2 3 4 5 6
do 
ln -s ${RESTART_DIR}/${YYYY}${MM}${DD}.${HH}0000.sfc_data.tile${tt}.nc .
done

ln -s ${EXEC_DIR}/calcfSCF.exe .

./calcfSCF.exe

