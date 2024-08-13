#!/bin/bash

ulimit -S -s unlimited


YYYY=2023
MM=05
DD=01
HH=12

RES=48

EXEC_DIR=../exec/bin/
RESTART_DIR=/scratch2/BMC/gsienkf/Clara.Draper/DA_test_cases/land-IMSproc/restarts/

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

TSTUB=oro_C${RES}

cat >> fscf.nml << EOF
 &fSCF_nml
  source=1,
  idim=$RES, 
  jdim=$RES,
  jdate=$JDATE,
  otype=${TSTUB},
  yyyymmddhh=${YYYY}${MM}${DD}.${HH},
  imsformat=1,
  imsversion=1.3,
  imsres=4km,
  IMS_obs_path="${SNOW_OBS_DIR}/snow_ice_cover/IMS/ascii/${YYYY}/",
  IMS_ind_path="${SNOW_OBS_DIR}/snow_ice_cover/IMS/index_files/",
  /
EOF


# stage restarts
for tt in 1 2 3 4 5 6
do 
ln -s ${RESTART_DIR}/${YYYY}${MM}${DD}.${HH}0000.sfc_data.tile${tt}.nc .
done

ln -s ${EXEC_DIR}/calcfSCF.exe .

./calcfSCF.exe

