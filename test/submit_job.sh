#! /bin/sh -l

IY=2019
IM=12
ID=15
IH=18
JDAY=350

RES=48

EXEC_DIR=/scratch2/BMC/gsienkf/Clara.Draper/gerrit-hera/snowDA/IMSobsproc_Mike/exec/

############# NO CHANGES NEEDED BELOW THIS LINE
export LD_LIBRARY_PATH=/apps/hdf5/1.10.5/intel/18.0.5.274/lib:/apps/nco/4.7.0/intel/18.0.3.051/lib:/apps/netcdf/4.7.4/intel/18.0.5/lib:/apps/pnetcdf/1.10.0/intel/16.1.150/impi/5.1.2.150/lib:/apps/wgrib2/2.0.8/intel/18.0.3.222/lib:/apps/intel/compilers_and_libraries_2018/linux/mpi/intel64/lib::/apps/slurm/default/lib:/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/compiler/lib/intel64:/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/ipp/lib/intel64:/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/compiler/lib/intel64_lin:/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/mkl/lib/intel64_lin:/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/tbb/lib/intel64/gcc4.7:/apps/intel/parallel_studio_xe_2018.4.057/debugger_2018/libipt/intel64/lib:/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/daal/lib/intel64_lin:/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/daal/../tbb/lib/intel64_lin/gcc4.4:$LD_LIBRARY_PATH

export SNOW_OBS_DIR=/scratch1/NCEPDEV/stmp2/Michael.Barlage/DA/IMSdata

rm fims.nml

DATE_STRING=$IY$JDAY

cat >> fims.nml << EOF
 &fIMS_nml
  idim=$RES, jdim=$RES,
  date_str=$DATE_STRING,
  IMS_SNOWCOVER_PATH="${SNOW_OBS_DIR}/IMS/", 
  IMS_INDEX_PATH="/scratch2/NCEPDEV/land/Youlong.Xia/save/obs-product/ims/ims_index/"
  / 
EOF

${EXEC_DIR}/calcfIMS

