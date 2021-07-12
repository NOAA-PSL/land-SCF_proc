#! /bin/sh -l
#BATCH --job-name=global_cycle
#SBATCH --account=da-cpu ## change
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --tasks-per-node=6
#SBATCH --cpus-per-task=1
#SBATCH -t 00:30:00
#SBATCH -o global_cycle_run.log
#SBATCH -e global_cycle_run.log
#SBATCH --export=NONE
#SBATCH --comment=376d154725287bda4729a1014e6f2d89

IY=2019
IM=12
ID=15
IH=18

RES=96

EXEC_DIR=/scratch1/NCEPDEV/da/Youlong.Xia/psl_ClaraDraper/IMSobsproc/exec # change

############# NO CHANGES NEEDED BELOW THIS LINE
export LD_LIBRARY_PATH=/apps/hdf5/1.10.5/intel/18.0.5.274/lib:/apps/nco/4.7.0/intel/18.0.3.051/lib:/apps/netcdf/4.7.4/intel/18.0.5/lib:/apps/pnetcdf/1.10.0/intel/16.1.150/impi/5.1.2.150/lib:/apps/wgrib2/2.0.8/intel/18.0.3.222/lib:/apps/intel/compilers_and_libraries_2018/linux/mpi/intel64/lib::/apps/slurm/default/lib:/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/compiler/lib/intel64:/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/ipp/lib/intel64:/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/compiler/lib/intel64_lin:/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/mkl/lib/intel64_lin:/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/tbb/lib/intel64/gcc4.7:/apps/intel/parallel_studio_xe_2018.4.057/debugger_2018/libipt/intel64/lib:/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/daal/lib/intel64_lin:/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/daal/../tbb/lib/intel64_lin/gcc4.4:$LD_LIBRARY_PATH

export FIXDIR=/scratch1/NCEPDEV/global/glopara/fix/fix_fv3_gmted2010/C${RES}/
export SNOW_OBS_DIR=/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/

rm fims.nml

cat >> fims.nml << EOF
 &fIMS_nml
  idim=$RES, jdim=$RES,
  date_str=$IY$IM$ID$IH,
  IMS_SNOWCOVER_PATH="${SNOW_OBS_DIR}/IMS/", 
  IMS_INDEXES_PATH="${SNOW_OBS_DIR}/IMS_INDEXES/"
  / 
EOF

# link inputs 
for tt in 1 2 3 4 5 6
do
if [[ ! -e C$RES.vegetation_type.tile${tt}.nc  ]]; then 
ln -s $FIXDIR/fix_sfc/C$RES.vegetation_type.tile${tt}.nc .
fi 
done 


srun '--export=ALL' -n 6 ${EXEC_DIR}/calcfIMS

