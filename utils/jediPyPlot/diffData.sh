#! /bin/sh

module load nco/4.9.3

for tt in 1 2 3 4 5 6
do
ncdiff test1.tile$tt.nc test2.tile$tt.nc diff.tile$tt.nc
done


