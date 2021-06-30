Test case for processing the IMS observations onto UFS model grid.

To run, 
in submit_job.sh  change 1) EXEC_DIR to point to directory with your executable. 
                         2) account to one you have access to 

>sbatch submit_job.sh

output:
IMSfSCA.tile*.nc for each tile. 

output should match files in  ./output/ 
check_diffs.sh  will compare them 

It should be noted that output result is for C96 so that the user can check if its build is correc.

Clara Draper and Tseganeh Gichamo with inputs from Youlong Xia
June 30, 2021

