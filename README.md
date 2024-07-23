Code for processing IMS or VIIRS (based on user selection) input files on the UFS model grid. 

Inputs: IMS ascii file, IMS index file (generated offline).
        IMS index file is  model resolution specific.
        VIIRS h5 file, VIIRS index file (generated offline).
        VIIRS index file is model resolution specific.

Output: file with i) IMS(VIIRS)-derived snow cover fraction over land on UFS model grid, and ii) snow depth derived form the IMS(VIIRS) snow cover fraction, using an inversion of the noah model snow depletion curve.
        

To compile:
1.ln -s $your_GDAS_App_installation*
2 ./build.sh
On hera, can use:
*/scratch2/NCEPDEV/land/data/DA/GDASApp

See ./test for test case (hera only)

Clara Draper, Tseganeh Gichamo, with input from Youlong Xia. June, 2021.
Yuan Xue: add VIIRS processing capability and merge with original IMS'. July 2024.


