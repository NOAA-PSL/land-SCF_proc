Code for processing IMS input ascii files on the UFS model grid.

Inputs: IMS ascii file, IMS index file (generated offline).
        IMS index file is  model resolution specific.

Output: file with i) IMS-derived snow cover fraction over land on UFS model grid, and ii) snow depth derived form the IMS snow cover fraction, using an inversion of the noah model snow depletion curve.
        

To compile: 
ln -s $your_GDAS_App_installation**
./build.sh

On hera, can use: 

**/scratch2/NCEPDEV/land/data/DA/GDASApp

See ./test for test case (hera only).

Clara Draper, Tseganeh Gichamo, with input from Youlong Xia. June, 2021.


