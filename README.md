Code for processing IMS input ascii files on the UFS model grid. 

Inputs: IMS ascii file, IMS index file (generated offline).
        IMS index file is  model resolution specific.

Output: file with i) IMS-derived snow cover fraction over land on UFS model grid, and ii) snow depth derived form the IMS snow cover fraction, using an inversion of the noah model snow depletion curve.
        

To compile: 
./build.sh

Test case on hera: /scratch2/BMC/gsienkf/Clara.Draper/DA_test_cases/snow/fIMS
(2024 Feb 19: has not been maintained, may not work)


Clara Draper, Tseganeh Gichamo, with input from Youlong Xia. June, 2021.


