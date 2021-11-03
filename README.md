Code for processing IMS input ascii files on the UFS model grid. 

Inputs: IMS ascii file, IMS index file (generated offline).
        IMS index file is  model resolution specific.

Output: file with i) IMS-derived snow cover fraction over land on UFS model grid, and ii) snow depth derived form the IMS snow cover fraction, using an inversion of the noah model snow depletion curve.
        

To compile: 
./build.sh

To run, see example in test/

Clara Draper. June, 2021.


