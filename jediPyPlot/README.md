This is a python plot script to plot FV3-tile fields. 

On Hera machine (module use and load):

module use -a /scratch2/NCEPDEV/marineda/Jong.Kim/save/modulefiles/
module load anaconda/3.15.1

Command line:

python plot_jedi_tiles.py -i inpout -t type -v variableName -l level -g geographicalData

Examples:

python plot_jedi_tiles.py -i ../test/output/IMSfSCA. -t landSurface -v IMSfsca -l 1 -g /scratch1/NCEPDEV/global/glopara/fix/fix_fv3_gmted2010/C96/C96_oro_data

python plot_surface_jedi_tiles.py -i ../test/output/IMSfSCA. -v IMSfsca -g /scratch1/NCEPDEV/global/glopara/fix/fix_fv3_gmted2010/C96/C96_oro_data

If people want to make a difference plot, the following steps need to be done:

module load NCO
ncdiff test1.tile1.nc test2.tile1.nc diff.tile1.nc (do all 6 tile files manually)

And then used either python plot_jedi_tiles.py or python plot_surface_jedi_tiles.py diff.tile1.nc to make plots.

An example shell script diffData.sh can be used more easily to make the data difference.

Youlong Xia, Tseganeh Gichamo, and Clara Draper
June 30, 2021


