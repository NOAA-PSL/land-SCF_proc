#!/usr/bin/env python3
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import numpy as np
import argparse
import glob
import os

def plot_world_map(lons, lats, data, metadata, plotpath):
    # plot generic world map
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))
    ax.add_feature(cfeature.GSHHSFeature(scale='auto'))
    ax.set_extent([-180, 180, -90, 90])
    vmin = 0.0 #vmin = np.nanmin(data)
    vmax = np.nanmax(data)
    cmap = 'viridis'
    cbarlabel = '%s' % (metadata['var'])
    plttitle = 'JEDI FV3 variable %s by %s' % (metadata['var'],os.environ['LOGNAME'])
    for t in range(0,6):
        cs = ax.pcolormesh(lons[...,t], lats[...,t], data[...,t],vmin=vmin,vmax=vmax,cmap=cmap)
    cb = plt.colorbar(cs, orientation='horizontal', shrink=0.5, pad=.04)
    cb.set_label(cbarlabel, fontsize=12)
    plt.title(plttitle)
    plt.savefig(plotpath)
    plt.close('all')

def read_var(datapath, geopath, varname):
    # read in data from 6 tile files for global DA
    tilefiles = glob.glob(datapath+'*')
    geotilefiles = glob.glob(geopath+'*')
    # get resolution to create the empty output array
    tmpdata = nc.Dataset(geotilefiles[0],'r')
    tmplat = tmpdata.variables['geolat'][:]
    tmpdata.close()
    arrayshape = tmplat.shape + (6,)
    dataout = np.empty(arrayshape)
    lonout = np.empty(arrayshape)
    latout = np.empty(arrayshape)
    for t in range(1,7):
        tilestr = 'tile%d.nc' % (t)
        datafile = [i for i in tilefiles if tilestr in i][0]
        geofile = [i for i in geotilefiles if tilestr in i][0]
        geonc = nc.Dataset(geofile)
        datanc = nc.Dataset(datafile)
        lats = geonc.variables['geolat'][:]
        lons = geonc.variables['geolon'][:]
        data = datanc.variables[varname][0, ...]
        if t in [3,4,5]:
            lats = np.rot90(lats)
            lons = np.rot90(lons)
            data = np.rot90(data)
        latout[:,:,t-1] = lats
        lonout[:,:,t-1] = lons
        dataout[:,:,t-1] = data
        geonc.close()
        datanc.close()
    return dataout, lonout, latout


def gen_figure(inpath, geopath, outpath, varname):
    # read the files to get the 2D array to plot
    data, lons, lats = read_var(inpath, geopath, varname)
    plotpath = outpath+'/%s.png' % (varname)
    metadata = {'var': varname,
                }
    plot_world_map(lons, lats, data, metadata, plotpath)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--output', help="path to output directory", default="./")
    ap.add_argument('-i', '--input', help="path to prefix of input files (no tileN.nc)", required=True)
    ap.add_argument('-v', '--variable', help="variable name to plot", required=True)
    ap.add_argument('-g', '--geoin', help="path to prefix of input files with geolat/geolon", required=True)
    MyArgs = ap.parse_args()
    gen_figure(MyArgs.input, MyArgs.geoin, MyArgs.output, MyArgs.variable)
