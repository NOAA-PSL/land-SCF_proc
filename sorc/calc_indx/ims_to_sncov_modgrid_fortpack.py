# Calls the function: imsgridmap.select_imscells_within_tile_grid
# to generate files with ims grid indices the fall within a model (FV3) grid cell
# imsgridmap compiled with  python -m numpy.f2py -c -m imsgridmap imsgridmap.f90
# copy the files .pyd and .dll(from imsgridmap/.libs) to dir accessible by this source. The file names are similar to:
#       imsgridmap.cp38-win_amd64.pyd
#       libimsgridm.RKR3EIZJMKB7DQUFN5BDV4SB3QQQQSX3.gfortran-win_amd64.dll

# inputs (all them in netcdf file):
#     1. fv3 grid file with coordinates of the bourndaries of the model grid
#     2. oro_files wit the lat/lon coordinates of the target fv3 file
#     3. one ims file, to read the ims coordinates from
#  output: for each tile, a netcdf file with the list of IMS grid indices


from netCDF4 import Dataset
from my_mapping_functions import *

try:
    from osgeo import gdal, osr, ogr
except:
    import gdal, osr, ogr
import os
from helperFunctions import callSubprocess_ps
import imsgridmap   # compiled with  python -m numpy.f2py -c -m imsgridmap imsgridmap.f90

res = 96    #48   #768 #
max_dx = 50

n_lon = res
n_lat = res
n_til = 6
ny = 2*n_lat + 1
nx = 2*n_lon + 1
nout = n_lat * n_lon

num_sub = 3840  #42600  # 1920   #7800

n_lat_ims = 2176
n_lon_ims = 8703

working_dir = "C:/jedi-vm/vagrant_data/C" + str(res) + '/'
latlon_dir = "C:/jedi-vm/vagrant_data/C" + str(res) + '/'
os.chdir(working_dir)

lat_grid = np.full((n_til, n_lat, n_lon), np.nan)
lon_grid = np.full((n_til, n_lat, n_lon), np.nan)
y_grid = np.full((n_til, ny, nx), np.nan)
x_grid = np.full((n_til, ny, nx), np.nan)
area_grid = np.full((n_til, ny-1, nx-1), np.nan)
data_grid_ims = np.full((n_lat_ims, n_lon_ims), np.nan)

lat_grid_ims_1 = np.full((n_lat_ims), np.nan)
lon_grid_ims_ = np.full((n_lon_ims), np.nan)
lon_grid_ims_1 = np.full((n_lon_ims), np.nan)
lat_grid_ims = np.full((n_lat_ims, n_lon_ims), np.nan)
lon_grid_ims = np.full((n_lat_ims, n_lon_ims), np.nan)

# Read coords at the boundaris/faces of the grid cell
for t in np.arange(n_til):
    ncid = Dataset(latlon_dir+'C'+str(res)+'_grid.tile'+str(t+1)+'.nc')
    y_grid[t,:,:] = ncid['y'][:]
    x_grid[t,:,:] = ncid['x'][:]
    area_grid[t,:,:] = ncid['area'][:]
    ncid.close()
print("min area = ", np.min(area_grid), " max area = ", np.max(area_grid))
# lat/lon coordinates of the destination grid
for t in np.arange(n_til):
    ncid = Dataset(latlon_dir+'C'+str(res)+'_oro_data.tile'+str(t+1)+'.nc')
    lat_grid[t,:,:] = ncid['geolat'][:]
    lon_grid[t,:,:] = ncid['geolon'][:]
    ncid.close()

# cmdString = "gdalwarp -t_srs EPSG:4326  -te -180.0 0.0 180.0 90.0 -r near  -overwrite -of netCDF ims2019319_4km_GIS_v1.3.tif wgs_ims2019319_4km_GIS_v1.3.tif.nc"
# #"gdalwarp -t_srs EPSG:4326 -te -180.0 0.0 180.0 90.0 -tr 0.036 0.036 -r near  -overwrite -of netCDF -co \"FORMAT=NC4\" ims2019319_4km_GIS_v1.3.tif wgs_ims2019319_4km_GIS_v1.3.tif.nc"
# callSubprocess_ps(cmdString)
#Windows
#for %i in (ims*_4km_GIS_v1.3.tif) do (gdalwarp -t_srs EPSG:4326 -te -180.0 0.0 180.0 90.0 -tr 0.036 0.036 -r near  -overwrite -of netCDF -co "FORMAT=NC4" %i wgs_%i.nc)
input_nc = 'wgs_ims2019319_4km_GIS_v1.3.tif.nc'
ncOut = Dataset(input_nc, "r+")  # , format='NETCDF4')
lat_grid_ims_1 = ncOut['lat'][:]
lon_grid_ims_ = ncOut['lon'][:]
# data_grid_ims[:,:] = ncOut['Band1'][:,:]
ncOut.close()
# print(lat_grid_ims_1)
# print(lon_grid_ims_1)
# print(data_grid_ims)
# exit(2)

# print(lon_grid_ims_)
# lon_grid_ims_[lon_grid_ims_ < 0.] = lon_grid_ims_[lon_grid_ims_ < 0.] + 360.
# print(lon_grid_ims_)
for ilon in range(n_lon_ims):
    if lon_grid_ims_[ilon] < 0.:
        lon_grid_ims_1[ilon] = lon_grid_ims_[ilon] + 360.
    else:
        lon_grid_ims_1[ilon] = lon_grid_ims_[ilon]
#lon_grid_ims_1[lon_grid_ims_1 < 0.] = lon_grid_ims_1[lon_grid_ims_1 < 0.] + 360.

for til in [0, 1, 2, 3, 4, 5]: #range(n_til):
    file_in = 'sfc_data.tile' + str(til + 1) + '.nc'
    file_out = "C"+str(res)+".IMS.Indices.tile" + str(til + 1) + '.nc'
    cmdString = "nccopy -V xaxis_1,yaxis_1,zaxis_1,Time " + file_in + " " + file_out
    callSubprocess_ps(cmdString)
    print("tile = "+str(til + 1))
    lat_min = np.amin(lat_grid[til]); lat_max = np.amax(lat_grid[til])
    lon_min = np.amin(lon_grid[til]); lon_max = np.amax(lon_grid[til])
    print("min/max lat ", lat_min, lat_max)
    print("min/max lon ", lon_min, lon_max)
    #continue
    print("min/max lat ", np.nanmin(y_grid[til]), np.nanmax(y_grid[til]))
    print("min/max lon ", np.nanmin(x_grid[til]), np.nanmax(x_grid[til]))

    for ilon in range(n_lon_ims):
        lat_grid_ims[:, ilon] = lat_grid_ims_1[:]
    for ilat in range(n_lat_ims):
        lon_grid_ims[ilat, :] = lon_grid_ims_1[:]

    data_grid = imsgridmap.select_imscells_within_tile_grid(max_dx, lat_grid_ims.flatten(), lon_grid_ims.flatten(),
                                                   y_grid[til], x_grid[til], n_lat, n_lon, num_sub) #ny, nx,
    print(data_grid.shape)
    ncOut = Dataset(file_out, "r+")  # , format='NETCDF4')
    ncOut.createDimension("Indices", num_sub)  # n_tim)
    ncOut.createDimension("lat", n_lat)  # n_tim)
    ncOut.createDimension("lon", n_lon)  # n_tim)
    ncOut.createVariable("IMS_Indices", np.int,
                         ("lat", "lon", "Indices",))  # double sncovr(Time, yaxis_1, xaxis_1) ;
    attDict = {'long_name': 'Indices of corr. IMS', 'units': "-"}
    ncOut.variables['IMS_Indices'].setncatts(attDict)
    ncOut.variables['IMS_Indices'][:, :, :] = data_grid[:]
    # delete shleg
    ncOut.close()

