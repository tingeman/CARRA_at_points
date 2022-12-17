import pathlib
import numpy as np
from netCDF4 import Dataset
import pandas as pd


datapath = pathlib.Path('../3h')
coords_file = pathlib.Path('../site_coords/TIN_Greenland_ccordinates_for_CARRA.csv')
data_freq = '3H'

netcdf_path = pathlib.Path('../netcdf')


# Ensure path exist
netcdf_path.mkdir(parents=True, exist_ok=True)


# Get sites information
df = pd.read_csv(coords_file)
num_sites = len(df)
site_names = list(df['name'])


# Get timerange of data in npy files
files = list(datapath.glob('*.npy'))
year_list = []
for fname in files:
    year_list.append(int(fname.stem))
year_list = sorted(year_list)

all_times = pd.date_range(start='{0:4.0f}-01-01'.format(year_list[0]),
                          end='{0:4.0f}-01-01'.format(year_list[-1]+1),
                          freq=data_freq, inclusive='right')
# Using "inclusive='right'" excludes the first time stamp of the series, but includes the last.
# So effectively the series starts at 03:00h and ends with last stamp at 00:00h



for site_num, site in enumerate(site_names):
    print('Processing CARRA data for site: {0}'.format(site))

    # Getting grid coordinates and elevations
    lat_lon_elev_path = pathlib.Path('../TIN_lat_lon_elev/lat_lon_elev_{0}.npy'.format(site))
    with open(lat_lon_elev_path, 'rb') as f:
        lat = np.load(f)
        lon = np.load(f)
        elev = np.load(f)

    print('Generating netcdf-file...')

    # set up netcdf data storage
    out_file = netcdf_path / 'CARRA_data_{0}.nc'.format(site)
    rootgrp = Dataset(out_file, "w", format="NETCDF4")

    # Create dimensions
    idy = rootgrp.createDimension("idy", lat.shape[0])
    idx = rootgrp.createDimension("idx", lat.shape[1])
    time = rootgrp.createDimension("time", len(all_times))

    # Create variables
    times = rootgrp.createVariable("time","f8",("time",))
    idys = rootgrp.createVariable("idy","i4",("idy",))
    idxs = rootgrp.createVariable("idx","i4",("idx",))

    t2m = rootgrp.createVariable("t2m","f4",("time","idy","idx",))
    t2m.units = "K"
    rf = rootgrp.createVariable("rf","f4",("time","idy","idx",))
    rf.units = "mm"
    tp = rootgrp.createVariable("tp","f4",("time","idy","idx",))
    tp.units = "mm"

    latitudes = rootgrp.createVariable("latitude","f4",("idy","idx",))
    longitudes = rootgrp.createVariable("longitude","f4",("idy","idx",))
    elevations = rootgrp.createVariable("elevation","f4",("idy","idx",))

    # populate data
    longitudes[:] = lon
    latitudes[:] = lat
    elevations[:] = elev

    idys[:] = np.arange(lat.shape[0])
    idxs[:] = np.arange(lat.shape[1])

    times[:] = all_times


    for fname in files:
        print('Processing data-file: {0}'.format(fname))
        
        yr = int(fname.stem)

        # read data file for this particular year
        with open(fname, 'rb') as f:
            rf_3h_5dcube=np.load(f)
            tp_3h_5dcube=np.load(f)
            t2m_3h_5dcube=np.load(f)

        # Generate timestamps
        times = pd.date_range(start='{0:4.0f}-01-01'.format(yr),
                        end='{0:4.0f}-01-01'.format(yr+1),
                        freq=data_freq, inclusive='right')

        # Find indices of where to place the data in the full timeseries
        idt = np.logical_and(all_times >= times[0], all_times <= times[-1])

        # loop over all grid cells, and add full year of data to netcdf storage
        for iy in range(lat.shape[1]):
            for ix in range(lat.shape[1]):
                t2m[idt, iy, ix] = t2m_3h_5dcube[site_num, :, :, iy, ix].T.reshape(-1)
                rf[idt, iy, ix] = rf_3h_5dcube[site_num, :, :, iy, ix].T.reshape(-1)
                tp[idt, iy, ix] = tp_3h_5dcube[site_num, :, :, iy, ix].T.reshape(-1)

    rootgrp.close()
    print('Finished processing site!')
    print('')


# data = Dataset('CARRA_data_sisimiut_2.nc')

# dat_t2m = data['t2m'][:,0,0]
# plt.plot(all_times, dat_t2m)

# dat_rf = data['rf'][:,0,0]
# plt.plot(all_times, dat_rf)

# dat_sf = data['tp'][:,0,0] - data['rf'][:,0,0]
# plt.plot(all_times, dat_sf)




