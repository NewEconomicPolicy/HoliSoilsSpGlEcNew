"""
#-------------------------------------------------------------------------------
# Name:        getClimGenFns.py
# Purpose:     additional functions for getClimGenNC.py
# Author:      s03mm5
# Created:     08/02/2018
# Copyright:   (c) s03mm5 2015
# Licence:     <your licence>
# -------------------------------------------------------------------------------
"""
__prog__ = 'getClimGenFns.py'
__author__ = 's03mm5'

from numpy.ma.core import MaskedConstant, MaskError
from netCDF4 import Dataset
from warnings import filterwarnings
from PyQt5.QtWidgets import QApplication

GRANULARITY = 120

ERROR_STR = '*** Error *** '
WARNING = '*** Warning *** '

def join_hist_fut_to_all_wthr(climgen, pettmp_hist, pettmp_fut):
    """
    join historic and future weather
     """
    fut_yr_strt = climgen.fut_wthr_set_defn['year_start']
    hist_yr_end = climgen.hist_wthr_set_defn['year_end']
    overlap_yrs = hist_yr_end - fut_yr_strt + 1
    if overlap_yrs < 0:
        mess = ERROR_STR + 'historic and future weather datasets are not contiguous or do not overlap'
        print(mess + '  - cannot proceed')
        QApplication.processEvents()
        return

    indx_hist_abbrev = overlap_yrs*12

    pettmp_all = {}
    for metric in pettmp_fut:
        pettmp_all[metric] = {}


    for gran_coord in pettmp_hist['precipitation']:

        if gran_coord not in pettmp_fut['precipitation']:
            lat, lon = pettmp_hist['lat_lons'][gran_coord]
            mess = WARNING + 'granular coordinate {} with lat: {}\tlong: {}'.format(gran_coord, lat, lon)
            print(mess + ' not present in future weather')
            QApplication.processEvents()
            continue

        for metric in pettmp_hist:
            if metric == 'lat_lons':
                continue

            if pettmp_hist[metric][gran_coord] is None:
                continue

            hist_seg = pettmp_hist[metric][gran_coord][:-indx_hist_abbrev]
            fut_seg = pettmp_fut[metric][gran_coord][:]

            pettmp_all[metric][gran_coord] = hist_seg + fut_seg

    return pettmp_all

def fetch_wthr_dset_overlap(wthr_set1, wthr_set2):
    """
    return overlap of two weather sets
    """
    lon_ll = max(wthr_set1['lon_ll'], wthr_set2['lon_ll'])
    lat_ll = max(wthr_set1['lat_ll'], wthr_set2['lat_ll'])
    lon_ur = min(wthr_set1['lon_ur'], wthr_set2['lon_ur'])
    lat_ur = min(wthr_set1['lat_ur'], wthr_set2['lat_ur'])

    return (lon_ll, lat_ll, lon_ur, lat_ur)

def genLocalGrid(wthr_set, bbox_wthr, bbox_aoi):
    """
    return the weather indices for the area which encloses the supplied bounding box
    this function does not alter the ClimGenNC (self) object
    """
    # junk = seterr(all='ignore') # switch off warning messages

    lon_aoi_ll, lat_aoi_ll, lon_aoi_ur, lat_aoi_ur = bbox_aoi
    lon_wthr_ll, lat_wthr_ll, lon_wthr_ur, lat_wthr_ur = bbox_wthr

    # determine bounds for climate grid which will enclose the greatest area of the bounding box 
    # ==========================================================================================
    lon_ll = max(lon_wthr_ll, lon_aoi_ll)
    lat_ll = max(lat_wthr_ll, lat_aoi_ll)
    lon_ur = min(lon_wthr_ur, lon_aoi_ur)
    lat_ur = min(lat_wthr_ur, lat_aoi_ur)

    lat_ur_indx, lon_ur_indx = get_wthr_nc_coords(wthr_set, lat_ur, lon_ur)
    lat_ll_indx, lon_ll_indx = get_wthr_nc_coords(wthr_set, lat_ll, lon_ll)

    """
    istep_lat = _coord_indices(lat_ll_indx,lat_ur_indx)
    lat_indices = [indx for indx in range(lat_ll_indx, lat_ur_indx, istep_lat)]

    istep_lon = _coord_indices(lon_ll_indx, lon_ur_indx)
    lon_indices = [indx for indx in range(lon_ll_indx, lon_ur_indx, istep_lon)]
    """

    lat_indx_min, lat_indx_max = _coord_order(lat_ll_indx,lat_ur_indx)
    lon_indx_min, lon_indx_max = _coord_order(lon_ll_indx,lon_ur_indx)

    return lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max  # aoi_indices

def _coord_order(ll_indx, ur_indx):
    """
    C
    """
    if ur_indx > ll_indx:
        max_indx = ur_indx
        min_indx = ll_indx
    else:
        max_indx = ll_indx
        min_indx = ur_indx

    return (min_indx, max_indx)

def _coord_indices(ll_indx, ur_indx):
    """
    C
    """
    if ur_indx > ll_indx:
        istep = 1
    else:
        istep = -1

    return istep

def _apply_start_year_correction(sim_strt_yr, hist_dset_defn, pettmp):
    """
    assume monthly datasets
    check and, if necessary, correct situation where sim_strt_yr is before historic dataset start year
    """
    repeat_period = (hist_dset_defn['year_start'] - sim_strt_yr)*12
    if repeat_period <= 0:
        return pettmp

    new_pettmp = {}
    for metric in pettmp:
        new_pettmp[metric] = pettmp[metric][0:repeat_period] + pettmp[metric]

    return new_pettmp

def _fetch_wthrset_indices(wthr_set_defn, sim_strt_yr, sim_end_yr):
    """
    get indices for simulation years for monthly weather set
    """
    wthr_yr_strt = wthr_set_defn['year_start']
    wthr_yr_end = wthr_set_defn['year_end']

    # simulation end year is before start of this dataset - nothing to do
    # ===================================================================
    if sim_end_yr < wthr_yr_strt:
        return 3*[None]

    # simulation start year is after end of this dataset - nothing to do
    # ===================================================================
    if wthr_yr_strt > sim_end_yr:
        return 3 * [None]

    indx_strt = max(0, (sim_strt_yr - wthr_yr_strt)*12)

    if sim_end_yr >= wthr_yr_end:

        # simulation end year is in future and beyond this dataset end year
        # =================================================================
        indx_end = -1
        next_strt_yr = wthr_yr_end + 1
    else:
        # simulation end year is before this dataset end year
        # ===================================================
        indx_end = (sim_end_yr - wthr_yr_end)*12
        next_strt_yr = -1

    return indx_strt, indx_end, next_strt_yr

def open_wthr_NC_sets(climgen):
    """
    C
    """
    hist_wthr_dsets = {}
    fut_wthr_dsets = {}

    for metric, ds_fname in zip(list(['precip', 'tas']), list(['ds_precip', 'ds_tas'])):
        hist_wthr_dsets[metric] = Dataset(climgen.hist_wthr_set_defn[ds_fname])
        fut_wthr_dsets[metric] = Dataset(climgen.fut_wthr_set_defn[ds_fname])

    return hist_wthr_dsets, fut_wthr_dsets

def fetch_WrldClim_data(lgr, lat, lon, climgen, nc_dsets, lat_indx, lon_indx, hist_flag=False):
    """
    C
    """
    filterwarnings("error")

    pettmp = {}
    for metric in list(['precip', 'tas']):
        if hist_flag:
            varname = climgen.hist_wthr_set_defn[metric]
            try:
                slice = nc_dsets[metric].variables[varname][:, lat_indx, lon_indx]
            except BaseException as err:
                print(ERROR_STR + str(err))
                return None
        else:
            varname = climgen.fut_wthr_set_defn[metric]
            slice = nc_dsets[metric].variables[varname][:, lat_indx, lon_indx]

        # test to see if cell data is valid, if not then this location is probably sea
        # =============================================================================
        if type(slice[0]) is MaskedConstant:
            pettmp = None
            mess = 'No data at lat: {} {}\tlon: {} {}\thist_flag: {}'.format(lat, lat_indx, lon, lon_indx, hist_flag)
            lgr.info(mess)
            print(mess)
        else:
            pettmp[metric] = [float(val) for val in slice]

    return pettmp

def get_wthr_nc_coords(dset_defn, latitude, longitude):
    """
    C
    """
    lon_frst = dset_defn['lon_frst']
    lat_frst = dset_defn['lat_frst']
    resol_lat = dset_defn['resol_lat']
    resol_lon = dset_defn['resol_lon']
    max_lat_indx = len(dset_defn['latitudes']) - 1
    max_lon_indx = len(dset_defn['longitudes']) - 1

    lat_indx = int(round((latitude - lat_frst)/resol_lat))
    lon_indx = int(round((longitude - lon_frst)/resol_lon))

    if lat_indx < 0 or lat_indx > max_lat_indx:
        print(WARNING + 'latitude index {} out of bounds for latitude {}\tmax indx: {}'
                                                            .format(lat_indx, round(latitude, 4), max_lat_indx))
        return -1, -1

    if lon_indx < 0 or lon_indx > max_lon_indx:
        print(WARNING + 'longitude index {} out of bounds for longitude {}\tmax indx: {}'
                                                            .format(lon_indx, round(longitude, 4), max_lon_indx))
        return -1, -1

    return lat_indx, lon_indx

def check_clim_nc_limits(form, bbox_aoi = None, wthr_rsrce = 'CRU') :
    """
    this function checks that the specified bounding box lies within extent of the requested weather dataset
    """
    func_name =  __prog__ + ' check_clim_nc_limits'

    limits_ok_flag = True
    '''
    if hasattr(form, 'combo10w'):
        wthr_rsrce = form.combo10w.currentText()
    
    if wthr_rsrce == 'NASA' or wthr_rsrce == 'CRU':
        return limits_ok_flag
    '''
    lon_ll_aoi = float(form.w_ll_lon.text())
    lat_ll_aoi = float(form.w_ll_lat.text())
    lon_ur_aoi = float(form.w_ur_lon.text())
    lat_ur_aoi = float(form.w_ur_lat.text())

    wthr_rsrce = wthr_rsrce + '_hist'      # was + '_Day'
    lat_ur_dset = form.wthr_sets[wthr_rsrce]['lat_ur']
    lon_ur_dset = form.wthr_sets[wthr_rsrce]['lon_ur']
    lat_ll_dset = form.wthr_sets[wthr_rsrce]['lat_ll']
    lon_ll_dset = form.wthr_sets[wthr_rsrce]['lon_ll']

    # similar functionality in lu_extract_fns.py in LU_extract project
    # ================================================================
    if (lon_ll_dset < lon_ll_aoi and lon_ur_dset > lon_ur_aoi) and \
                    (lat_ll_dset < lat_ll_aoi and lat_ur_dset > lat_ur_aoi):
        print('AOI lies within ' + wthr_rsrce + ' weather datasets')
    else:
        print('AOI lies outwith ' + wthr_rsrce + ' weather datasets - LL long/lat: {} {}\tUR long/lat: {} {}'
              .format(lon_ll_dset, lat_ll_dset, lon_ur_dset, lat_ur_dset))
        limits_ok_flag = False

    return limits_ok_flag
