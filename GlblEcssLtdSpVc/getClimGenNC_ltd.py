# -------------------------------------------------------------------------------
# Name:        getClimGenNC.py
# Purpose:     read netCDF files comprising ClimGen data
# Author:      s03mm5
# Created:     08/12/2015
# Copyright:   (c) s03mm5 2015
# Licence:     <your licence>
# -------------------------------------------------------------------------------
#
__prog__ = 'getClimGenNC.py'
__author__ = 's03mm5'

# Version history
# ---------------
#
from os.path import normpath, isfile, join, split, isdir
from os import makedirs
from calendar import monthrange
from netCDF4 import Dataset
from math import ceil
from numpy import arange, seterr, ma
from time import time
import warnings
import csv
from pandas import read_csv
from PyQt5.QtWidgets import QApplication

from prepare_ecosse_files import update_progress
from getClimGenFns import fetch_days_per_month
from thornthwaite import thornthwaite

null_value = -9999
set_spacer_len = 12
MAX_CELLS = 9999999

NUMSECSDAY = 3600*24

GRANULARITY = 120
weather_resource_permitted = list(['CRU', 'EObs', 'EWEMBI'])

WARNING = '*** Warning *** '
ERROR_STR = '*** Error *** '

def _consistency_check(pettmp, varnams_mapped):
    """
    make sure for a each key if one metric is zero length then all other metrics for that key are also blank
    TODO: this function only works for two metrics and is unpythonic!
    """
    metric_list = list(varnams_mapped.values())
    metric0 = metric_list[0]
    metric1 = metric_list[1]
    for key in pettmp[metric0]:
        len_key0 = len(pettmp[metric0][key])

        if len_key0 == 0:
            pettmp[metric1][key] = []

        len_key1 = len(pettmp[metric1][key])
        if len_key1 == 0:
            pettmp[metric0][key] = []

    return pettmp

def _check_list_for_none(metric_list):
    """
    if a None is found then return an empty list
    """
    for indx, val in enumerate(metric_list):
        if val is None:
            return []

    return metric_list

def _input_txt_line_layout(data, comment):
    """
    C
    """
    spacer_len = max(set_spacer_len - len(data), 2)
    spacer = ' ' * spacer_len
    return '{}{}# {}\n'.format(data, spacer, comment)

class ClimGenNC(object,):
    """
    C
    """
    def __init__(self, form, fut_clim_scen=None):

        func_name =  __prog__ +  ' ClimGenNC __init__'

        if hasattr(form, 'w_mnthly'):
            if form.w_mnthly.isChecked():  # monthly timestep
                sim_mnthly_flag = True
            else:
                sim_mnthly_flag = False  # daily timestep
        else:
            sim_mnthly_flag = True

            # determine user choices
        # ======================
        if hasattr(form, 'combo10w'):
            wthr_rsrce = form.combo10w.currentText()
            if fut_clim_scen is None:
                fut_clim_scen = form.combo10.currentText()
            hist_start_year = int(form.combo09s.currentText())
            hist_end_year = int(form.combo09e.currentText())
            ave_weather_flag = form.w_ave_weather.isChecked()
            sim_start_year = int(form.combo11s.currentText())   # might need these later
            sim_end_year = int(form.combo11e.currentText())
        else:
            wthr_rsrce = form.weather_resource
            fut_clim_scen = form.scenario
            hist_start_year = form.hist_strt_year
            hist_end_year = form.hist_end_year
            ave_weather_flag = form.ave_weather_flag
            sim_start_year = form.sim_strt_year
            sim_end_year = form.sim_end_year

        self.weather_resource = wthr_rsrce
        if hasattr(form, 'sims_dir'):
            wthr_out_dir = join(split(form.sims_dir)[0], 'Wthr', fut_clim_scen)
            if not isdir(wthr_out_dir):
                makedirs(wthr_out_dir)
        else:
            wthr_out_dir = None

        self.wthr_out_dir = wthr_out_dir
        self.coords_lookup = None
        self.sim_mnthly_flag = sim_mnthly_flag

        # African Monsoon Multidisciplinary Analysis (AMMA) 2050 datasets
        # ===============================================================
        if wthr_rsrce in form.amma_2050_allowed_gcms:
            wthr_set_key = wthr_rsrce + '_' + fut_clim_scen
            if wthr_set_key not in form.weather_sets:
                print('key {} not in weather sets in function {} - cannot continue'.format(wthr_set_key, func_name))
                return
            hist_wthr_set = form.weather_sets[wthr_rsrce + '_historical']
            fut_wthr_set = form.weather_sets[wthr_set_key]
            lat = 'lat'
            lon = 'lon'
        elif wthr_rsrce== 'HARMONIE':
            hist_wthr_set = form.weather_sets['HARMONIE_V2']
            fut_wthr_set = form.weather_sets['HARMONIE_V2']
            lat = 'lat'
            lon = 'lon'
        elif wthr_rsrce == 'EObs':
            hist_wthr_set = form.weather_sets['EObs_Mnth']
            fut_wthr_set = form.weather_sets['EObs_Mnth']
            lat = 'latitude'
            lon = 'longitude'
        elif wthr_rsrce == 'CRU':
            wthr_set_key = 'ClimGen_' + fut_clim_scen
            if wthr_set_key not in form.weather_sets:
                print('key {} not in weather sets in function {} - cannot continue'.format(wthr_set_key, func_name))
                return
            hist_wthr_set = form.weather_sets['CRU_hist']
            fut_wthr_set = form.weather_sets[wthr_set_key]
            lat = 'latitude'
            lon = 'longitude'
        elif wthr_rsrce == 'EFISCEN-ISIMIP':
            wthr_set_key = 'EFISCEN-ISIMIP_' + fut_clim_scen
            if wthr_set_key not in form.weather_sets:
                print('key {} not in weather sets in function {} - cannot continue'.format(wthr_set_key, func_name))
                return
            hist_wthr_set = form.weather_sets['CRU_hist']
            fut_wthr_set = form.weather_sets[wthr_set_key]
            lat = 'lat'
            lon = 'lon'
        elif wthr_rsrce == 'NCAR_CCSM4':
            hist_wthr_set = form.weather_sets['NCAR_CCSM4']
            fut_wthr_set = form.weather_sets['NCAR_CCSM4']
            lat = 'lat'
            lon = 'lon'
        else:
            print('weather resource ' + wthr_rsrce + ' not recognised in ' + func_name + ' - cannot continue')
            return

        # ===============================================================
        self.fut_wthr_set_defn = fut_wthr_set
        self.hist_wthr_set_defn = hist_wthr_set

        # make sure start and end years are within dataset limits
        # =======================================================
        hist_start_year = max(hist_wthr_set['year_start'], hist_start_year)
        hist_end_year = min(hist_wthr_set['year_end'], hist_end_year)

        self.ave_weather_flag = ave_weather_flag
        num_hist_years = hist_end_year - hist_start_year + 1
        self.num_hist_years = num_hist_years
        self.hist_start_year = hist_start_year
        self.hist_end_year = hist_end_year
        self.months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        self.fut_clim_scen = fut_clim_scen

        # future data
        # ===========
        self.fut_precip_fname = fut_wthr_set['ds_precip']
        self.fut_tas_fname = fut_wthr_set['ds_tas']
        self.resol_fut_lon = fut_wthr_set['resol_lon']
        self.resol_fut_lat = fut_wthr_set['resol_lat']

        # granularity
        # ===========
        self.lon = lon
        self.lat = lat
        self.lon_min = fut_wthr_set['lon_ll']
        self.lat_min = fut_wthr_set['lat_ll']
        self.longitudes = fut_wthr_set['longitudes']
        self.latitudes = fut_wthr_set['latitudes']
        self.pettmp = {}        # dictionary whose keys will reference the climate grid, pt_grid
        self.lgr = form.lgr

        # past (monthly) data
        # ===================
        self.hist_precip_fname = hist_wthr_set['ds_precip']
        self.hist_tas_fname = hist_wthr_set['ds_tas']
        self.longitudes_hist = hist_wthr_set['longitudes']
        self.latitudes_hist = hist_wthr_set['latitudes']

        # New stanza to facilitate option when user selects "use average weather"
        # =======================================================================
        if num_hist_years == 1:
            self.met_ave_file = 'met' + str(hist_start_year) + 'a.txt'
        else:
            self.met_ave_file = 'met' + str(hist_start_year) + '_' + str(hist_end_year) + 'a.txt'

        self.num_ave_wthr_years = num_hist_years    # only years when there is historic data will be taken into account

        self.sim_start_year = sim_start_year
        self.sim_end_year = sim_end_year
        self.num_sim_years = sim_end_year - sim_start_year + 1

        self.ave_file = 'met' + str(hist_start_year) + '_' + str(hist_end_year) + 'a.txt'

        self.ave_sim_text_fn = 'met{}_to_{}_ave.txt'.format(sim_start_year, sim_end_year)
        self.ave_sim_met_fn = 'met' + str(sim_start_year) + '_' + str(sim_end_year) + 'a.txt'

    def readCoordsLookup(self):
        """
        C
        """
        lookup_flag = False
        if isdir(self.wthr_out_dir):
            lookup_fn = join(self.wthr_out_dir, 'coords_lookup.csv')
            if isfile(lookup_fn):
                df = read_csv(lookup_fn, sep='\t')
                coords_lookup = df.to_dict('list')
                lookup_flag = True

        if lookup_flag:
            self.coords_lookup = coords_lookup
            print('coords lookup file {} has {} records'.format(lookup_fn, len(coords_lookup)))
        else:
            print(WARNING + 'no coords lookup file found - cannot proceed')

        QApplication.processEvents()
        return lookup_flag

    def genLocalGrid(self, bbox, hwsd, snglPntFlag = False):
        """
        return the weather indices for the area which encloses the supplied bounding box
        this function does not alter the ClimGenNC (self) object
        """
        # junk = seterr(all='ignore') # switch off warning messages

        bb_lon_min, bb_lat_min, bb_lon_max, bb_lat_max =  bbox
        if snglPntFlag:
            bb_lon_max = bb_lon_min
            bb_lat_max = bb_lat_min

        # determine bounds for climate grid which will enclose the supplied bounding box
        # ==============================================================================
        resol_fut_lat = self.resol_fut_lat   # negative for future CRU data
        lat_indices_fut = []
        lat_indices_hist = []
        clim_lat_min = self.lat_min
        num_lats = ceil( abs((bb_lat_max - clim_lat_min)/resol_fut_lat) )
        lat_max = round(abs(num_lats*resol_fut_lat) + clim_lat_min, 8)   # rounding introduced for NCAR_CCSM4
        lat_indices_fut.append(self.latitudes.index(lat_max))
        lat_indices_hist.append(self.latitudes_hist.index(lat_max))

        num_lats = int( abs((bb_lat_min - clim_lat_min)/resol_fut_lat) )
        lat_min = round(abs(num_lats*resol_fut_lat) + clim_lat_min, 8)   # rounding introduced for NCAR_CCSM4
        lat_indices_fut.append(self.latitudes.index(lat_min))
        lat_indices_hist.append(self.latitudes_hist.index(lat_min))

        # longitudes
        # ==========
        lon_indices_fut = []
        lon_indices_hist = []
        resol_fut_lon = self.resol_fut_lon
        clim_lon_min = self.lon_min
        num_lons = ceil((bb_lon_max - clim_lon_min)/resol_fut_lon)
        lon_max = round(num_lons*resol_fut_lon + clim_lon_min, 8)
        lon_indices_fut.append(self.longitudes.index(lon_max))
        lon_indices_hist.append(self.longitudes_hist.index(lon_max))

        num_lons = int((bb_lon_min - clim_lon_min)/resol_fut_lon)
        lon_min = round(num_lons*resol_fut_lon + clim_lon_min, 8)
        lon_indices_fut.append(self.longitudes.index(lon_min))
        lon_indices_hist.append(self.longitudes_hist.index(lon_min))
        """
        # generate ClimGen grid    NB need to add one when slicing!!!
        # =====================    ==================================
        alons = arange(lon_min, lon_max, resol_fut_lon)
        alats = arange(lat_min, lat_max, resol_fut_lat)
        nlats = len(alats)
        nlons = len(alons)
        
        granlons = arange(nlons)
        for ic, lon in enumerate(alons):
            granlons[ic] = (180.0 + lon)*hwsd.granularity
        granlons.sort()

        granlats = arange(nlats)
        for ic, lat in enumerate(alats):
            granlats[ic] = (90.0 - lat)*hwsd.granularity
        granlats.sort()
        """
        # must be in correct order
        # ========================
        lat_indices_hist.sort()
        lon_indices_hist.sort()
        lat_indices_fut.sort()
        lon_indices_fut.sort()

        aoi_indices_fut = lat_indices_fut + lon_indices_fut
        aoi_indices_hist = lat_indices_hist + lon_indices_hist
        return aoi_indices_fut, aoi_indices_hist

    def fetch_isimip_NC_data(self, aoi_indices, dset_strt_yr, nmnths, max_cells=MAX_CELLS):
        """
        fetch precipitation - units: kg m-2 s-1, and temperature - units: Kelvin, for a given lat/long AOI
            for ISIMAP precipitation - convert kg m-2 s-1 to mm month-1
            1 kg water = 1000 mm-3      1 m-2 = 1 million mm-2 so 1 kg m-2 = 1000 / 1000000  = 0.001 mm
            so to convert to mm day-1 = 0.001 * NUMSECSDAY = 86.4
        """
        cnvrt_isimip_pr = 1.0 * NUMSECSDAY

        # warnings.simplefilter('default')
        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices

        # read all times
        # ==============
        days_per_month = fetch_days_per_month(dset_strt_yr, nmnths)
        varnams_mapped = {'pr':'precipitation','tas':'temperature'}
        varnams = sorted(varnams_mapped.keys())

        icount = 0
        pettmp = {}
        pettmp['lat_lons'] = {}
        for varname, fname in zip(varnams, list([self.fut_precip_fname, self.fut_tas_fname])):
            num_key_masked = 0
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            nc_dset = Dataset(fname, mode='r')

            # collect readings for all time values
            # ====================================
            slice = nc_dset.variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]

            if ma.is_masked(slice):
                slice_is_masked_flag = True
                self.lgr.info('Future slice is masked')
            else:
                slice_is_masked_flag = False

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                lat = self.latitudes[lat_indx]
                gran_lat = round((90.0 - lat)*GRANULARITY)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    lon = self.longitudes[lon_indx]
                    gran_lon = round((180.0 + lon)*GRANULARITY)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = null_value
                    if slice_is_masked_flag :
                        val = slice[0, ilat, ilon]
                        if val is ma.masked:
                            self.lgr.info('val is ma.masked for key ' + key)
                            pettmp[varnam_map][key] = None
                            num_key_masked += 1

                    # add data for this coordinate
                    # ============================
                    if pettmp[varnam_map][key] == null_value:
                        if varname == 'tas':
                            pettmp[varnam_map][key] = [round(float(val) - 273.15, 2) for val in slice[:, ilat, ilon]]
                        else:
                            # precipitation - convert kg m-2 s-1 to mm month-1
                            # ================================================
                            # pettmp[varnam_map][key] =  = [round(val * ndays * cnvrt_isimip_pr, 1) for val, ndays in
                            #                                                zip(slice[:, ilat, ilon], days_per_month)]
                            pettmp[varnam_map][key] = []
                            for val, ndays in zip(slice[:, ilat, ilon], days_per_month):
                                # val_mm = round(val * ndays * cnvrt_isimip_pr, 3)
                                val_mm = round(val * cnvrt_isimip_pr, 2)
                                pettmp[varnam_map][key].append(val_mm)

                    pettmp['lat_lons'][key] = [lat, lon]
                    icount += 1
                    if icount >= max_cells:
                        break

                if icount >= max_cells:
                    break

            nc_dset.close()     # close netCDF file
            if num_key_masked > 0:
                print('# masked weather keys: {}'.format(num_key_masked))
                QApplication.processEvents()

        return pettmp

    def fetch_ewembi_NC_data(self, aoi_indices, num_band, future_flag = True):
        """
        get precipitation or temperature data for a given variable and lat/long index for all times
        CRU uses NETCDF4 format
        """
        func_name = __prog__ +  ' fetch_fut_future_NC_data'
        warnings.simplefilter('default')

        num_key_masked = 0
        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}

        weather_resource = self.weather_resource
        if future_flag:
            precip_fname = self.fut_precip_fname
            tas_fname = self.fut_tas_fname
            start_year = self.sim_start_year
        else:
            precip_fname = self.hist_precip_fname
            tas_fname = self.hist_tas_fname
            start_year = self.hist_start_year

        varnams_mapped = {'pr': 'precipitation', 'tas': 'temperature'}
        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([precip_fname, tas_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            nc_dset = Dataset(fname)

            # collect readings for all time values
            # ====================================
            slice = nc_dset.variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]

            if ma.isMaskedArray(slice):
                slice_is_masked_flag = True
                self.lgr.info('Future slice is a masked array in band {}'.format(num_band))
            else:
                slice_is_masked_flag = False

            # generate days per month
            # ======================
            if varname == 'pr':
                days_per_month = []
                nmonths = len(nc_dset.variables[varname])
                for year in range(start_year, start_year + int(nmonths/12) + 1):
                    for imnth in range(12):
                        dummy, ndays = monthrange(year, imnth + 1)
                        days_per_month.append(ndays)

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                gran_lat = round((90.0 - self.latitudes[lat_indx])*GRANULARITY)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    gran_lon = round((180.0 + self.longitudes[lon_indx])*GRANULARITY)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = null_value
                    if slice_is_masked_flag :
                        val = slice[0,ilat,ilon]
                        if val is ma.masked:
                            self.lgr.info('val is ma.masked for key ' + key)
                            pettmp[varnam_map][key] = None
                            num_key_masked += 1

                    # add data for this coordinate
                    # ============================
                    if pettmp[varnam_map][key] == null_value:
                        if varname == 'pr':
                            pettmp[varnam_map][key] = [round(val*ndays*NUMSECSDAY, 1) for val, ndays in
                                                                                zip(slice[:,ilat,ilon], days_per_month)]
                        elif varname == 'tas':

                            # risk of masked constant
                            # =======================
                            try:
                                pettmp[varnam_map][key] = [round(float(val) - 273.15, 1) for val in slice[:,ilat,ilon]]
                            except UserWarning as e:
                                print('{}\n\tLatitude:{}\tLongitude {}'
                                      .format(e, self.latitudes[lat_indx], self.longitudes[lon_indx]))

            # close netCDF file
            nc_dset.close()
            print('# masked weather keys: {}'.format(num_key_masked))

        return pettmp

    def fetch_eobs_NC_data(self, aoi_indices, num_band, future_flag = True):
        """
        get precipitation or temperature data for a given variable and lat/long index for all times
        EObs uses NETCDF format
        """
        func_name = __prog__ +  ' fetch_eobs_NC_data'
        warnings.simplefilter('default')

        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}
        if future_flag == True:
            precip_fname = self.fut_precip_fname
            temper_fname = self.fut_tas_fname
        else:
            precip_fname = self.hist_precip_fname
            temper_fname = self.hist_tas_fname

        # process future climate
        # ======================
        varnams_mapped = {'rr':'precipitation','tg':'temperature'}
        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([precip_fname, temper_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            nc_dset = Dataset(fname, mode='r')

            # collect readings for all time values
            # ====================================
            try:
                slice = nc_dset.variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]
            except RuntimeWarning as e:
                print(e)

            if ma.is_masked(slice):
                slice_is_masked_flag = True
                self.lgr.info('Future slice is masked in band {}'.format(num_band))
            else:
                slice_is_masked_flag = False

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                lat = self.latitudes[lat_indx]
                gran_lat = round((90.0 - lat)*GRANULARITY)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    long = self.longitudes[lon_indx]
                    gran_lon = round((180.0 + long)*GRANULARITY)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = null_value
                    if slice_is_masked_flag :
                        val = slice[0,ilat,ilon]
                        if val is ma.masked:
                            self.lgr.info('val is ma.masked for key {} lat {} long {}'.format(key, lat, long))
                            pettmp[varnam_map][key] = []

                    # add data for this coordinate
                    # ============================
                    if pettmp[varnam_map][key] == null_value:
                        pettmp[varnam_map][key] = _check_list_for_none(slice[:,ilat,ilon].tolist())

            # close netCDF file
            nc_dset.close()

        pettmp = _consistency_check(pettmp, varnams_mapped)
        return pettmp

    def fetch_ncar_ccsm4_NC_data(self, aoi_indices, num_band, future_flag = True):
        """
        get precipitation or temperature data for a given variable and lat/long index for all times
        CORDEX uses NETCDF3_64BIT format
        """
        func_name = __prog__ +  ' fetch_ncar_ccsm4_NC_data'
        warnings.simplefilter('default')

        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}
        if future_flag == True:
            precip_fname = self.fut_precip_fname
            temper_fname = self.fut_tas_fname
            start_year   = self.sim_start_year
        else:
            precip_fname = self.hist_precip_fname
            temper_fname = self.hist_tas_fname
            start_year   = self.hist_start_year

        # process future climate
        # ======================
        varnams_mapped = {'pr':'precipitation','tas':'temperature'}
        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([precip_fname, temper_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            nc_dset = Dataset(fname, mode='r')

            # collect readings for all time values
            # ====================================
            try:
                slice = nc_dset.variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]
            except RuntimeWarning as e:
                print(e)

            if ma.is_masked(slice):
                slice_is_masked_flag = True
                print('Future slice is masked in band {}'.format(num_band))
            else:
                slice_is_masked_flag = False

            # generate days per month
            # ======================
            if varname == 'pr':
                days_per_month = []
                nmonths = len(nc_dset.variables[varname])
                for year in range(start_year, start_year + int(nmonths/12) + 1):
                    for imnth in range(12):
                        dummy, ndays = monthrange(year, imnth + 1)
                        days_per_month.append(ndays)

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                gran_lat = round((90.0 - self.latitudes[lat_indx])*GRANULARITY)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    gran_lon = round((180.0 + self.longitudes[lon_indx])*GRANULARITY)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = null_value
                    if slice_is_masked_flag :
                        val = slice[0,ilat,ilon]
                        if val is ma.masked:
                            self.lgr.info('val is ma.masked for key ' + key)
                            pettmp[varnam_map][key] = None

                    # add data for this coordinate
                    # ============================
                    if pettmp[varnam_map][key] == null_value:
                        if varname == 'pr':
                            pettmp[varnam_map][key] = [round(val*ndays*NUMSECSDAY, 1) for val in slice[:,ilat,ilon]]
                        elif varname == 'tas':
                            pettmp[varnam_map][key] = [round(val - 273.15, 1) for val in slice[:,ilat,ilon]]

            # close netCDF file
            nc_dset.close()

        return pettmp

    def fetch_harmonie_NC_data(self, aoi_indices, num_band, future_flag = True):
        """
        get precipitation or temperature data for a given variable and lat/long index for all times
        """
        func_name = __prog__ +  ' fetch_harmonie_NC_data'
        warnings.simplefilter('default')

        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}
        if future_flag == True:
            precip_fname = self.fut_precip_fname
            temper_fname = self.fut_tas_fname
            start_year   = self.sim_start_year
        else:
            precip_fname = self.hist_precip_fname
            temper_fname = self.hist_tas_fname
            start_year   = self.hist_start_year

        # process future climate
        # ======================
        varnams_mapped = {'Precipalign':'precipitation','Tairalign':'temperature'}
        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([precip_fname, temper_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            nc_dset = Dataset(fname, mode='r')

            # collect readings for all time values
            # ====================================
            try:
                slice = nc_dset.variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]
            except RuntimeWarning as e:
                print(e)

            if ma.is_masked(slice):
                slice_is_masked_flag = True
                print('Future slice is masked in band {}'.format(num_band))
            else:
                slice_is_masked_flag = False

            # generate days per month
            # ======================
            if varname == 'Precipalign':
                days_per_month = []
                nmonths = len(nc_dset.variables[varname])
                for year in range(start_year, start_year + int(nmonths/12) + 1):
                    for imnth in range(12):
                        dummy, ndays = monthrange(year, imnth + 1)
                        days_per_month.append(ndays)

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                gran_lat = round((90.0 - self.latitudes[lat_indx])*GRANULARITY)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    gran_lon = round((180.0 + self.longitudes[lon_indx])*GRANULARITY)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = null_value
                    if slice_is_masked_flag :
                        val = slice[0,ilat,ilon]
                        if val is ma.masked:
                            self.lgr.info('val is ma.masked for key ' + key)
                            pettmp[varnam_map][key] = None

                    # add data for this coordinate
                    # ============================
                    if pettmp[varnam_map][key] == null_value:
                        if varname == 'Precipalign':
                            pettmp[varnam_map][key] = [round(val, 2) for val in slice[:,ilat,ilon]]
                        elif varname == 'Tairalign':
                            pettmp[varnam_map][key] = [round(val - 273.15, 1) for val in slice[:,ilat,ilon]]

            # close netCDF file
            nc_dset.close()

        return pettmp

    def fetch_cru_future_NC_data(self, aoi_indices, num_band, fut_start_indx=0):
        """
        get precipitation or temperature data for a given variable and lat/long index for all times
        CRU uses NETCDF4 format
        """
        warnings.simplefilter('default')

        num_key_masked = 0
        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}

        # process future climate
        # ======================
        varnams_mapped = {'precipitation':'precipitation','temperature':'temperature'}

        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([self.fut_precip_fname, self.fut_tas_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            nc_dset = Dataset(fname)

            # collect readings for all time values
            # ====================================
            slice = nc_dset.variables[varname][lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1, :]

            if ma.is_masked(slice):
                slice_is_masked_flag = True
                self.lgr.info('Future slice is masked in band {}'.format(num_band))
            else:
                slice_is_masked_flag = False

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                gran_lat = round((90.0 - self.latitudes[lat_indx])*GRANULARITY)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    gran_lon = round((180.0 + self.longitudes[lon_indx])*GRANULARITY)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = null_value
                    if slice_is_masked_flag :
                        val = slice[ilat,ilon,0]
                        if val is ma.masked:
                            self.lgr.info('val is ma.masked for key ' + key)
                            pettmp[varnam_map][key] = None
                            num_key_masked += 1

                    # add data for this coordinate
                    # ============================
                    if pettmp[varnam_map][key] == null_value:
                        # remove overlap with historic data - for CRU data only
                        record = [round(val, 1) for val in slice[ilat,ilon,:]]
                        pettmp[varnam_map][key] = record[fut_start_indx:]

            # close netCDF file
            nc_dset.close()
            if num_key_masked > 0:
                print('# masked weather keys: {}'.format(num_key_masked))

        return pettmp

    def fetch_cru_historic_NC_data(self, aoi_indices, num_band, max_cells=MAX_CELLS):
        """
        get precipitation or temperature data for a given variable and lat/long index for all times
        CRU uses NETCDF4 format
        """
        warnings.simplefilter('default')

        num_key_masked = 0
        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}
        pettmp['lat_lons'] = {}
        nlats = lat_indx_max - lat_indx_min + 1

        # process historic climate
        # ========================
        varnams_mapped = {'pre': 'precipitation', 'tmp': 'temperature'}

        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([self.hist_precip_fname, self.hist_tas_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            nc_dset = Dataset(fname)

            # for user feedback
            # =================
            last_time = time()
            start_time = None
            completed_lats = 0
            skipped = 0
            warn_count = 0

            # collect readings for all time values
            # ====================================
            slice = nc_dset.variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]

            if ma.is_masked(slice):
                slice_is_masked_flag = True
                self.lgr.info('Historic weather slice is masked in band {}'.format(num_band))
            else:
                slice_is_masked_flag = False

            # reform slice
            # ============
            icount = 0
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                lat = self.latitudes_hist[lat_indx]
                gran_lat = round((90.0 - lat)*GRANULARITY)

                last_time = update_progress(last_time, start_time, completed_lats, nlats, skipped, warn_count)
                QApplication.processEvents()

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    lon = self.longitudes_hist[lon_indx]
                    gran_lon = round((180.0 + lon)*GRANULARITY)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = null_value
                    if slice_is_masked_flag :
                        val = slice[0, ilat, ilon]
                        if val is ma.masked:
                            self.lgr.info('val is ma.masked for key ' + key)
                            pettmp[varnam_map][key] = None
                            num_key_masked += 1

                    # add data for this coordinate
                    # ============================
                    if pettmp[varnam_map][key] == null_value:
                        pettmp[varnam_map][key] = [round(val, 1) for val in slice[:, ilat,ilon]]

                    pettmp['lat_lons'][key] = [lat, lon]
                    icount += 1
                    if icount >= max_cells:
                        break

                completed_lats += 1
                if icount >= max_cells:
                    break

            nc_dset.close()     # close netCDF file
            if num_key_masked > 0:
                print('# masked weather keys: {}'.format(num_key_masked))
                QApplication.processEvents()

        return pettmp

    def create_FutureAverages(self, clim_dir, lat_inp, gran_coord, site, lta_precip, lta_tmean):
        """
        use prexisting metyyyys.txt files to generate a text file of average weather which will subsequently
        be included in the input.txt file
        also create a climate file for each of the simulation years based on average weather from the CRU year range
        """
        func_name = ' create_FutureAverages'
        full_func_name = __prog__ + func_name

        sim_start_year = self.sim_start_year
        sim_end_year = self.sim_end_year
        num_sim_yrs = sim_end_year - sim_start_year + 1
        months = self.months

        # skip if already exists
        ave_sim_text_fn = join(normpath(clim_dir), self.ave_sim_text_fn)
        ave_sim_met_fn = join(normpath(clim_dir), self.ave_sim_met_fn)
        if isfile(ave_sim_text_fn) and isfile(ave_sim_met_fn):
            mess = 'Files:\n\t' + ave_sim_met_fn + ' and ' + ave_sim_text_fn + '\n\talready exist - will overwrite'
            self.lgr.info(WARNING + mess)

        # read precipitation and temperature
        sim_precip = {}
        sim_tmean = {}
        for month in months:
            sim_precip[month] = 0.0
            sim_tmean[month] = 0.0

        for year in range(sim_start_year, sim_end_year + 1):
            fname = 'met{0}s.txt'.format(year)
            met_fpath = join(clim_dir, fname)

            if not isfile(met_fpath):
                print('File ' + met_fpath + ' does not exist - will abandon average weather creation')
                return -1

            with open(met_fpath, 'r', newline='') as fpmet:
                lines = fpmet.readlines()

            for line, month in zip(lines, months):
                tlst = line.split('\t')
                sim_precip[month] += float(tlst[1])
                sim_tmean[month] += float(tlst[3].rstrip('\r\n'))

        # write stanza for input.txt file consisting of long term average climate
        # =======================================================================
        output = []
        num_sim_years = self.num_sim_years
        for month in self.months:
            ave_precip = sim_precip[month]/num_sim_years
            output.append(_input_txt_line_layout('{}'.format(round(ave_precip, 1)),
                                                '{} long term average monthly precipitation [mm]'.format(month)))

        for month in self.months:
            ave_tmean = sim_tmean[month]/num_sim_years
            output.append(_input_txt_line_layout('{}'.format(round(ave_tmean, 2)),
                                                '{} long term average monthly temperature [degC]'.format(month)))

        # write text file of average simulated weather which will subsequently be included in the input.txt file
        # ======================================================================================================
        try:
            fhand = open(ave_sim_text_fn, 'w')
        except IOError:
            raise IOError(ERROR_STR + 'Unable to open file: ' + ave_sim_text_fn)
        else:
            fhand.writelines(output)
            fhand.close()

        self.lgr.info('Successfully wrote average weather file {} in function {}'.format(ave_sim_text_fn, func_name))

        # write long term average climate file
        # ====================================
        ave_precip = [round(sim_precip[month]/num_sim_years, 1) for month in months]
        ave_tmean = [round(sim_tmean[month]/num_sim_years, 1) for month in months]

        # pet
        if max(ave_tmean) > 0.0:
            pet = thornthwaite(ave_tmean, lat_inp, year)
        else:
            pet = [0.0]*12
            mess = WARNING + 'monthly average temperatures are below zero in ' + full_func_name
            mess += ' for latitude: {}\tgranular coord: {}'.format(lat_inp, gran_coord)
            print(mess)

        pot_evapotrans = [round(p, 1) for p in pet]

        # write file
        output = []
        for tstep, mean_temp in enumerate(ave_tmean):
            output.append([tstep+1, ave_precip[tstep], pot_evapotrans[tstep], mean_temp])

        with open(ave_sim_met_fn, 'w', newline='') as fpout:
            writer = csv.writer(fpout, delimiter='\t')
            writer.writerows(output)
            fpout.close()

        # note float conversion from float32 otherwise rounding does not work as expected
        # ===============================================================================
        if site is not None:
            lta = {'pet': [], 'precip': lta_precip, 'tas': lta_tmean}
            lta['pet'] = thornthwaite(lta['tas'], lat_inp, year)

            site.lta_pet = [round(float(pet), 1) for pet in lta['pet']]
            site.lta_precip = [round(float(precip), 1) for precip in lta['precip']]
            site.lta_tmean = [round(float(tmean), 1) for tmean in lta['tas']]

        self.lgr.info('Successfully wrote average weather file {} in function {}'.format(ave_sim_met_fn, func_name))

        return 0