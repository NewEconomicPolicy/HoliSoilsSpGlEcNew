# -------------------------------------------------------------------------------
# Name:        litter_and_efiscen_fns.py
# Purpose:     script to create objects describing NC data sets
# Author:      Mike Martin
# Created:     31/05/2020
# Licence:     <your licence>
# -------------------------------------------------------------------------------

__prog__ = 'litter_and_efiscen_fns.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os.path import normpath, exists
from numpy import ma, zeros
from netCDF4 import Dataset
from PyQt5.QtWidgets import QApplication

ERROR_STR = '*** Error *** '
WARN_STR = '*** Warning *** '

EFISCEN_CARBON_VAR = 'TOTAL_BM_LITTER_c'
EFISCEN_MANDAT_VARS = ('lat', 'lon', EFISCEN_CARBON_VAR)       # EFISCEN | European Forest Institute

def fetch_nc_litter(form, fname):
    """
    currently permit only a single cell
    """
    if not exists(fname):
        mess = WARN_STR + 'ORCHIDEE or EFISCEN NetCDF litter file name ' + fname
        if fname.isspace() or fname == '':
            mess += ' must not be blank'
        else:
            mess += '  does not exist'

        print(mess)
        return None

    if check_efiscen_dset(fname):
        form.w_create_files.setEnabled(True)
    else:
        form.w_create_files.setEnabled(False)
        form.w_nc_extnt.setText('')
        form.w_ave_val.setText('')
        return None

    carbon_var = form.combo08.currentText()
    lttr_defn = EfiscenSet(fname)
    form.litter_defn = lttr_defn
    form.w_nc_extnt.setText(lttr_defn.nc_extnt)

    # reset PFT combo widget
    # ======================
    form.w_combo_pfts.clear()

    form.combo08.setCurrentText(EFISCEN_CARBON_VAR)
    form.combo08.setEnabled(False)
    form.w_var_desc.setText(form.carbon_vars[carbon_var])
    form.w_combo_pfts.addItem('total biomass litter C')

    # report average value
    # ====================
    pft_name = form.w_combo_pfts.currentText()
    if pft_name == '':
        mess = 'No data for ' + carbon_var
        form.w_ave_val.setText(mess)
    else:
        pft_key = '00'
        ave_val = form.litter_defn.aves[pft_key]
        mess = 'average value: ' + str(round(float(ave_val), 2))

    form.w_ave_val.setText(mess)

    return None

def check_efiscen_dset(nc_fname):
    """
    C
    """
    nc_fname = normpath(nc_fname)
    nc_dset = Dataset(nc_fname)

    # check for EFISCEN dataset
    # =========================
    efiscen_flag = True
    vars_not_prsnt = []
    for var in EFISCEN_MANDAT_VARS:
        if var not in nc_dset.variables:
            vars_not_prsnt.append(var)

    if len(vars_not_prsnt) > 0 or len(nc_dset.variables[EFISCEN_CARBON_VAR].shape) != 3:
        efiscen_flag = False

    nc_dset.close()

    if efiscen_flag:
        print('Dataset ' + nc_fname + ' is identified as EFISCEN')
    else:
        print(ERROR_STR + 'invalid dataset, could not be identified as EFISCEN')

    return efiscen_flag

class EfiscenSet(object, ):
    """
    Create object from an EFISCEN NC file
    """
    def __init__(self, nc_fname):
        """
        assumption is that dataset has been pre-checked using check_efiscen_dset function
        """
        lat_var = 'lat'
        lon_var = 'lon'

        nc_fname = normpath(nc_fname)

        nc_dset = Dataset(nc_fname)
        print('\nReading EFISCEN file ' + nc_fname)

        lats = nc_dset.variables[lat_var][:]
        nlats = len(lats)
        lons = nc_dset.variables[lon_var][:]
        nlons = len(lons)

        dset_type = 'EFISCEN'
        nyears = -999

        # record var names
        # ================
        exclude_vars = list([lat_var, lon_var, 'time'])
        start_year = None
        end_year = None
        var_names = []
        for var in nc_dset.variables:
            if var not in exclude_vars:
                var_names.append(var)

            # stanza to get start of time series
            # ==================================
            if var == 'time':
                time_var = nc_dset.variables[var]
                nyears = len(time_var)
                start_year = time_var[0]
                end_year = time_var[-1]

        lat_frst = round(float(lats[0]), 4)
        lon_frst = round(float(lons[0]), 4)
        lat_last = round(float(lats[-1]), 4)
        lon_last = round(float(lons[-1]), 4)

        # required for bounding box
        # =========================
        if lat_last > lat_frst:
            lat_ll = lat_frst
            lat_ur = lat_last
        else:
            lat_ll = lat_last
            lat_ur = lat_frst

        if lon_last > lon_frst:
            lon_ll = lon_frst
            lon_ur = lon_last
        else:
            lon_ll = lon_last
            lon_ur = lon_frst

        self.dset_type = dset_type
        self.lat_frst = lat_frst
        self.lon_frst = lon_frst
        self.lat_last = lat_last
        self.lon_last = lon_last

        self.lat_var = lat_var
        self.lon_var = lon_var
        self.bbox = lon_ll, lat_ll, lon_ur, lat_ur

        self.nc_fname = nc_fname
        self.var_names = var_names
        self.nc_dset = None

        # resolutions
        # ===========
        self.resol_lon = round((lons[-1] - lons[0])/(nlons - 1), 4)
        self.resol_lat = round((lats[-1] - lats[0])/(nlats - 1), 4)
        self.max_lat_indx = nlats - 1
        self.max_lon_indx = nlons - 1

        #
        self.lats = list(lats)
        self.lons = list(lons)

        extent_lats = 'N latitudes: {}   extent: {} {}\t'.format(nlats, lat_frst, lat_last)
        extent_lons = 'N longitudes: {}   extent: {} {}\t'.format(nlons, lon_frst, lon_last)
        time_span = 'Start year: {}   End year: {}\t'.format(start_year, end_year)
        grid_resol_str = 'grid resolution: {}'.format(self.resol_lat)
        self.nc_extnt = extent_lats + extent_lons + time_span + grid_resol_str

        # Create a boolean table of cells with and without data
        # =====================================================
        vals, aves, lookup_table = fetch_efiscen_vals(nc_dset, nlats, nlons)

        nc_dset.close()

        self.lookup_table = lookup_table
        self.vals = vals
        self.aves = aves
        self.nyears = nyears
        self.start_year = start_year
        self.end_year = end_year

    def get_efiscen_nc_data(self, pft_key, lat, long, baseline_flag):
        """
        retrieve data on condition that lat, long are within bounds and data is present
        """
        wthn_bnds = True

        lat_indx = int(round((lat - self.lat_frst)/self.resol_lat))
        lon_indx = int(round((long - self.lon_frst)/self.resol_lon))

        if lat_indx < 0 or lat_indx > self.max_lat_indx:
            wthn_bnds = False
            print(WARN_STR + 'latitude index {} out of bounds for latitude {}\tmax indx: {}'.format(lat_indx,
                                                                            round(lat, 4), self.max_lat_indx))

        if lon_indx < 0 or lon_indx > self.max_lon_indx:
            wthn_bnds = False
            print(WARN_STR + 'longitude index {} out of bounds for longitude {}\tmax indx: {}'.format(lon_indx,
                                                                            round(long, 4), self.max_lon_indx))

        if wthn_bnds and self.lookup_table[lat_indx][lon_indx]:
            plnt_inpts = {'yrs': [], 'pis': []}
            strt_yr = self.start_year
            nyears = self.nyears
            if pft_key is None:
                vals = nyears * [0]
            else:
                slice = self.vals[pft_key][lat_indx, lon_indx, :]
                is_masked = ma.is_masked(slice)
                if is_masked:
                    vals = ma.getdata(slice)
                else:
                    vals = list(slice)

            plnt_inpts['yrs'] = [yr for yr in range(strt_yr, strt_yr + nyears)]
            if baseline_flag:
                plnt_inpts['pis'] = nyears*[0]
            else:
                plnt_inpts['pis'] = [val for val in vals]
        else:
            plnt_inpts = None

        return plnt_inpts

def resize_yrs_pi(sim_strt_yr, sim_end_yr, yrs_pi):
    """
    patch to enable adjust yrs_pi to correspond to user specified simulation period
    """
    if yrs_pi is None:
        return None

    yr_frst = yrs_pi['yrs'][0]
    yr_last = yrs_pi['yrs'][-1]

    pi_frst = yrs_pi['pis'][0]
    pi_last = yrs_pi['pis'][-1]

    # trap error
    # ==========
    frst_yrs = []
    frst_pis = []
    if sim_strt_yr < yr_frst:
        frst_yrs = list(range(sim_strt_yr, yr_frst))
        frst_pis = [pi_frst for yr in range(sim_strt_yr, yr_frst)]

    last_yrs = []
    last_pis = []

    if sim_end_yr > yr_last:
        last_yrs = list(range(yr_last + 1, sim_end_yr))
        last_pis = [pi_last for yr in range(yr_last + 1, sim_end_yr)]

    sim_yrs = frst_yrs + yrs_pi['yrs'] + last_yrs
    sim_pis = frst_pis + yrs_pi['pis'] + last_pis

    new_yrs_pi = {'yrs': sim_yrs, 'pis': sim_pis}

    return new_yrs_pi

def fetch_efiscen_vals(nc_dset, nlats, nlons):
    """
    Read plant inputs from EFISCEN dataset and create a lookup table
    """
    vals = {}
    aves = {}
    pft_key = '00'
    lookup_table = zeros((nlats, nlons), dtype=bool)

    tmp_vals = nc_dset.variables[EFISCEN_CARBON_VAR][:, :, :]
    vals[pft_key] = tmp_vals
    aves[pft_key] = tmp_vals.mean()

    for lat_indx in range(nlats):
        for lon_indx in range(nlons):

            val_mean = tmp_vals[lat_indx, lon_indx, :].mean()
            if val_mean > 0.0:
                lookup_table[lat_indx][lon_indx] = True
            else:
                lookup_table[lat_indx][lon_indx] = False

    return vals, aves, lookup_table
