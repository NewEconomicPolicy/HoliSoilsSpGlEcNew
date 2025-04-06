"""
#-------------------------------------------------------------------------------
# Name:        glbl_ecsse_low_level_fns.py
# Purpose:     consist of low level functions invoked by high level module
# Author:      Mike Martin
# Created:     11/12/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#
"""

__prog__ = 'glbl_ecsse_low_level_fns.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

from os.path import exists, join, split, isfile
from json import load as json_load
from pandas import read_csv
from pandas.errors import EmptyDataError
from sys import stdout
from time import time
from netCDF4 import Dataset
from pandas import DataFrame, Series
from glob import glob
from numpy.ma.core import MaskError
from locale import LC_ALL, setlocale, format_string
import numpy
from PyQt5.QtWidgets import QApplication

GRANULARITY = 120
SOIL_DIR = 'soil_metrics'

WARN_STR = '*** Warning *** '
ERROR_STR = '*** Error *** '

def fetch_coord_nearest_xy(coords_lookup, y_point, x_point):
    """
    C
    """
    y_array = numpy.array(coords_lookup['Lat'])
    x_array = numpy.array(coords_lookup['Lon'])

    distance = (y_array - y_point)**2 + (x_array - x_point)**2
    id = numpy.where(distance==distance.min())
    indx = id[0].item()

    gran_coord = coords_lookup['gran_coord'][indx]

    return gran_coord

def check_cell_within_csv(hwsd_mu_globals, lat, lon):

    gran_lat = round((90.0 - lat) * GRANULARITY)
    gran_lon = round((180.0 + lon) * GRANULARITY)
    aoi_chunk = hwsd_mu_globals.loc[(hwsd_mu_globals['gran_lat'] == gran_lat)
                                    & (hwsd_mu_globals['gran_lon'] == gran_lon)]
    if len(aoi_chunk) == 0:
        inside = False
    else:
        inside = True

    return inside

def set_region_study(form):

    study = form.w_study.text()
    study = study.replace(' ', '_')

    region = form.w_combo00a.currentText()
    region_indx = form.w_combo00a.currentIndex()
    crop_name = form.w_combo00b.currentText()

    region_study = (region + '_' + crop_name).replace(' ', '_') + '_' + study   # obsolete
    form.setup['region_study'] = study
    form.setup['region_wthr_dir'] = form.regions['Wthr dir'][region_indx]

    return

def check_rotation_json_fname(form):
    """
    Crop Rotation file comprises a start year, crop names and corresponding lookup codes from the CROP_SUN.DAT file
    for example:
    {
       "CropRotation": {
            "crops": { "Spring Wheat": 5, "Setaside": 12, "Linseed":19 },
            "start_year": 120
            }
    }
    """
    rota_json_fname = form.w_lbl16.text()

    if not exists(rota_json_fname):
        return 'Crop rotation json file does not exist'
    try:
        with open(rota_json_fname, 'r') as fcultiv:
            rota_content = json_load(fcultiv)
            print('Read crop rotation json input file ' + rota_json_fname)

    except (OSError, IOError) as err:
            print(err)
            return 'Could not read crop rotation json file'

    # check file contents
    # ===================
    rota_key = 'CropRotation'
    try:
        crop_pattern = rota_content[rota_key]['crops']
        start_year = rota_content[rota_key]['start_year']
        rota_pattern = rota_content
    except KeyError as bad_key:
        rota_pattern = None
        dummy, short_fname = split(rota_json_fname)
        mess = 'Key {} not recognised in {}'.format(bad_key, short_fname)

    # verify crops
    # ===========
    if rota_pattern is not None:
        mess = 'Crop rotation starts at year {} and comprises: '.format(start_year)
        for crop in crop_pattern:
            if crop in form.crop_defns:
                mess += crop + ', '
            else:
                mess = crop + ' not found in CROP_SUN.DAT file'
                rota_pattern = None
                break

    if rota_pattern is None:
        form.rota_pattern = None
    else:
        form.rota_pattern = rota_pattern['CropRotation']

    return mess.rstrip(', ')

def check_cultiv_json_fname(form):
    """
    Cultivation file comprises three variables for each year: cultivation, vigour and residues incorporated
    see page 53 on ECOSSE manual: cultivation is integer 0 to 3
                                 vigour determines proportion of humus released to biomass float 0 to 1
                                 residues incorporated can be 0 = no, 1 = yes
    Example:
        {
           "Cultivation": {
            "FrstYr": [3, 0.5, 0],
            "ChngYr": [0, 0.5, 1]
          }
        }
    """
    form.cultiv_pattern = None
    cultiv_json_fname = form.w_lbl13.text()

    if not exists(cultiv_json_fname):
        return 'Cultivation json file does not exist'

    try:
        with open(cultiv_json_fname, 'r') as fcultiv:
            cultiv_content = json_load(fcultiv)
            print('Read cultivation json input file ' + cultiv_json_fname)

    except (OSError, IOError) as err:
            print(err)
            return 'Could not read cultivation json file'

    # check file contents
    # ===================
    cultiv_key = 'Cultivation'
    if cultiv_key in cultiv_content.keys():
        cultiv_pattern = cultiv_content[cultiv_key]
        mess = 'Cultivation input file is valid'
    else:
        mess = WARN_STR + 'Key ' + cultiv_key + ' must be present in Cultivation file'
        cultiv_pattern = None

    form.cultiv_pattern = cultiv_pattern

    return mess

def Cell_hwsd_data_frame(lggr, hwsd):
    """
    set values within the HWSD grid to zero if they fall outside the shapefile polygon
    Argument description:
       hwsd:   HWSD object comprising a grid of mu_globals
    """
    func_name =  __prog__ + '  dump_AOI'

    total_number_cells = hwsd.nlats*hwsd.nlons
    data_frame = DataFrame()

    lggr.info('Function: {}\ttotal number of cells: {}'.format(func_name, total_number_cells))

    # now take each soil coord (integers)
    # move north to south and west to east (decreasing latitude and increasing longitude)
    # build list of mu_globals with their coordinates (granular and lat lon)
    num_zeros_before = 0    # from HWSD
    num_zeros_added = 0     # counter for points outside of main boundary
    num_points_ok = 0   # points inside
    num_boundary_bad = 0   # boundary violations

    gran_lats = []
    gran_lons = []
    mu_globals = []
    latitudes = []
    longitudes = []
    land_uses = []
    for iyhws, irow in enumerate(range(hwsd.nrow1, hwsd.nrow2 + 1)):

        for ixhws, icol in enumerate(range(hwsd.ncol1, hwsd.ncol2 + 1)):

            # skip if mu_global is 0 (usually sea)
            # ====================================
            mu_global = hwsd.rows[iyhws][ixhws]
            if mu_global == 0:
                num_zeros_before += 1
            elif mu_global in hwsd.bad_muglobals:
                hwsd.rows[iyhws][ixhws] = 0
                num_zeros_added += 1
            else:
                latitude = 90.0 - irow/hwsd.granularity
                longitude = icol/hwsd.granularity - 180.0

                # build lists - TODO
                # ==================
                gran_lats.append(irow)
                gran_lons.append(icol)
                mu_globals.append(mu_global)
                latitudes.append(latitude)
                longitudes.append(longitude)
                land_uses.append(-999)
                num_points_ok += 1

    # create data frame
    # =================
    data_frame['gran_lat']  = Series(gran_lats)
    data_frame['gran_lon']  = Series(gran_lons)
    data_frame['mu_global'] = Series(mu_globals)
    data_frame['latitude']  = Series(latitudes)
    data_frame['longitude'] = Series(longitudes)
    data_frame['land_use']  = Series(land_uses)

    return data_frame

def _bbox_locate(cntry_bboxes, lat, lon):
    """
    locate lat lon from dictionary of countries or provinces and their respective bounding boxes
    """
    cntry_list = []
    for cntry in cntry_bboxes:
        iso, lon_ll, lat_ll, lon_ur, lat_ur = cntry_bboxes[cntry][:5]
        if (lat >= lat_ll and lat <= lat_ur) and (lon >= lon_ll and lon <= lon_ur):
            cntry_list.append(cntry)

    ncntrys = len(cntry_list)
    if ncntrys == 0:
        cntry_list.append(None)
    elif ncntrys >= 2:
        for cntry in ['Russia', 'United States', 'United States Minor Outlying Islands']:
            try:
                cntry_list.remove(cntry)
            except ValueError:
                pass

    return cntry_list[0]

def _major_states_lookup(lggr, glbl_n_inpts, prvnc_bboxes, cntry, lat, lon):
    """
    locate province using lat/lon bounding boxes
    """
    glbl_amt = None
    mess = 'Major country: ' + cntry + '\tLat/lon: {} {}'.format(lat, lon)
    found_flag = False

    for prvnce in prvnc_bboxes[cntry]:
        lon_ll, lat_ll, lon_ur, lat_ur, area = prvnc_bboxes[cntry][prvnce]
        if (lat >= lat_ll and lat <= lat_ur) and (lon >= lon_ll and lon <= lon_ur):

            # identify province
            # =================
            mess += '\tcould not locate province: ' + prvnce + ' in lookup table'
            for rec_prvnc in glbl_n_inpts[cntry]:
                if prvnce == rec_prvnc[0]:
                    glbl_amt = rec_prvnc[2]
                    found_flag = True
                    break

        if found_flag:
            break

    return glbl_amt, found_flag, mess

def _fetch_glbl_amnt(lggr, n_inpts_obj, glbl_n_flag, lat, lon):
    """
    n_inpts_obj - object consisting of all components of global N application
    identify country from netCDF countries file with fall back using bounding boxes
    once a coordinate has been assigned a country or province then locate entry in N application from Excel file
    """
    glbl_amt = None
    if not glbl_n_flag:
        return glbl_amt

    lat_indx, lon_indx, ret_code = n_inpts_obj.cntries_defn.get_nc_coords(lat, lon)
    try:
        cntry_indx = int(n_inpts_obj.cntries_defn.nc_dset['countries'][lat_indx, lon_indx])
    except MaskError as err:
        mess = WARN_STR + 'lat/indx: {} {} long/indx: {} {} {}'.format(lat, lat_indx, lon, lon_indx, str(err))
        lggr.info(mess)

        # locate country using bounding boxes
        # ===================================
        cntry = _bbox_locate(n_inpts_obj.cntry_bboxes, lat, lon)
    else:
        try:
            cntry = n_inpts_obj.cntry_dict[cntry_indx]
        except KeyError as err:
            print('\n' + WARN_STR + 'No country found for lat/long: {} {} {}'.format(lat, lon, str(err)))
            return glbl_amt

    glbl_n_inpts = n_inpts_obj.glbl_n_inpts
    major_states = n_inpts_obj.major_states
    prvnc_bboxes = n_inpts_obj.prvnc_bboxes

    # set defaults in case of failure
    # ===============================
    found_flag = False
    mess = 'Country ' + cntry + ' in soils NC file not found in Excel World sheet'

    for rec in glbl_n_inpts['World']:
        if cntry == rec[0]:
            if cntry in major_states:
                glbl_amt, found_flag, mess = _major_states_lookup(lggr, glbl_n_inpts, prvnc_bboxes, cntry, lat, lon)
                break
            else:
                glbl_amt = rec[2]   # mean N Application
                found_flag = True
                break

    if not found_flag:
        lggr.info(mess)
        print('\n' + mess)

    return glbl_amt

def check_run_mask(mask_defn, lon_ll_indx, lat_ll_indx, lon_ur_indx, lat_ur_indx):
    """

    """
    setlocale(LC_ALL, '')

    lat_indx_min = lat_ll_indx
    lat_indx_max = lat_ur_indx + 1
    lon_indx_min = lon_ll_indx
    lon_indx_max = lon_ur_indx + 1

    nc_dset = Dataset(mask_defn.nc_fname, mode='r')

    # run_mask = nc_dset.variables['cropmask'][lat_indx_min:lat_indx_max, lon_indx_min:lon_indx_max]

    # adjust lat_ur_indx so that main loop starts at beginning of crop area
    # =====================================================================
    for lat_indx in range(lat_ur_indx, lat_ll_indx + 1, -1):
        arr = nc_dset.variables['cropmask'][lat_indx, lon_ll_indx:lon_ur_indx]
        num_non_zeros = len(arr.nonzero()[0])
        if num_non_zeros > 0:
            lat_ur_indx = lat_indx
            break

    # sum total number of active cells
    # ================================
    ngrow_cells = 0
    for lat_indx in range(lat_ur_indx, lat_ll_indx + 1, -1):
        arr = nc_dset.variables['cropmask'][lat_indx, lon_ll_indx:lon_ur_indx]
        num_non_zeros = len(arr.nonzero()[0])
        ngrow_cells += num_non_zeros

    start_lat = nc_dset.variables['lat'][lat_ur_indx]
    nc_dset.close()

    ntotal = (lon_indx_max - lon_indx_min)*(lat_indx_max - lat_indx_min)

    mess = 'Start of crop found at latitude: {}'.format(start_lat)
    ngrow_cells_str = format_string('%d', ngrow_cells, grouping=True)
    ntotal_str = format_string('%d', ntotal, grouping=True)
    mess += '\t\tN cells for crop growing: ' + ngrow_cells_str + '\tfrom total of ' + ntotal_str + ' cells (includes sea)'
    print(mess)

    return lat_ur_indx, ngrow_cells

def generate_cells(form):
    """

    """
    func_name =  __prog__ + ' generate_cells'

    project_path = form.setup['proj_loc']
    crop = form.w_combo00b.currentText()
    crop_mask_path = join(project_path, 'cropmasks',crop.lower())
    req_resol_deg = str(form.req_resol_deg)
    lon_ll, lat_ll, lon_ur, lat_ur = form.bbox
    mask_fname = glob(crop_mask_path + '\\*' + req_resol_deg + '*.nc')[0]
    nc_dset = Dataset(mask_fname, mode='r')
    # slice = nc_dset.variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]
    nc_dset.close()

    return

def fetch_soil_metrics(form):
    """
    verify and load soil CSV data
    """
    soil_dir = join(split(form.sims_dir)[0], SOIL_DIR)
    soil_mtrcs_nm =  'HWSD_recs' + '_' + '{:0=2d}'.format(form.req_resol_upscale)
    fname = join(soil_dir, soil_mtrcs_nm + '.csv')

    soil_df = None
    if isfile(fname):
        try:
            soil_df = read_csv(fname, sep=',', index_col=0)
        except EmptyDataError as err:
            print(WARN_STR + 'soil CSV file ', fname + ' has no soil records')
            soil_df = None
        else:
            print('soil CSV file ', fname + ' has ' + str(soil_df.shape[0] - 1) + ' soil records')
    else:
        print(WARN_STR + 'soil CSV file ' + fname  + ' does not exist, cannot proceed')

    QApplication.processEvents()

    return soil_df

def update_progress(last_time, ncompleted, nskipped, ntotal_grow, ngrowing, nno_grow, hwsd = None):
    """
    Update progress bar
    """
    new_time = time()
    if new_time - last_time > 5:
        # used to be: Estimated remaining
        mess = '\rCompleted: {:=6d} Growing cells: {:=6d} No grow cells: {:=6d}'.format(ncompleted, ngrowing, nno_grow)

        if hwsd is None:
            bad_muglobals = []
        else:
            bad_muglobals = hwsd.bad_muglobals
        mess += ' Skipped: {:=5d} Bad mu globals: {:=5d}'.format(nskipped, len(bad_muglobals))
        mess += ' Remaining: {:=6d}'.format(ntotal_grow - ncompleted)
        stdout.flush()
        stdout.write(mess)
        last_time = new_time

    return last_time

def update_wthr_progress(last_time, ncompleted, nskipped=0, ntotal_grow=0, ngrowing=0, nno_grow=0, region='Europe'):
    """
    Update progress bar
    """
    new_time = time()
    if new_time - last_time > 5:
        mess = '\rCompleted: {:=6d} Growing cells: {:=6d} No grow cells: {:=6d}'.format(ncompleted, ngrowing, nno_grow)
        mess += ' Skipped: {:=5d} Region: {:15s}'.format(nskipped, region)
        mess += ' Remaining: {:=6d}'.format(ntotal_grow - ncompleted)
        stdout.flush()
        stdout.write(mess)
        last_time = new_time

    return last_time

def update_avemet_progress(last_time, wthr_rsrce, scnr, region, nwrote):
    """
    Update progress bar
    """
    new_time = time()
    if new_time - last_time > 5:
        mess = '\rWrote: {:=6d} Avemet files\t'.format(nwrote)
        mess += '\twthr_rsrce: {}\tscnr:{}\tregion: {}'.format(wthr_rsrce, scnr, region)
        stdout.flush()
        stdout.write(mess)
        last_time = new_time

    return last_time