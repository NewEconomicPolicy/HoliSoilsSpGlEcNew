
# -------------------------------------------------------------------------------
# Name:
# Purpose:     consist of high level functions invoked by main GUI
# Author:      Mike Martin
# Created:     11/12/2015
# Licence:     <your licence>
# Description:
#
# ------------------------------------------------------------------------------
#
__prog__ = 'glec_new_high_level_sp.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

from time import time
from operator import itemgetter
from copy import copy
from netCDF4 import Dataset
from PyQt5.QtWidgets import QApplication

import make_ltd_data_files
from getClimGenNC_ltd import ClimGenNC
import hwsd_bil

from hwsd_mu_globals_fns import gen_grid_cells_for_band
from litter_and_orchidee_fns import resize_yrs_pi
from prepare_ecosse_files_ss import update_progress, make_ecosse_file
from glbl_ecss_cmmn_cmpntsGUI import calculate_grid_cell
from mngmnt_fns_and_class import ManagementSet, check_mask_location, get_hilda_land_uses
from glbl_ecsse_low_level_fns_sv import fetch_coord_nearest_xy, fetch_soil_metrics
from wthr_generation_fns import create_wthr_averages

WARN_STR = '*** Warning *** '
ERROR_STR = '*** Error *** '

def generate_sims(form):
    """
    called from GUI
    """
    soil_df = fetch_soil_metrics(form)
    if soil_df is None:
        return

    # weather choice
    # ==============
    num_meta_cells = len(soil_df)
    study = form.w_study.text()
    max_cells = int(form.w_max_sims.text())
    weather_resource = form.combo10w.currentText()

    # mask
    # ====
    if form.mask_fn is None or form.w_hilda_lus['all'].isChecked():
        mask_defn = None
        land_uses = None
    else:
        mask_defn = ManagementSet(form.mask_fn, 'cropmasks')
        land_uses = get_hilda_land_uses(form.w_hilda_lus)
        mask_defn.nc_dset = Dataset(mask_defn.nc_fname)

    # pft
    # ===
    pft_name = form.w_combo_pfts.currentText()
    if pft_name == '':
        pft_key = None
    else:
        pft_key = '00'

    # =========
    climgen = ClimGenNC(form)  # create climate object
    if not (climgen.readCoordsLookup()):
        return

    # instantiate a soil grid object and add requested grid resolution attributes to the form object
    # ==============================================================================================
    '''
    hwsd = hwsd_bil.HWSD_bil(form.lgr, form.hwsd_dir)
    calculate_grid_cell(form, hwsd.granularity)
    bbox = form.bbox

    # create grid of mu_globals based on bounding box
    # ===============================================
    print('\nRetrieving soil data for study ' + study)
    QApplication.processEvents()

    nvals_read = hwsd.read_bbox_hwsd_mu_globals(bbox, form.hwsd_mu_globals, form.req_resol_upscale)

    # retrieve dictionary consisting of mu_globals (keys) and number of occurences (values)
    # =====================================================================================
    mu_globals = hwsd.get_mu_globals_dict()
    if mu_globals is None:
        print('No soil records for AOI: {}\n'.format(bbox))
        return
   '''

    # main loop
    # =========
    landuse_yes, landuse_no, skipped, no_wthr, no_yrs_pi, completed, warn_count, no_pis = 8*[0]
    last_time = time()
    start_time = time()
    for site_index, row in soil_df.iterrows():

        # generate sets of Ecosse files for each site where each site has one or more soils
        # each soil can have one or more dominant soils
        # =======================================================================
        area = -999
        site_rec = [int(row['HWSD_Y']), int(row['HWSD_X']),
                    row['Latitude'], row['Longitude'], area, {int(row['MU_GLOBAL']): 1.0}]

        # land use mask
        # =============
        if mask_defn is not None:
            if check_mask_location(mask_defn, site_rec, land_uses, form.req_resol_deg):
                landuse_yes += 1
            else:
                landuse_no += 1
                skipped += 1
                continue

        gran_lat, gran_lon, lat, long, area, mu_globals_props = site_rec
        wthr_gran_coord = fetch_coord_nearest_xy(climgen.coords_lookup, lat, long)
        if wthr_gran_coord is None:
            mess = WARN_STR + 'No weather data for lat/lon: {}/{}\t'.format(lat, long)
            mess += 'granular lat/lon: {}/{}'.format(gran_lat, gran_lon)
            form.lgr.info(mess)
            no_wthr += 1
            continue

        lta_wthr_recs = create_wthr_averages(form.lgr, climgen, lat, wthr_gran_coord, 'simulation', text_flag=True)
        create_wthr_averages(form.lgr, climgen, lat, wthr_gran_coord, 'simulation', text_flag=False)

        yrs_pi = form.litter_defn.get_efiscen_nc_data(pft_key, lat, long, form.w_baseline.isChecked())
        if yrs_pi is None:
            no_yrs_pi += 1
            continue

        yrs_pi = resize_yrs_pi(climgen.sim_start_year, climgen.sim_end_year, yrs_pi)

        # create limited data object
        # ==========================
        ltd_data = make_ltd_data_files.MakeLtdDataFiles(form, climgen, yrs_pi)
        soil_list = [[123456.0, 1.1, 5.3, 23, 36, 41, 123456.0, 1.1, 5.3, 23, 36, 41, 100.0]]
        make_ecosse_file(form, climgen, ltd_data, site_rec, study, lta_wthr_recs, wthr_gran_coord, soil_list)
        completed += 1
        if completed >= max_cells:
            break

        last_time = update_progress(last_time, start_time, completed, num_meta_cells, skipped, warn_count)
        QApplication.processEvents()

    # close plant input NC dataset
    # ============================
    if mask_defn is not None:
        mask_defn.nc_dset.close()

    mess = '\n\tforest yes: {} no: {}\tno weather: {}\t'.format(landuse_yes, landuse_no, no_wthr)
    mess += 'no plant inputs: {}\tcompleted: {}'.format(no_pis + no_yrs_pi, completed)
    print(mess);
    QApplication.processEvents()

    print('')  # spacer
    return mess

    print('Finished processing after {} bands of latitude extents'.format(num_band))
    #      ======================================================
    QApplication.processEvents()

    for ichan in range(len(form.fstudy)):
        form.fstudy[ichan].close()

    return
