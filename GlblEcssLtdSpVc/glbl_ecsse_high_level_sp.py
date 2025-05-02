
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
__prog__ = 'glbl_ecsse_high_level_sp.py'
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

def _generate_ecosse_files(form, climgen, mask_defn, num_band):
    """
    Main loop for generating ECOSSE outputs
    """
    study = form.study
    print('Gathering soil data for study {}...'.format(study))

    # instantiate a soil grid object
    hwsd = hwsd_bil.HWSD_bil(form.lgr, form.hwsd_dir)

    # add requested grid resolution attributes to the form object
    calculate_grid_cell(form, hwsd.granularity)
    bbox = form.bbox

    # create grid of mu_globals based on bounding box
    # ===============================================
    print('\nRetrieving soil data for band ' + str(num_band))
    QApplication.processEvents()

    nvals_read = hwsd.read_bbox_hwsd_mu_globals(bbox, form.hwsd_mu_globals, form.req_resol_upscale)

    # retrieve dictionary consisting of mu_globals (keys) and number of occurences (values)
    # =====================================================================================
    mu_globals = hwsd.get_mu_globals_dict()
    if mu_globals is None:
        print('No soil records for AOI: {}\n'.format(bbox))
        return

    mess = 'Retrieved {} values of HWSD grid '.format(nvals_read)
    mess += 'consisting of {} rows and {} columns: '.format(hwsd.nlats, hwsd.nlons)
    mess += '\n\tnumber of unique mu_globals: {}'.format(len(mu_globals))
    form.lgr.info(mess)
    print(mess)
    QApplication.processEvents()

    # for each grid point in the band return a quintuples of integer Lat/Lon 30 seconds coords, Lat/Lons and mu_global
    hwsd.bad_muglobals = form.hwsd_mu_globals.bad_mu_globals
    aoi_res, bbox = gen_grid_cells_for_band(hwsd, form.req_resol_upscale)
    if form.w_use_high_cover.isChecked():
        aoi_res = _simplify_aoi(aoi_res)

    lon_ll_aoi, lat_ll_aoi, lon_ur_aoi, lat_ur_aoi = bbox
    num_meta_cells = len(aoi_res)

    mess = 'Band aoi LL lon/lat: {} {}\t'.format(lon_ll_aoi, lat_ll_aoi)
    mess += 'UR lon/lat: {} {}\t# meta cells: {}'.format(lon_ur_aoi, lat_ur_aoi, num_meta_cells)
    print(mess)
    QApplication.processEvents()

    if num_meta_cells == 0:
        mess = 'No aoi_res recs therefore unable to create simulation files... \n'
        print(mess); form.lgr.info(mess)
        return

    # 4.5 = estimated mean number of dominant soils per cell
    # ======================================================
    est_num_sims = 0
    for site_rec in aoi_res:
        est_num_sims += len(site_rec[-1])
    est_num_sims = int(est_num_sims * 4.5)
    mess = 'Generated {} Area of Interest grid cell records '.format(num_meta_cells)
    mess += 'for band {}\twhich will result in an estimated {} simulations'.format(num_band, est_num_sims)
    form.lgr.info(mess)
    print(mess)
    QApplication.processEvents()

    print('Creating simulation files for band {}...'.format(num_band))
    QApplication.processEvents()
    #      =========================================

    # open land use NC dataset
    # ========================
    if mask_defn is not None:
        mask_defn.nc_dset = Dataset(mask_defn.nc_fname)

    pft_name = form.w_combo_pfts.currentText()
    if pft_name == '':
        pft_key = None
    else:
        pft_key = '00'

    last_time = time()
    start_time = time()
    completed = 0
    skipped = 0
    landuse_yes = 0
    landuse_no = 0
    warn_count = 0
    no_yrs_pi = 0
    no_pis = 0
    no_wthr = 0

    if mask_defn is not None:
        land_uses = get_hilda_land_uses(form.w_hilda_lus)

    # generate sets of Ecosse files for each site where each site has one or more soils
    # each soil can have one or more dominant soils
    # =======================================================================
    max_cells = int(form.w_max_sims.text())
    for site_indx, site_rec in enumerate(aoi_res):

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
        make_ecosse_file(form, climgen, ltd_data, site_rec, study, lta_wthr_recs, wthr_gran_coord)
        completed += 1
        if completed >= max_cells:
            break

        last_time = update_progress(last_time, start_time, completed, num_meta_cells, skipped, warn_count)
        QApplication.processEvents()

    # close plant input NC dataset
    # ============================
    if mask_defn is not None:
        mask_defn.nc_dset.close()

    mess = '\nBand: {}\tforest yes: {} no: {}\tno weather: {}\t'.format(num_band, landuse_yes, landuse_no, no_wthr)
    mess += 'no plant inputs: {}\tcompleted: {}'.format(no_pis + no_yrs_pi, completed)
    print(mess); QApplication.processEvents()

    print('')   # spacer
    return mess

def generate_banded_sims(form):
    """
    called from GUI
    """
    ret_code = fetch_soil_metrics(form)
    if ret_code is None:
        return

    if form.hwsd_mu_globals is None:
        print('Undetermined HWSD aoi - please select a valid HSWD csv file')
        return

    if form.w_use_dom_soil.isChecked():
        use_dom_soil_flag = True
    else:
        use_dom_soil_flag = False

    # make sure bounding box is correctly set
    # =======================================
    lon_ll = form.litter_defn.lon_frst
    lat_ll = form.litter_defn.lat_frst
    lon_ur = form.litter_defn.lon_last
    lat_ur_lttr = form.litter_defn.lat_last
    form.bbox = list([lon_ll, lat_ll, lon_ur, lat_ur_lttr])

    # lat_ll_aoi is the floor i.e. least latitude, of the HWSD aoi which marks the end of the banding loop
    # ====================================================================================================
    lat_ll_aoi = form.hwsd_mu_globals.lat_ll_aoi
    lon_ll_aoi = form.hwsd_mu_globals.lon_ll_aoi
    lat_ur_aoi = form.hwsd_mu_globals.lat_ur_aoi
    lon_ur_aoi = form.hwsd_mu_globals.lon_ur_aoi
    bbox_aoi = list([lon_ll_aoi,lat_ll_aoi,lon_ur_aoi,lat_ur_aoi])

    # check overlap - study too far to west or east or too far south or north of AOI file
    # ===================================================================================
    if (lon_ur < lon_ll_aoi) or (lon_ll > lon_ur_aoi) or  (lat_ur_lttr < lat_ll_aoi) or (lat_ll > lat_ur_aoi):
        print(ERROR_STR + 'Study bounding box and HWSD CSV file do not overlap - no simulations are possible')
        return

    # weather choice
    # ==============
    weather_resource = form.combo10w.currentText()

    # mask
    # ====
    if form.mask_fn is None or form.w_hilda_lus['all'].isChecked():
        mask_defn = None
    else:
        mask_defn = ManagementSet(form.mask_fn, 'cropmasks')

    # ============================ for each PFT end =====================================
    # print('Study bounding box and HWSD CSV file overlap')
    #        ============================================

    # extract required values from the HWSD database and simplify if requested
    # ========================================================================
    hwsd = hwsd_bil.HWSD_bil(form.lgr, form.hwsd_dir)

    # TODO: patch to be sorted
    # ========================
    mu_global_pairs = {}
    for mu_global in form.hwsd_mu_globals.mu_global_list:
        mu_global_pairs[mu_global] = None

    soil_recs = hwsd.get_soil_recs(mu_global_pairs)  # list is already sorted

    # TODO: patch to be sorted
    # ========================
    for mu_global in hwsd.bad_muglobals:
        del soil_recs[mu_global]

    form.hwsd_mu_globals.soil_recs = simplify_soil_recs(soil_recs, use_dom_soil_flag)
    form.hwsd_mu_globals.bad_mu_globals = [0] + hwsd.bad_muglobals
    del hwsd; del soil_recs
    lat_ur_aoi = form.hwsd_mu_globals.lat_ur_aoi

    # main banding loop
    # =================
    climgen = ClimGenNC(form)   # create climate object
    if not(climgen.readCoordsLookup()):
        return

    lat_step = 0.5
    lat_ur = copy(lat_ur_lttr)
    nbands = int((lat_ur - lat_ll)/lat_step) + 1
    start_at_band = int(form.w_strt_band.text())
    end_at_band = int(form.w_end_band.text())
    nbands_to_prcss = end_at_band - start_at_band + 1
    print('Total # bands: {}\tstarting and ending at bands {} and {}'.format(nbands, start_at_band, end_at_band))
    QApplication.processEvents()

    num_band = -999
    band_reports = []
    for iband in range(nbands):
        lat_ll_new = lat_ur - lat_step
        num_band = iband + 1
        if lat_ll_new > lat_ur_aoi:
            mess = 'Skipping simulations at band {} since new band latitude floor '.format(num_band)
            print(mess + '{} exceeds AOI upper latitude {}'.format(round(lat_ll_new,6), round(lat_ur_aoi, 6)))

        elif num_band < start_at_band:
            print('Skipping out of area band {} of {} with latitude extent of min: {}\tmax: {}'
                                            .format(num_band, nbands, round(lat_ll_new, 6), round(lat_ur, 6)))
        elif num_band > end_at_band:
            print('Exiting from processing after {} bands'.format(num_band - 1))
            break
        else:
            form.bbox = list([lon_ll, lat_ll_new, lon_ur, lat_ur])

            print('\nProcessing band {} of {} with latitude extent of min: {}\tmax: {}'
                                            .format(num_band, nbands_to_prcss, round(lat_ll_new,6), round(lat_ur, 6)))
            QApplication.processEvents()

            report = _generate_ecosse_files(form, climgen, mask_defn, num_band)  # does actual work
            band_reports.append(report)

        # check to see if the last band is completed
        # ==========================================
        if lat_ll_aoi > lat_ll_new:
            break

        lat_ur = lat_ll_new

    print('Finished processing after {} bands of latitude extents'.format(num_band))
    #      ======================================================
    QApplication.processEvents()
    form.band_reports = band_reports
    for report in band_reports:
        form.lgr.info(report)

    for ichan in range(len(form.fstudy)):
        form.fstudy[ichan].close()

    return

# ===============================================================
#
def simplify_soil_recs(soil_recs, use_dom_soil_flag):
    """
    compress soil records if duplicates are present
    simplify soil records if requested
    each mu_global points to a group of soils
    a soil group can have up to ten soils
    """
    func_name = __prog__ + ' _simplify_soil_recs'

    num_raw = 0         # total number of sub-soils
    num_compress = 0    # total number of sub-soils after compressions

    new_soil_recs = {}
    for mu_global in soil_recs:

        # no processing necessary
        # =======================
        num_sub_soils = len(soil_recs[mu_global])
        num_raw += num_sub_soils
        if num_sub_soils == 1:
            num_compress += 1
            new_soil_recs[mu_global] = soil_recs[mu_global]
            continue

        # check each soil for duplicates
        # ==============================
        new_soil_group = []
        soil_group = sorted(soil_recs[mu_global])

        # skip empty groups
        # =================
        if len(soil_group) == 0:
            continue

        first_soil = soil_group[0]
        metrics1 = first_soil[:-1]
        share1   = first_soil[-1]
        for soil in soil_group[1:]:
            metrics2 = soil[:-1]
            share2 =   soil[-1]
            if metrics1 == metrics2:
                share1 += share2
            else:
                new_soil_group.append(metrics1 + [share1])
                metrics1 = metrics2
                share1 = share2

        new_soil_group.append(metrics1 + [share1])
        num_sub_soils = len(new_soil_group)
        num_compress += num_sub_soils
        if num_sub_soils == 1:
            new_soil_recs[mu_global] = new_soil_group
            continue

        if use_dom_soil_flag:
            # assign 100% to the first entry of sorted list
            # =============================================
            dom_soil = copy(sorted(new_soil_group, reverse = True, key=itemgetter(-1))[0])
            dom_soil[-1] = 100.0
            new_soil_recs[mu_global] = list([dom_soil])

    mess = 'Leaving {}\trecords in: {} out: {}'.format(func_name, len(soil_recs),len(new_soil_recs))
    print(mess + '\tnum raw sub-soils: {}\tafter compression: {}'.format(num_raw, num_compress))
    return new_soil_recs

def _simplify_aoi(aoi_res):
    """
    simplify AOI records
    """
    aoi_res_new = []
    j = 0
    for site_rec in aoi_res:
        content = site_rec[-1]
        npairs = len(content)
        if npairs == 0:
            print('No soil information for AOI cell at lat/long: {} {} - will skip'
                                                                .format(round(site_rec[2], 4), round(site_rec[3], 4)))
        elif npairs == 1:
            aoi_res_new.append(site_rec)
        else:
            site_rec_list = list(site_rec)  # convert tuple to a list so we can edit last element
            new_content = sorted(content.items(), reverse = True, key = itemgetter(1))  # sort content so we can pick up most dominant mu_global
            total_proportion = sum(content.values())    # add up proportions
            site_rec_list[-1] = {new_content[0][0]: total_proportion}       # create a new single mu global with summed proportions

            aoi_res_new.append(tuple(site_rec_list)) # convert list to tuple

    return aoi_res_new
