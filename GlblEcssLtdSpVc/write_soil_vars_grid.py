"""
#-------------------------------------------------------------------------------
# Name:        write_soil_vars_grid.py
# Purpose:     consist of high level functions invoked by main GUI
# Author:      Mike Martin
# Created:     11/12/2015
# Licence:     <your licence>
# Description:
#
#-------------------------------------------------------------------------------
#
"""
__prog__ = 'write_soil_vars_grid.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

import time
import csv
from pandas import read_csv
from os.path import join, split, isdir, isfile
from os import mkdir
from PyQt5.QtWidgets import QApplication

from glbl_ecsse_high_level_sp import simplify_soil_recs, _simplify_aoi
import hwsd_bil

from hwsd_mu_globals_fns import gen_grid_cells_for_band
from prepare_ecosse_files import update_progress
from glbl_ecss_cmmn_cmpntsGUI import calculate_grid_cell

SOIL_HDRS = list(['UID', 'HWSD_Y', 'HWSD_X', 'Latitude', 'Longitude', 'MU_GLOBAL', 'SHARE',
                  'T_OC', 'T_BULK_DENSITY', 'T_PH_H2O', 'T_SAND', 'T_SILT', 'T_CLAY',
                  'S_OC', 'S_BULK_DENSITY', 'S_PH_H2O', 'S_SAND', 'S_SILT', 'S_CLAY'])
SOIL_DIR = 'soil_metrics'

WARN_STR = '*** Warning *** '

def fetch_soil_metrics(form):
    """
    Create doil data frame
    """
    soil_dir = join(split(form.sims_dir)[0], SOIL_DIR)
    fname = join(soil_dir, 'HWSD_recs' + '.csv')
    
    if isfile(fname):
        soil_df = read_csv(fname, sep=',')
    else:
        print(WARN_STR + 'soil CSV file does not exist, cannot proceed')
        return False

    return True

class SoilCsvOutputs(object):
    """
    Class to write CSV files of soil data
    """
    def __init__(self, form):

        self.lgr = form.lgr
        self.hdrs = SOIL_HDRS
        soil_dir = join(split(form.sims_dir)[0], SOIL_DIR)
        if not isdir(soil_dir):
            mkdir(soil_dir)

        self.soil_dir = soil_dir
        self.study = form.w_study.text()
        self.req_resol_upscale = form.req_resol_upscale

        self.output_fobj = None
        self.writer = None

    def create_soil_file(self):
        """
        Create empty results file
        """
        req_resol_str =  '{:0=2d}'.format(self.req_resol_upscale)
        fname = join(self.soil_dir, 'HWSD_recs_' + req_resol_str + '.csv')
        try:
            self.output_fobj = open(fname, 'w', newline='')
        except (OSError, IOError) as err:
            err_mess = 'Unable to open output file. {} {}'.format(fname, err)
            self.lgr.critical(err_mess)
            print(err_mess)
            QApplication.processEvents()

        self.writer = csv.writer(self.output_fobj, delimiter=',')
        self.writer.writerow(self.hdrs)
        return

def _write_to_soil_file(form, soil_csv, num_band):
    """
    Main loop for generating soil data outputs
    """
    print('Gathering soil data for band {}'.format(num_band))

    # instantiate a soil grid and climate objects
    hwsd = hwsd_bil.HWSD_bil(form.lgr, form.hwsd_dir)

    # add requested grid resolution attributes to the form object
    calculate_grid_cell(form, hwsd.granularity)
    bbox = form.bbox

    # create grid of mu_globals based on bounding box
    # ===============================================
    nvals_read = hwsd.read_bbox_hwsd_mu_globals(bbox, form.hwsd_mu_globals, form.req_resol_upscale)

    # retrieve dictionary consisting of mu_globals (keys) and number of occurences (values)
    # =====================================================================================
    print('\nRetrieving soil data for band ' + str(num_band))
    QApplication.processEvents()

    mu_globals = hwsd.get_mu_globals_dict()
    if mu_globals is None:
        print('No soil records for AOI: {}\n'.format(bbox))
        return

    mess = 'Retrieved {} values  of HWSD grid consisting of {} rows and {} columns: ' \
          '\n\tnumber of unique mu_globals: {}'.format(nvals_read, hwsd.nlats, hwsd.nlons, len(mu_globals))
    form.lgr.info(mess)

    # for each grid point in the band) return a quintuples of integer Lat/Lon 30 seconds coords, Lat/Lons and mu_global
    hwsd.bad_muglobals = form.hwsd_mu_globals.bad_mu_globals
    aoi_res, bbox = gen_grid_cells_for_band(hwsd, form.req_resol_upscale)
    if form.w_use_high_cover.isChecked():
        aoi_res = _simplify_aoi(aoi_res)

    lon_ll_aoi, lat_ll_aoi, lon_ur_aoi, lat_ur_aoi = bbox
    num_meta_cells = len(aoi_res)
    print('Band aoi LL lat/lon: {} {}\tUR lat/lon: {} {}\t# meta cells: {}'
                                    .format(lat_ll_aoi, lon_ll_aoi, lat_ur_aoi, lon_ur_aoi, num_meta_cells))
    QApplication.processEvents()

    if num_meta_cells == 0:
        mess = 'No aoi_res recs therefore unable to create simulation files... \n'
        print(mess)
        form.lgr.info(mess)

        return

    mess = 'Generated {} Area of Interest grid cells for band {} '.format(num_meta_cells, num_band)
    form.lgr.info(mess)
    print(mess)

    print('Writing soil data for band {}...'.format(num_band))
    QApplication.processEvents()

    # open land use NC dataset
    # ========================
    last_time = time.time()
    start_time = time.time()
    completed = 0
    skipped = 0
    warn_count = 0

    # each soil can have one or more dominant soils
    # =======================================================================
    for site_rec in aoi_res:
        gran_lat, gran_lon, lat, long, area, mu_globals_props = site_rec
        gran_coord = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

        for pair in mu_globals_props.items():
            mu_global, proportion = pair

            soil_list = form.hwsd_mu_globals.soil_recs[mu_global]

            # soil_c, bulk_dens, ph = soil[:3]
            # ===============================
            for soil in soil_list:
                share = soil[-1]
                if len(soil) == 13:
                    sub_soil = soil[6:-1]
                else:
                    sub_soil = 6*[-999]

                t_oc, t_bulk_density, t_ph_h2o, t_sand, t_silt, t_clay = soil[:6]
                s_oc, s_bulk_density, s_ph_h2o, s_sand, s_silt, s_clay = sub_soil

                out_rec = list([gran_coord, gran_lat, gran_lon, round(lat, 4), round(long, 4), mu_global, share])
                out_rec += list([t_oc, t_bulk_density, t_ph_h2o, t_sand, t_silt, t_clay])
                out_rec += list([s_oc, s_bulk_density, s_ph_h2o, s_sand, s_silt, s_clay])
                soil_csv.writer.writerow(out_rec)

        completed += 1
        last_time = update_progress(last_time, start_time, completed, num_meta_cells, skipped, warn_count)
        QApplication.processEvents()

    mess = '\nBand: {}\tskipped: {}\tcompleted: {}'.format(num_band, skipped, completed)
    print(mess)
    QApplication.processEvents()

    print('')   # spacer
    return completed

def generate_all_soil_metrics(form, max_cells=100000000):
    """
    called from GUI
    """
    if hasattr(form, 'w_max_sims'):
        max_cells = int(form.w_max_sims.text())
    elif hasattr(form, 'w_max_cells'):
        max_cells = int(form.w_max_cells.text())

    if form.hwsd_mu_globals is None:
        print(WARN_STR + 'Undetermined HWSD aoi - please select a valid HSWD csv file')
        QApplication.applicationVersion()
        return

    if form.w_use_dom_soil.isChecked():
        use_dom_soil_flag = True
    else:
        use_dom_soil_flag = False

    # lat_ll is the floor i.e. least latitude, of the HWSD aoi which marks the end of the banding loop
    # ====================================================================================================
    lat_ll = form.hwsd_mu_globals.lat_ll_aoi
    lon_ll = form.hwsd_mu_globals.lon_ll_aoi
    lat_ur = form.hwsd_mu_globals.lat_ur_aoi
    lon_ur = form.hwsd_mu_globals.lon_ur_aoi
    bbox = list([lon_ll, lat_ll, lon_ur, lat_ur])
    form.bbox = bbox

    # ============================ for each PFT end =====================================
    # print('Study bounding box and HWSD CSV file overlap')
    #        ============================================
    start_at_band = form.start_at_band   # from setup file, generally zero
    print('Starting at band {}'.format(start_at_band))

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
        del(soil_recs[mu_global])

    form.hwsd_mu_globals.soil_recs = simplify_soil_recs(soil_recs, use_dom_soil_flag)
    form.hwsd_mu_globals.bad_mu_globals = [0] + hwsd.bad_muglobals
    del hwsd, soil_recs

    # Create empty soil files
    # =======================
    soil_csv = SoilCsvOutputs(form)
    soil_csv.create_soil_file()

    # main banding loop
    # =================
    ncompleted = 0
    lat_step = 0.5
    nsteps = int((lat_ur-lat_ll)/lat_step) + 1
    for isec in range(nsteps):
        lat_ll_new = lat_ur - lat_step
        num_band = isec + 1

        # if the latitude floor of the band has not reached the ceiling of the HWSD aoi then skip this band
        # =================================================================================================
        str_lat_extent = ' with latitude extent of min: {}\tmax: {}\n'.format(round(lat_ll_new, 6), round(lat_ur, 6))
        if lat_ll_new > form.hwsd_mu_globals.lat_ur_aoi or num_band < start_at_band:
            print('Skipping out of area band {} of {}'.format(num_band, nsteps)  + str_lat_extent)
            continue

        print('\nProcessing band {} of {}'.format(num_band, nsteps)  + str_lat_extent)
        QApplication.processEvents()

        form.bbox = list([lon_ll, lat_ll_new, lon_ur, lat_ur])
        completed = _write_to_soil_file(form, soil_csv, num_band)  # does actual work
        ncompleted += 1
        if ncompleted >= max_cells:
            print('Finished processing after generating {} cells'.format(ncompleted))
            QApplication.processEvents()
            break

        # check to see if the last band is completed
        # ==========================================
        if lat_ll > lat_ll_new or num_band == nsteps:
            print('Finished processing after {} bands of latitude extents'.format(num_band))
            for ichan in range(len(form.fstudy)):
                form.fstudy[ichan].close()
            break

        lat_ur = lat_ll_new

    # ==============
    soil_csv.output_fobj.close()    # close CSV file
    print('\nFinished soil metric writing')
    QApplication.processEvents()

    return
