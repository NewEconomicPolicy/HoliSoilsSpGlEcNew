"""
#-------------------------------------------------------------------------------
# Name:
# Purpose:     consist of high level functions invoked by main GUI
# Author:      Mike Martin
# Created:     11/12/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#
"""
__prog__ = 'wthr_generation_fns'
__version__ = '0.0.1'
__author__ = 's03mm5'

from time import time
import csv
from os.path import join, normpath, isdir, split, isfile
from os import listdir, walk, makedirs
from pandas import Series, read_excel, DataFrame
from PyQt5.QtWidgets import QApplication

from getClimGenNC_ltd import ClimGenNC
from getClimGenFns_ss import (genLocalGrid, fetch_wthr_dset_overlap, join_hist_fut_to_all_wthr)
from glbl_ecsse_low_level_fns_sv import update_wthr_progress, update_avemet_progress
from prepare_ecosse_low_level import fetch_long_term_ave_wthr_recs, make_met_files
from hwsd_soil_class import _gran_coords_from_lat_lon as gran_coords_from_lat_lon
from weather_datasets import write_csv_wthr_file

from thornthwaite import thornthwaite

ERROR_STR = '*** Error *** '
WARN_STR = '*** Warning *** '
QUICK_FLAG = False       # forces break from loops after max cells reached in first GCM and SSP

GRANULARITY = 120
NEXPCTD_MET_FILES = 302
MAX_BANDS = 500000
MAX_CELLS = 9999999
LTA_RECS_FN = 'lta_ave.txt'

SPACER_LEN = 12

def make_wthr_coords_lookup(form):
    """
    C
    """
    func_name = __prog__ + ' ClimGenNC __init__'

    if not hasattr(form, 'combo10'):
        print(ERROR_STR + 'This function ' + func_name + ' requires attribute combo10')
        QApplication.processEvents()
        return -1

    fut_clim_scen = form.combo10.currentText()
    wthr_out_dir = join(split(form.sims_dir)[0], 'Wthr', fut_clim_scen)
    if not isdir(wthr_out_dir):
        makedirs(wthr_out_dir)

    for directory, subdirs_raw, files in walk(wthr_out_dir):
        num_sims = len(subdirs_raw)
        break

    if num_sims == 0:
        print(WARN_STR + 'no sub-directories under path ' + wthr_out_dir)
        QApplication.processEvents()
        return

    #
    # ===============================================================
    recs = []
    for gran_coord in subdirs_raw:
        if gran_coord.find('_') == -1:
            print(WARN_STR + 'non compliant directory found in weather directory ' + join(wthr_out_dir, gran_coord))
            QApplication.processEvents()
        else:
            gran_lat, gran_lon = gran_coord.split('_')
            lat = 90.0 - int(gran_lat) / GRANULARITY
            lon = (int(gran_lon) / GRANULARITY) - 180.0
            recs.append([gran_coord, lat, lon])

    # write coords file
    # =================
    coords_fn = join(wthr_out_dir, 'coords_lookup.csv')
    df = DataFrame(recs, columns=['gran_coord', 'Lat', 'Lon'])
    df.to_csv(coords_fn, sep='\t', index=False)
    mess = 'Wrote coordinates lookup file: ' + coords_fn
    print(mess)
    QApplication.processEvents()

    return

def generate_all_weather(form):
    """
    C
    """
    lat_ll_aoi = form.hwsd_mu_globals.lat_ll_aoi
    lon_ll_aoi = form.hwsd_mu_globals.lon_ll_aoi
    lat_ur_aoi = form.hwsd_mu_globals.lat_ur_aoi
    lon_ur_aoi = form.hwsd_mu_globals.lon_ur_aoi
    bbox_aoi = list([lon_ll_aoi, lat_ll_aoi, lon_ur_aoi, lat_ur_aoi])

    max_cells = int(form.w_max_sims.text())

    resol_deg = 0.5
    resol_d2 = resol_deg/2

    sims_dir = form.sims_dir
    study = form.w_study.text()

    wthr_set_nm = form.weather_set_linkages['EFISCEN-ISIMIP'][1]
    this_gcm = wthr_set_nm.split('_')[0]
    fut_wthr_set = form.weather_sets[wthr_set_nm]

    wthr_set_nm = form.weather_set_linkages['EFISCEN-ISIMIP'][0]
    hist_wthr_set = form.weather_sets[wthr_set_nm]

    bbox_wthr = fetch_wthr_dset_overlap(hist_wthr_set, fut_wthr_set)

    # generate weather dataset indices which enclose the AOI
    # ======================================================
    scnr = form.combo10.currentText()
    climgen = ClimGenNC(form, scnr)
    sim_start_year = climgen.sim_start_year
    sim_end_year = climgen.sim_end_year

    # development only
    # ================
    if max_cells <= 3:
        bbox_wthr = (12.0, 47.0, 14.0, 49.0)
        max_cells = MAX_CELLS

    aoi_indices_fut = genLocalGrid(fut_wthr_set, bbox_wthr, bbox_aoi)
    aoi_indices_hist = genLocalGrid(hist_wthr_set, bbox_wthr, bbox_aoi)

    # for each GCM and SSP dataset group e.g. UKESM1-0-LL 585
    # =======================================================
    print('')    
    last_time = time()
    num_band = -999
    write_csv_wthr_flag = False

    print('Getting historic weather data from weather set: ' + hist_wthr_set['ds_precip'])
    QApplication.processEvents()

    pettmp_hist = climgen.fetch_cru_historic_NC_data(aoi_indices_hist, num_band, max_cells)
    if pettmp_hist is None:
        print('\nHistorical data retrieval failed from weather set: ' + 'CRU' + '\tScenario: ' + scnr)
        QApplication.processEvents()
        return -1

    print('\nGetting future data from weather set: ' + this_gcm + '\tScenario: ' + scnr)
    QApplication.processEvents()

    #      =============================
    dset_strt_yr = climgen.fut_wthr_set_defn['year_start']
    dset_end_yr = climgen.fut_wthr_set_defn['year_end']
    nmnths = (dset_end_yr - dset_strt_yr + 1) * 12
    pettmp_fut = climgen.fetch_isimip_NC_data(aoi_indices_fut, dset_strt_yr, nmnths, max_cells)
    if pettmp_fut is None:
        print('\nFuture data retrieval failed from weather set: ' + this_gcm + '\tScenario: ' + scnr)
        QApplication.processEvents()

    keys_hist = list(pettmp_hist['precipitation'].keys())
    keys_fut = list(pettmp_fut['precipitation'].keys())
    keys_hist, keys_fut = _check_and_sync_keys(keys_fut, keys_hist)
    
    pettmp_all = join_hist_fut_to_all_wthr(climgen, pettmp_hist, pettmp_fut)

    # create weather
    # ==============
    nwrttn = 0    
    site_obj = MakeSiteObj(form, climgen)
    for gran_coord in keys_hist:
        
        if gran_coord in pettmp_fut['precipitation']:
            lat, lon = pettmp_hist['lat_lons'][gran_coord]
            clim_dir = make_wthr_files(site_obj, lat, gran_coord, climgen, pettmp_hist, pettmp_all)
            if write_csv_wthr_flag:
                write_csv_wthr_file(form.lgr, study, this_gcm, scnr, lat, lon, sim_start_year, sim_end_year,
                        pettmp_fut['precipitation'][gran_coord], pettmp_fut['temperature'][gran_coord], clim_dir)
            nwrttn += 1
            if nwrttn >= max_cells:
                print('\nFinished checking after {} cells completed'.format(nwrttn))
                QApplication.processEvents()
                break

        last_time = update_wthr_progress(last_time, nwrttn)

    # close NC files
    # ==============
    """
    for metric in list(['precip', 'tas']):
        hist_wthr_dsets[metric].close()
        fut_wthr_dsets[metric].close()
    """

    # write coords file
    # =================
    make_wthr_coords_lookup(form)
    mess = 'Completed weather set: ' + this_gcm + '\tScenario: ' + scnr + '\n'   
    print(mess)

    print('Finished weather generation - total number of sets written: {}'.format(nwrttn))

    return

def _check_and_sync_keys(keys_fut, keys_hist):
    """
    check future first, then historic
    """
    not_in_fut = []
    new_keys_hist = []
    for key in keys_hist:
        if key in keys_fut:
            new_keys_hist.append(key)
        else:
            not_in_fut.append(key)

    not_in_hist = []
    new_keys_fut = []
    for key in keys_fut:
        if key in keys_hist:
            new_keys_fut.append(key)
        else:
            not_in_hist.append(key)

    # ======= report =========
    if len(not_in_fut) > 0:
        print('keys not in future: {}'.format(not_in_fut))
        
    if len(not_in_hist) > 0:
        print('keys not in history: {}'.format(not_in_hist))
    
    QApplication.processEvents()

    return new_keys_hist, new_keys_fut

def make_avemet_file(clim_dir, lta_precip, lta_pet, lta_tmean):
    """
    will be copied
    """
    avemet_dat = join(clim_dir, 'AVEMET.DAT')
    with open(avemet_dat, 'w') as fobj:
        for imnth, (precip, pet, tmean) in enumerate(zip(lta_precip, lta_pet, lta_tmean)):
            fobj.write('{} {} {} {}\n'.format(imnth + 1, precip, pet, tmean))

    return

def make_wthr_files(site, lat, gran_coord, climgen, pettmp_hist, pettmp_all):
    """
    generate ECOSSE historic and future weather data
    """
    clim_dir = normpath(join(site.wthr_prj_dir, gran_coord))

    if pettmp_hist is None:
        return

    gran_lon = gran_coord.split('_')[1]
    lon = (int(gran_lon) / GRANULARITY) - 180.0
    mess = 'granular coord {} with lat/lon: {} {}\t'.format(gran_coord, lat, lon)

    if gran_coord not in pettmp_hist['precipitation']:
        print(WARN_STR + mess + 'not in historic weather')
        QApplication.processEvents()
        return

    if gran_coord not in pettmp_all['precipitation']:
        print(WARN_STR + + mess + 'not in simulation weather')
        QApplication.processEvents()
        return

    if not isdir(clim_dir):
        makedirs(clim_dir)  # only create if weather data all present

    # calculate historic average weather
    # ==================================
    '''
    pettmp_hist_site = {'precip': pettmp_hist['precipitation'][gran_coord],
                        'tas': pettmp_hist['temperature'][gran_coord]}
    hist_lta_precip, hist_lta_tmean, hist_weather_recs = fetch_long_term_ave_wthr_recs(climgen, pettmp_hist_site)
    '''
    # write a single set of met files for all simulations for this grid cell
    # ======================================================================
    pettmp_all_site = {'precip': pettmp_all['precipitation'][gran_coord], 'tas': pettmp_all['temperature'][gran_coord]}

    year_start = climgen.hist_wthr_set_defn['year_start']
    met_fnames = make_met_files(clim_dir, lat, climgen, pettmp_all_site, year_start)  # all weather
    nmet_fns = len(met_fnames)

    # create additional weather related files from already existing met files
    # =======================================================================
    '''
    irc = climgen.create_FutureAverages(clim_dir, lat, gran_coord, site, hist_lta_precip, hist_lta_tmean)
    if irc == 0:
        lta_ave_fn = _make_lta_file(site, clim_dir)
    '''
    return clim_dir

def fetch_hist_lta_from_lat_lon(sims_dir, climgen, lat, lon):
    """
    check existence of weather cell
    """
    read_lta_flag = True
    integrity_flag, hist_lta_recs, met_fnames = _check_wthr_cell_exstnc(sims_dir, climgen, lat, lon, read_lta_flag)

    return integrity_flag, hist_lta_recs, met_fnames

def _check_wthr_cell_exstnc(sims_dir, climgen, lat, lon, read_lta_flag=False):
    """
    check existence and integrity of weather cell
    allowable criteria are 1) a full set of weather files, namely 300 met files e.g. met2014s.txt, lta_ave.txt and AVEMET.DAT
                           2) an empty directory
    """
    integrity_flag = False
    hist_lta_recs = None
    met_fnames = None
    gran_lat, gran_lon = gran_coords_from_lat_lon(lat, lon)
    gran_coord = '{0:0=5g}_{1:0=5g}'.format(gran_lat, gran_lon)
    clim_dir = normpath(join(sims_dir, climgen.region_wthr_dir, gran_coord))
    if isdir(clim_dir):
        fns = listdir(clim_dir)
        nfiles = len(fns)
        if nfiles == 0 or nfiles >= 302:
            if nfiles == 0:
                integrity_flag = True
                hist_lta_recs, met_fnames = None, None
            else:
                if LTA_RECS_FN in fns:
                    if read_lta_flag:
                        lta_ave_fn = join(clim_dir, LTA_RECS_FN)
                        hist_lta_recs = []
                        with open(lta_ave_fn, 'r') as fave:
                            for line in fave:
                                line = line.rstrip()  # strip out all tailing whitespace
                                hist_lta_recs.append(line)

                    integrity_flag = True
                    met_fnames = fns[2:]

    return integrity_flag, hist_lta_recs, met_fnames

def write_avemet_files(form):
    """
    traverse each GCM and SSP dataset group e.g. UKESM1-0-LL 585
    """
    print('')
    max_cells = int(form.w_max_cells.text())
    sims_dir = form.setup['sims_dir']

    nwrote = 0
    for wthr_set in form.weather_set_linkages['EFISCEN-ISIMIP']:
        wthr_rsrce, scnr = wthr_set.split('_')
        if scnr == 'hist':  # mod
            continue

        # for each region
        # ===============
        for irow, region in enumerate(form.regions['Region']):
            lon_ll, lon_ur, lat_ll, lat_ur, wthr_dir_abbrv = form.regions.iloc[irow][1:]

            # main traversal loop
            # ===================
            region_wthr_dir = wthr_dir_abbrv + wthr_rsrce + '_' + scnr
            clim_dir = normpath(join(sims_dir, region_wthr_dir))

            mess = '\nProcessing weather set: ' + wthr_rsrce + '\tScenario: ' + scnr + '\tRegion: ' + region
            mess += '\t\tabbrev: ' + wthr_dir_abbrv + '\n\tclim_dir: ' + clim_dir
            print(mess)

            if not isdir(clim_dir):
                print(clim_dir + ' *** does not exist ***')
                break

            # step through each directory comprising ECOSSE met files
            # =======================================================
            last_time = time()
            nwrote = 0
            for drctry, subdirs, files in walk(clim_dir):
                nfiles = len(files)
                nsubdirs = len(subdirs)
                if nsubdirs > 0:  # first directory has four scenarios e.g.  Y:\GlblEcssOutputsSv\EcosseSims\AfUKESM1-0-LL_126
                    continue

                # there should be 300 met files plus lta_ave.txt and AVEMET.DAT
                # =============================================================
                last_time = update_avemet_progress(last_time, wthr_rsrce, scnr, region, nwrote)
                if nfiles >= NEXPCTD_MET_FILES:
                    continue

                # if lta_ave.txt is not present then something is wrong
                # =====================================================
                if LTA_RECS_FN in files:
                    lta_ave_fn = join(drctry, LTA_RECS_FN)
                    with open(lta_ave_fn, 'r') as flta_ave:

                        lta_recs = flta_ave.readlines()
                        vals = [float(rec.split('#')[0]) for rec in lta_recs]
                        lta_precip, lta_tmean = vals[:12], vals[12:]

                        gran_coord = split(drctry)[1]
                        gran_lat = int(gran_coord.split('_')[0])
                        cell_lat = 90.0 - gran_lat / GRANULARITY
                        lta_pet = thornthwaite(lta_tmean, cell_lat)

                        make_avemet_file(drctry, lta_precip, lta_pet, lta_tmean)
                        nwrote += 1
                else:
                    print(WARN_STR + LTA_RECS_FN + ' file should be present in ' + drctry)

            if nwrote >= max_cells:
                print('\nFinished checking having written {} AVEMET.DAT files'.format(nwrote))
                break

            print('Completed Region: ' + region)

        print('Completed weather set: ' + wthr_rsrce + '\tScenario: ' + scnr + '\n')

    print('Finished AVEMET creation - checked: {} cells'.format(nwrote))
    return

def _make_lta_file(site, clim_dir):
    """
    write long term average climate section of site.txt file
    """
    lines = []
    lta_precip, lta_tmean = site.lta_precip, site.lta_tmean
    if lta_precip is None or lta_tmean is None:
        return

    for precip, month in zip(lta_precip, site.months):
        lines.append(_make_line('{}'.format(precip), '{} long term average monthly precipitation [mm]'.format(month)))

    for tmean, month in zip(lta_tmean, site.months):
        lines.append(_make_line('{}'.format(tmean), '{} long term average monthly temperature [mm]'.format(month)))

    lta_ave_fn = join(clim_dir, LTA_RECS_FN)
    with open(lta_ave_fn, 'w') as fhand:
        fhand.writelines(lines)

    # will be copied
    # ==============
    make_avemet_file(clim_dir, site.lta_precip, site.lta_pet, site.lta_tmean)

    return lta_ave_fn

def _make_line(data, comment):
    """

    """
    spacer_len = max(SPACER_LEN - len(data), 2)
    spacer = ' ' * spacer_len

    return '{}{}# {}\n'.format(data, spacer, comment)

class MakeSiteObj(object,):
    """
    C
    """
    def __init__(self, form, climgen):

        func_name =  __prog__ +  ' ClimGenNC __init__'

        self.wthr_prj_dir = climgen.wthr_out_dir
        self.months = climgen.months

# ==========================
def create_wthr_averages(lggr, climgen, lat_inp, gran_coord, period, text_flag):
    """
    use prexisting metyyyys.txt files to generate a text file of average weather which will subsequently
    be included in the input.txt file
    also create a climate file which is the average of the year range
    """
    func_name = ' create_lta_averages'
    full_func_name = __prog__ + func_name

    clim_dir = climgen.wthr_out_dir
    months = climgen.months
    output_recs = None

    if period == 'historic':
        yr_strt = climgen.hist_start_year
        yr_end = climgen.hist_end_year
    else:
        yr_strt = climgen.sim_start_year
        yr_end = climgen.sim_end_year

    ave_text_fn = 'met{}_to_{}_ave.txt'.format(yr_strt, yr_end)
    ave_met_fn = 'met' + str(yr_strt) + '_' + str(yr_end) + 'a.txt'

    num_period_yrs = yr_end - yr_strt + 1

    # skip if already exists
    ave_long_text_fn = join(clim_dir, gran_coord, ave_text_fn)
    ave_long_met_fn = join(clim_dir, gran_coord, ave_met_fn)
    if isfile(ave_long_met_fn) and isfile(ave_long_text_fn):
        mess = 'Files:\n\t' + ave_long_met_fn + ' and ' + ave_long_text_fn + '\n\talready exist - will overwrite'
        lggr.info(WARN_STR + mess)

    # read precipitation and temperature
    sim_precip = {}
    sim_tmean = {}
    for month in months:
        sim_precip[month] = 0.0
        sim_tmean[month] = 0.0

    for year in range(yr_strt, yr_end + 1):
        fname = 'met{0}s.txt'.format(year)
        met_fpath = join(clim_dir, gran_coord, fname)

        if not isfile(met_fpath):
            print('File ' + met_fpath + ' does not exist - will abandon average weather creation')
            QApplication.processEvents()
            return -1

        with open(met_fpath, 'r', newline='') as fpmet:
            lines = fpmet.readlines()

        for line, month in zip(lines, months):
            tlst = line.split('\t')
            sim_precip[month] += float(tlst[1])
            sim_tmean[month] += float(tlst[3].rstrip('\r\n'))

    # write stanza for input.txt file consisting of long term average climate
    # =======================================================================
    if text_flag:
        output_recs = []
        for month in months:
            ave_precip = sim_precip[month]/num_period_yrs
            output_recs.append(_input_txt_line_layout('{}'.format(round(ave_precip, 1)),
                                                '{} long term average monthly precipitation [mm]'.format(month)))

        for month in months:
            ave_tmean = sim_tmean[month]/num_period_yrs
            output_recs.append(_input_txt_line_layout('{}'.format(round(ave_tmean, 2)),
                                                '{} long term average monthly temperature [degC]'.format(month)))

        # write text file of average simulated weather which will subsequently be included in the input.txt file
        # ======================================================================================================
        try:
            fhand = open(ave_long_text_fn, 'w')
        except IOError:
            raise IOError(ERROR_STR + 'Unable to open file: ' + ave_text_fn)
        else:
            fhand.writelines(output_recs)
            fhand.close()

        lggr.info('Successfully wrote average weather file {} in function {}'.format(ave_text_fn, func_name))

    # write long term average climate file
    # ====================================
    ave_precip = [round(sim_precip[month]/num_period_yrs, 1) for month in months]
    ave_tmean = [round(sim_tmean[month]/num_period_yrs, 1) for month in months]

    # pet
    if max(ave_tmean) > 0.0:
        pet = thornthwaite(ave_tmean, lat_inp, year)
    else:
        pet = [0.0]*12
        mess = WARN_STR + 'monthly average temperatures are below zero in ' + full_func_name
        mess += ' for latitude: {}\tgranular coord: {}'.format(lat_inp, gran_coord)
        print(mess)

    pot_evapotrans = [round(p, 1) for p in pet]

    # write file
    output = []
    for tstep, mean_temp in enumerate(ave_tmean):
        output.append([tstep+1, ave_precip[tstep], pot_evapotrans[tstep], mean_temp])

    with open(ave_long_met_fn, 'w', newline='') as fpout:
        writer = csv.writer(fpout, delimiter='\t')
        writer.writerows(output)
        fpout.close()

    lggr.info('Successfully wrote average weather file {} in function {}'.format(ave_met_fn, func_name))

    return output_recs

def _input_txt_line_layout(data, comment):
    """
    C
    """
    set_spacer_len = 12
    spacer_len = max(set_spacer_len - len(data), 2)
    spacer = ' ' * spacer_len
    return '{}{}# {}\n'.format(data, spacer, comment)
