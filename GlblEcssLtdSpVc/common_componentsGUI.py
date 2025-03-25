"""
#-------------------------------------------------------------------------------
# Name:        common_componentsGUI.p
# Purpose:     consist of high level functions invoked by main GUI
# Author:      Mike Martin
# Created:     11/12/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#
"""

__prog__ = 'common_componentsGUI.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

# Version history
# ---------------
#
from os.path import normpath, isfile

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QLabel, QLineEdit, QComboBox, QPushButton, QCheckBox, QRadioButton, QButtonGroup)

from initialise_funcs import write_study_definition_file, read_config_file, write_config_file

WDGT_SIZE_40 = 40
WDGT_SIZE_60 = 60
WDGT_SIZE_80 = 80
WDGT_SIZE_100 = 100
WDGT_SIZE_120 = 120
WDGT_SIZE_200 = 200

RESOLUTIONS = {120:'30"', 30:'2\'', 20:'3\'', 10:'6\'', 8:'7\' 30"', 6:'10\'', 4:'15\'', 3:'20\'', 2:'30\''}
LU_DEFNS = {'lu_type' : ['Arable','Forestry','Miscanthus','Grassland','Semi-natural', 'SRC', 'Rapeseed', 'Sugar cane'],
                   'abbrev': ['ara',   'for',      'mis',      'gra',      'nat',     'src', 'rps',      'sgc'],
                        'ilu':[1,        3,          5,          2,          4,          6,     7,          7]}

HILDA_LANDUSES = ['cropland', 'pasture', 'other', 'forest', 'grassland', 'all']

# run modes
# =========
SPATIAL = 1
XLS_FILE = 2

# ========================================

def _chck_box_inpt_choices(form, grid, irow):
    """

    """
    irow += 1

    form.w_hilda_lus = {}
    for icol, lu in enumerate(HILDA_LANDUSES):
        w_hilda_lu = QCheckBox(lu.title())
        helpText = ''
        w_hilda_lu.setToolTip(helpText)
        if lu == 'all':
            w_hilda_lu.clicked.connect(form.adjustLuChckBoxes)

        grid.addWidget(w_hilda_lu, irow, icol)
        form.w_hilda_lus[lu] = w_hilda_lu

    return irow

def commonSection(form, grid, irow):
    """
    C
    """

    # default to EWEMBI
    # =================
    # hist_syears, hist_eyears, fut_syears, fut_eyears, scenarios = get_weather_parms(form, 'CRU')
    equimodeDflt = '9.5'
    form.depths = list([30,100]) # soil depths

    luTypes = {}; lu_type_abbrevs = {}
    for lu_type, abbrev, ilu in zip(LU_DEFNS['lu_type'], LU_DEFNS['abbrev'], LU_DEFNS['ilu']):
        luTypes[lu_type] = ilu
        lu_type_abbrevs[lu_type] = abbrev

    form.land_use_types = luTypes
    form.lu_type_abbrevs = lu_type_abbrevs

    # resources
    # =========
    irow += 1
    lbl10w = QLabel('Weather resource:')
    lbl10w.setAlignment(Qt.AlignRight)
    helpText = 'permissable weather dataset resources include CRU, Euro-CORDEX - see: http://www.euro-cordex.net, MERA and EObs'
    lbl10w.setToolTip(helpText)
    grid.addWidget(lbl10w, irow, 0)

    combo10w = QComboBox()
    for weather_resource in form.weather_resources_generic:
        combo10w.addItem(weather_resource)
    combo10w.setFixedWidth(WDGT_SIZE_120)
    form.combo10w = combo10w
    grid.addWidget(combo10w, irow, 1)

    # line 9: scenarios
    # =================
    lbl10 = QLabel('Climate Scenario:')
    lbl10.setAlignment(Qt.AlignRight)
    helpText = 'Ecosse requires future average monthly precipitation and temperature derived from climate models.\n' \
        + 'The data used here is ClimGen v1.02 created on 16.10.08 developed by the Climatic Research Unit\n' \
        + ' and the Tyndall Centre. See: http://www.cru.uea.ac.uk/~timo/climgen/'

    lbl10.setToolTip(helpText)
    grid.addWidget(lbl10, irow, 2)

    # weather scenarios are populated when the config file is read or weather resource is changed
    # ===========================================================================================
    combo10 = QComboBox()
    combo10.setFixedWidth(WDGT_SIZE_80)
    grid.addWidget(combo10, irow, 3)
    form.combo10 = combo10

    # equilibrium mode
    # ================
    lbl12 = QLabel('Equilibrium mode:')
    lbl12.setAlignment(Qt.AlignRight)
    helpText = 'mode of equilibrium run, generally OK with 9.5'
    lbl12.setToolTip(helpText)
    grid.addWidget(lbl12, irow, 4)

    w_equimode = QLineEdit()
    w_equimode.setText(equimodeDflt)
    w_equimode.setFixedWidth(WDGT_SIZE_60)
    grid.addWidget(w_equimode, irow, 5)
    form.w_equimode = w_equimode

    # Historic
    # ========
    irow += 1
    lbl09s = QLabel('Historic start year:')
    lbl09s.setAlignment(Qt.AlignRight)
    helpText = 'Ecosse requires long term average monthly precipitation and temperature\n' \
            + 'which is derived from datasets managed by Climatic Research Unit (CRU).\n' \
            + ' See: http://www.cru.uea.ac.uk/about-cru'
    lbl09s.setToolTip(helpText)
    grid.addWidget(lbl09s, irow, 0)

    combo09s = QComboBox()
    combo09s.setFixedWidth(WDGT_SIZE_60)
    grid.addWidget(combo09s, irow, 1)
    form.combo09s = combo09s

    lbl09e = QLabel('End year:')
    lbl09e.setAlignment(Qt.AlignRight)
    grid.addWidget(lbl09e, irow, 2)

    combo09e = QComboBox()
    combo09e.setFixedWidth(WDGT_SIZE_60)
    grid.addWidget(combo09e, irow, 3)
    form.combo09e = combo09e

    # Simulation years    
    # ================
    irow += 1
    lbl11s = QLabel('Simulation start year:')
    helpText = 'Simulation start and end years determine the number of growing seasons to simulate\n' \
            + 'CRU and CORDEX resources run to 2100 whereas EObs resource runs to 2017'
    lbl11s.setToolTip(helpText)
    lbl11s.setAlignment(Qt.AlignRight)
    grid.addWidget(lbl11s, irow, 0)

    combo11s = QComboBox()
    combo11s.setFixedWidth(WDGT_SIZE_60)
    grid.addWidget(combo11s, irow, 1)
    form.combo11s = combo11s

    lbl11e = QLabel('End year:')
    lbl11e.setAlignment(Qt.AlignRight)
    grid.addWidget(lbl11e, irow, 2)

    combo11e = QComboBox()
    combo11e.setFixedWidth(WDGT_SIZE_60)
    grid.addWidget(combo11e, irow, 3)
    form.combo11e = combo11e
    
    w_ave_weather = QCheckBox('Use average weather')
    helpText = 'Select this option to use average weather, from the CRU year range, for\n' \
               ' the climate file for each of the simulation years'
    w_ave_weather.setToolTip(helpText)
    grid.addWidget(w_ave_weather, irow, 4, 1, 2)
    form.w_ave_weather = w_ave_weather

    irow += 1
    grid.addWidget(QLabel(''), irow, 2)     # spacer

    # ======
    irow = _chck_box_inpt_choices(form, grid, irow)

    # ======
    irow += 1
    grid.addWidget(QLabel(''), irow, 2)  # spacer

    return irow

def save_clicked(form):
        """
        write last GUI selections
        """
        #
        write_config_file(form)
        write_study_definition_file(form)

        return

def exit_clicked(form, write_config_flag = True):

    # write last GUI selections
    if write_config_flag:
        write_config_file(form)
        write_study_definition_file(form)

    # close various files
    if hasattr(form, 'fobjs'):
        for key in form.fobjs:
            form.fobjs[key].close()

    # close logging
    try:
        form.lgr.handlers[0].close()
    except AttributeError:
        pass

    form.close()

    return

def changeConfigFile(form):

    # identify and read the new configuration file
    new_study = form.combo00s.currentText()
    new_config = 'global_ecosse_config_hwsd_' + new_study
    config_file = normpath(form.config_dir + '/' + new_config + '.txt')

    if isfile(config_file):
        form.config_file = config_file
        read_config_file(form)
        form.study = new_study
        form.w_study.setText(new_study)
    else:
        print('Could not locate ' + config_file)

    return

def studyTextChanged(form):

        # replace spaces with underscores and rebuild study list
        # ======================================================
        study = form.w_study.text()

        form.w_study.setText(study.replace(' ','_'))

        return
