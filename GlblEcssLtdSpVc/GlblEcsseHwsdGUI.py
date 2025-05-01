# -------------------------------------------------------------------------------
# Name:
# Purpose:     Creates a GUI with five adminstrative levels plus country
# Author:      Mike Martin
# Created:     11/12/2015
# Licence:     <your licence>
# -------------------------------------------------------------------------------

__prog__ = 'GlblEcsseHwsdGUI.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

import sys
from os.path import normpath, join
from os import system, getcwd
import subprocess
from time import time
from pandas import DataFrame

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap, QFont
from PyQt5.QtWidgets import (QLabel, QWidget, QApplication, QHBoxLayout, QVBoxLayout, QGridLayout, QLineEdit,
                             QComboBox, QPushButton, QCheckBox, QFileDialog, QTextEdit, QMessageBox)

from common_componentsGUI import (exit_clicked, commonSection, changeConfigFile, studyTextChanged, save_clicked)
from glbl_ecss_cmmn_cmpntsGUI import calculate_grid_cell, grid_resolutions, glblecss_limit_sims

from glbl_ecsse_high_level_sp import generate_banded_sims
from glec_new_high_level_sp import generate_sims
from generate_soil_vars_grid import generate_soil_outputs
from generate_soil_vars_nc import make_soil_nc_outputs
from write_soil_vars_grid import generate_all_soil_metrics
from glbl_ecsse_low_level_fns_sv import fetch_soil_metrics

from weather_datasets import change_weather_resource
from initialise_funcs import read_config_file
from initialise_common_funcs import initiation, build_and_display_studies, write_runsites_config_file
from litter_and_orchidee_fns import fetch_nc_litter
from wthr_generation_fns import generate_all_weather, make_wthr_coords_lookup
from set_up_logging import OutLog

STD_BTN_SIZE_120 = 120
STD_BTN_SIZE_80 = 80
STD_FLD_SIZE_180 = 180
STD_FLD_SIZE_200 = 200
STD_FLD_SIZE_250 = 250

CARBON_VARS = {'TOTAL_BM_LITTER_c': 'total conversion of biomass to litter',
                                                            'TOTAL_LITTER_SOIL_c': 'total litter and soil carbon'}
ERROR_STR = '*** Error *** '
WARN_STR = '*** Warning *** '

class Form(QWidget):
    """
    C
    """
    def __init__(self, parent=None):

        super(Form, self).__init__(parent)

        self.version = 'HWSD_grid'
        initiation(self, '_consol')
        font = QFont(self.font())
        font.setPointSize(font.pointSize() + 2)
        self.setFont(font)

        # The layout is done with the QGridLayout
        grid = QGridLayout()
        grid.setSpacing(10)  # set spacing between widgets

        # line 0
        # ======
        irow = 0
        lbl00 = QLabel('Study:')
        lbl00.setAlignment(Qt.AlignRight)
        grid.addWidget(lbl00, irow, 0)

        w_study = QLineEdit()
        w_study.setFixedWidth(STD_FLD_SIZE_180)
        grid.addWidget(w_study, irow, 1, 1, 2)
        self.w_study = w_study

        lbl00s = QLabel('studies:')
        lbl00s.setAlignment(Qt.AlignRight)
        helpText = 'list of studies'
        lbl00s.setToolTip(helpText)
        grid.addWidget(lbl00s, irow, 3)

        combo00s = QComboBox()
        for study in self.studies:
            combo00s.addItem(study)
        grid.addWidget(combo00s, irow, 4, 1, 2)
        combo00s.currentIndexChanged[str].connect(self.changeConfigFile)
        self.combo00s = combo00s

        # =====================
        irow += 1
        lbl03a = QLabel('Study extents:')
        lbl03a.setAlignment(Qt.AlignRight)
        grid.addWidget(lbl03a, irow, 0)

        '''
        self.w_bbox = QLabel()
        grid.addWidget(self.w_bbox, irow, 1, 1, 5)
        '''

        # =====================
        w_nc_lttr_pshb = QPushButton("NetCDF file of plant litter")
        helpText = 'Select NetCDF file of plant inputs'
        w_nc_lttr_pshb.setToolTip(helpText)
        w_nc_lttr_pshb.setEnabled(True)
        grid.addWidget(w_nc_lttr_pshb, irow, 0)
        w_nc_lttr_pshb.clicked.connect(self.fetchPiNcFile)

        w_nc_lttr_fn = QLabel()
        grid.addWidget(w_nc_lttr_fn, irow, 1, 1, 4)
        self.w_nc_lttr_fn = w_nc_lttr_fn

        irow += 1
        w_nc_extnt = QLabel('')  # number of lat lons
        grid.addWidget(w_nc_extnt, irow, 1, 1, 4)
        self.w_nc_extnt = w_nc_extnt

        # ============= select ochidee var =====================
        irow += 1
        lbl08 = QLabel('Carbon variable:')
        lbl08.setAlignment(Qt.AlignRight)
        helpText = ''
        lbl08.setToolTip(helpText)
        grid.addWidget(lbl08, irow, 0)

        combo08 = QComboBox()
        for carbon_var in CARBON_VARS:
            combo08.addItem(carbon_var)
        combo08.setFixedWidth(STD_FLD_SIZE_180)
        grid.addWidget(combo08, irow, 1)
        combo08.currentIndexChanged[str].connect(self.changeCarbonVar)
        self.combo08 = combo08

        w_var_desc = QLabel('')  # variable description
        grid.addWidget(w_var_desc, irow, 2, 1, 3)
        self.w_var_desc = w_var_desc

        self.carbon_vars = CARBON_VARS

        # ======== PFTs ==========
        irow += 1
        w_lbl_pfts = QLabel('Plant functional types:')
        w_lbl_pfts.setAlignment(Qt.AlignRight)
        grid.addWidget(w_lbl_pfts, irow, 0)

        w_combo_pfts = QComboBox()
        w_combo_pfts.setFixedWidth(STD_FLD_SIZE_250)
        w_combo_pfts.addItem('dummy')
        grid.addWidget(w_combo_pfts, irow, 1, 1, 2)
        w_combo_pfts.currentIndexChanged[str].connect(self.changePlntFncType)
        self.w_combo_pfts = w_combo_pfts

        w_ave_val = QLabel('')
        grid.addWidget(w_ave_val, irow, 3, 1, 2)
        self.w_ave_val = w_ave_val

        irow += 1
        grid.addWidget(QLabel(''), irow, 2)  # spacer

        # soil switches
        # =============
        irow += 1
        lbl04 = QLabel('Options:')
        lbl04.setAlignment(Qt.AlignRight)
        grid.addWidget(lbl04, irow, 0)

        w_use_dom_soil = QCheckBox('Use most dominant soil')
        helpText = 'Each HWSD grid cell can have up to 10 soils. Select this option to use most dominant soil and\n' \
                   ' discard all others. The the most dominant soil is defined as having the highest percentage coverage ' \
                   ' of all the soils for that grid cell'
        w_use_dom_soil.setToolTip(helpText)
        # grid.addWidget(w_use_dom_soil, irow, 1, 1, 2)
        grid.addWidget(w_use_dom_soil, irow, 1)
        self.w_use_dom_soil = w_use_dom_soil

        w_use_high_cover = QCheckBox('Use highest coverage soil')
        helpText = 'Each meta-cell has one or more HWSD mu global keys with each key associated with a coverage expressed \n' \
                   ' as a proportion of the area of the meta cell. Select this option to use the mu global with the highest coverage,\n' \
                   ' discard the others and aggregate their coverages to the selected mu global'
        w_use_high_cover.setToolTip(helpText)
        # grid.addWidget(w_use_high_cover, irow, 3, 1, 2)
        grid.addWidget(w_use_high_cover, irow, 2)
        self.w_use_high_cover = w_use_high_cover

        w_baseline = QCheckBox('Create baseline')
        helpText = 'set PIs to zero'
        w_baseline.setToolTip(helpText)
        # grid.addWidget(w_baseline, irow, 3, 1, 2)
        grid.addWidget(w_baseline, irow, 3)
        self.w_baseline = w_baseline

        # AOI bounding box detail
        # =======================
        irow += 1
        w_lbl07 = QLabel('AOI bounding box:')
        helpText = 'Select NetCDF file of plant inputs'
        w_lbl07.setToolTip(helpText)
        w_lbl07.setAlignment(Qt.AlignRight)
        grid.addWidget(w_lbl07, irow, 0)

        w_hwsd_bbox = QLabel('')
        w_hwsd_bbox.setToolTip(helpText)
        w_hwsd_bbox.setAlignment(Qt.AlignLeft)
        self.w_hwsd_bbox = w_hwsd_bbox
        grid.addWidget(self.w_hwsd_bbox, irow, 1, 1, 5)

        # irow += 1
        # grid.addWidget(QLabel(''), irow, 2)  # spacer

        # create weather and grid resolution
        # ==================================
        irow = commonSection(self, grid, irow)
        irow = grid_resolutions(self, grid, irow)
        irow += 1
        grid.addWidget(QLabel(''), irow, 2)  # spacer

        # command line
        # ============
        irow += 1
        w_create_files = QPushButton("Create sim files")
        helpText = 'Generate ECOSSE simulation file sets corresponding to ordered HWSD global mapping unit set in CSV file'
        w_create_files.setToolTip(helpText)
        # w_create_files.setEnabled(False)
        w_create_files.setFixedWidth(STD_BTN_SIZE_120)
        grid.addWidget(w_create_files, irow, 0, )
        w_create_files.clicked.connect(self.createSimsClicked)
        self.w_create_files = w_create_files

        w_auto_spec = QCheckBox('Auto run Ecosse')
        helpText = 'Select this option to automatically run Ecosse'
        w_auto_spec.setToolTip(helpText)
        grid.addWidget(w_auto_spec, irow, 1)
        self.w_auto_spec = w_auto_spec

        w_run_ecosse = QPushButton('Run Ecosse')
        helpText = 'Select this option to create a configuration file for the spec.py script and run it.\n' \
                   + 'The spec.py script runs the ECOSSE programme'
        w_run_ecosse.setToolTip(helpText)
        w_run_ecosse.setFixedWidth(STD_BTN_SIZE_80)
        w_run_ecosse.clicked.connect(self.runEcosseClicked)
        grid.addWidget(w_run_ecosse, irow, 2)
        self.w_run_ecosse = w_run_ecosse

        w_save = QPushButton("Save")
        helpText = 'Save configuration and study definition files'
        w_save.setToolTip(helpText)
        w_save.setFixedWidth(STD_BTN_SIZE_80)
        grid.addWidget(w_save, irow, 3)
        w_save.clicked.connect(self.saveClicked)

        w_cancel = QPushButton("Cancel")
        helpText = 'Leaves GUI without saving configuration and study definition files'
        w_cancel.setToolTip(helpText)
        w_cancel.setFixedWidth(STD_BTN_SIZE_80)
        grid.addWidget(w_cancel, irow, 4)
        w_cancel.clicked.connect(self.cancelClicked)

        w_exit = QPushButton("Exit", self)
        grid.addWidget(w_exit, irow, 5)
        w_exit.setFixedWidth(STD_BTN_SIZE_80)
        w_exit.clicked.connect(self.exitClicked)

        # =============================================
        irow += 1
        irow = glblecss_limit_sims(self, grid, irow)

        # ================= row 3 ============================
        irow += 1
        icol = 0
        w_wthr_only = QPushButton('Create weather')
        helpText = 'Generate weather only'
        w_wthr_only.setToolTip(helpText)
        w_wthr_only.setEnabled(True)
        grid.addWidget(w_wthr_only, irow, icol)
        w_wthr_only.clicked.connect(self.gnrtWthrClicked)
        self.w_wthr_only = w_wthr_only

        icol += 1
        w_wthr_lookup = QPushButton("Make wthr lookup")
        helpText = 'Make weather coords lookup file'
        w_wthr_lookup.setToolTip(helpText)
        w_wthr_lookup.setFixedWidth(STD_BTN_SIZE_120)
        grid.addWidget(w_wthr_lookup, irow, icol)
        w_wthr_lookup.clicked.connect(self.makeWthrLookupClicked)
        self.w_wthr_lookup = w_wthr_lookup

        icol += 1
        w_soil_outpts = QPushButton("Make soil files")
        helpText = 'Generate CSV data of soil carbon (Dominant), pH and bulk density for the HoliSoils project'
        w_soil_outpts.setToolTip(helpText)
        w_soil_outpts.setFixedWidth(STD_BTN_SIZE_120)
        grid.addWidget(w_soil_outpts, irow, icol)
        w_soil_outpts.clicked.connect(self.genSoilOutptsClicked)
        self.w_soil_outpts = w_soil_outpts

        icol += 1
        w_soil_nc = QPushButton("Make soil NC")
        helpText = 'Generate NetCDF file of soil carbon (Dominant), pH and bulk density for both layers'
        w_soil_nc.setToolTip(helpText)
        w_soil_nc.setFixedWidth(STD_BTN_SIZE_120)
        grid.addWidget(w_soil_nc, irow, icol)
        w_soil_nc.clicked.connect(self.genSoilNcClicked)
        self.w_soil_nc = w_soil_nc

        icol += 1
        w_clear = QPushButton("Del sims")
        helpText = 'Remove all simulation files for this study'
        w_clear.setToolTip(helpText)
        w_clear.setFixedWidth(STD_BTN_SIZE_80)
        w_clear.setEnabled(False)
        grid.addWidget(w_clear, irow, icol)
        w_clear.clicked.connect(self.cleanSimsClicked)

        icol += 1
        w_clear = QPushButton("Clear window", self)
        helpText = 'Clear reporting window'
        w_clear.setToolTip(helpText)
        w_clear.setFixedWidth(STD_BTN_SIZE_120)
        w_clear.clicked.connect(self.clearReporting)
        grid.addWidget(w_clear, irow, icol)

        # ================= row 4 ============================
        irow += 1
        icol = 0
        w_new_sims = QPushButton("Create new sims")
        helpText = 'Generate ECOSSE simulation file sets corresponding to ordered HWSD global mapping unit set in CSV file'
        w_new_sims.setToolTip(helpText)
        # w_new_sims.setEnabled(False)
        w_new_sims.setFixedWidth(STD_BTN_SIZE_120)
        grid.addWidget(w_new_sims, irow, 0, )
        w_new_sims.clicked.connect(self.createNewSims)
        self.w_new_sims = w_new_sims

        icol += 2
        w_soil_all = QPushButton("Make soil CSV")
        helpText = 'Generate CSV data of soil carbon (Dominant) for all metrics'
        w_soil_all.setToolTip(helpText)
        w_soil_all.setFixedWidth(STD_BTN_SIZE_120)
        grid.addWidget(w_soil_all, irow, icol)
        w_soil_all.clicked.connect(lambda: self.genSoilOutptsClicked(True))
        self.w_soil_all = w_soil_all

        icol += 2
        w_check_soil = QPushButton("Check soil CSV")
        helpText = 'Generate CSV data of soil carbon (Dominant) for all metrics'
        w_check_soil.setToolTip(helpText)
        w_check_soil.setFixedWidth(STD_BTN_SIZE_120)
        grid.addWidget(w_check_soil, irow, icol)
        w_check_soil.clicked.connect(self.checkSoilCsv)
        self.w_check_soil = w_check_soil

        # LH vertical box consists of png image
        # =====================================
        lh_vbox = QVBoxLayout()

        lbl20 = QLabel()
        lbl20.setPixmap(QPixmap(self.fname_png))
        lbl20.setScaledContents(True)
        lh_vbox.addWidget(lbl20)

        # add grid consisting of combo boxes, labels and buttons to RH vertical box
        # =========================================================================
        rh_vbox = QVBoxLayout()
        rh_vbox.addLayout(grid)

        # add reporting
        # =============
        bot_hbox = QHBoxLayout()
        w_report = QTextEdit()
        w_report.verticalScrollBar().minimum()
        w_report.setMinimumHeight(250)
        w_report.setMinimumWidth(1000)
        w_report.setStyleSheet('font: bold 10.5pt Courier')  # big jump to 11pt
        bot_hbox.addWidget(w_report, 1)
        self.w_report = w_report

        sys.stdout = OutLog(self.w_report, sys.stdout)
        # sys.stderr = OutLog(self.w_report, sys.stderr, QColor(255, 0, 0))

        # add LH and RH vertical boxes to main horizontal box
        # ===================================================
        main_hbox = QHBoxLayout()
        main_hbox.setSpacing(10)
        main_hbox.addLayout(lh_vbox)
        main_hbox.addLayout(rh_vbox, stretch=1)

        # feed horizontal boxes into the window
        # =====================================
        outer_layout = QVBoxLayout()
        outer_layout.addLayout(main_hbox)
        outer_layout.addLayout(bot_hbox)
        self.setLayout(outer_layout)

        # posx, posy, width, height
        self.setGeometry(200, 100, 690, 250)
        self.setWindowTitle('Global Ecosse Holisoils spatial variation - uses EFISCEN NetCDF plant inputs')

        # reads and set values from last run
        # ==================================
        read_config_file(self)
        if len(self.weather_set_linkages) == 0:
            self.w_wthr_only.setEnabled(False)
        else:
            self.w_wthr_only.setEnabled(True)

        self.combo10w.currentIndexChanged[str].connect(self.weatherResourceChanged)

    def createNewSims(self):
        """

        """
        calculate_grid_cell(self)
        generate_sims(self)
        return

    def checkSoilCsv(self):
        """

        """
        fetch_soil_metrics(self)
        return

    def makeWthrLookupClicked(self):
        """
        Make weather coords lookup file
        """
        make_wthr_coords_lookup(self)

        return
    def gnrtWthrClicked(self):
        """
        generate weather for all regions, scenarios and GCMs
        """
        generate_all_weather(self)

        return

    def changeCarbonVar(self):
        """

        """
        carbon_var = self.combo08.currentText()
        self.w_var_desc.setText(self.carbon_vars[carbon_var])
        fname = self.w_nc_lttr_fn.text()
        fetch_nc_litter(self, fname)

    def viewRunReport(self):
        """
        C
        """
        if self.band_reports is None:
            print(WARN_STR + 'Nothing to report')
            QApplication.processEvents()
            return

        notepad_flag = True
        dictr = {}
        for nline, line in enumerate(self.band_reports):
            if line is None:
                continue

            atoms = line.split()
            if nline == 0:
                headers = [atoms[0][0:-1], atoms[2] + ' yes', atoms[2] + ' no', atoms[8][0:-1], 'no PIs',
                                                                                                    atoms[-2][0:-1]]
                dictr = {field: [] for field in headers}

            dictr['Band'].append(atoms[1])
            dictr['forest yes'].append(atoms[4])
            dictr['forest no'].append(atoms[6])
            dictr['weather'].append(atoms[9])
            dictr['no PIs'].append(atoms[13])
            dictr['completed'].append(atoms[15])

        if len(dictr) == 0:
            print(WARN_STR + 'Nothing to report')
            QApplication.processEvents()
            return

        dictr_df = DataFrame(dictr)
        if notepad_flag:
            scrtch_file = join(self.settings['log_dir'], 'run_report')
            dictr_df.to_csv(scrtch_file, index=False, sep='\t')
            ret_code = subprocess.run(['notepad.exe', scrtch_file])
        else:
            mess_content = dictr_df.to_string(index=False, justify='center', col_space=10)
            w_mess_box = QMessageBox()
            w_mess_box.setWindowTitle("Banded simulations report")
            w_mess_box.setText(mess_content)
            w_mess_box.setStandardButtons(QMessageBox.Cancel)
            ret_code = w_mess_box.exec()
        return

    def genSoilOutptsClicked(self, all_metrics_flag=False):
        """

        """
        if all_metrics_flag:
            generate_all_soil_metrics(self)
        else:
            generate_soil_outputs(self)

    def genSoilNcClicked(self):
        """

        """
        make_soil_nc_outputs(self)

    def changePlntFncType(self):
        """
        C
        """
        pft_name = self.w_combo_pfts.currentText()
        if pft_name != '':
            pft_key = '00'
            ave_val = self.litter_defn.aves[pft_key]
            mess = 'average value: ' + str(round(float(ave_val), 2))
            self.w_ave_val.setText(mess)

    def cleanSimsClicked(self):
        """
        C
        """
        print('under construction')

    def clearReporting(self):
        """
        C
        """
        self.w_report.clear()

    def adjustLuChckBoxes(self):
        """
        C
        """
        for lu in self.w_hilda_lus:
            if lu == 'all':
                continue
            else:
                if self.w_hilda_lus['all'].isChecked():
                    self.w_hilda_lus[lu].setEnabled(False)
                else:
                    self.w_hilda_lus[lu].setEnabled(True)
        return

    def weatherResourceChanged(self):
        """
        C
        """
        change_weather_resource(self)

    def resolutionChanged(self):
        """
        C
        """
        granularity = 120
        calculate_grid_cell(self, granularity)

    def studyTextChanged(self):
        """
        C
        """
        studyTextChanged(self)

    def createSimsClicked(self):
        """
        C
        """
        study = self.w_study.text()
        if study == '':
            print('study cannot be blank')
            return

        # check for spaces
        # ================
        if study.find(' ') >= 0:
            print('*** study name must not have spaces ***')
            return

        self.study = study

        generate_banded_sims(self)

        # run further steps...
        if self.w_auto_spec.isChecked():
            self.runEcosseClicked()

    def runEcosseClicked(self):
        """
        C
        """
        # components of the command string have been checked at startup
        # =============================================================
        if write_runsites_config_file(self):
            # run the make simulations script
            # ===============================
            print('Working dir: ' + getcwd())
            start_time = time()
            cmd_str = self.python_exe + ' ' + self.runsites_py + ' ' + self.runsites_config_file
            system(cmd_str)
            end_time = time()
            print('Time taken: {}'.format(round(end_time - start_time)))

    def saveClicked(self):
        """
        C
        """
        func_name = __prog__ + ' saveClicked'

        # check for spaces
        # ================
        study = self.w_study.text()
        if study == '':
            print('study cannot be blank')
        else:
            if study.find(' ') >= 0:
                print('*** study name must not have spaces ***')
            else:
                save_clicked(self)
                build_and_display_studies(self)

    def cancelClicked(self):
        """
        C
        """
        func_name = __prog__ + ' cancelClicked'

        exit_clicked(self, write_config_flag=False)

    def exitClicked(self):
        """
        exit cleanly
        """
        # check for spaces
        # ================
        study = self.w_study.text()
        if study == '':
            print('study cannot be blank')
        else:
            if study.find(' ') >= 0:
                print('*** study name must not have spaces ***')
            else:
                exit_clicked(self)

    def changeConfigFile(self):
        """
        permits change of configuration file
        """
        changeConfigFile(self)

    def fetchPiNcFile(self):
        """
        Select NetCDF file of plant inputs
        if user cancels then fname is returned as an empty string
        """
        fname = self.w_nc_lttr_fn.text()
        fname, dummy = QFileDialog.getOpenFileName(self, 'Select NetCDF of plant inputs', fname, 'NetCDF file (*.nc)')
        if fname != '':
            fname = normpath(fname)
            self.w_nc_lttr_fn.setText(fname)
            fetch_nc_litter(self, fname)

def main():
    """
    C
    """
    app = QApplication(sys.argv)  # create QApplication object
    form = Form()  # instantiate form
    # display the GUI and start the event loop if we're not running batch mode
    form.show()  # paint form
    sys.exit(app.exec_())  # start event loop

if __name__ == '__main__':
    main()
