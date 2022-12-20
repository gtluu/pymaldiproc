import os
from functools import reduce
import numpy as np
import pandas as pd
from uuid import uuid4
from pyteomics import mzml as pyt
from scipy.signal import savgol_filter, find_peaks, find_peaks_cwt, peak_widths
from BaselineRemoval import BaselineRemoval
from pyMSpec.smoothing import sg_smooth, apodization, rebin, fast_change, median
from pyMSpec.normalisation import tic, rms, mad, sqrt


# Barebones class to hold basic MALDI spectrum data.
class MALDISpectrum(object):
    def __init__(self, pyteomics_dict, software):
        self.name = None
        self.uuid = str(uuid4())
        self.replicate = '1'  # only looks at replicate for MSConvert data. needs to be manually edited for timsTOF
        self.spectrum_id = None
        self.chip = None
        self.spot = None
        self.ms_level = None
        self.raw_mz_array = None
        self.raw_intensity_array = None
        self.preprocessed_mz_array = None
        self.preprocessed_intensity_array = None
        self.peak_picked_mz_array = None
        self.peak_picked_intensity_array = None
        self.data_processing = {}

        self.parse_pyteomics_dict(pyteomics_dict, software)

    def parse_pyteomics_dict(self, pyteomics_dict, software):
        # metadata
        if 'pwiz' in software:
            spectrum_id = pyteomics_dict['spectrum title'].split(' ')[0]
            spectrum_filename = pyteomics_dict['spectrum title'].split(' ')[1]
            self.name = spectrum_id.split('.')[0]
            self.replicate = spectrum_id.split('.')[1]
            self.chip = spectrum_filename.split('\\')[-4].split('_')[0]
            self.spot = spectrum_filename.split('\\')[-4].split('_')[1]
        elif software == 'psims-writer':  # assume psims usage == TIMSCONVERT
            self.name = str(pyteomics_dict['spectrum title']) + '_' + str(pyteomics_dict['index'])
            self.spot = pyteomics_dict['maldi spot identifier']

        self.spectrum_id = self.name + '|' + self.replicate + '|' + self.uuid

        # spectra
        self.ms_level = pyteomics_dict['ms level']
        self.raw_mz_array = pyteomics_dict['m/z array']
        self.raw_intensity_array = pyteomics_dict['intensity array']

    # get most up to date mz_array (raw or preprocessed)
    def get_mz_array(self):
        if self.preprocessed_mz_array is None and self.peak_picked_mz_array is None:
            return self.raw_mz_array
        elif self.preprocessed_mz_array is not None and self.peak_picked_mz_array is None:
            return self.preprocessed_mz_array
        elif self.peak_picked_mz_array is not None:
            return self.peak_picked_mz_array

    # get most up to date intensity_array (raw or preprocessed)
    def get_intensity_array(self):
        if self.preprocessed_intensity_array is None and self.peak_picked_intensity_array is None:
            return self.raw_intensity_array
        elif self.preprocessed_intensity_array is not None and self.peak_picked_intensity_array is None:
            return self.preprocessed_intensity_array
        elif self.peak_picked_intensity_array is not None:
            return self.peak_picked_intensity_array

    # get peak picked (centroided) feature list
    def get_peak_list(self):
        if self.peak_picked_mz_array is not None and self.peak_picked_intensity_array is not None:
            return pd.DataFrame(data={'mz': self.peak_picked_mz_array,
                                      'intensity': self.peak_picked_intensity_array})
        else:
            return pd.DataFrame(data={'mz': [0], 'intensity': [0]})
