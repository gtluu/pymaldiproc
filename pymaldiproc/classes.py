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
    def __init__(self, pyteomics_dict):
        self.name = None
        self.uuid = str(uuid4())
        self.replicate = '1'
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

        self.parse_pyteomics_dict(pyteomics_dict)

    def parse_pyteomics_dict(self, pyteomics_dict):
        # metadata
        # currently only supports MSConvert style metadata
        if 'spectrum title' in pyteomics_dict.keys():
            spectrum_id = pyteomics_dict['spectrum title'].split(' ')[0]
            spectrum_filename = pyteomics_dict['spectrum title'].split(' ')[1]
            self.name = spectrum_id.split('.')[0]
            self.replicate = spectrum_id.split('.')[1]
            self.chip = spectrum_filename.split('\\')[-4].split('_')[0]
            self.spot = spectrum_filename.split('\\')[-4].split('_')[1]
        # placeholder else clause
        # TODO: update TIMSCONVERT to include spectrum title
        else:
            self.name = pyteomics_dict['pymaldiproc_spectra_count']
            self.replicate = '1'

        self.spectrum_id = self.name + '|' + self.replicate + '|' + self.uuid

        # spectra
        self.ms_level = pyteomics_dict['ms level']
        self.raw_mz_array = pyteomics_dict['m/z array']
        self.raw_intensity_array = pyteomics_dict['intensity array']

    # get most up to date mz_array (raw or preprocessed)
    def get_mz_array(self):
        if self.preprocessed_mz_array is None:
            return self.raw_mz_array
        elif self.preprocessed_mz_array is not None:
            return self.preprocessed_mz_array

    # get most up to date intensity_array (raw or preprocessed)
    def get_intensity_array(self):
        if self.preprocessed_intensity_array is None:
            return self.raw_intensity_array
        elif self.preprocessed_intensity_array is not None:
            return self.preprocessed_intensity_array

    # get peak picked (centroided) feature list
    def get_peak_list(self):
        if self.peak_picked_mz_array is not None and self.peak_picked_intensity_array is not None:
            return pd.DataFrame(data={'mz': self.peak_picked_mz_array,
                                      'intensity': self.peak_picked_intensity_array})
        else:
            return pd.DataFrame(data={'mz': [0], 'intensity': [0]})


class ConsensusPeakList(object):
    def __init__(self, consensus_peak_list_dict):
        self.id = consensus_peak_list_dict['id']
        self.ms_level = consensus_peak_list_dict['ms_level']
        self.consensus_mz_array = consensus_peak_list_dict['mz_array']
        self.consensus_intensity_array = consensus_peak_list_dict['intensity_array']
        self.data_processing = consensus_peak_list_dict['data_processing']


class MALDIDataset(object):
    def __init__(self, input_path):
        self.spectra = []
        #self.spectra = {}
        self.feature_matrix = None
        self.cos_distance_matrix = None
        self.cos_distance_matrix_shape = None

        self.import_mzml(input_path)
        #self.spectra_ids = [spectrum.spectrum_id for spectrum in self.spectra]
        #self.spectra_ids = [spectrum_id + '|' + spectrum_replicate
        #                    for spectrum_id, value in self.spectra.items()
        #                    for spectrum_replicate, value in self.spectra[spectrum_id].items()]

    def import_mzml(self, input_path):
        # find mzML files
        if input_path.endswith('.mzML'):
            input_files = [input_path]
        elif not input_path.endswith('.mzML') and os.path.isdir(input_path):
            input_files = [os.path.join(dirpath, filename) for dirpath, dirnames, filenames in os.walk(input_path)
                           for filename in filenames if filename.endswith('.mzML')]

        # read in data with pyteomics
        spectra_count = 1
        for mzml_filename in input_files:
            mzml_data = list(pyt.read(mzml_filename))
            for scan_dict in mzml_data:
                scan_dict['pymaldiproc_spectra_count'] = str(spectra_count)
                spectrum = MALDISpectrum(scan_dict)
                self.spectra.append(spectrum)
                #self.spectra[spectrum.id] = {spectrum.replicate: spectrum}
                spectra_count += 1

    def get_feature_matrix(self, missing_value_imputation=True):
        spectra_dfs_peak_picked = []
        spectra_dfs_preprocessed = []
        for spectrum in self.spectra:
            spectra_dfs_peak_picked.append(pd.DataFrame(data={'mz': spectrum.peak_picked_mz_array,
                                                              spectrum.spectrum_id: spectrum.peak_picked_intensity_array}))
            spectra_dfs_preprocessed.append(pd.DataFrame(data={'mz': spectrum.preprocessed_mz_array,
                                                               spectrum.spectrum_id: spectrum.preprocessed_intensity_array}))
        self.feature_matrix = reduce(lambda x, y: pd.merge(x, y, how='outer', on='mz'), spectra_dfs_peak_picked).sort_values(by='mz')
        if missing_value_imputation:
            ref_matrix = reduce(lambda x, y: pd.merge(x, y, how='outer', on='mz'), spectra_dfs_preprocessed).sort_values(by='mz')
            for colname in self.feature_matrix.columns:
                if colname != 'mz':
                    tmp_df = pd.merge(self.feature_matrix[['mz', colname]],
                                      ref_matrix[['mz', colname]],
                                      how='left',
                                      on='mz')
                    #self.feature_matrix[colname].fillna(tmp_df.drop('mz', axis=1).mean(axis=1), inplace=True)
                    self.feature_matrix[colname] = tmp_df.drop('mz', axis=1).mean(axis=1).values

    def get_cos_distance_matrix(self):
        pass