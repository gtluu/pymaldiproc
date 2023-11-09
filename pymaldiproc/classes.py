import os
import pandas as pd
from uuid import uuid4
from pyTDFSDK.classes import TsfSpectrum, TdfSpectrum


# Barebones class to hold basic MALDI spectrum data.
class OpenMALDISpectrum(object):
    def __init__(self, pyteomics_dict, filename):
        self.source = filename
        self.name = None
        self.uuid = str(uuid4())
        self.spectrum_id = None
        self.coord = None
        self.ms_level = None
        self.centroided = None
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
        if self.source.lower().endswith('mzml'):
            self.name = str(os.path.splitext(os.path.split(self.source)[-1])[0]) + '_' + str(pyteomics_dict['index'])
        elif self.source.lower().endswith('mzxml'):
            self.name = str(os.path.splitext(os.path.split(self.source)[-1])[0]) + '_' + str(pyteomics_dict['num'])
        else:
            self.name = ''

        try:
            self.coord = pyteomics_dict['maldi spot identifier']
            self.spectrum_id = self.name + '|' + self.coord + '|' + self.uuid
        except KeyError:
            self.spectrum_id = self.name + '|' + self.uuid

        if 'centroid spectrum' in pyteomics_dict and 'profile spectrum' not in pyteomics_dict:
            self.centroided = True
        elif 'centroid spectrum' not in pyteomics_dict and 'profile spectrum' in pyteomics_dict:
            self.centroided = False

        # spectra
        if self.source.lower().endswith('mzml'):
            self.ms_level = pyteomics_dict['ms level']
        elif self.source.lower().endswith('mzxml'):
            self.ms_level = pyteomics_dict['msLevel']
        self.raw_mz_array = pyteomics_dict['m/z array']
        self.raw_intensity_array = pyteomics_dict['intensity array']

    # get most up to date mz_array (raw or preprocessed)
    def get_mz_array(self, array_type=None):
        if array_type is None:
            if self.preprocessed_mz_array is None and self.peak_picked_mz_array is None:
                return self.raw_mz_array
            elif self.preprocessed_mz_array is not None and self.peak_picked_mz_array is None:
                return self.preprocessed_mz_array
            elif self.peak_picked_mz_array is not None:
                return self.peak_picked_mz_array
        elif array_type == 'raw':
            return self.raw_mz_array
        elif array_type == 'preprocessed':
            return self.preprocessed_mz_array
        elif array_type == 'peak_picked':
            return self.peak_picked_mz_array

    # get most up to date intensity_array (raw or preprocessed)
    def get_intensity_array(self, array_type=None):
        if array_type is None:
            if self.preprocessed_intensity_array is None and self.peak_picked_intensity_array is None:
                return self.raw_intensity_array
            elif self.preprocessed_intensity_array is not None and self.peak_picked_intensity_array is None:
                return self.preprocessed_intensity_array
            elif self.peak_picked_intensity_array is not None:
                return self.peak_picked_intensity_array
        elif array_type == 'raw':
            return self.raw_intensity_array
        elif array_type == 'preprocessed':
            return self.preprocessed_intensity_array
        elif array_type == 'peak_picked':
            return self.peak_picked_intensity_array

    # get peak picked (centroided) feature list
    def get_peak_list(self):
        if self.peak_picked_mz_array is not None and self.peak_picked_intensity_array is not None:
            return pd.DataFrame(data={'mz': self.peak_picked_mz_array,
                                      'intensity': self.peak_picked_intensity_array})
        else:
            return pd.DataFrame(data={'mz': [0], 'intensity': [0]})


class PMPTsfSpectrum(TsfSpectrum):
    def __init__(self, tsf_data, frame: int, mode: str, profile_bins=0, encoding=64):
        super().__init__(tsf_data, frame, mode, profile_bins, encoding)

        self.source = tsf_data.source_file
        self.name = None
        self.uuid = str(uuid4())
        self.spectrum_id = None
        self.raw_mz_array = self.mz_array
        self.raw_intensity_array = self.intensity_array
        self.preprocessed_mz_array = None
        self.preprocessed_intensity_array = None
        self.peak_picked_mz_array = None
        self.peak_picked_intensity_array = None
        self.data_processing = {}

        self.get_spectrum_metadata()

    def get_spectrum_metadata(self):
        self.name = str(os.path.splitext(os.path.split(self.source)[-1])[0]) + '_' + str(self.frame)
        self.spectrum_id = self.name + '|' + self.coord + '|' + self.uuid

    # get most up to date mz_array (raw or preprocessed)
    def get_mz_array(self, array_type=None):
        if array_type is None:
            if self.preprocessed_mz_array is None and self.peak_picked_mz_array is None:
                return self.raw_mz_array
            elif self.preprocessed_mz_array is not None and self.peak_picked_mz_array is None:
                return self.preprocessed_mz_array
            elif self.peak_picked_mz_array is not None:
                return self.peak_picked_mz_array
        elif array_type == 'raw':
            return self.raw_mz_array
        elif array_type == 'preprocessed':
            return self.preprocessed_mz_array
        elif array_type == 'peak_picked':
            return self.peak_picked_mz_array

    # get most up to date intensity_array (raw or preprocessed)
    def get_intensity_array(self, array_type=None):
        if array_type is None:
            if self.preprocessed_intensity_array is None and self.peak_picked_intensity_array is None:
                return self.raw_intensity_array
            elif self.preprocessed_intensity_array is not None and self.peak_picked_intensity_array is None:
                return self.preprocessed_intensity_array
            elif self.peak_picked_intensity_array is not None:
                return self.peak_picked_intensity_array
        elif array_type == 'raw':
            return self.raw_intensity_array
        elif array_type == 'preprocessed':
            return self.preprocessed_intensity_array
        elif array_type == 'peak_picked':
            return self.peak_picked_intensity_array

    # get peak picked (centroided) feature list
    def get_peak_list(self):
        if self.peak_picked_mz_array is not None and self.peak_picked_intensity_array is not None:
            return pd.DataFrame(data={'mz': self.peak_picked_mz_array,
                                      'intensity': self.peak_picked_intensity_array})
        else:
            return pd.DataFrame(data={'mz': [0], 'intensity': [0]})


class PMPTdfSpectrum(TdfSpectrum):
    def __init__(self, tdf_data, frame: int, mode: str, precursor=0, diapasef_window=None, profile_bins=0, encoding=64,
                 exclude_mobility=False):
        super().__init__(tdf_data, frame, mode, precursor, diapasef_window, profile_bins, encoding, exclude_mobility)

        self.source = tdf_data.source_file
        self.name = None
        self.uuid = str(uuid4())
        self.spectrum_id = None
        self.raw_mz_array = self.mz_array
        self.raw_intensity_array = self.intensity_array
        self.preprocessed_mz_array = None
        self.preprocessed_intensity_array = None
        self.peak_picked_mz_array = None
        self.peak_picked_intensity_array = None
        self.data_processing = {}

        self.get_spectrum_metadata()

    def get_spectrum_metadata(self):
        self.name = str(os.path.splitext(os.path.split(self.source)[-1])[0]) + '_' + str(self.frame)
        self.spectrum_id = self.name + '|' + self.coord + '|' + self.uuid

    # get most up to date mz_array (raw or preprocessed)
    def get_mz_array(self, array_type=None):
        if array_type is None:
            if self.preprocessed_mz_array is None and self.peak_picked_mz_array is None:
                return self.raw_mz_array
            elif self.preprocessed_mz_array is not None and self.peak_picked_mz_array is None:
                return self.preprocessed_mz_array
            elif self.peak_picked_mz_array is not None:
                return self.peak_picked_mz_array
        elif array_type == 'raw':
            return self.raw_mz_array
        elif array_type == 'preprocessed':
            return self.preprocessed_mz_array
        elif array_type == 'peak_picked':
            return self.peak_picked_mz_array

    # get most up to date intensity_array (raw or preprocessed)
    def get_intensity_array(self, array_type=None):
        if array_type is None:
            if self.preprocessed_intensity_array is None and self.peak_picked_intensity_array is None:
                return self.raw_intensity_array
            elif self.preprocessed_intensity_array is not None and self.peak_picked_intensity_array is None:
                return self.preprocessed_intensity_array
            elif self.peak_picked_intensity_array is not None:
                return self.peak_picked_intensity_array
        elif array_type == 'raw':
            return self.raw_intensity_array
        elif array_type == 'preprocessed':
            return self.preprocessed_intensity_array
        elif array_type == 'peak_picked':
            return self.peak_picked_intensity_array

    # get peak picked (centroided) feature list
    def get_peak_list(self):
        if self.peak_picked_mz_array is not None and self.peak_picked_intensity_array is not None:
            return pd.DataFrame(data={'mz': self.peak_picked_mz_array,
                                      'intensity': self.peak_picked_intensity_array})
        else:
            return pd.DataFrame(data={'mz': [0], 'intensity': [0]})