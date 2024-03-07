import os
import pandas as pd
from uuid import uuid4
from pyTDFSDK.classes import TsfSpectrum, TdfSpectrum


# Barebones class to hold basic MALDI spectrum data.
# uses maldi data to get the m/z, intensity arrays and peak picked list

# -----changes-------
# use get() instead of lower().endswith() to grab strings -> got rid of try except
# simplify class functions : more efficent with less if-else
# initialized raw_mz_array and raw_intensity_array as an empty list 

class OpenMALDISpectrum(object):    # for mzml files

    def __init__(self, pyteomics_dict: dict, filename: str):
        self.source = filename  #store filename
        self.name = None    
        self.uuid = str(uuid4())    # generate identifier for spectrum instance
        self.spectrum_id = None
        self.coord = None
        self.ms_level = None
        self.centroided = None
        self.raw_mz_array = []  # empty list
        self.raw_intensity_array = []   # empty list
        self.preprocessed_mz_array = None
        self.preprocessed_intensity_array = None
        self.peak_picked_mz_array = None
        self.peak_picked_intensity_array = None
        self.data_processing = {}

        self.parse_pyteomics_dict(pyteomics_dict)   # parse to fill dict 

    # use dict to parse metadata and spectra
    def parse_pyteomics_dict(self, pyteomics_dict: dict):
        # metadata
        base_name = os.path.splitext(os.path.basename(self.source))[0]
        index_or_num = pyteomics_dict.get('index', pyteomics_dict.get('num', ''))
        self.name = f"{base_name}_{index_or_num}"   #set spectrum name based on file

        # dont need try except bc get() doesnt error if empty, returns None
        self.coord = pyteomics_dict.get('maldi spot identifier', '')
        self.spectrum_id = f"{self.name}|{self.coord}|{self.uuid}"

        # spectra
        self.centroided = 'centroid spectrum' in pyteomics_dict and 'profile spectrum' not in pyteomics_dict    # check for peak 
        self.ms_level = pyteomics_dict.get('ms level', pyteomics_dict.get('msLevel', self.ms_level))    # get ms lvls from dict
        # fill arrays from dict
        self.raw_mz_array = pyteomics_dict.get('m/z array', []) 
        self.raw_intensity_array = pyteomics_dict.get('intensity array', [])

    def get_mz_array(self, array_type: str = None) -> list:
    # return array if it exists and is not None.
        if array_type is not None:
            specific_array = getattr(self, f"{array_type}_mz_array", None)
            if specific_array is not None:
                return specific_array
            else:
                return [] 
        # return most updated array if no specific type given
        for array_attr in ['peak_picked_mz_array', 'preprocessed_mz_array', 'raw_mz_array']:
            array = getattr(self, array_attr, None)
            if array is not None:
                return array
        return []        # empty list if arrays are NONE or if they didnt fall in the above criteria


    # get most up to date intensity_array (raw or preprocessed)
    def get_intensity_array(self, array_type=None):
    # get rid of usage of logical operators on array
        if array_type is None:
            if self.peak_picked_intensity_array is not None:
                return self.peak_picked_intensity_array
            elif self.preprocessed_intensity_array is not None:
                return self.preprocessed_intensity_array
            else:
                return self.raw_intensity_array
        elif array_type == 'raw':
            return self.raw_intensity_array
        elif array_type == 'preprocessed':
            return self.preprocessed_intensity_array
        elif array_type == 'peak_picked':
            return self.peak_picked_intensity_array


    # get peak picked (centroided) feature list
    # if not available -> still return placeholder dataframe
    def get_peak_list(self) -> pd.DataFrame:
        if self.peak_picked_mz_array is not None and self.peak_picked_intensity_array is not None:
            return pd.DataFrame({'mz': self.peak_picked_mz_array, 'intensity': self.peak_picked_intensity_array})
        return pd.DataFrame({'mz': [0], 'intensity': [0]})  # if no data -> placeholder
    

class PMPTsfSpectrum(TsfSpectrum):  # for tsf files
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
    def get_mz_array(self, array_type: str = None) -> list:
        if array_type is None:
            return self.peak_picked_mz_array or self.preprocessed_mz_array or self.raw_mz_array
        return getattr(self, f"{array_type}_mz_array", [])

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