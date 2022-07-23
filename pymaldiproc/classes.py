import os
from pyteomics import mzml as pyt


# Barebones class to hold basic MALDI spectrum data.
class MALDISpectrum(object):
    def __init__(self, pyteomics_dict):
        self.id = None
        self.replicate = None
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
        spectrum_id = pyteomics_dict['spectrum title'].split(' ')[0]
        spectrum_filename = pyteomics_dict['spectrum title'].split(' ')[1]
        self.id = spectrum_id.split('.')[0]
        self.replicate = spectrum_id.split('.')[1]
        self.chip = spectrum_filename.split('\\')[-4].split('_')[0]
        self.spot = spectrum_filename.split('\\')[-4].split('_')[1]

        # spectra
        self.ms_level = pyteomics_dict['ms level']
        self.raw_mz_array = pyteomics_dict['m/z array']
        self.raw_intensity_array = pyteomics_dict['intensity array']


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
        self.consensus_peak_lists = []

        self.import_mzml(input_path)

    def import_mzml(self, input_path):
        # find mzML files
        if input_path.endswith('.mzML'):
            input_files = [input_path]
        elif not input_path.endswith('.mzML') and os.path.isdir(input_path):
            input_files = [os.path.join(dirpath, filename) for dirpath, dirnames, filenames in os.walk(input_path)
                           for filename in filenames if filename.endswith('.mzML')]

        # read in data with pyteomics
        for mzml_filename in input_files:
            mzml_data = list(pyt.read(mzml_filename))
            for scan_dict in mzml_data:
                self.spectra.append(MALDISpectrum(scan_dict))

    #
    def generate_consensus_peak_lists(self, method):
        pass
