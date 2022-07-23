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

    def mz_array(self, mode):
        if mode == 'raw':
            return self.raw_mz_array
        elif mode == 'preprocessed':
            return self.preprocessed_mz_array
        elif mode == 'peak_picked':
            return self.peak_picked_mz_array

    def intensity_array(self, mode):
        if mode == 'raw':
            return self.raw_intensity_array
        elif mode == 'preprocessed':
            return self.preprocessed_intensity_array
        elif mode == 'peak_picked':
            return self.peak_picked_intensity_array

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
