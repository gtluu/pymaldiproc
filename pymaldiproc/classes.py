# Barebones class to hold basic MALDI spectrum data.
class MALDISpectrum(object):
    def __init__(self, pyteomics_dict):
        self.id = None
        self.replicate = None
        self.chip = None
        self.spot = None
        self.ms_level = None
        self.centroid = None
        self.mz_array = None
        self.intensity_array = None
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
        if 'positive scan' in pyteomics_dict.keys() and 'centroid scan' not in pyteomics_dict.keys():
            self.centroid = False
        elif 'positive scan' not in pyteomics_dict.keys() and 'centroid scan' in pyteomics_dict.keys():
            self.centroid = False
        self.mz_array = pyteomics_dict['m/z array']
        self.intensity_array = pyteomics_dict['intensity array']
