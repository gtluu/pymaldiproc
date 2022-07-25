import os
import numpy as np
import pandas as pd
from pyteomics import mzml as pyt
from scipy.signal import savgol_filter, find_peaks, find_peaks_cwt, peak_widths
from BaselineRemoval import BaselineRemoval
from pyMSpec.smoothing import sg_smooth, apodization, rebin, fast_change, median
from pyMSpec.normalisation import tic, rms, mad, sqrt


# Barebones class to hold basic MALDI spectrum data.
class MALDISpectrum(object):
    def __init__(self, pyteomics_dict):
        self.id = None
        self.replicate = 1
        self.chip = None
        self.spot = None
        self.ms_level = None
        self.mz_array = self.__get_mz_array()
        self.intensity_array = self.__get_intensity_array()
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
            self.id = spectrum_id.split('.')[0]
            self.replicate = spectrum_id.split('.')[1]
            self.chip = spectrum_filename.split('\\')[-4].split('_')[0]
            self.spot = spectrum_filename.split('\\')[-4].split('_')[1]
        # placeholder else clause
        # TODO: update TIMSCONVERT to include spectrum title
        else:
            self.id = pyteomics_dict['pymaldiproc_spectra_count']
            self.replicate = '1'

        # spectra
        self.ms_level = pyteomics_dict['ms level']
        self.raw_mz_array = pyteomics_dict['m/z array']
        self.raw_intensity_array = pyteomics_dict['intensity array']

    # get most up to date mz_array (raw or preprocessed)
    def __get_mz_array(self):
        if self.preprocessed_mz_array is None:
            return self.raw_mz_array
        elif self.preprocessed_mz_array is not None:
            return self.preprocessed_mz_array

    # get most up to date intensity_array (raw or preprocessed)
    def __get_intensity_array(self):
        if self.preprocessed_intensity_array is None:
            return self.raw_intensity_array
        elif self.preprocessed_intensity_array is not None:
            return self.preprocessed_intensity_array

    def trim_spectrum(self, lower_mass_range, upper_mass_range):
        spectrum_df = pd.DataFrame(data={'mz': self.mz_array, 'intensity': self.intensity_array})
        spectrum_df = spectrum_df[~(spectrum_df['mz'] <= lower_mass_range) & ~(spectrum_df['mz'] >= upper_mass_range)]
        self.preprocessed_mz_array = spectrum_df['mz'].values
        self.preprocessed_intensity_array = spectrum_df['intensity'].values
        self.data_processing = {'lower mass range': lower_mass_range,
                                'upper mass range': upper_mass_range}

    def transform_intensity(self, method='sqrt'):
        if method == 'sqrt':
            self.preprocessed_intensity_array = np.sqrt(self.intensity_array)
        elif method == 'log':
            self.preprocessed_intensity_array = np.log(self.intensity_array)
        elif method == 'log10':
            self.preprocessed_intensity_array = np.log10(self.intensity_array)
        elif method == 'log2':
            self.preprocessed_intensity_array = np.log2(self.intensity_array)
        else:
            raise Exception('Method must be "sqrt", "log", "log2", or "log10"')
        self.data_processing['intensity transformation'] = {'method': method}

    # TODO: implement MovingAverage method
    def smooth_baseline(self, method='SavitzkyGolay', window_length=20, polyorder=3, delta_mz=0.2, diff_thresh=0.01):
        self.data_processing['baseline smoothing'] = {'method': method}
        if method == 'SavitzkyGolay':
            self.preprocessed_intensity_array = savgol_filter(self.intensity_array,
                                                              window_length=window_length,
                                                              polyorder=polyorder)
            self.data_processing['baseline_smoothing']['window length'] = window_length
            self.data_processing['baseline smoothing']['polyorder'] = polyorder
        elif method == 'MovingAverage':
            pass
        elif method == 'apodization':
            self.preprocessed_mz_array, self.preprocessed_intensity_array = apodization(self.mz_array,
                                                                                        self.intensity_array,
                                                                                        w_size=window_length)
            self.data_processing['baseline smoothing']['window length'] = window_length
        elif method == 'rebin':
            self.preprocessed_mz_array, self.preprocessed_intensity_array = rebin(self.mz_array,
                                                                                  self.intensity_array,
                                                                                  delta_mz=delta_mz)
            self.data_processing['baseline smoothing']['delta m/z'] = delta_mz
        elif method == 'fast_change':
            self.preprocessed_mz_array, self.preprocessed_intensity_array = fast_change(self.mz_array,
                                                                                        self.intensity_array,
                                                                                        diff_thresh=diff_thresh)
            self.data_processing['baseline smoothing']['difference threshold'] = diff_thresh
        elif method == 'median':
            self.preprocessed_mz_array, self.preprocessed_intensity_array = median(self.mz_array,
                                                                                   self.intensity_array,
                                                                                   w_size=window_length)
            self.data_processing['baseline smoothing']['window length'] = window_length
        else:
            #raise Exception('Method must be "SavitzkyGolay", "MovingAverage", "apodization", "rebin", "fast_change", or "median"')
            raise Exception('Method must be "SavitzkyGolay", "apodization", "rebin", "fast_change", or "median"')

    # TODO: implement other algorithms
    def remove_baseline(self, method='ZhangFit'):
        if method == 'ZhangFit':
            self.preprocessed_intensity_array = BaselineRemoval(self.intensity_array).ZhangFit()
        elif method == 'SNIP':
            pass
        elif method == 'TopHat':
            pass
        elif method == 'convexhull':
            pass
        elif method == 'median':
            pass
        elif method == 'modpoly':
            pass
        elif method == 'imodpoly':
            pass
        else:
            raise Exception('Method must be "snip", "tophat", "convexhull", "median", "ZhangFit", "modpoly", or "imodpoly"')
        self.data_processing['baseline removal'] = {'method': method}

    # TODO: implement PQN and median from MALDIquant
    def normalize_intensity(self, method='tic'):
        if method == 'tic':
            self.preprocessed_intensity_array = tic(self.mz_array, self.intensity_array)
        elif method == 'rms':
            self.preprocessed_intensity_array = rms(self.mz_array, self.intensity_array)
        elif method == 'mad':
            self.preprocessed_intensity_array = mad(self.mz_array, self.intensity_array)
        elif method == 'sqrt':
            self.preprocessed_intensity_array = sqrt(self.mz_array, self.intensity_array)
        else:
            raise Exception('Method must be "tic", "rms", "mad", or "sqrt"')
        self.data_processing['intensity normalization'] = {'method': method}

    def estimate_peak_widths(self):
        peak_indices, peak_properties = find_peaks(self.intensity_array)
        widths = peak_widths(self.intensity_array, peak_indices)
        return widths[0]

    def peak_picking(self, method='cwt', widths=None, snr=3):
        self.data_processing['peak picking'] = {'method': method}
        if method == 'locmax':
            peak_indices, peak_properties = find_peaks(self.intensity_array)
            self.peak_picked_mz_array = self.mz_array[peak_indices]
            self.peak_picked_intensity_array = self.mz_array[peak_indices]
        elif method == 'cwt':
            if widths is None:
                estimated_widths = self.estimate_peak_widths()
                widths_start = np.min(estimated_widths)
                widths_stop = 2 * np.mean(estimated_widths)
                widths_step = (widths_stop - widths_start) / 10
                widths = np.arange(widths_start, widths_stop, widths_step)
            peak_indices = find_peaks_cwt(self.intensity_array, widths, min_snr=snr)
            self.peak_picked_mz_array = self.mz_array[peak_indices]
            self.peak_picked_intensity_array = self.intensity_array[peak_indices]
            self.data_processing['peak picking']['lower peak width'] = np.min(widths)
            self.data_processing['peak picking']['upper peak width'] = np.max(widths)


class ConsensusPeakList(object):
    def __init__(self, consensus_peak_list_dict):
        self.id = consensus_peak_list_dict['id']
        self.ms_level = consensus_peak_list_dict['ms_level']
        self.consensus_mz_array = consensus_peak_list_dict['mz_array']
        self.consensus_intensity_array = consensus_peak_list_dict['intensity_array']
        self.data_processing = consensus_peak_list_dict['data_processing']


class MALDIDataset(object):
    def __init__(self, input_path):
        #self.spectra = []
        self.spectra = {}
        self.consensus_peak_lists = []

        self.import_mzml(input_path)
        self.spectra_ids = [spectrum_id + '|' + spectrum_replicate
                            for spectrum_id, value in self.spectra.items()
                            for spectrum_replicate, value in self.spectra[spectrum_id].items()]

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
                #self.spectra.append(spectrum)
                self.spectra[spectrum.id] = {spectrum.replicate: spectrum}
                spectra_count += 1

    def apply_preprocessing(self, preprocessing_steps):
        pass

    #
    def generate_consensus_peak_lists(self, method):
        pass
