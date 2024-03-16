import os
import copy
import gc
import numpy as np
import pandas as pd
from uuid import uuid4
from pyTDFSDK.classes import TsfSpectrum, TdfSpectrum
from functools import reduce
from scipy.signal import savgol_filter, find_peaks, find_peaks_cwt, peak_widths
from scipy.stats import median_abs_deviation
from pyMSpec.smoothing import sg_smooth, apodization, rebin, fast_change, median
from pybaselines.smooth import snip, noise_median
from pybaselines.morphological import tophat
from BaselineRemoval import BaselineRemoval
from pyMSpec.normalisation import tic, rms, mad, sqrt
from icoshift import icoshift
import seaborn as sns
import matplotlib.pyplot as plt


class PMPMethods(object):
    # get peak picked (centroided) feature list
    def get_peak_list(self):
        if self.peak_picked_mz_array is None and self.peak_picked_intensity_array is None:
            # If no peak list found, peak picking with default settings is used.
            self.peak_picking()
        return pd.DataFrame(data={'mz': self.peak_picked_mz_array,
                                  'intensity': self.peak_picked_intensity_array})

    def trim_spectrum(self, lower_mass_range, upper_mass_range):
        indices = np.where((self.preprocessed_mz_array >= lower_mass_range) &
                           (self.preprocessed_mz_array <= upper_mass_range))[0]
        self.preprocessed_mz_array = copy.deepcopy(self.preprocessed_mz_array[indices])
        self.preprocessed_intensity_array = copy.deepcopy(self.preprocessed_intensity_array[indices])
        self.data_processing['spectrum trimming'] = {'lower mass range': lower_mass_range,
                                                     'upper mass range': upper_mass_range}
        gc.collect()

    def transform_intensity(self, method='sqrt'):
        # check method
        if method not in ['sqrt', 'log', 'log2', 'log10']:
            raise Exception('Method must be "sqrt", "log", "log2", or "log10"')
        else:
            if method == 'sqrt':
                self.preprocessed_intensity_array = np.sqrt(copy.deepcopy(self.preprocessed_intensity_array))
            elif method == 'log':
                self.preprocessed_intensity_array = np.log(copy.deepcopy(self.preprocessed_intensity_array))
            elif method == 'log2':
                self.preprocessed_intensity_array = np.log2(copy.deepcopy(self.preprocessed_intensity_array))
            elif method == 'log10':
                self.preprocessed_intensity_array = np.log10(copy.deepcopy(self.preprocessed_intensity_array))
            self.data_processing['intensity transformation'] = {'method': method}
            gc.collect()

    def smooth_baseline(self, method='SavitzkyGolay', window_length=21, polyorder=3, delta_mz=0.2, diff_thresh=0.01):
        # check method
        if method not in ['SavitzkyGolay', 'apodization', 'rebin', 'fast_change', 'median']:
            raise Exception('Method must be "SavitzkyGolay", "apodization", "rebin", "fast_change", "median"')
        else:
            self.data_processing['baseline smoothing'] = {'method': method}
            if method == 'SavitzkyGolay':
                self.preprocessed_intensity_array = savgol_filter(copy.deepcopy(self.preprocessed_intensity_array),
                                                                  window_length=window_length,
                                                                  polyorder=polyorder)
                self.data_processing['baseline smoothing']['window length'] = window_length
                self.data_processing['baseline smoothing']['polyorder'] = polyorder
            elif method == 'apodization':
                self.preprocessed_mz_array, self.preprocessed_intensity_array = apodization(
                    self.preprocessed_mz_array,
                    copy.deepcopy(self.preprocessed_intensity_array),
                    w_size=window_length)
                self.data_processing['baseline smoothing']['window length'] = window_length
            elif method == 'rebin':
                self.preprocessed_mz_array, self.preprocessed_intensity_array = rebin(
                    self.preprocessed_mz_array,
                    copy.deepcopy(self.preprocessed_intensity_array),
                    delta_mz=delta_mz)
                self.data_processing['baseline smoothing']['delta m/z'] = delta_mz
            elif method == 'fast_change':
                self.preprocessed_mz_array, self.preprocessed_intensity_array = fast_change(
                    self.preprocessed_mz_array,
                    copy.deepcopy(self.preprocessed_intensity_array),
                    diff_thresh=diff_thresh)
                self.data_processing['baseline smoothing']['difference threshold'] = diff_thresh
            elif method == 'median':
                self.preprocessed_mz_array, self.preprocessed_intensity_array = median(
                    self.preprocessed_mz_array,
                    copy.deepcopy(self.preprocessed_intensity_array),
                    w_size=window_length)
                self.data_processing['baseline smoothing']['window length'] = window_length
            gc.collect()

    def remove_baseline(self, method='SNIP', min_half_window=1, max_half_window=100, decreasing=True,
                        smooth_half_window=None, filter_order=2, sigma=None, increment=1, max_hits=1,
                        window_tol=0.000001, lambda_=100, porder=1, repetition=None, degree=2, gradient=0.001):
        # check method
        if method not in ['SNIP', 'TopHat', 'Median', 'ZhangFit', 'ModPoly', 'IModPoly']:
            raise Exception('Method must be "SNIP", "TopHat", "Median", "ZhangFit", "ModPoly", or "IModPoly"')
        else:
            self.data_processing['baseline removal'] = {'method': method}
            if method == 'SNIP':
                baseline = snip(data=copy.deepcopy(self.preprocessed_intensity_array),
                                max_half_window=max_half_window,
                                decreasing=decreasing,
                                smooth_half_window=smooth_half_window,
                                filter_order=filter_order)[0]
                self.preprocessed_intensity_array = copy.deepcopy(self.preprocessed_intensity_array) - baseline
                self.data_processing['baseline removal']['max half window'] = max_half_window
                self.data_processing['baseline removal']['decreasing'] = decreasing
                self.data_processing['baseline removal']['smooth half window'] = smooth_half_window
                self.data_processing['baseline removal']['filter order'] = filter_order
            elif method == 'TopHat':
                baseline = tophat(data=copy.deepcopy(self.preprocessed_intensity_array),
                                  half_window=max_half_window,
                                  increment=increment,
                                  max_hits=max_hits,
                                  window_tol=window_tol,
                                  max_half_window=max_half_window,
                                  min_half_window=min_half_window)[0]
                self.preprocessed_intensity_array = copy.deepcopy(self.preprocessed_intensity_array) - baseline
                self.data_processing['baseline removal']['half window'] = max_half_window
                self.data_processing['baseline removal']['increment'] = increment
                self.data_processing['baseline removal']['max hits'] = max_hits
                self.data_processing['baseline removal']['window tolerance'] = window_tol
                self.data_processing['baseline removal']['max half window'] = max_half_window
                self.data_processing['baseline removal']['min half window'] = min_half_window
            elif method == 'Median':
                baseline = noise_median(data=copy.deepcopy(self.preprocessed_intensity_array),
                                        half_window=max_half_window,
                                        smooth_half_window=smooth_half_window,
                                        sigma=sigma)[0]
                self.preprocessed_intensity_array = copy.deepcopy(self.preprocessed_intensity_array) - baseline
                self.data_processing['baseline removal']['half window'] = max_half_window
                self.data_processing['baseline removal']['smooth half window'] = smooth_half_window
                self.data_processing['baseline removal']['sigma'] = sigma
            elif method == 'ZhangFit':
                if repetition is None:
                    repetition = 15
                self.preprocessed_intensity_array = BaselineRemoval(copy.deepcopy(
                    self.preprocessed_intensity_array)).ZhangFit(
                    lambda_=lambda_,
                    porder=porder,
                    repitition=repetition)
                self.data_processing['baseline removal']['lambda'] = lambda_
                self.data_processing['baseline removal']['porder'] = porder
                self.data_processing['baseline removal']['repetition'] = repetition
            elif method == 'ModPoly':
                if repetition is None:
                    repetition = 100
                self.preprocessed_intensity_array = BaselineRemoval(copy.deepcopy(
                    self.preprocessed_intensity_array)).ModPoly(
                    degree=degree,
                    repitition=repetition,
                    gradient=gradient)
                self.data_processing['baseline removal']['degree'] = degree
                self.data_processing['baseline removal']['repetition'] = repetition
                self.data_processing['baseline removal']['gradient'] = gradient
            elif method == 'IModPoly':
                if repetition is None:
                    repetition = 100
                self.preprocessed_intensity_array = BaselineRemoval(copy.deepcopy(
                    self.preprocessed_intensity_array)).IModPoly(
                    degree=degree,
                    repitition=repetition,
                    gradient=gradient)
                self.data_processing['baseline removal']['degree'] = degree
                self.data_processing['baseline removal']['repetition'] = repetition
                self.data_processing['baseline removal']['gradient'] = gradient
            gc.collect()

    def normalize_intensity(self, method='tic'):
        # check method
        if method not in ['tic', 'rms', 'mad', 'sqrt']:
            raise Exception('Method must be "tic", "rms", "mad", or "sqrt"')
        else:
            if method == 'tic':
                self.preprocessed_intensity_array = tic(self.preprocessed_mz_array,
                                                        copy.deepcopy(self.preprocessed_intensity_array))
            elif method == 'rms':
                self.preprocessed_intensity_array = rms(self.preprocessed_mz_array,
                                                        copy.deepcopy(self.preprocessed_intensity_array))
            elif method == 'mad':
                self.preprocessed_intensity_array = mad(self.preprocessed_mz_array,
                                                        copy.deepcopy(self.preprocessed_intensity_array))
            elif method == 'sqrt':
                self.preprocessed_intensity_array = sqrt(self.preprocessed_mz_array,
                                                         copy.deepcopy(self.preprocessed_intensity_array))
            self.data_processing['intensity normalization'] = {'method': method}
            gc.collect()

    def bin_spectrum(self, n_bins, lower_mass_range, upper_mass_range):
        bins = np.linspace(lower_mass_range, upper_mass_range, n_bins, dtype=np.float64)
        unique_indices, inverse_indices = np.unique(np.digitize(copy.deepcopy(self.preprocessed_mz_array),
                                                                bins),
                                                    return_inverse=True)
        bin_counts = np.bincount(inverse_indices)
        np.place(bin_counts, bin_counts < 1, [1])
        self.preprocessed_mz_array = np.bincount(inverse_indices,
                                                 weights=copy.deepcopy(self.preprocessed_mz_array)) / bin_counts
        self.preprocessed_intensity_array = np.bincount(inverse_indices,
                                                        weights=copy.deepcopy(self.preprocessed_intensity_array))

        self.data_processing['spectra binning'] = {'lower mass range': lower_mass_range,
                                                   'upper mass range': upper_mass_range,
                                                   'number of bins': n_bins}
        gc.collect()

    def estimate_peak_widths(self):
        peak_indices, peak_properties = find_peaks(copy.deepcopy(self.preprocessed_intensity_array))
        widths = peak_widths(copy.deepcopy(self.preprocessed_intensity_array), peak_indices)
        return widths[0]

    def peak_picking(self, method='cwt', widths=None, snr=3):
        # check method
        if method not in ['locmax', 'cwt']:
            raise Exception('Method must be "locmax" or "cwt"')

        self.data_processing['peak picking'] = {'method': method}
        if method == 'locmax':
            #print(np.mean(self.preprocessed_intensity_array))
            #print(np.median(self.preprocessed_intensity_array))
            #print(median_abs_deviation(self.preprocessed_intensity_array))
            peak_indices, peak_properties = find_peaks(copy.deepcopy(self.preprocessed_intensity_array),
                                                       height=np.mean(self.preprocessed_intensity_array) * snr)
            self.peak_picking_indices = copy.deepcopy(peak_indices)
            self.peak_picked_mz_array = copy.deepcopy(self.preprocessed_mz_array)[peak_indices]
            self.peak_picked_intensity_array = copy.deepcopy(self.preprocessed_intensity_array)[peak_indices]
        elif method == 'cwt':
            # estimate peak widths if necessary
            if widths is None:
                estimated_widths = self.estimate_peak_widths()
                widths_start = np.min(estimated_widths)
                widths_stop = 2 * np.mean(estimated_widths)
                widths_step = ((2 * np.mean(estimated_widths)) - np.min(estimated_widths)) / 10
                widths = np.arange(widths_start, widths_stop, widths_step)
            peak_indices = find_peaks_cwt(copy.deepcopy(self.preprocessed_intensity_array), widths, min_snr=snr)
            self.peak_picking_indices = copy.deepcopy(peak_indices)
            self.peak_picked_mz_array = copy.deepcopy(self.preprocessed_mz_array)[peak_indices]
            self.peak_picked_intensity_array = copy.deepcopy(self.preprocessed_intensity_array)[peak_indices]
            self.data_processing['peak picking']['lower peak width'] = np.min(widths)
            self.data_processing['peak picking']['upper peak width'] = np.max(widths)

        gc.collect()

    def undo_all_processing(self):
        self.preprocessed_mz_array = None
        self.preprocessed_intensity_array = None
        self.peak_picked_mz_array = None
        self.peak_picked_intensity_array = None
        self.peak_picking_indices = None
        self.data_processing = None
        gc.collect()
        self.preprocessed_mz_array = copy.deepcopy(self.mz_array)
        self.preprocessed_intensity_array = copy.deepcopy(self.intensity_array)
        self.data_processing = {}
        gc.collect()

    def plot_spectrum(self):
        spectrum_df = pd.DataFrame({'m/z': copy.deepcopy(self.preprocessed_mz_array),
                                    'Intensity': copy.deepcopy(self.preprocessed_intensity_array)})
        fig = sns.lineplot(data=spectrum_df, x='m/z', y='Intensity')
        plt.show()


# Barebones class to hold basic MALDI spectrum data.
class OpenMALDISpectrum(PMPMethods):
    def __init__(self, pyteomics_dict, filename):
        self.source = filename
        self.name = None
        self.uuid = str(uuid4())
        self.spectrum_id = None
        self.coord = None
        self.ms_level = None
        self.centroided = None
        self.mz_array = None
        self.intensity_array = None
        self.preprocessed_mz_array = None
        self.preprocessed_intensity_array = None
        self.peak_picked_mz_array = None
        self.peak_picked_intensity_array = None
        self.peak_picking_indices = None
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
        self.mz_array = pyteomics_dict['m/z array']
        self.intensity_array = pyteomics_dict['intensity array']
        self.preprocessed_mz_array = pyteomics_dict['m/z array']
        self.preprocessed_intensity_array = pyteomics_dict['intensity array']


class PMPTsfSpectrum(TsfSpectrum, PMPMethods):
    def __init__(self, tsf_data, frame: int, mode: str, profile_bins=0, encoding=64):
        super().__init__(tsf_data, frame, mode, profile_bins, encoding)

        self.source = tsf_data.source_file
        self.name = None
        self.uuid = str(uuid4())
        self.spectrum_id = None
        self.preprocessed_mz_array = copy.deepcopy(self.mz_array)
        self.preprocessed_intensity_array = copy.deepcopy(self.intensity_array)
        self.peak_picked_mz_array = None
        self.peak_picked_intensity_array = None
        self.peak_picking_indices = None
        self.data_processing = {}

        self.get_spectrum_metadata()

    def get_spectrum_metadata(self):
        self.name = str(os.path.splitext(os.path.split(self.source)[-1])[0]) + '_' + str(self.frame)
        self.spectrum_id = self.name + '|' + self.coord + '|' + self.uuid


class PMPTdfSpectrum(TdfSpectrum, PMPMethods):
    def __init__(self, tdf_data, frame: int, mode: str, precursor=0, diapasef_window=None, profile_bins=0, encoding=64,
                 exclude_mobility=False):
        super().__init__(tdf_data, frame, mode, precursor, diapasef_window, profile_bins, encoding, exclude_mobility)

        self.source = tdf_data.source_file
        self.name = None
        self.uuid = str(uuid4())
        self.spectrum_id = None
        self.preprocessed_mz_array = copy.deepcopy(self.mz_array)
        self.preprocessed_intensity_array = copy.deepcopy(self.intensity_array)
        self.peak_picked_mz_array = None
        self.peak_picked_intensity_array = None
        self.peak_picking_indices = None
        self.data_processing = {}

        self.get_spectrum_metadata()

    def get_spectrum_metadata(self):
        self.name = str(os.path.splitext(os.path.split(self.source)[-1])[0]) + '_' + str(self.frame)
        self.spectrum_id = self.name + '|' + self.coord + '|' + self.uuid
