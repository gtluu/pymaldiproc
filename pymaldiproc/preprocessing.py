import numpy as np
import pandas as pd
from scipy.signal import savgol_filter, find_peaks, find_peaks_cwt, peak_widths
from pyMSpec.smoothing import sg_smooth, apodization, rebin, fast_change, median
from BaselineRemoval import BaselineRemoval
from pyMSpec.normalisation import tic, rms, mad, sqrt


def trim_spectra(list_of_spectra, lower_mass_range, upper_mass_range):
    #trimmed_spectra = []
    for spectrum in list_of_spectra:
        spectrum_df = pd.DataFrame(data={'mz': spectrum.raw_mz_array, 'intensity': spectrum.raw_intensity_array})
        spectrum_df = spectrum_df[~(spectrum_df['mz'] <= lower_mass_range) & ~(spectrum_df['mz'] >= upper_mass_range)]
        spectrum.preprocessed_mz_array = spectrum_df['mz'].values
        spectrum.preprocessed_intensity_array = spectrum_df['intensity'].values
        spectrum.data_processing['spectrum trimming'] = {'lower mass range': lower_mass_range,
                                                         'upper mass range': upper_mass_range}
        #trimmed_spectra.append(spectrum)
    return list_of_spectra


def transform_intensity(list_of_spectra, method='sqrt'):
    # check method
    if method not in ['sqrt', 'log', 'log2', 'log10']:
        raise Exception('Method must be "sqrt", "log", "log2", or "log10"')

    # transform intensity
    for spectrum in list_of_spectra:
        if method == 'sqrt':
            spectrum.preprocessed_intensity_array = np.sqrt(spectrum.preprocessed_intensity_array)
        elif method == 'log':
            spectrum.preprocessed_intensity_array = np.log(spectrum.preprocessed_intensity_array)
        elif method == 'log2':
            spectrum.preprocessed_intensity_array = np.log2(spectrum.preprocessed_intensity_array)
        elif method == 'log10':
            spectrum.preprocessed_intensity_array = np.log10(spectrum.preprocessed_intensity_array)
        spectrum.data_processing['intensity transformation'] = {'method': method}

    return list_of_spectra


# TODO: implement MovingAverage method
def smooth_baseline(list_of_spectra, method='SavitzkyGolay', window_length=20, polyorder=3, delta_mz=0.2,
                    diff_thresh=0.01):
    # check method
    #if method not in ['SavitzkyGolay', 'MovingAverage', 'apodization', 'rebin', 'fast_change', 'median']:
    #    raise Exception('Method must be "SavitzkyGolay", "MovingAverage", "apodization", "rebin", "fast_change", '
    #                    '"median"')

    if method not in ['SavitzkyGolay', 'apodization', 'rebin', 'fast_change', 'median']:
        raise Exception('Method must be "SavitzkyGolay", "apodization", "rebin", "fast_change", "median"')

    # smooth baseline
    for spectrum in list_of_spectra:
        spectrum.data_processing['baseline smoothing'] = {'method': method}
        if method == 'SavitzkyGolay':
            spectrum.preprocessed_intensity_array = savgol_filter(spectrum.preprocessed_intensity_array,
                                                                  window_length=window_length,
                                                                  polyorder=polyorder)
            spectrum.data_processing['baseline smoothing']['window length'] = window_length
            spectrum.data_processing['baseline smoothing']['polyorder'] = polyorder
        elif method == 'MovingAverage':
            pass
        elif method == 'apodization':
            spectrum.preprocessed_mz_array, spectrum.preprocessed_intensity_array = apodization(spectrum.preprocessed_mz_array,
                                                                                                spectrum.preprocessed_intensity_array,
                                                                                                w_size=window_length)
            spectrum.data_processing['baseline smoothing']['window length'] = window_length
        elif method == 'rebin':
            spectrum.preprocessed_mz_array, spectrum.preprocessed_intensity_array = rebin(spectrum.preprocessed_mz_array,
                                                                                          spectrum.preprocessed_intensity_array,
                                                                                          delta_mz=delta_mz)
            spectrum.data_processing['baseline smoothing']['delta m/z'] = delta_mz
        elif method == 'fast_change':
            spectrum.preprocessed_mz_array, spectrum.preprocessed_intensity_array = fast_change(spectrum.preprocessed_mz_array,
                                                                                                spectrum.preprocessed_intensity_array,
                                                                                                diff_thresh=diff_thresh)
            spectrum.data_processing['baseline smoothing']['difference threshold'] = diff_thresh
        elif method == 'median':
            spectrum.preprocessed_mz_array, spectrum.preprocessed_intensity_array = median(spectrum.preprocessed_mz_array,
                                                                                           spectrum.preprocessed_intensity_array,
                                                                                           w_size=window_length)
            spectrum.data_processing['baseline smoothing']['window length'] = window_length

    return list_of_spectra


# TODO: implement other algorithms
def remove_baseline(list_of_spectra, method='ZhangFit'):
    # check method
    #if method not in ['snip', 'tophat', 'convexhull', 'median', 'ZhangFit', 'modpoly', 'imodpoly']:
    #    raise Exception('Method must be "snip", "tophat", "convexhull", "median", "ZhangFit", "modpoly", or "imodpoly"')

    # remove baseline
    for spectrum in list_of_spectra:
        if method == 'ZhangFit':
            spectrum.preprocessed_intensity_array = BaselineRemoval(spectrum.preprocessed_intensity_array).ZhangFit()
        spectrum.data_processing['baseline removal'] = {'method': method}

    return list_of_spectra


# TODO: implement PQN and median from MALDIquant
def normalize_intensity(list_of_spectra, method='tic'):
    # check method
    if method not in ['tic', 'rms', 'mad', 'sqrt']:
        raise Exception('Method must be "tic", "rms", "mad", or "sqrt"')

    # normalize intensity
    for spectrum in list_of_spectra:
        if method == 'tic':
            spectrum.preprocessed_intensity_array = tic(spectrum.preprocessed_mz_array,
                                                        spectrum.preprocessed_intensity_array)
        elif method == 'rms':
            spectrum.preprocessed_intensity_array = rms(spectrum.preprocessed_mz_array,
                                                        spectrum.preprocessed_intensity_array)
        elif method == 'mad':
            spectrum.preprocessed_intensity_array = mad(spectrum.preprocessed_mz_array,
                                                        spectrum.preprocessed_intensity_array)
        elif method == 'sqrt':
            spectrum.preprocessed_intensity_array = sqrt(spectrum.preprocessed_mz_array,
                                                         spectrum.preprocessed_intensity_array)
        spectrum.data_processing['intensity normalization'] = {'method': method}

    return list_of_spectra


def estimate_peak_widths(intensity_array):
    peak_indices, peak_properties = find_peaks(intensity_array)
    widths = peak_widths(intensity_array, peak_indices)
    return widths[0]


def peak_picking(list_of_spectra, method='cwt', widths=None, snr=3):
    # check method
    if method not in ['locmax', 'cwt']:
        raise Exception('Method must be "locmax" or "cwt"')

    # peak picking
    for spectrum in list_of_spectra:
        spectrum.data_processing['peak picking'] = {'method': method}
        if method == 'locmax':
            peak_indices, peak_properties = find_peaks(spectrum.preprocessed_intensity_array)
            spectrum.peak_picked_mz_array = spectrum.preprocessed_mz_array[peak_indices]
            spectrum.peak_picked_intensity_array = spectrum.preprocessed_intensity_array[peak_indices]
        elif method == 'cwt':
            # estimate peak widths if necessary
            if widths is None:
                estimated_widths = estimate_peak_widths(spectrum.preprocessed_intensity_array)
                widths_start = np.min(estimated_widths)
                widths_stop = 2 * np.mean(estimated_widths)
                widths_step = ((2 * np.mean(estimated_widths)) - np.min(estimated_widths)) / 10
                widths = np.arange(widths_start, widths_stop, widths_step)
            peak_indices = find_peaks_cwt(spectrum.preprocessed_intensity_array, widths, min_snr=snr)
            spectrum.peak_picked_mz_array = spectrum.preprocessed_mz_array[peak_indices]
            spectrum.peak_picked_intensity_array = spectrum.preprocessed_intensity_array[peak_indices]
            spectrum.data_processing['peak picking']['lower peak width'] = np.min(widths)
            spectrum.data_processing['peak picking']['upper peak width'] = np.max(widths)

    return list_of_spectra
