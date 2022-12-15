import numpy as np
import pandas as pd
from functools import reduce
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


# TODO: add peak alignment functionality


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


def bin_peaks(list_of_spectra, lower_mass_range, upper_mass_range, num_bins=0, bin_width=0):
    if num_bins == 0 and bin_width == 0:
        raise Exception('Number of bins or bin width to use for binning must be specified')
    elif num_bins != 0 and bin_width != 0:
        raise Exception('Only specify number of bins or bin width. Do not specify both')
    elif num_bins < 0 or bin_width < 0:
        raise Exception('Number of bins/bin width must be > 0')
    elif num_bins > 0:
        for spectrum in list_of_spectra:
            bins = np.linspace(lower_mass_range, upper_mass_range, num_bins, dtype=np.float64)
            unique_indices, inverse_indices = np.unique(np.digitize(spectrum.peak_picked_mz_array, bins), return_inverse=True)
            bin_counts = np.bincount(inverse_indices)
            np.place(bin_counts, bin_counts < 1, [1])
            spectrum.peak_picked_mz_array = np.bincount(inverse_indices, weights=spectrum.peak_picked_mz_array) / bin_counts
            spectrum.peak_picked_intensity_array = np.bincount(inverse_indices, weights=spectrum.peak_picked_intensity_array)
    elif bin_width > 0:
        for spectrum in list_of_spectra:
            bins = np.arange(lower_mass_range, upper_mass_range, bin_width, dtype=np.float64)
            unique_indices, inverse_indices = np.unique(np.digitize(spectrum.peak_picked_mz_array, bins), return_inverse=True)
            bin_counts = np.bincount(inverse_indices)
            np.place(bin_counts, bin_counts < 1, [1])
            spectrum.peak_picked_mz_array = np.bincount(inverse_indices, weights=spectrum.peak_picked_mz_array) / bin_counts
            spectrum.peak_picked_intensity_array = np.bincount(inverse_indices, weights=spectrum.peak_picked_intensity_array)

    return list_of_spectra


def get_feature_matrix(list_of_spectra, missing_value_imputation=True):
    # get a consensus m/z array
    peak_picked_mz_arrays = [spectrum.peak_picked_mz_array for spectrum in list_of_spectra]
    peak_picked_consensus = pd.DataFrame(data={'mz': np.unique(np.concatenate(peak_picked_mz_arrays))}).sort_values(by='mz')
    preprocessed_mz_arrays = [spectrum.preprocessed_mz_array for spectrum in list_of_spectra]
    preprocessed_concensus = pd.DataFrame(data={'mz': np.unique(np.concatenate(preprocessed_mz_arrays))}).sort_values(by='mz')

    spectra_dfs_peak_picked = [peak_picked_consensus]
    spectra_dfs_preprocessed = [preprocessed_concensus]
    for spectrum in list_of_spectra:
        spectra_dfs_peak_picked.append(pd.DataFrame(data={'mz': spectrum.peak_picked_mz_array,
                                                          spectrum.spectrum_id: spectrum.peak_picked_intensity_array}))
        spectra_dfs_preprocessed.append(pd.DataFrame(data={'mz': spectrum.preprocessed_mz_array,
                                                           spectrum.spectrum_id: spectrum.preprocessed_intensity_array}))
    #feature_matrix = reduce(lambda x, y: pd.merge(x, y, how='outer', on='mz'), spectra_dfs_peak_picked).sort_values(by='mz')
    feature_matrix = reduce(lambda x, y: pd.merge_asof(x, y, on='mz', tolerance=0.5, direction='nearest'), spectra_dfs_peak_picked).sort_values(by='mz')
    if missing_value_imputation:
        #ref_matrix = reduce(lambda x, y: pd.merge(x, y, how='outer', on='mz'), spectra_dfs_preprocessed).sort_values(by='mz')
        ref_matrix = reduce(lambda x, y: pd.merge_asof(x, y, on='mz', tolerance=0.5, direction='nearest'), spectra_dfs_preprocessed).sort_values(by='mz')
        for colname in feature_matrix.columns:
            if colname != 'mz':
                #tmp_df = pd.merge(feature_matrix[['mz', colname]],
                #                  ref_matrix[['mz', colname]],
                #                  how='left',
                #                  on='mz')
                tmp_df = pd.merge_asof(feature_matrix[['mz', colname]],
                                       ref_matrix[['mz', colname]],
                                       on='mz',
                                       tolerance=0.5,
                                       direction='nearest')
                #feature_matrix[colname].fillna(tmp_df.drop('mz', axis=1).mean(axis=1), inplace=True)
                feature_matrix[colname] = tmp_df.drop('mz', axis=1).mean(axis=1).values
    feature_matrix = feature_matrix.fillna(0)
    return feature_matrix


def export_feature_list(feature_matrix, output):
    feature_matrix.to_csv(output, index=False)


def get_cos_distance_matrix():
    pass
