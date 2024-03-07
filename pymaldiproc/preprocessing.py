import numpy as np
import pandas as pd
from functools import reduce
from scipy.signal import savgol_filter, find_peaks, find_peaks_cwt, peak_widths
from pyMSpec.smoothing import sg_smooth, apodization, rebin, fast_change, median
from pybaselines.smooth import snip, noise_median
from pybaselines.morphological import tophat
from BaselineRemoval import BaselineRemoval
from pyMSpec.normalisation import tic, rms, mad, sqrt
from icoshift import icoshift


def trim_spectra(list_of_spectra, lower_mass_range, upper_mass_range):
    print('Trimming spectra')
    for spectrum in list_of_spectra:
        spectrum_df = pd.DataFrame(data={'mz': spectrum.get_mz_array(), 'intensity': spectrum.get_intensity_array()})
        spectrum_df = spectrum_df[~(spectrum_df['mz'] <= lower_mass_range) & ~(spectrum_df['mz'] >= upper_mass_range)]
        spectrum.preprocessed_mz_array = spectrum_df['mz'].values
        spectrum.preprocessed_intensity_array = spectrum_df['intensity'].values
        spectrum.data_processing['spectrum trimming'] = {'lower mass range': lower_mass_range,
                                                         'upper mass range': upper_mass_range}
    return list_of_spectra


def transform_intensity(list_of_spectra, method='sqrt'):
    print('Performing intensity transformation')
    # check method
    if method not in ['sqrt', 'log', 'log2', 'log10']:
        raise Exception('Method must be "sqrt", "log", "log2", or "log10"')

    # transform intensity
    for spectrum in list_of_spectra:
        if method == 'sqrt':
            spectrum.preprocessed_intensity_array = np.sqrt(spectrum.get_intensity_array())
        elif method == 'log':   #ln
            spectrum.preprocessed_intensity_array = np.log(spectrum.get_intensity_array())
        elif method == 'log2':
            spectrum.preprocessed_intensity_array = np.log2(spectrum.get_intensity_array())
        elif method == 'log10':
            spectrum.preprocessed_intensity_array = np.log10(spectrum.get_intensity_array())
        spectrum.data_processing['intensity transformation'] = {'method': method}
    return list_of_spectra


def smooth_baseline(list_of_spectra, method='SavitzkyGolay', window_length=21, polyorder=3, delta_mz=0.2,
                    diff_thresh=0.01):
    print('Smoothing baseline')
    # check method
    if method not in ['SavitzkyGolay', 'apodization', 'rebin', 'fast_change', 'median']:
        raise Exception('Method must be "SavitzkyGolay", "apodization", "rebin", "fast_change", "median"')

    # smooth baseline
    for spectrum in list_of_spectra:
        spectrum.data_processing['baseline smoothing'] = {'method': method}
        if method == 'SavitzkyGolay':
            spectrum.preprocessed_intensity_array = savgol_filter(spectrum.get_intensity_array(),
                                                                  window_length=window_length,
                                                                  polyorder=polyorder)
            spectrum.data_processing['baseline smoothing']['window length'] = window_length
            spectrum.data_processing['baseline smoothing']['polyorder'] = polyorder
        elif method == 'apodization':
            spectrum.preprocessed_mz_array, spectrum.preprocessed_intensity_array = apodization(spectrum.get_mz_array(),
                                                                                                spectrum.get_intensity_array(),
                                                                                                w_size=window_length)
            spectrum.data_processing['baseline smoothing']['window length'] = window_length
        elif method == 'rebin':
            spectrum.preprocessed_mz_array, spectrum.preprocessed_intensity_array = rebin(spectrum.get_mz_array(),
                                                                                          spectrum.get_intensity_array(),
                                                                                          delta_mz=delta_mz)
            spectrum.data_processing['baseline smoothing']['delta m/z'] = delta_mz
        elif method == 'fast_change':
            spectrum.preprocessed_mz_array, spectrum.preprocessed_intensity_array = fast_change(spectrum.get_mz_array(),
                                                                                                spectrum.get_intensity_array(),
                                                                                                diff_thresh=diff_thresh)
            spectrum.data_processing['baseline smoothing']['difference threshold'] = diff_thresh
        elif method == 'median':
            spectrum.preprocessed_mz_array, spectrum.preprocessed_intensity_array = median(spectrum.get_mz_array(),
                                                                                           spectrum.get_intensity_array(),
                                                                                           w_size=window_length)
            spectrum.data_processing['baseline smoothing']['window length'] = window_length

    return list_of_spectra


def remove_baseline(list_of_spectra, method='SNIP', min_half_window=1, max_half_window=100, decreasing=True,
                    smooth_half_window=None, filter_order=2, sigma=None, increment=1, max_hits=1, window_tol=0.000001,
                    lambda_=100, porder=1, repetition=None, degree=2, gradient=0.001):
    print('Removing baseline')
    # check method
    if method not in ['SNIP', 'TopHat', 'Median','ZhangFit', 'ModPoly', 'IModPoly']:
        raise Exception('Method must be "SNIP", "TopHat", "Median", "ZhangFit", "ModPoly", or "IModPoly"')

    # remove baseline
    for spectrum in list_of_spectra:
        spectrum.data_processing['baseline removal'] = {'method': method}
        if method == 'SNIP':
            baseline = snip(data=spectrum.get_intensity_array(),
                            max_half_window=max_half_window,
                            decreasing=decreasing,
                            smooth_half_window=smooth_half_window,
                            filter_order=filter_order)[0]
            spectrum.preprocessed_intensity_array = spectrum.get_intensity_array() - baseline
            spectrum.data_processing['baseline removal']['max half window'] = max_half_window
            spectrum.data_processing['baseline removal']['decreasing'] = decreasing
            spectrum.data_processing['baseline removal']['smooth half window'] = smooth_half_window
            spectrum.data_processing['baseline removal']['filter order'] = filter_order
        elif method == 'TopHat':
            baseline = tophat(data=spectrum.get_intensity_array(),
                              half_window=max_half_window,
                              increment=increment,
                              max_hits=max_hits,
                              window_tol=window_tol,
                              max_half_window=max_half_window,
                              min_half_window=min_half_window)[0]
            spectrum.preprocessed_intensity_array = spectrum.get_intensity_array() - baseline
            spectrum.data_processing['baseline removal']['half window'] = max_half_window
            spectrum.data_processing['baseline removal']['increment'] = increment
            spectrum.data_processing['baseline removal']['max hits'] = max_hits
            spectrum.data_processing['baseline removal']['window tolerance'] = window_tol
            spectrum.data_processing['baseline removal']['max half window'] = max_half_window
            spectrum.data_processing['baseline removal']['min half window'] = min_half_window
        elif method == 'Median':
            baseline = noise_median(data=spectrum.get_intensity_array(),
                                    half_window=max_half_window,
                                    smooth_half_window=smooth_half_window,
                                    sigma=sigma)[0]
            spectrum.preprocessed_intensity_array = spectrum.get_intensity_array() - baseline
            spectrum.data_processing['baseline removal']['half window'] = max_half_window
            spectrum.data_processing['baseline removal']['smooth half window'] = smooth_half_window
            spectrum.data_processing['baseline removal']['sigma'] = sigma
        elif method == 'ZhangFit':
            if repetition is None:
                repetition = 15
            spectrum.preprocessed_intensity_array = BaselineRemoval(spectrum.get_intensity_array()).ZhangFit(lambda_=lambda_,
                                                                                                             porder=porder,
                                                                                                             repitition=repetition)
            spectrum.data_processing['baseline removal']['lambda'] = lambda_
            spectrum.data_processing['baseline removal']['porder'] = porder
            spectrum.data_processing['baseline removal']['repetition'] = repetition
        elif method == 'ModPoly':
            if repetition is None:
                repetition = 100
            spectrum.preprocessed_intensity_array = BaselineRemoval(spectrum.get_intensity_array()).ModPoly(degree=degree,
                                                                                                            repitition=repetition,
                                                                                                            gradient=gradient)
            spectrum.data_processing['baseline removal']['degree'] = degree
            spectrum.data_processing['baseline removal']['repetition'] = repetition
            spectrum.data_processing['baseline removal']['gradient'] = gradient
        elif method == 'IModPoly':
            if repetition is None:
                repetition = 100
            spectrum.preprocessed_intensity_array = BaselineRemoval(spectrum.get_intensity_array()).IModPoly(degree=degree,
                                                                                                             repitition=repetition,
                                                                                                             gradient=gradient)
            spectrum.data_processing['baseline removal']['degree'] = degree
            spectrum.data_processing['baseline removal']['repetition'] = repetition
            spectrum.data_processing['baseline removal']['gradient'] = gradient

    return list_of_spectra


def normalize_intensity(list_of_spectra, method='tic'):
    print('Normalizing intensity')
    # check method
    if method not in ['tic', 'rms', 'mad', 'sqrt']:
        raise Exception('Method must be "tic", "rms", "mad", or "sqrt"')

    # normalize intensity
    for spectrum in list_of_spectra:
        if method == 'tic':
            spectrum.preprocessed_intensity_array = tic(spectrum.get_mz_array(),
                                                        spectrum.get_intensity_array())
        elif method == 'rms':
            spectrum.preprocessed_intensity_array = rms(spectrum.get_mz_array(),
                                                        spectrum.get_intensity_array())
        elif method == 'mad':
            spectrum.preprocessed_intensity_array = mad(spectrum.get_mz_array(),
                                                        spectrum.get_intensity_array())
        elif method == 'sqrt':
            spectrum.preprocessed_intensity_array = sqrt(spectrum.get_mz_array(),
                                                         spectrum.get_intensity_array())
        spectrum.data_processing['intensity normalization'] = {'method': method}

    return list_of_spectra


def bin_spectra(list_of_spectra, n_bins, lower_mass_range, upper_mass_range, array_type='preprocessed'):
    print('Binning spectra')
    bins = np.linspace(lower_mass_range, upper_mass_range, n_bins, dtype=np.float64)

    for spectrum in list_of_spectra:
        unique_indices, inverse_indices = np.unique(np.digitize(spectrum.get_mz_array(), bins), return_inverse=True)
        bin_counts = np.bincount(inverse_indices)
        np.place(bin_counts, bin_counts < 1, [1])
        if array_type == 'raw':
            spectrum.raw_mz_array = np.bincount(inverse_indices,
                                                weights=spectrum.raw_mz_array) / bin_counts
            spectrum.raw_intensity_array = np.bincount(inverse_indices,
                                                       weights=spectrum.raw_intensity_array)
        elif array_type == 'preprocessed':
            spectrum.preprocessed_mz_array = np.bincount(inverse_indices,
                                                         weights=spectrum.preprocessed_mz_array) / bin_counts
            spectrum.preprocessed_intensity_array = np.bincount(inverse_indices,
                                                                weights=spectrum.preprocessed_intensity_array)
        elif array_type == 'peak_picked':
            spectrum.peak_picked_mz_array = np.bincount(inverse_indices,
                                                        weights=spectrum.peak_picked_mz_array) / bin_counts
            spectrum.peak_picked_intensity_array = np.bincount(inverse_indices,
                                                               weights=spectrum.peak_picked_intensity_array) / bin_counts

        spectrum.data_processing['spectra binning'] = {'lower mass range': lower_mass_range,
                                                       'upper mass range': upper_mass_range,
                                                       'number of bins': n_bins}

    return list_of_spectra


def align_spectra(list_of_spectra, method='average', inter='whole', n='f', scale=None, coshift_preprocessing=False,
                  coshift_preprocessing_max_shift=None, fill_with_previous=True, average2_multiplier=3):
    print('Aligning spectra')
    # check method
    if method not in ['average', 'median', 'max', 'average2']:
        raise Exception('Method must be "average", "median", "max", or "average2"')

    # intensity arrays must be binned/on same m/z axis prior to alignment
    intensity_arrays = [spectrum.get_intensity_array() for spectrum in list_of_spectra]
    intensity_array_2d = np.stack(intensity_arrays)
    xcs, ints, ind, target = icoshift(method, intensity_array_2d, inter=inter, n=n, scale=scale,
                                      coshift_preprocessing=coshift_preprocessing,
                                      coshift_preprocessing_max_shift=coshift_preprocessing_max_shift,
                                      fill_with_previous=fill_with_previous,
                                      average2_multiplier=average2_multiplier)

    for spectrum, row in zip(list_of_spectra, xcs):
        spectrum.preprocessed_intensity_array = row
        spectrum.data_processing['spectra alignment'] = {'method': method}
        spectrum.data_processing['spectra alignment']['inter'] = inter
        spectrum.data_processing['spectra alignment']['n'] = n
        spectrum.data_processing['spectra alignment']['scale'] = scale
        spectrum.data_processing['spectra alignment']['coshift_preprocessing'] = coshift_preprocessing
        spectrum.data_processing['spectra alignment']['coshift_preprocessing_max_shift'] = coshift_preprocessing_max_shift
        spectrum.data_processing['spectra alignment']['fill_with_previous'] = fill_with_previous
        spectrum.data_processing['spectra alignment']['average2_multiplier'] = average2_multiplier

    return list_of_spectra


def estimate_peak_widths(intensity_array):
    print('Estimating peak widths')
    peak_indices, peak_properties = find_peaks(intensity_array)
    widths = peak_widths(intensity_array, peak_indices)
    return widths[0]


def peak_picking(list_of_spectra, method='cwt', widths=None, snr=3):
    print('Performing peak picking')
    # check method
    if method not in ['locmax', 'cwt']:
        raise Exception('Method must be "locmax" or "cwt"')

    # peak picking
    for spectrum in list_of_spectra:
        spectrum.data_processing['peak picking'] = {'method': method}
        if method == 'locmax':
            peak_indices, peak_properties = find_peaks(spectrum.get_intensity_array())
            spectrum.peak_picked_mz_array = spectrum.get_mz_array()[peak_indices]
            spectrum.peak_picked_intensity_array = spectrum.get_intensity_array()[peak_indices]
        elif method == 'cwt':
            # estimate peak widths if necessary
            if widths is None:
                estimated_widths = estimate_peak_widths(spectrum.get_intensity_array())
                widths_start = np.min(estimated_widths)
                widths_stop = 2 * np.mean(estimated_widths)
                widths_step = ((2 * np.mean(estimated_widths)) - np.min(estimated_widths)) / 10
                widths = np.arange(widths_start, widths_stop, widths_step)
            peak_indices = find_peaks_cwt(spectrum.get_intensity_array(), widths, min_snr=snr)
            spectrum.peak_picked_mz_array = spectrum.get_mz_array()[peak_indices]
            spectrum.peak_picked_intensity_array = spectrum.get_intensity_array()[peak_indices]
            spectrum.data_processing['peak picking']['lower peak width'] = np.min(widths)
            spectrum.data_processing['peak picking']['upper peak width'] = np.max(widths)

    return list_of_spectra


def get_feature_matrix(list_of_spectra, tolerance=0.05, decimals=4, missing_value_imputation=True):
    print('Creating feature intensity matrix')
    # get a consensus m/z array
    peak_picked_mz_arrays = [spectrum.peak_picked_mz_array for spectrum in list_of_spectra]
    peak_picked_consensus = pd.DataFrame(data={'mz': np.unique(np.concatenate(peak_picked_mz_arrays))}).sort_values(by='mz')
    preprocessed_mz_arrays = [spectrum.get_mz_array() for spectrum in list_of_spectra]
    preprocessed_concensus = pd.DataFrame(data={'mz': np.unique(np.concatenate(preprocessed_mz_arrays))}).sort_values(by='mz')

    spectra_dfs_peak_picked = [peak_picked_consensus]
    spectra_dfs_preprocessed = [preprocessed_concensus]
    for spectrum in list_of_spectra:
        spectra_dfs_peak_picked.append(pd.DataFrame(data={'mz': spectrum.peak_picked_mz_array,
                                                          spectrum.spectrum_id: spectrum.get_intensity_array()}))
        spectra_dfs_preprocessed.append(pd.DataFrame(data={'mz': spectrum.get_mz_array(),
                                                           spectrum.spectrum_id: spectrum.get_intensity_array()}))
    feature_matrix = reduce(lambda x, y: pd.merge_asof(x, y, on='mz', tolerance=tolerance, direction='nearest'), spectra_dfs_peak_picked).sort_values(by='mz')
    if missing_value_imputation:
        ref_matrix = reduce(lambda x, y: pd.merge_asof(x, y, on='mz', tolerance=tolerance, direction='nearest'), spectra_dfs_preprocessed).sort_values(by='mz')
        for colname in feature_matrix.columns:
            if colname != 'mz':
                tmp_df = pd.merge_asof(feature_matrix[['mz', colname]],
                                       ref_matrix[['mz', colname]],
                                       on='mz',
                                       tolerance=tolerance,
                                       direction='nearest')
                feature_matrix[colname] = tmp_df.drop('mz', axis=1).mean(axis=1).values
    feature_matrix = feature_matrix.fillna(0)
    feature_matrix = feature_matrix.round({'mz': decimals})
    feature_matrix = feature_matrix.groupby(['mz'], as_index=False).aggregate(sum)
    return feature_matrix


def export_feature_list(feature_matrix, output):
    print('Exporting feature matrix to ' + output)
    feature_matrix.to_csv(output, index=False)
