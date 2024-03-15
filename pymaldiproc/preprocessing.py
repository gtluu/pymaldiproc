import copy
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


def align_spectra(list_of_spectra, method='average', inter='whole', n='f', scale=None, coshift_preprocessing=False,
                  coshift_preprocessing_max_shift=None, fill_with_previous=True, average2_multiplier=3):
    print('Aligning spectra')
    # check method
    if method not in ['average', 'median', 'max', 'average2']:
        raise Exception('Method must be "average", "median", "max", or "average2"')

    # intensity arrays must be binned/on same m/z axis prior to alignment
    intensity_arrays = [copy.deepcopy(spectrum.preprocessed_intensity_array) for spectrum in list_of_spectra]
    intensity_array_2d = np.stack(intensity_arrays)
    xcs, ints, ind, target = icoshift(method, intensity_array_2d, inter=inter, n=n, scale=scale,
                                      coshift_preprocessing=coshift_preprocessing,
                                      coshift_preprocessing_max_shift=coshift_preprocessing_max_shift,
                                      fill_with_previous=fill_with_previous,
                                      average2_multiplier=average2_multiplier)

    for spectrum, row in zip(list_of_spectra, xcs):
        spectrum.preprocessed_intensity_array = copy.deepcopy(row)
        spectrum.data_processing['spectra alignment'] = {'method': method}
        spectrum.data_processing['spectra alignment']['inter'] = inter
        spectrum.data_processing['spectra alignment']['n'] = n
        spectrum.data_processing['spectra alignment']['scale'] = scale
        spectrum.data_processing['spectra alignment']['coshift_preprocessing'] = coshift_preprocessing
        spectrum.data_processing['spectra alignment']['coshift_preprocessing_max_shift'] = coshift_preprocessing_max_shift
        spectrum.data_processing['spectra alignment']['fill_with_previous'] = fill_with_previous
        spectrum.data_processing['spectra alignment']['average2_multiplier'] = average2_multiplier

    return list_of_spectra


def get_feature_matrix(list_of_spectra, tolerance=0.05, decimals=4, missing_value_imputation=True):
    print('Creating feature intensity matrix')
    # get a consensus m/z array
    peak_picked_mz_arrays = [copy.deepcopy(spectrum.peak_picked_mz_array) for spectrum in list_of_spectra]
    peak_picked_consensus = pd.DataFrame(data={'mz': np.unique(np.concatenate(peak_picked_mz_arrays))}).sort_values(by='mz')
    preprocessed_mz_arrays = [copy.deepcopy(spectrum.get_mz_array()) for spectrum in list_of_spectra]
    preprocessed_concensus = pd.DataFrame(data={'mz': np.unique(np.concatenate(preprocessed_mz_arrays))}).sort_values(by='mz')

    spectra_dfs_peak_picked = [peak_picked_consensus]
    spectra_dfs_preprocessed = [preprocessed_concensus]
    for spectrum in list_of_spectra:
        spectra_dfs_peak_picked.append(pd.DataFrame(data={'mz': spectrum.peak_picked_mz_array,
                                                          spectrum.spectrum_id: spectrum.get_intensity_array()}))
        spectra_dfs_preprocessed.append(pd.DataFrame(data={'mz': spectrum.get_mz_array(),
                                                           spectrum.spectrum_id: spectrum.get_intensity_array()}))
    feature_matrix = reduce(lambda x, y: pd.merge_asof(x,
                                                       y,
                                                       on='mz',
                                                       tolerance=tolerance,
                                                       direction='nearest'),
                            spectra_dfs_peak_picked).sort_values(by='mz')
    if missing_value_imputation:
        ref_matrix = reduce(lambda x, y: pd.merge_asof(x,
                                                       y,
                                                       on='mz',
                                                       tolerance=tolerance,
                                                       direction='nearest'),
                            spectra_dfs_preprocessed).sort_values(by='mz')
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


def export_feature_matrix(feature_matrix, output):
    print('Exporting feature matrix to ' + output)
    feature_matrix.to_csv(output, index=False)
