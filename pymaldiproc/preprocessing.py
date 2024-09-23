import copy
import numpy as np
import pandas as pd
from functools import reduce


def get_feature_matrix(list_of_spectra, tolerance=0.05, decimals=4, missing_value_imputation=True):
    """
    Create a feature matrix from two or more preprocessed spectra in which peak picking has been performed.

    :param list_of_spectra: List of spectrum objects.
    :type list_of_spectra: list[pymaldiproc.classes.OpenMALDISpectrum|pymaldiproc.classes.PMPTsfSpectrum|pymaldiproc.classes.PMP2DTdfSpectrum]
    :param tolerance: Tolerance in Da used to consider whether features from different spectra are the same feature.
    :type tolerance: float
    :param decimals: Number of decimal places to report for m/z values.
    :type decimals: int
    :param missing_value_imputation: Whether to perform missing value imputation for peak picked features missing a
        corresponding intensity value. Missing values are obtained from corresponding preprocessed intensity arrays.
    :type missing_value_imputation: bool
    :return: Feature matrix with intensity values from each spectrum.
    :rtype: pandas.DataFrame
    """
    print('Creating feature intensity matrix')
    # get a consensus m/z array
    peak_picked_mz_arrays = [copy.deepcopy(spectrum.peak_picked_mz_array) for spectrum in list_of_spectra]
    peak_picked_consensus = pd.DataFrame(data={'mz': np.unique(np.concatenate(peak_picked_mz_arrays))}).sort_values(by='mz')
    preprocessed_mz_arrays = [copy.deepcopy(spectrum.preprocessed_mz_array) for spectrum in list_of_spectra]
    preprocessed_concensus = pd.DataFrame(data={'mz': np.unique(np.concatenate(preprocessed_mz_arrays))}).sort_values(by='mz')

    spectra_dfs_peak_picked = [peak_picked_consensus]
    spectra_dfs_preprocessed = [preprocessed_concensus]
    for spectrum in list_of_spectra:
        spectra_dfs_peak_picked.append(pd.DataFrame(data={'mz': spectrum.peak_picked_mz_array,
                                                          spectrum.spectrum_id: spectrum.peak_picked_intensity_array}))
        spectra_dfs_preprocessed.append(pd.DataFrame(data={'mz': spectrum.preprocessed_mz_array,
                                                           spectrum.spectrum_id: spectrum.preprocessed_intensity_array}))
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
    feature_matrix = feature_matrix.groupby(feature_matrix.columns[1:].values.tolist()).mean()
    feature_matrix = feature_matrix.sort_values(by='mz').round(decimals).reset_index()
    feature_matrix = feature_matrix.groupby(['mz'], as_index=False).aggregate(sum)
    return feature_matrix


def export_feature_matrix(feature_matrix, output):
    """
    Export a feature matrix to a CSV file.

    :param feature_matrix: Feature matrix with intensity values from each spectrum.
    :type feature_matrix: pandas.DataFrame
    :param output: Output CSV file path.
    :type output: str
    """
    print('Exporting feature matrix to ' + output)
    feature_matrix.to_csv(output, index=False)
