import copy
import gc
import numpy as np
import pandas as pd
from functools import reduce
from msalign import Aligner
import matplotlib.pyplot as plt


def align_spectra(list_of_spectra, lower_mass_range, upper_mass_range, method='pchip', ref_index=0, ref_peaks=None,
                  ref_peaks_weights=None, num_ref_peaks=10, n_bins=0, width=10.0, ratio=2.5, resolution=100,
                  iterations=5, grid_steps=20, shift_range=(-100, 100), align_by_index=True, only_shift=False,
                  return_shifts=True, plot=False):
    """
    Align two or more similar spectra to account for mass shifts across replicates. If not already done, spectrum
    binning and peak picking using the local maxima algorithm and default parameters will be applied.

    :param list_of_spectra: List of spectrum objects.
    :type list_of_spectra: list[pymaldiproc.classes.OpenMALDISpectrum|pymaldiproc.classes.PMPTsfSpectrum|pymaldiproc.classes.PMP2DTdfSpectrum]
    :param lower_mass_range: Mass in daltons to use for the lower mass range.
    :type lower_mass_range: int
    :param upper_mass_range: Mass in Daltons to use for the upper mass range.
    :type upper_mass_range: int
    :param method: Interpolation method to be used for alignment: "pchip" or "cubic".
    :type method: str
    :param ref_index: Index of the spectrum (in which the first spectrum would be index 0) in the provided list of
        spectra to be used as the reference spectrum to align to.
    :type ref_index: int
    :param ref_peaks: List of reference peaks to be aligned to from the reference spectrum. These peaks will act as
        (re)calibration points in each spectrum and should be able to be detected in all spectra for best results.
    :type ref_peaks: list|numpy.array
    :param ref_peaks_weights: List of weights associated with the list of reference peaks. Must be the same length of
        ref_peaks.
    :type ref_peaks_weights: list|numpy.array
    :param num_ref_peaks: If no reference peaks are provided, reference peaks up to this maximum number of reference
        peaks will be automatically detected and used for spectra alignment.
    :type num_ref_peaks: int
    :param n_bins: Number of bins to use for spectra binning if binning has not already been performed. If 0, the number
        of bins will automatically be calculated using a bin width of 0.05 Da.
    :type n_bins: int
    :param width: Width of the gaussian peak in separation units.
    :type width: float
    :param ratio: Scaling value that determines the size of the window around every alignment peak. The synthetic
        signal is compared to the input signal within these regions.
    :type ratio: float
    :param resolution: Resolution of the peaks.
    :type resolution: int
    :param iterations: Number of iterations. Increasing this value will (slightly) slow down the function but will
        improve performance.
    :param grid_steps: Number of steps to be used in the grid search.
    :type grid_steps: int
    :param shift_range: Maximum allowed shifts.
    :type shift_range: list|numpy.array
    :param align_by_index: Decide whether alignment should be done based on index rather than the preprocessed m/z
        array.
    :type align_by_index: bool
    :param only_shift: Determines if the signal should be shifted (True) or rescaled (False).
    :type only_shift: bool
    :param return_shifts: Decide whether shift parameter "shift_opt" should also be returned.
    :type return_shifts: bool
    :param plot: Plot the overlaid spectra using matplotlib
    :type plot: bool
    :return: List of spectrum objects.
    :rtype: list[pymaldiproc.classes.OpenMALDISpectrum|pymaldiproc.classes.PMPTsfSpectrum|pymaldiproc.classes.PMP2DTdfSpectrum]
    """
    for spectrum in list_of_spectra:
        # bin spectra if not already done
        if 'spectra binning' not in spectrum:
            if n_bins == 0:
                n_bins = int((upper_mass_range - lower_mass_range) / 0.05)
            spectrum.bin_spectrum(n_bins=n_bins,
                                  lower_mass_range=lower_mass_range,
                                  upper_mass_range=upper_mass_range)
        # pick peaks if not already done
        if spectrum.peak_picked_mz_array is None or \
                spectrum.peak_picked_intensity_array is None or \
                spectrum.peak_picking_indices is None:
            spectrum.peak_picking(method='locmax')

    # get 2D array of intensities
    intensity_array_2d = np.stack([copy.deepcopy(spectrum.preprocessed_intensity_array)
                                   for spectrum in list_of_spectra])
    # use first spectrum in the list as a reference spectrum if no reference peaks are provided by the user
    ref_mz_array = copy.deepcopy(list_of_spectra[ref_index].preprocessed_mz_array)
    if ref_peaks is None:
        ref_peaks = pd.DataFrame({'m/z': copy.deepcopy(list_of_spectra[ref_index].peak_picked_mz_array),
                                  'Intensity': copy.deepcopy(list_of_spectra[ref_index].peak_picked_intensity_array)})
        subset_step = int((upper_mass_range - lower_mass_range) / num_ref_peaks)
        ref_peaks_indices = []
        for i in range(lower_mass_range, upper_mass_range, subset_step):
            if i == lower_mass_range:
                pass
            else:
                subset_peaks = ref_peaks[((i - subset_step) <= ref_peaks['m/z']) & (ref_peaks['m/z'] < i)]
                if not subset_peaks.empty:
                    ref_peaks_indices.append(np.where(ref_peaks['m/z'] ==
                                                      subset_peaks.sort_values(by='Intensity',
                                                                               ascending=False)['m/z'].values[0])[0][0])
        ref_peaks = copy.deepcopy(list_of_spectra[ref_index].peak_picked_mz_array)[ref_peaks_indices]
    # align all intensity arrays based on reference m/z array and reference points (i.e. calibration points)
    aligner = Aligner(x=ref_mz_array,
                      array=intensity_array_2d,
                      method=method,
                      peaks=ref_peaks,
                      weights=ref_peaks_weights,
                      width=width,
                      ratio=ratio,
                      resolution=resolution,
                      iterations=iterations,
                      grid_steps=grid_steps,
                      shift_range=shift_range,
                      return_shifts=return_shifts,
                      align_by_index=align_by_index,
                      only_shift=only_shift)
    aligner.run()
    aligned_intensity_array_2d, shifts_out = aligner.apply()
    # write aligned intensity arrays to spectra and remove previous peak picking results
    for i in range(0, len(list_of_spectra)):
        list_of_spectra[i].preprocessed_intensity_array = copy.deepcopy(aligned_intensity_array_2d[i,])
        list_of_spectra[i].peak_picked_mz_array = None
        list_of_spectra[i].peak_picked_intensity_array = None
        list_of_spectra[i].peak_picking_indices = None
        del list_of_spectra[i].data_processing['peak picking']
        gc.collect()

    # plot
    def overlay_plot(ax, x, array, peak):
        """Generate overlay plot, showing each signal and the alignment peak(s)"""
        for i, y in enumerate(array):
            y = (y / y.max()) + (i * 0.2)
            ax.plot(x, y, lw=3)
        ax.axes.get_yaxis().set_visible(False)
        ax.set_xlabel("Index", fontsize=18)
        ax.set_xlim((x[0], x[-1]))
        ax.vlines(peak, *ax.get_ylim())
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 10))
    overlay_plot(ax[0], ref_mz_array, intensity_array_2d, ref_peaks)
    overlay_plot(ax[1], ref_mz_array, aligned_intensity_array_2d, ref_peaks)
    plt.show()

    return list_of_spectra


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
