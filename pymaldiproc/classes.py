import os
import copy
import gc
import pandas as pd
from uuid import uuid4
from pyTDFSDK.classes import TsfSpectrum, TdfSpectrum
from scipy.signal import savgol_filter, find_peaks, find_peaks_cwt, peak_widths
from pyMSpec.smoothing import apodization, rebin, fast_change, median
from pybaselines.smooth import snip, noise_median
from pybaselines.morphological import tophat
from BaselineRemoval import BaselineRemoval
from pyMSpec.normalisation import tic, rms, mad, sqrt
import seaborn as sns
import matplotlib.pyplot as plt
from pyopenms import *


class PMPMethods(object):
    """
    Class containing preprocessing methods for MALDI-TOF and MALDI-qTOF mass spectra.
    """
    def get_peak_list(self):
        """
        Obtain a peak picked (centroided) peak list. Performs peak picking if no peak list is found.

        :return: Pandas DataFrame containing m/z and intensity values.
        :rtype: pandas.DataFrame
        """
        if self.peak_picked_mz_array is None and self.peak_picked_intensity_array is None:
            # If no peak list found, peak picking with default settings is used.
            self.peak_picking()
        return pd.DataFrame(data={'mz': self.peak_picked_mz_array,
                                  'intensity': self.peak_picked_intensity_array})

    def trim_spectrum(self, lower_mass_range, upper_mass_range):
        """
        Trim the mass spectrum to only include features between the user specified lower and upper mass ranges
        (inclusive). m/z and intensity arrays are stored in PMPMethods.preprocessed_mz_array and
        PMP.preprocessed_intensity_array, respectively.

        :param lower_mass_range: Mass in daltons to use for the lower mass range.
        :type lower_mass_range: int
        :param upper_mass_range: Mass in Daltons to use for the upper mass range.
        :type upper_mass_range: int
        """
        indices = np.where((self.preprocessed_mz_array >= lower_mass_range) &
                           (self.preprocessed_mz_array <= upper_mass_range))[0]
        self.preprocessed_mz_array = copy.deepcopy(self.preprocessed_mz_array[indices])
        self.preprocessed_intensity_array = copy.deepcopy(self.preprocessed_intensity_array[indices])
        self.data_processing['spectrum trimming'] = {'lower mass range': lower_mass_range,
                                                     'upper mass range': upper_mass_range}
        gc.collect()

    def transform_intensity(self, method='sqrt'):
        """
        Apply a mathematical transformation to intensity values in the intensity array.

        :param method: Method to use for intensity transformation. Either square root ('sqrt'), natural log ('log'),
            log2 ('log2'), or log10 ('log10') transformation.
        :type method: str
        """
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
        """
        Perform baseline smoothing. Documentation for each algorithm is taken from their respective documentation pages.

        :param method: Method to use for baseline smoothing. Either Savitzky Golay ('SavitzkyGolay'), apodization
            ('apodization'), rebin ('rebin'), fast change ('fast_change'), or median ('median').
        :type method: str
        :param window_length: The length of the filter window (i.e. number of coefficients).
        :type window_length: int
        :param polyorder: The order of the polynomial used to fit the samples. Must be less than window_length.
        :type polyorder: int
        :param delta_mz: New m/z dimension bin width.
        :type delta_mz: float
        :param diff_thresh: Numeric change to remove.
        :type diff_thresh: float
        """
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
        """
        Perform baseline removal. Documentation for each algorithm is taken from their respective documentation pages.

        :param method: Method to use for baseline removal. Either statistics-sensitive non-linear iterative
            peak-clipping ('SNIP'), TopHat ('TopHat'), median ('Median'), ZhangFit ('ZhangFit'), modified polynomial
            fit ('ModPoly'), or improved modified polynomial fit ('IModPoly').
        :type method: str
        :param min_half_window: The minimum half window size used for morphological operations.
        :type min_half_window: int
        :param max_half_window: The maximum number of iterations/maximum half window size used for morphological
            operations. Should be (w-1)/2 where w is the index-based width of feature or peak.
        :type max_half_window: int
        :param decreasing: If False, will iterate through window sizes from 1 to max_half_window. If True, will reverse
            the order and iterate from max_half_window to 1 (gives smoother baseline).
        :type decreasing: bool
        :param smooth_half_window: The half window to use for smoothing the data. If greater than 0, will perform a
            moving average smooth on the data for each window to give better results for noisy data.
        :type smooth_half_window: int | None
        :param filter_order: If the measured data has a more complicated baseline consisting of other elements such as
            Compton edges, thena  higher filter_order should be selected.
        :type filter_order: int
        :param sigma: The standard deviation of the smoothing Gaussian kernal. If None, uses
            (2 * smooth_half_window + 1) / 6.
        :type sigma: float | None
        :param increment: The step size for iterating half windows.
        :type increment: int
        :param max_hits: The number of consecutive half windows that must produce the same morphological opening before
            accepting the half window as the optimum value.
        :type max_hits: int
        :param window_tol: The tolerance value for considering two morphological openings as equivalent.
        :type window_tol: float
        :param lambda_: Affects smoothness of the resulting background. The larger the lambda, the smoother the
            background.
        :type lambda_: int
        :param porder: Adaptive iteratively reweighted penalized least squares for baseline fitting.
        :type porder: int
        :param repetition: How many iterations to run.
        :type repetition: int | None
        :param degree: Polynomial degree.
        :type degree: int
        :param gradient: Gradient for polynomial loss. Measures incremental gain over each iteration. If gain in any
            iteration is less than this, further improvement will stop.
        :type gradient: float
        """
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
        """
        Apply normalization to intensity values in the intensity array.

        :param method: Method to use for normalizaton. Either total ion count ('tic'), root mean squared ('rms'),
            median absolute deviation ('mad'), or square root ('sqrt').
        :type method: str
        """
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
        """
        Bin spectrum to N bins between a user specified lower and upper mass range (inclusive).

        :param n_bins: Number of bins to use.
        :type n_bins: int
        :param lower_mass_range: Mass in daltons to use for the lower mass range.
        :type lower_mass_range: int
        :param upper_mass_range: Mass in Daltons to use for the upper mass range.
        :type upper_mass_range: int
        """
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
        """
        Estimate peaks widths using a first pass at peak picking via scipy.signal.peak_widths.

        :return: The widths for each peak in samples.
        :rtype: numpy.ndarray
        """
        peak_indices, peak_properties = find_peaks(copy.deepcopy(self.preprocessed_intensity_array))
        widths = peak_widths(copy.deepcopy(self.preprocessed_intensity_array), peak_indices)
        return widths[0]

    def peak_picking(self, method='cwt', widths=None, snr=3, deisotope=False, fragment_tolerance=0.05,
                     fragment_unit_ppm=False, min_charge=1, max_charge=1, keep_only_deisotoped=False, min_isopeaks=2,
                     max_isopeaks=10, make_single_charged=True, annotate_charge=False, annotate_iso_peak_count=False,
                     use_decreasing_model=True, start_intensity_check=1, add_up_intensity=False):
        """
        Perform peak picking on mass spectrum. Documentation for various parameters is taken from their respective
        documentation pages.

        :param method: Method to use for peak picking. Either local maxima ('locmax') or continuous wavelet
            transformation ('cwt').
        :type method: str
        :param widths: Required width of peaks in samples. If using 'cwt' method, used for calculating the CWT matrix.
            Range should cover the expected width of peaks of interest.
        :type widths: float | None
        :param snr: Minimum signal-to-noise ratio required to consider peak.
        :type snr: int
        :param deisotope: Whether to perform deisotoping/ion deconvolution. Deisotoping performed using
            pyopenms.Deisotoper.
        :type deisotope: bool
        :param fragment_tolerance: The tolerance used to match isotopic peaks.
        :type fragment_tolerance: float
        :param fragment_unit_ppm: Whether ppm or m/z is used as tolerance.
        :type fragment_unit_ppm: bool
        :param min_charge: The minimum charge considered.
        :type min_charge: int
        :param max_charge: The maximum charge considered.
        :type max_charge: int
        :param keep_only_deisotoped: If True, only monoisotopic peaks of fragments with isotopic pattern are retained.
        :type keep_only_deisotoped: bool
        :param min_isopeaks: The minimum number of isotopic peaks (at least 2) required for an isotopic cluster.
        :type min_isopeaks: int
        :param max_isopeaks: The maximum number of isotopic peaks (at least 2) required for an isotopic cluster.
        :type max_isopeaks: int
        :param make_single_charged: Whether to convert deisotoped monoisotopic peak to single charge.
        :type make_single_charged: bool
        :param annotate_charge: Whether to annotate the charge to the peaks in pyopenms.MSSpectrum.IntegerDataArray:
            'charge'.
        :type annotate_charge: bool
        :param annotate_iso_peak_count: Whether to annotate the number of isotopic peaks in a pattern for each
            monoisotopic peak in pyopenms.MSSpectrum.IntegerDataArray: 'iso_peak_count'.
        :type annotate_iso_peak_count: bool
        :param use_decreasing_model: Whether to use a simple averagine model that expects heavier isotopes to have less
            intensity. If False, no intensity checks are applied.
        :type use_decreasing_model: bool
        :param start_intensity_check: Number of the isotopic peak from which the decreasing model should be applied.
            <= 1 will force the monoisotopic peak to be most intense. 2 will allow the monoisotopic peak to be less
            intense than the 2nd peak. 3 will allow the monoisotopic peak and the 2nd peak to be less intense than the
            3rd, etc. A number higher than max_isopeaks will effectively disable use_decreasing_model completely.
        :type start_intensity_check: int
        :param add_up_intensity: Whether to sum up the total intensity of each isotopic pattern into the intensity of
            the reported monoisotopic peak.
        :type add_up_intensity: bool
        """
        # check method
        if method not in ['locmax', 'cwt']:
            raise Exception('Method must be "locmax" or "cwt"')

        self.data_processing['peak picking'] = {'method': method}
        if method == 'locmax':
            peak_indices, peak_properties = find_peaks(copy.deepcopy(self.preprocessed_intensity_array),
                                                       height=np.mean(self.preprocessed_intensity_array) * snr)
            self.peak_picking_indices = copy.deepcopy(peak_indices)
            self.peak_picked_mz_array = copy.deepcopy(self.preprocessed_mz_array)[peak_indices]
            self.peak_picked_intensity_array = copy.deepcopy(self.preprocessed_intensity_array)[peak_indices]
            self.data_processing['peak picking']['signal to noise ratio'] = snr
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
            self.data_processing['peak picking']['signal to noise ratio'] = snr
            self.data_processing['peak picking']['lower peak width'] = np.min(widths)
            self.data_processing['peak picking']['upper peak width'] = np.max(widths)

        if deisotope:
            poms_spec = MSSpectrum()
            poms_spec.set_peaks((copy.deepcopy(self.peak_picked_mz_array),
                                 copy.deepcopy(self.peak_picked_intensity_array)))
            Deisotoper.deisotopeAndSingleCharge(poms_spec,
                                                fragment_tolerance,
                                                fragment_unit_ppm,
                                                min_charge,
                                                max_charge,
                                                keep_only_deisotoped,
                                                min_isopeaks,
                                                max_isopeaks,
                                                make_single_charged,
                                                annotate_charge,
                                                annotate_iso_peak_count,
                                                use_decreasing_model,
                                                start_intensity_check,
                                                add_up_intensity)
            poms_spec_mz, poms_spec_intensity = poms_spec.get_peaks()
            self.peak_picking_indices = np.searchsorted(copy.deepcopy(self.preprocessed_mz_array),
                                                        copy.deepcopy(poms_spec_mz))
            self.peak_picked_mz_array = copy.deepcopy(poms_spec_mz)
            self.peak_picked_intensity_array = copy.deepcopy(poms_spec_intensity)

        gc.collect()

    def undo_all_processing(self):
        """
        Undo all processing that has been performed up to this and clear PMPMethods.data_processing dict.
        """
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
        """
        Plot and show a basic spectrum into a basic seaborn.lineplot.
        """
        spectrum_df = pd.DataFrame({'m/z': copy.deepcopy(self.preprocessed_mz_array),
                                    'Intensity': copy.deepcopy(self.preprocessed_intensity_array)})
        fig = sns.lineplot(data=spectrum_df, x='m/z', y='Intensity')
        plt.show()


class OpenMALDISpectrum(PMPMethods):
    """
    Class for parsing and storing spectrum metadata and data arrays from open format MALDI-TOF and MALDI-qTOF mass
    spectra in *.mzML and *.mzXML files. Preprocessing methods are inherited from pymaldiproc.classes.PMPMethods.

    :param pyteomics_dict: Dictionary containing spectrum metadata obtained from pyteomics.mzml.read or
        pyteomics.mzxml.read
    :type pyteomics_dict: dict
    :param filename: Path to the source file.
    :type filename: str
    """
    def __init__(self, pyteomics_dict, filename):
        """
        Constructor Method
        """
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
        """
        Parse metadata and spectral data from dictionary obtained from pyteomics.mzml.read or pyteomics.mzxml.read.

        :param pyteomics_dict: Dictionary containing spectrum metadata obtained from pyteomics.mzml.read or
            pyteomics.mzxml.read
        :type pyteomics_dict: dict
        """
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

    def info(self):
        """
        Get basic spectrum metadata.

        :return: Dictionary containing spectrum name, coordinate, MS level, mode, source file, and UUID.
        :rtype: dict
        """
        return {'Name': self.name,
                'Coordinate': self.coord,
                'MS Level': self.ms_level,
                'Mode': 'Centroid' if self.centroided else 'Profile',
                'Source File': self.source,
                'UUID': self.uuid}


class PMPTsfSpectrum(TsfSpectrum, PMPMethods):
    """
    Class for parsing and storing spectrum metadata and data arrays from open format MALDI-TOF and MALDI-qTOF mass
    spectra in Bruker TSF files. Data parsing methods are inherited from pyTDFSDK.classes.TsfSpectrum. Preprocessing
    methods are inherited from pymaldiproc.classes.PMPMethods.

    :param tsf_data: TsfData object containing metadata from analysis.tsf database.
    :type tsf_data: pyTDFSDK.classes.TsfData
    :param frame: ID of the frame of interest.
    :type frame: int
    :param mode: Data array mode, either "profile", "centroid", or "raw".
    :type mode: str
    :param profile_bins: Number of bins to bin spectrum to.
    :type profile_bins: int
    :param encoding: Encoding bit mode, either "64" or "32"
    :type encoding: int
    """
    def __init__(self, tsf_data, frame: int, mode: str, profile_bins=0, encoding=64):
        """
        Constructor Method
        """
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
        """
        Obtain the name of a spectrum and generate an ID based on a spectrum's name, position on the target plate, and
        ID.
        """
        self.name = str(os.path.splitext(os.path.split(self.source)[-1])[0]) + '_' + str(self.frame)
        self.spectrum_id = self.name + '|' + self.coord + '|' + self.uuid

    def info(self):
        """
        Get basic spectrum metadata.

        :return: Dictionary containing spectrum name, coordinate, MS level, mode, source file, and UUID.
        :rtype: dict
        """
        return {'Name': self.name,
                'Coordinate': self.coord,
                'MS Level': self.ms_level,
                'Mode': 'Centroid' if self.centroided else 'Profile',
                'Source File': self.source,
                'UUID': self.uuid}


class PMPTdfSpectrum(TdfSpectrum, PMPMethods):
    """
    Class for parsing and storing spectrum metadata and data arrays from open format MALDI-TIMS-TOF and MALDI-TIMS-qTOF
    mass spectra in Bruker TDF files. Data parsing methods are inherited from pyTDFSDK.classes.TdfSpectrum.
    Preprocessing methods are inherited from pymaldiproc.classes.PMPMethods.

    :param tdf_data: TdfData object containing metadata from analysis.tdf database.
    :type tdf_data: pyTDFSDK.classes.TdfData
    :param frame: ID of the frame of interest.
    :type frame: int
    :param mode: Data array mode, either "profile", "centroid", or "raw".
    :type mode: str
    :param precursor: ID of the precursor of interest for ddaPASEF data. If specified, overrides the frame ID during
        data parsing.
    :type precursor: int
    :param diapasef_window: Dictionary containing a row of metadata from the
        pyTDFSDK.classesTdfData.analysis['DiaFrameMsMsWindows'] table required for parsing diaPASEF data.
    :type diapasef_window: dict
    :param profile_bins: Number of bins to bin spectrum to.
    :type profile_bins: int
    :param encoding: Encoding bit mode, either "64" or "32"
    :type encoding: int
    """
    def __init__(self, tdf_data, frame: int, mode: str, precursor=0, diapasef_window=None, profile_bins=0, encoding=64,
                 exclude_mobility=False):
        """
        Construcor Method
        """
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
        """
        Obtain the name of a spectrum and generate an ID based on a spectrum's name, position on the target plate, and
        ID.
        """
        self.name = str(os.path.splitext(os.path.split(self.source)[-1])[0]) + '_' + str(self.frame)
        self.spectrum_id = self.name + '|' + self.coord + '|' + self.uuid

    def info(self):
        """
        Get basic spectrum metadata.

        :return: Dictionary containing spectrum name, coordinate, MS level, mode, source file, and UUID.
        :rtype: dict
        """
        return {'Name': self.name,
                'Coordinate': self.coord,
                'MS Level': self.ms_level,
                'Mode': 'Centroid' if self.centroided else 'Profile',
                'Source File': self.source,
                'UUID': self.uuid}
