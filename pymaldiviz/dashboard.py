import os
import copy
import gc
import pandas as pd
from pymaldiproc.data_import import import_mzml, import_timstof_raw_data
from pymaldiproc.classes import PMP3DTdfSpectrum
from pymaldiviz.layout import get_dashboard_layout
from pymaldiviz.util import (get_preprocessing_params, get_spectrum, toggle_rebin_style, toggle_apodization_style,
                             toggle_cwt_style, toggle_snip_style, toggle_locmax_style, toggle_tophat_style,
                             toggle_smoothing_median_style, toggle_modpoly_style, toggle_imodpoly_style,
                             toggle_zhangfit_style, toggle_deisotope_on_style, toggle_deisotope_off_style,
                             toggle_fast_change_style, toggle_savitzky_golay_style, toggle_removal_median_style,
                             cleanup_file_system_backend, get_peakmap)
from pymaldiviz.tmpdir import FILE_SYSTEM_BACKEND
from dash import State, callback_context, no_update
from dash_extensions.enrich import (Input, Output, DashProxy, MultiplexerTransform, Serverside,
                                    ServersideOutputTransform, FileSystemBackend)
import dash_bootstrap_components as dbc
from plotly_resampler import FigureResampler
import tkinter
from tkinter.filedialog import askopenfilenames, askdirectory, asksaveasfilename

# will be a dictionary of MALDISpectrum objects used for the spectrum plot
INDEXED_DATA = {}

# Use DashProxy instead of Dash to allow for multiple callbacks to the same plot
app = DashProxy(prevent_initial_callbacks=True,
                transforms=[MultiplexerTransform(),
                            ServersideOutputTransform(backends=[FileSystemBackend(cache_dir=FILE_SYSTEM_BACKEND)])],
                external_stylesheets=[dbc.themes.BOOTSTRAP])
app.layout = get_dashboard_layout(get_preprocessing_params())


@app.callback([Output('spectrum_id', 'options', allow_duplicate=True),
               Output('spectrum_id', 'value', allow_duplicate=True)],
              [Input('upload_mzml', 'n_clicks'),
               Input('upload_d', 'n_clicks')])
def upload_data(n_clicks_mzml, n_clicks_d):
    """
    Dash callback for upload data buttons. Obtains path to Bruker .d directory or mzML file and loads all spectra to a
    global variable INDEXED_DATA.

    :param n_clicks_mzml: Input signal if the upload_mzml button is clicked.
    :param n_clicks_d: Input signal if the upload_d button is clicked.
    :return: Tuple of options and value attributes for the dcc.Dropdown object.
    """
    global INDEXED_DATA
    changed_id = [i['prop_id'] for i in callback_context.triggered][0]
    if changed_id == 'upload_mzml.n_clicks':
        main_tk_window = tkinter.Tk()
        main_tk_window.attributes('-topmost', True, '-alpha', 0)
        filenames = askopenfilenames(filetypes=[('mzML Files', '*.mzML')])
        main_tk_window.destroy()
        for filename in filenames:
            data = import_mzml(filename)
            for spectrum in data:
                INDEXED_DATA[spectrum.spectrum_id] = spectrum
    elif changed_id == 'upload_d.n_clicks':
        main_tk_window = tkinter.Tk()
        main_tk_window.attributes('-topmost', True, '-alpha', 0)
        dirname = askdirectory(mustexist=True)
        main_tk_window.destroy()
        if dirname.endswith('.d'):
            data = import_timstof_raw_data(dirname, mode='profile', exclude_mobility=False)
            for spectrum in data:
                INDEXED_DATA[spectrum.spectrum_id] = spectrum
    options = [{'label': i, 'value': i} for i in INDEXED_DATA.keys()]
    value = [i for i in INDEXED_DATA.keys()]
    return options, value


@app.callback([Output('spectrum_plot', 'figure'),
               Output('store_plot', 'data'),
               Output('trim_spectrum', 'disabled'),
               Output('transform_intensity', 'disabled'),
               Output('smooth_baseline', 'disabled'),
               Output('remove_baseline', 'disabled'),
               Output('normalize_intensity', 'disabled'),
               Output('bin_spectrum', 'disabled'),
               Output('toggle_log_intensity', 'disabled')],
              Input('spectrum_id', 'value'))
def plot_spectrum(value):
    """
    Dash callback to plot the spectrum selected from the spectrum_id dropdown using plotly.express and
    plotly_resampler.FigureResampler. If a pymaldiproc.classes.PMP3DTdfSpectrum is detected, a 3D peakmap will be
    plotted instead.

    :param value: Input signal spectrum_id used as the key in INDEXED_DATA.
    :return: Tuple of spectrum figure as a plotly.express.line plot, data store for plotly_resampler, and whether or
        not preprocessing buttons are to be disabled.
    """
    global INDEXED_DATA
    if isinstance(INDEXED_DATA[value], PMP3DTdfSpectrum):
        fig = get_peakmap(INDEXED_DATA[value])
        return fig, None, True, True, True, True, True, True, False
    else:
        fig = get_spectrum(INDEXED_DATA[value])
        cleanup_file_system_backend(FILE_SYSTEM_BACKEND)
        return fig, Serverside(fig), False, False, False, False, False, False, True


@app.callback([Output('spectrum_plot', 'figure'),
               Output('store_use_log_intensity', 'data')],
              Input('toggle_log_intensity', 'n_clicks'),
              [State('spectrum_id', 'value'),
               State('store_use_log_intensity', 'data')])
def toggle_peakmap_log_intensity(n_clicks, value, use_log_intensity):
    global INDEXED_DATA
    changed_id = [i['prop_id'] for i in callback_context.triggered][0]
    if changed_id == 'toggle_log_intensity.n_clicks':
        if isinstance(INDEXED_DATA[value], PMP3DTdfSpectrum):
            use_log_intensity = not use_log_intensity
            fig = get_peakmap(INDEXED_DATA[value], use_log_intensity=use_log_intensity)
            return fig, use_log_intensity
        else:
            return no_update


@app.callback([Output('edit_processing_parameters_modal', 'is_open'),
               Output('store_preprocessing_params', 'data')],
              [Input('edit_preprocessing_parameters', 'n_clicks'),
               Input('edit_processing_parameters_save', 'n_clicks'),
               Input('edit_processing_parameters_cancel', 'n_clicks'),
               Input('trim_spectrum_lower_mass_range_value', 'value'),
               Input('trim_spectrum_upper_mass_range_value', 'value'),
               Input('transform_intensity_method', 'value'),
               Input('smooth_baseline_method', 'value'),
               Input('smooth_baseline_window_length_value', 'value'),
               Input('smooth_baseline_polyorder_value', 'value'),
               Input('smooth_baseline_delta_mz_value', 'value'),
               Input('smooth_baseline_diff_thresh_value', 'value'),
               Input('remove_baseline_method', 'value'),
               Input('remove_baseline_min_half_window_value', 'value'),
               Input('remove_baseline_max_half_window_value', 'value'),
               Input('remove_baseline_decreasing', 'value'),
               Input('remove_baseline_smooth_half_window_value', 'value'),
               Input('remove_baseline_filter_order_value', 'value'),
               Input('remove_baseline_sigma_value', 'value'),
               Input('remove_baseline_increment_value', 'value'),
               Input('remove_baseline_max_hits_value', 'value'),
               Input('remove_baseline_window_tol_value', 'value'),
               Input('remove_baseline_lambda__value', 'value'),
               Input('remove_baseline_porder_value', 'value'),
               Input('remove_baseline_repetition_value', 'value'),
               Input('remove_baseline_degree_value', 'value'),
               Input('remove_baseline_gradient_value', 'value'),
               Input('normalize_intensity_method', 'value'),
               Input('bin_spectrum_n_bins_value', 'value'),
               Input('bin_spectrum_lower_mass_range_value', 'value'),
               Input('bin_spectrum_upper_mass_range_value', 'value'),
               Input('peak_picking_method', 'value'),
               Input('peak_picking_snr_value', 'value'),
               Input('peak_picking_widths_value', 'value'),
               Input('peak_picking_deisotope', 'value'),
               Input('peak_picking_deisotope_fragment_tolerance_value', 'value'),
               Input('peak_picking_deisotope_fragment_unit_ppm', 'value'),
               Input('peak_picking_deisotope_min_charge_value', 'value'),
               Input('peak_picking_deisotope_max_charge_value', 'value'),
               Input('peak_picking_deisotope_keep_only_deisotoped', 'value'),
               Input('peak_picking_deisotope_min_isopeaks_value', 'value'),
               Input('peak_picking_deisotope_max_isopeaks_value', 'value'),
               Input('peak_picking_deisotope_make_single_charged', 'value'),
               Input('peak_picking_deisotope_annotate_charge', 'value'),
               Input('peak_picking_deisotope_annotate_iso_peak_count', 'value'),
               Input('peak_picking_deisotope_use_decreasing_model', 'value'),
               Input('peak_picking_deisotope_start_intensity_check_value', 'value'),
               Input('peak_picking_deisotope_add_up_intensity', 'value'),
               Input('peak_picking_3d_min_distance_value', 'value'),
               Input('peak_picking_3d_noise_value', 'value'),
               Input('peak_picking_3d_snr_value', 'value'),
               Input('peak_picking_3d_exclude_border_value', 'value')],
              [State('store_preprocessing_params', 'data'),
               State('edit_processing_parameters_modal', 'is_open')])
def toggle_edit_preprocessing_parameters_modal(n_clicks_button,
                                               n_clicks_save,
                                               n_clicks_cancel,
                                               trim_spectrum_lower_mass_range,
                                               trim_spectrum_upper_mass_range,
                                               transform_intensity_method,
                                               smooth_baseline_method,
                                               smooth_baseline_window_length,
                                               smooth_baseline_polyorder,
                                               smooth_baseline_delta_mz,
                                               smooth_baseline_diff_thresh,
                                               remove_baseline_method,
                                               remove_baseline_min_half_window,
                                               remove_baseline_max_half_window,
                                               remove_baseline_decreasing,
                                               remove_baseline_smooth_half_window,
                                               remove_baseline_filter_order,
                                               remove_baseline_sigma,
                                               remove_baseline_increment,
                                               remove_baseline_max_hits,
                                               remove_baseline_window_tol,
                                               remove_baseline_lambda_,
                                               remove_baseline_porder,
                                               remove_baseline_repetition,
                                               remove_baseline_degree,
                                               remove_baseline_gradient,
                                               normalize_intensity_method,
                                               bin_spectrum_n_bins,
                                               bin_spectrum_lower_mass_range,
                                               bin_spectrum_upper_mass_range,
                                               peak_picking_method,
                                               peak_picking_snr,
                                               peak_picking_widths,
                                               peak_picking_deisotope,
                                               peak_picking_fragment_tolerance,
                                               peak_picking_fragment_unit_ppm,
                                               peak_picking_min_charge,
                                               peak_picking_max_charge,
                                               peak_picking_keep_only_deisotoped,
                                               peak_picking_min_isopeaks,
                                               peak_picking_max_isopeaks,
                                               peak_picking_make_single_charged,
                                               peak_picking_annotate_charge,
                                               peak_picking_annotate_iso_peak_count,
                                               peak_picking_use_decreasing_model,
                                               peak_picking_start_intensity_check,
                                               peak_picking_add_up_intensity,
                                               peak_picking_3d_min_distance,
                                               peak_picking_3d_noise,
                                               peak_picking_3d_snr,
                                               peak_picking_3d_exclude_border,
                                               preprocessing_params,
                                               is_open):
    """
    Dash callback to toggle the preprocessing parameters modal window, populate the current preprocessing parameters
    saved in the global variable preprocessing_params, and save any modified preprocessing parameters to
    preprocessing_params if the Save button is clicked.

    :param n_clicks_button: Input signal if the edit_preprocessing_parameters button is clicked.
    :param n_clicks_save: Input signal if the edit_preprocessing_parameters_save button is clicked.
    :param n_clicks_cancel: Input signal if the edit_preprocessing_parameters_cancel button is clicked.
    :param trim_spectrum_lower_mass_range: Mass in daltons to use for the lower mass range during spectrum trimming.
    :param trim_spectrum_upper_mass_range: Mass in Daltons to use for the upper mass range during spectrum trimming.
    :param transform_intensity_method: Method to use for intensity transformation. Either square root ('sqrt'), natural
        log ('log'), log2 ('log2'), or log10 ('log10') transformation.
    :param smooth_baseline_method: Method to use for baseline smoothing. Either Savitzky Golay ('SavitzkyGolay'),
        apodization ('apodization'), rebin ('rebin'), fast change ('fast_change'), or median ('median').
    :param smooth_baseline_window_length: The length of the filter window (i.e. number of coefficients).
    :param smooth_baseline_polyorder: The order of the polynomial used to fit the samples. Must be less than
        window_length.
    :param smooth_baseline_delta_mz: New m/z dimension bin width.
    :param smooth_baseline_diff_thresh: Numeric change to remove.
    :param remove_baseline_method: Method to use for baseline removal. Either statistics-sensitive non-linear iterative
        peak-clipping ('SNIP'), TopHat ('TopHat'), median ('Median'), ZhangFit ('ZhangFit'), modified polynomial
        fit ('ModPoly'), or improved modified polynomial fit ('IModPoly').
    :param remove_baseline_min_half_window: The minimum half window size used for morphological operations.
    :param remove_baseline_max_half_window: The maximum number of iterations/maximum half window size used for
        morphological operations. Should be (w-1)/2 where w is the index-based width of feature or peak.
    :param remove_baseline_decreasing: If False, will iterate through window sizes from 1 to max_half_window. If True,
        will reverse the order and iterate from max_half_window to 1 (gives smoother baseline).
    :param remove_baseline_smooth_half_window: The half window to use for smoothing the data. If greater than 0, will
        perform a moving average smooth on the data for each window to give better results for noisy data.
    :param remove_baseline_filter_order: If the measured data has a more complicated baseline consisting of other
        elements such as Compton edges, thena  higher filter_order should be selected.
    :param remove_baseline_sigma: The standard deviation of the smoothing Gaussian kernal. If None, uses
        (2 * smooth_half_window + 1) / 6.
    :param remove_baseline_increment: The step size for iterating half windows.
    :param remove_baseline_max_hits: The number of consecutive half windows that must produce the same morphological
        opening before accepting the half window as the optimum value.
    :param remove_baseline_window_tol: The tolerance value for considering two morphological openings as equivalent.
    :param remove_baseline_lambda_: Affects smoothness of the resulting background. The larger the lambda, the smoother
        the background.
    :param remove_baseline_porder: Adaptive iteratively reweighted penalized least squares for baseline fitting.
    :param remove_baseline_repetition: How many iterations to run.
    :param remove_baseline_degree: Polynomial degree.
    :param remove_baseline_gradient: Gradient for polynomial loss. Measures incremental gain over each iteration. If
        gain in any iteration is less than this, further improvement will stop.
    :param normalize_intensity_method: Method to use for normalizaton. Either total ion count ('tic'), root mean
        squared ('rms'), median absolute deviation ('mad'), or square root ('sqrt').
    :param bin_spectrum_n_bins: Number of bins to use.
    :param bin_spectrum_lower_mass_range: Mass in daltons to use for the lower mass range during spectrum binning.
    :param bin_spectrum_upper_mass_range: Mass in Daltons to use for the upper mass range during spectrum binning.
    :param peak_picking_method: Method to use for peak picking. Either local maxima ('locmax') or continuous wavelet
        transformation ('cwt').
    :param peak_picking_snr: Minimum signal-to-noise ratio required to consider peak.
    :param peak_picking_widths: Required width of peaks in samples. If using 'cwt' method, used for calculating the CWT
        matrix. Range should cover the expected width of peaks of interest.
    :param peak_picking_deisotope: Whether to perform deisotoping/ion deconvolution. Deisotoping performed using
        pyopenms.Deisotoper.
    :param peak_picking_fragment_tolerance: The tolerance used to match isotopic peaks.
    :param peak_picking_fragment_unit_ppm: Whether ppm or m/z is used as tolerance.
    :param peak_picking_min_charge: The minimum charge considered.
    :param peak_picking_max_charge: The maximum charge considered.
    :param peak_picking_keep_only_deisotoped: If True, only monoisotopic peaks of fragments with isotopic pattern are
        retained.
    :param peak_picking_min_isopeaks: The minimum number of isotopic peaks (at least 2) required for an isotopic
        cluster.
    :param peak_picking_max_isopeaks: The maximum number of isotopic peaks (at least 2) required for an isotopic
        cluster.
    :param peak_picking_make_single_charged: Whether to convert deisotoped monoisotopic peak to single charge.
    :param peak_picking_annotate_charge: Whether to annotate the charge to the peaks in
        pyopenms.MSSpectrum.IntegerDataArray: 'charge'.
    :param peak_picking_annotate_iso_peak_count: Whether to annotate the number of isotopic peaks in a pattern for each
        monoisotopic peak in pyopenms.MSSpectrum.IntegerDataArray: 'iso_peak_count'.
    :param peak_picking_use_decreasing_model: Whether to use a simple averagine model that expects heavier isotopes to
        have less intensity. If False, no intensity checks are applied.
    :param peak_picking_start_intensity_check: Number of the isotopic peak from which the decreasing model should be
        applied. <= 1 will force the monoisotopic peak to be most intense. 2 will allow the monoisotopic peak to be
        less intense than the 2nd peak. 3 will allow the monoisotopic peak and the 2nd peak to be less intense than the
        3rd, etc. A number higher than max_isopeaks will effectively disable use_decreasing_model completely.
    :param peak_picking_add_up_intensity: Whether to sum up the total intensity of each isotopic pattern into the
        intensity of the reported monoisotopic peak.
    :param peak_picking_3d_min_distance: The minimal allowed distance separating peaks. To find the maximum number of
        peaks, use min_distance=1.
    :param peak_picking_3d_noise: Absolute intensity value to be used as baseline noise level. If None, the median
        value of the intensity array is used as an estimate of the noise.
    :param peak_picking_3d_snr: Minimum signal-to-noise ratio required to consider peak.
    :param peak_picking_3d_exclude_border: If positive integer, excludes peaks from within n pixels of the border of
        the image. If tuple of non-negative ints, the length of the tuple must match the input array's dimensionality.
        Each element of the tuple will exclude peaks within n pixels of the border of the image along that dimension.
        If True, takes the min_distance parameter as value. If zero or False, peaks are identified regardless of their
        distance from the border.
    :param preprocessing_params: State signal containing data from store_preprocessing_params.
    :param is_open: State signal to determine whether the edit_preprocessing_parameters_modal modal window is open.
    :return: Output signal to determine whether the edit_preprocessing_parameters_modal modal window is open.
    """
    changed_id = [i['prop_id'] for i in callback_context.triggered][0]
    print(changed_id)
    if (changed_id == 'edit_preprocessing_parameters.n_clicks' or
            changed_id == 'edit_processing_parameters_save.n_clicks' or
            changed_id == 'edit_processing_parameters_cancel.n_clicks'):
        if changed_id == 'edit_processing_parameters_save.n_clicks':
            preprocessing_params['TRIM_SPECTRUM']['lower_mass_range'] = trim_spectrum_lower_mass_range
            preprocessing_params['TRIM_SPECTRUM']['upper_mass_range'] = trim_spectrum_upper_mass_range
            preprocessing_params['TRANSFORM_INTENSITY']['method'] = transform_intensity_method
            preprocessing_params['SMOOTH_BASELINE']['method'] = smooth_baseline_method
            preprocessing_params['SMOOTH_BASELINE']['window_length'] = smooth_baseline_window_length
            preprocessing_params['SMOOTH_BASELINE']['polyorder'] = smooth_baseline_polyorder
            preprocessing_params['SMOOTH_BASELINE']['delta_mz'] = smooth_baseline_delta_mz
            preprocessing_params['SMOOTH_BASELINE']['diff_thresh'] = smooth_baseline_diff_thresh
            preprocessing_params['REMOVE_BASELINE']['method'] = remove_baseline_method
            preprocessing_params['REMOVE_BASELINE']['min_half_window'] = remove_baseline_min_half_window
            preprocessing_params['REMOVE_BASELINE']['max_half_window'] = remove_baseline_max_half_window
            preprocessing_params['REMOVE_BASELINE']['decreasing'] = remove_baseline_decreasing
            preprocessing_params['REMOVE_BASELINE']['smooth_half_window'] = remove_baseline_smooth_half_window
            preprocessing_params['REMOVE_BASELINE']['filter_order'] = remove_baseline_filter_order
            preprocessing_params['REMOVE_BASELINE']['sigma'] = remove_baseline_sigma
            preprocessing_params['REMOVE_BASELINE']['increment'] = remove_baseline_increment
            preprocessing_params['REMOVE_BASELINE']['max_hits'] = remove_baseline_max_hits
            preprocessing_params['REMOVE_BASELINE']['window_tol'] = remove_baseline_window_tol
            preprocessing_params['REMOVE_BASELINE']['lambda_'] = remove_baseline_lambda_
            preprocessing_params['REMOVE_BASELINE']['porder'] = remove_baseline_porder
            preprocessing_params['REMOVE_BASELINE']['repetition'] = remove_baseline_repetition
            preprocessing_params['REMOVE_BASELINE']['degree'] = remove_baseline_degree
            preprocessing_params['REMOVE_BASELINE']['gradient'] = remove_baseline_gradient
            preprocessing_params['NORMALIZE_INTENSITY']['method'] = normalize_intensity_method
            preprocessing_params['BIN_SPECTRUM']['n_bins'] = bin_spectrum_n_bins
            preprocessing_params['BIN_SPECTRUM']['lower_mass_range'] = bin_spectrum_lower_mass_range
            preprocessing_params['BIN_SPECTRUM']['upper_mass_range'] = bin_spectrum_upper_mass_range
            preprocessing_params['PEAK_PICKING']['method'] = peak_picking_method
            preprocessing_params['PEAK_PICKING']['snr'] = peak_picking_snr
            preprocessing_params['PEAK_PICKING']['widths'] = peak_picking_widths
            preprocessing_params['PEAK_PICKING']['deisotope'] = peak_picking_deisotope
            preprocessing_params['PEAK_PICKING']['fragment_tolerance'] = peak_picking_fragment_tolerance
            preprocessing_params['PEAK_PICKING']['fragment_unit_ppm'] = peak_picking_fragment_unit_ppm
            preprocessing_params['PEAK_PICKING']['min_charge'] = peak_picking_min_charge
            preprocessing_params['PEAK_PICKING']['max_charge'] = peak_picking_max_charge
            preprocessing_params['PEAK_PICKING']['keep_only_deisotoped'] = peak_picking_keep_only_deisotoped
            preprocessing_params['PEAK_PICKING']['min_isopeaks'] = peak_picking_min_isopeaks
            preprocessing_params['PEAK_PICKING']['max_isopeaks'] = peak_picking_max_isopeaks
            preprocessing_params['PEAK_PICKING']['make_single_charged'] = peak_picking_make_single_charged
            preprocessing_params['PEAK_PICKING']['annotate_charge'] = peak_picking_annotate_charge
            preprocessing_params['PEAK_PICKING']['annotate_iso_peak_count'] = peak_picking_annotate_iso_peak_count
            preprocessing_params['PEAK_PICKING']['use_decreasing_model'] = peak_picking_use_decreasing_model
            preprocessing_params['PEAK_PICKING']['start_intensity_check'] = peak_picking_start_intensity_check
            preprocessing_params['PEAK_PICKING']['add_up_intensity'] = peak_picking_add_up_intensity
            preprocessing_params['PEAK_PICKING_3D']['min_distance'] = peak_picking_3d_min_distance
            preprocessing_params['PEAK_PICKING_3D']['noise'] = peak_picking_3d_noise
            preprocessing_params['PEAK_PICKING_3D']['snr'] = peak_picking_3d_snr
            preprocessing_params['PEAK_PICKING_3D']['exclude_border'] = peak_picking_3d_exclude_border
        return not is_open, preprocessing_params
    return is_open, preprocessing_params


@app.callback(Output('edit_processing_parameters_modal_saved', 'is_open'),
              [Input('edit_processing_parameters_save', 'n_clicks'),
               Input('edit_processing_parameters_modal_saved_close', 'n_clicks')],
              State('edit_processing_parameters_modal_saved', 'is_open'))
def toggle_edit_processing_parameters_saved_modal(n_clicks_save, n_clicks_close, is_open):
    """
    Dash callback to toggle the preprocessing parameters save confirmation message modal window.

    :param n_clicks_save: Input signal if the edit_preprocessing_parameters_save button is clicked.
    :param n_clicks_close: Input signal if the edit_preprocessing_parameters_close button is clicked.
    :param is_open: State signal to determine whether the edit_preprocessing_parameters_modal_saved modal window is
        open.
    :return: Output signal to determine whether the edit_preprocessing_parameters_modal_saved modal window is open.
    """
    if n_clicks_save or n_clicks_close:
        return not is_open
    return is_open


@app.callback([Output('smooth_baseline_window_length', 'style'),
               Output('smooth_baseline_polyorder', 'style'),
               Output('smooth_baseline_delta_mz', 'style'),
               Output('smooth_baseline_diff_thresh', 'style')],
              [Input('edit_preprocessing_parameters', 'n_clicks'),
               Input('smooth_baseline_method', 'value')])
def toggle_smooth_baseline_method_parameters(n_clicks, value):
    """
    Dash callback to toggle which baseline smoothing parameters are visible depending on the baseline smoothing method
    selected in the preprocessing parameters modal window.

    :param n_clicks: Input signal if the edit_preprocessing_parameters button is clicked.
    :param value: Input signal to obtain the currently selected baseline smoothing method.
    :return: List of dictionaries containing style template to show or hide parameters.
    """
    if value == 'SavitzkyGolay':
        return toggle_savitzky_golay_style()
    elif value == 'apodization':
        return toggle_apodization_style()
    elif value == 'rebin':
        return toggle_rebin_style()
    elif value == 'fast_change':
        return toggle_fast_change_style()
    elif value == 'median':
        return toggle_smoothing_median_style()


@app.callback([Output('remove_baseline_min_half_window', 'style'),
               Output('remove_baseline_max_half_window', 'style'),
               Output('remove_baseline_decreasing', 'style'),
               Output('remove_baseline_smooth_half_window', 'style'),
               Output('remove_baseline_filter_order', 'style'),
               Output('remove_baseline_sigma', 'style'),
               Output('remove_baseline_increment', 'style'),
               Output('remove_baseline_max_hits', 'style'),
               Output('remove_baseline_window_tol', 'style'),
               Output('remove_baseline_lambda_', 'style'),
               Output('remove_baseline_porder', 'style'),
               Output('remove_baseline_repetition', 'style'),
               Output('remove_baseline_degree', 'style'),
               Output('remove_baseline_gradient', 'style')],
              [Input('edit_preprocessing_parameters', 'n_clicks'),
               Input('remove_baseline_method', 'value')])
def toggle_remove_baseline_method_parameters(n_clicks, value):
    """
    Dash callback to toggle which baseline removal parameters are visible depending on the baseline removal method
    selected in the preprocessing parameters modal window.

    :param n_clicks: Input signal if the edit_preprocessing_parameters button is clicked.
    :param value: Input signal to obtain the currently selected baseline removal method.
    :return: List of dictionaries containing style template to show or hide parameters.
    """
    if value == 'SNIP':
        return toggle_snip_style()
    elif value == 'TopHat':
        return toggle_tophat_style()
    elif value == 'Median':
        return toggle_removal_median_style()
    elif value == 'ZhangFit':
        return toggle_zhangfit_style()
    elif value == 'ModPoly':
        return toggle_modpoly_style()
    elif value == 'IModPoly':
        return toggle_imodpoly_style()


@app.callback([Output('peak_picking_snr', 'style'),
               Output('peak_picking_widths', 'style')],
              [Input('edit_preprocessing_parameters', 'n_clicks'),
               Input('peak_picking_method', 'value')])
def toggle_peak_picking_method_parameters(n_clicks, value):
    """
    Dash callback to toggle which peak picking parameters are visible depending on the peak picking method selected in
    the preprocessing parameters modal window.

    :param n_clicks: Input signal if the edit_preprocessing_parameters button is clicked.
    :param value: Input signal to obtain the currently selected peak picking method.
    :return: List of dictionaries containing style template to show or hide parameters.
    """
    if value == 'locmax':
        return toggle_locmax_style()
    elif value == 'cwt':
        return toggle_cwt_style()


@app.callback([Output('peak_picking_deisotope_fragment_tolerance', 'style'),
               Output('peak_picking_deisotope_fragment_unit_ppm_label', 'style'),
               Output('peak_picking_deisotope_fragment_unit_ppm', 'style'),
               Output('peak_picking_deisotope_min_charge', 'style'),
               Output('peak_picking_deisotope_max_charge', 'style'),
               Output('peak_picking_deisotope_keep_only_deisotoped', 'style'),
               Output('peak_picking_deisotope_min_isopeaks', 'style'),
               Output('peak_picking_deisotope_max_isopeaks', 'style'),
               Output('peak_picking_deisotope_make_single_charged', 'style'),
               Output('peak_picking_deisotope_annotate_charge', 'style'),
               Output('peak_picking_deisotope_annotate_iso_peak_count', 'style'),
               Output('peak_picking_deisotope_use_decreasing_model', 'style'),
               Output('peak_picking_deisotope_start_intensity_check', 'style'),
               Output('peak_picking_deisotope_add_up_intensity', 'style')],
              [Input('edit_preprocessing_parameters', 'n_clicks'),
               Input('peak_picking_deisotope', 'value')])
def toggle_peak_picking_deisotope_parameters(n_clicks, value):
    """
    Dash callback to toggle whether deisotoping parameters are visible.

    :param n_clicks: Input signal if the edit_preprocessing_parameters button is clicked.
    :param value: Input signal to obtain the current status of whether deisotoping is enabled or disabled.
    :return: List of dictionaries containing style template to show or hide parameters.
    """
    if value:
        return toggle_deisotope_on_style()
    elif not value:
        return toggle_deisotope_off_style()


@app.callback([Output('spectrum_plot', 'figure'),
               Output('store_plot', 'data')],
              Input('trim_spectrum', 'n_clicks'),
              [State('store_preprocessing_params', 'data'),
               State('spectrum_id', 'value')])
def trim_spectrum_button(n_clicks, preprocessing_params, value):
    """
    Dash callback to apply spectrum trimming to the currently selected spectrum.

    :param n_clicks: Input signal if the trim_spectrum button is clicked.
    :param preprocessing_params: Input signal containing data from store_preprocessing_params.
    :param value: Input signal spectrum_id used as key in INDEXED_DATA.
    :return: Tuple of spectrum figure as a plotly.express.line plot and data store for plotly_resampler.
    """
    global INDEXED_DATA
    INDEXED_DATA[value].trim_spectrum(**preprocessing_params['TRIM_SPECTRUM'])
    fig = get_spectrum(INDEXED_DATA[value])
    cleanup_file_system_backend(FILE_SYSTEM_BACKEND)
    return fig, Serverside(fig)


@app.callback([Output('spectrum_plot', 'figure'),
               Output('store_plot', 'data')],
              Input('transform_intensity', 'n_clicks'),
              [State('store_preprocessing_params', 'data'),
               State('spectrum_id', 'value')])
def transform_intensity_button(n_clicks, preprocessing_params, value):
    """
    Dash callback to apply intensity transformation to the currently selected spectrum.

    :param n_clicks: Input signal if the transform_intensity button is clicked.
    :param preprocessing_params: Input signal containing data from store_preprocessing_params.
    :param value: Input signal spectrum_id used as key in INDEXED_DATA.
    :return: Tuple of spectrum figure as a plotly.express.line plot and data store for plotly_resampler.
    """
    global INDEXED_DATA
    INDEXED_DATA[value].transform_intensity(**preprocessing_params['TRANSFORM_INTENSITY'])
    fig = get_spectrum(INDEXED_DATA[value])
    cleanup_file_system_backend(FILE_SYSTEM_BACKEND)
    return fig, Serverside(fig)


@app.callback([Output('spectrum_plot', 'figure'),
               Output('store_plot', 'data')],
              Input('smooth_baseline', 'n_clicks'),
              [State('store_preprocessing_params', 'data'),
               State('spectrum_id', 'value')])
def smooth_baseline_button(n_clicks, preprocessing_params, value):
    """
    Dash callback to apply baseline smoothing to the currently selected spectrum.

    :param n_clicks: Input signal if the smooth_baseline button is clicked.
    :param preprocessing_params: Input signal containing data from store_preprocessing_params.
    :param value: Input signal spectrum_id used as key in INDEXED_DATA.
    :return: Tuple of spectrum figure as a plotly.express.line plot and data store for plotly_resampler.
    """
    global INDEXED_DATA
    INDEXED_DATA[value].smooth_baseline(**preprocessing_params['SMOOTH_BASELINE'])
    fig = get_spectrum(INDEXED_DATA[value])
    cleanup_file_system_backend(FILE_SYSTEM_BACKEND)
    return fig, Serverside(fig)


@app.callback([Output('spectrum_plot', 'figure'),
               Output('store_plot', 'data')],
              Input('remove_baseline', 'n_clicks'),
              [State('store_preprocessing_params', 'data'),
               State('spectrum_id', 'value')])
def remove_baseline_button(n_clicks, preprocessing_params, value):
    """
    Dash callback to apply baseline removal to the currently selected spectrum.

    :param n_clicks: Input signal if the remove_baseline button is clicked.
    :param preprocessing_params: Input signal containing data from store_preprocessing_params.
    :param value: Input signal spectrum_id used as key in INDEXED_DATA.
    :return: Tuple of spectrum figure as a plotly.express.line plot and data store for plotly_resampler.
    """
    global INDEXED_DATA
    INDEXED_DATA[value].remove_baseline(**preprocessing_params['REMOVE_BASELINE'],
                                        )
    fig = get_spectrum(INDEXED_DATA[value])
    cleanup_file_system_backend(FILE_SYSTEM_BACKEND)
    return fig, Serverside(fig)


@app.callback([Output('spectrum_plot', 'figure'),
               Output('store_plot', 'data')],
              Input('normalize_intensity', 'n_clicks'),
              [State('store_preprocessing_params', 'data'),
               State('spectrum_id', 'value')])
def normalize_intensity_button(n_clicks, preprocessing_params, value):
    """
    Dash callback to apply intensity normalization to the currently selected spectrum.

    :param n_clicks: Input signal if the normalize_intensity button is clicked.
    :param preprocessing_params: Input signal containing data from store_preprocessing_params.
    :param value: Input signal spectrum_id used as key in INDEXED_DATA.
    :return: Tuple of spectrum figure as a plotly.express.line plot and data store for plotly_resampler.
    """
    global INDEXED_DATA
    INDEXED_DATA[value].normalize_intensity(**preprocessing_params['NORMALIZE_INTENSITY'])
    fig = get_spectrum(INDEXED_DATA[value])
    cleanup_file_system_backend(FILE_SYSTEM_BACKEND)
    return fig, Serverside(fig)


@app.callback([Output('spectrum_plot', 'figure'),
               Output('store_plot', 'data')],
              Input('bin_spectrum', 'n_clicks'),
              [State('store_preprocessing_params', 'data'),
               State('spectrum_id', 'value')])
def bin_spectrum_button(n_clicks, preprocessing_params, value):
    """
    Dash callback to apply spectrum binning to the currently selected spectrum.

    :param n_clicks: Input signal if the bin_spectrum button is clicked.
    :param preprocessing_params: Input signal containing data from store_preprocessing_params.
    :param value: Input signal spectrum_id used as key in INDEXED_DATA.
    :return: Tuple of spectrum figure as a plotly.express.line plot and data store for plotly_resampler.
    """
    global INDEXED_DATA
    INDEXED_DATA[value].bin_spectrum(**preprocessing_params['BIN_SPECTRUM'])
    fig = get_spectrum(INDEXED_DATA[value])
    cleanup_file_system_backend(FILE_SYSTEM_BACKEND)
    return fig, Serverside(fig)


@app.callback([Output('spectrum_plot', 'figure'),
               Output('store_plot', 'data')],
              Input('peak_picking', 'n_clicks'),
              [State('store_preprocessing_params', 'data'),
               State('spectrum_id', 'value')])
def peak_picking_button(n_clicks, preprocessing_params, value):
    """
    Dash callback to apply peak picking to the currently selected spectrum.

    :param n_clicks: Input signal if the peak_picking button is clicked.
    :param preprocessing_params: Input signal containing data from store_preprocessing_params.
    :param value: Input signal spectrum_id used as key in INDEXED_DATA.
    :return: Tuple of spectrum figure as a plotly.express.line plot and data store for plotly_resampler.
    """
    global INDEXED_DATA
    if isinstance(INDEXED_DATA[value], PMP3DTdfSpectrum):
        INDEXED_DATA[value].peak_picking(**preprocessing_params['PEAK_PICKING_3D'])
        fig = get_peakmap(INDEXED_DATA[value], label_peaks=True)
        return fig, None
    else:
        INDEXED_DATA[value].peak_picking(**preprocessing_params['PEAK_PICKING'])
        fig = get_spectrum(INDEXED_DATA[value], label_peaks=True)
        cleanup_file_system_backend(FILE_SYSTEM_BACKEND)
        return fig, Serverside(fig)


@app.callback([Output('spectrum_plot', 'figure'),
               Output('store_plot', 'data')],
              Input('undo_peak_picking', 'n_clicks'),
              State('spectrum_id', 'value'))
def undo_peak_picking(n_clicks, value):
    """
    Dash callback to undo peak picking in the currently selected spectrum.

    :param n_clicks: Input signal if the undo_peak_picking button is clicked.
    :param value: Input signal spectrum_id used as key in INDEXED_DATA.
    :return: Tuple of spectrum figure as a plotly.express.line plot and data store for plotly_resampler.
    """
    global INDEXED_DATA
    if isinstance(INDEXED_DATA[value], PMP3DTdfSpectrum):
        INDEXED_DATA[value].peak_picked_mz_array = None
        INDEXED_DATA[value].peak_picked_intensity_array = None
        INDEXED_DATA[value].peak_picked_mobility_array = None
        INDEXED_DATA[value].peak_picking_indices = None
        del INDEXED_DATA[value].data_processing['peak picking']
        gc.collect()
        fig = get_peakmap(INDEXED_DATA[value])
        return fig, None
    else:
        INDEXED_DATA[value].peak_picked_mz_array = None
        INDEXED_DATA[value].peak_picked_intensity_array = None
        INDEXED_DATA[value].peak_picking_indices = None
        del INDEXED_DATA[value].data_processing['peak picking']
        gc.collect()
        fig = get_spectrum(INDEXED_DATA[value])
        cleanup_file_system_backend(FILE_SYSTEM_BACKEND)
        return fig, Serverside(fig)


@app.callback(Output('dummy', 'children'),
              Input('export_peak_list', 'n_clicks'),
              [State('store_preprocessing_params', 'data'),
               State('spectrum_id', 'value')])
def export_peak_list(n_clicks, preprocessing_params, value):
    """
    Dash callback to export the current peak list obtained from peak picking as a CSV file.

    :param n_clicks: Input signal if the export_peak_list button is clicked.
    :param preprocessing_params: Input signal containing data from store_preprocessing_params.
    :param value: Input signal spectrum_id used as key in INDEXED_DATA.
    :return: Empty list to dummy div.
    """
    global INDEXED_DATA
    if isinstance(INDEXED_DATA[value], PMP3DTdfSpectrum):
        if (INDEXED_DATA[value].peak_picked_mz_array is None and
                INDEXED_DATA[value].peak_picked_intensity_array is None and
                INDEXED_DATA[value].peak_picked_mobility_array is None):
            INDEXED_DATA[value].peak_picking(**preprocessing_params['PEAK_PICKING_3D'])
        spectrum_df = pd.DataFrame(data={'m/z': copy.deepcopy(INDEXED_DATA[value].peak_picked_mz_array),
                                         '1/K0': copy.deepcopy(INDEXED_DATA[value].peak_picked_mobility_array),
                                         'Intensity': copy.deepcopy(INDEXED_DATA[value].peak_picked_intensity_array)})
    else:
        if INDEXED_DATA[value].peak_picked_mz_array is None and INDEXED_DATA[value].peak_picked_intensity_array is None:
            INDEXED_DATA[value].peak_picking(**preprocessing_params['PEAK_PICKING'])
        spectrum_df = pd.DataFrame(data={'m/z': copy.deepcopy(INDEXED_DATA[value].peak_picked_mz_array),
                                         'Intensity': copy.deepcopy(INDEXED_DATA[value].peak_picked_intensity_array)})
    export_filetypes = [('Comma Separated Values', '*.csv')]
    if INDEXED_DATA[value].ms_level == 2:
        export_filetypes.append(('MASCOT Generic Format', '*.mgf'))
    main_tk_window = tkinter.Tk()
    main_tk_window.attributes('-topmost', True, '-alpha', 0)
    export_filename = asksaveasfilename(confirmoverwrite=True,
                                        filetypes=export_filetypes,
                                        defaultextension='csv')
    main_tk_window.destroy()
    if os.path.splitext(export_filename)[1] == '.csv':
        spectrum_df.to_csv(export_filename, index=False)
    elif os.path.splitext(export_filename)[1] == '.mgf':
        INDEXED_DATA[value].to_mgf(export_filename, peak_list=True)
    return []


@app.callback([Output('spectrum_plot', 'figure'),
               Output('store_plot', 'data')],
              Input('undo_preprocessing', 'n_clicks'),
              State('spectrum_id', 'value'))
def undo_preprocessing(n_clicks, value):
    """
    Dash callback to undo all preprocessing in the currently selected spectrum.

    :param n_clicks: Input signal if the undo_preprocessing button is clicked.
    :param value: Input signal spectrum_id used as key in INDEXED_DATA.
    :return: Tuple of spectrum figure as a plotly.express.line plot and data store for plotly_resampler.
    """
    global INDEXED_DATA
    if isinstance(INDEXED_DATA[value], PMP3DTdfSpectrum):
        INDEXED_DATA[value].undo_all_processing()
        fig = get_peakmap(INDEXED_DATA[value])
        return fig, None
    else:
        INDEXED_DATA[value].undo_all_processing()
        fig = get_spectrum(INDEXED_DATA[value])
        cleanup_file_system_backend(FILE_SYSTEM_BACKEND)
        return fig, Serverside(fig)


@app.callback(Output('about_modal', 'is_open'),
              [Input('about', 'n_clicks'),
               Input('about_close', 'n_clicks')],
              State('about_modal', 'is_open'))
def toggle_about_modal(n_clicks_about, n_clicks_about_close, is_open):
    """
    Dash callback to toggle the about modal window.

    :param n_clicks_about: Input signal if the about button is clicked.
    :param n_clicks_about_close: Input signal if the about_close button is clicked.
    :param is_open: State signal to determine whether the about_modal modal window is open.
    :return: Output signal to determine whether the about_modal modal window is open.
    """
    if n_clicks_about or n_clicks_about_close:
        return not is_open
    return is_open


@app.callback(Output('spectrum_plot', 'figure', allow_duplicate=True),
              Input('spectrum_plot', 'relayoutData'),
              State('store_plot', 'data'),
              prevent_initial_call=True,
              memoize=True)
def resample_spectrum(relayoutdata: dict, fig: FigureResampler):
    """
    Dash callback used for spectrum resampling to improve plotly figure performance.

    :param relayoutdata: Input signal with dictionary with spectrum_plot relayoutData.
    :param fig: State signal for data store for plotly_resampler.
    :return: Figure object used to update spectrum_plot figure.
    """
    if fig is None:
        return no_update
    return fig.construct_update_data_patch(relayoutdata)


if __name__ == '__main__':
    app.run_server(debug=False)
