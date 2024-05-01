import os
import copy
import configparser
import numpy as np
import pandas as pd
import plotly.express as px
from plotly_resampler import FigureResampler
from pymaldiviz.tmpdir import FILE_SYSTEM_BACKEND


SHOWN = {'margin': '10px',
         'display': 'flex'}
HIDDEN = {'margin': '10px',
          'display': 'none'}


def get_preprocessing_params():
    """
    Parse preprocessing parameters from the configuration file provided with pyMALDIproc.

    :return: Nest dictionaries containing preprocessing parameters for each preprocessing step.
    :rtype: dict
    """
    config = configparser.ConfigParser()
    config.read(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'etc', 'preprocessing.cfg'))
    trim_spectrum_params = {'lower_mass_range': int(config['trim_spectrum']['lower_mass_range']),
                            'upper_mass_range': int(config['trim_spectrum']['upper_mass_range'])}
    transform_intensity_params = {'method': config['transform_intensity']['method']}
    smooth_baseline_params = {'method': config['smooth_baseline']['method'],
                              'window_length': int(config['smooth_baseline']['window_length']),
                              'polyorder': int(config['smooth_baseline']['polyorder']),
                              'delta_mz': float(config['smooth_baseline']['delta_mz']),
                              'diff_thresh': float(config['smooth_baseline']['diff_thresh'])}
    remove_baseline_params = {'method': config['remove_baseline']['method'],
                              'min_half_window': int(config['remove_baseline']['min_half_window']),
                              'max_half_window': int(config['remove_baseline']['max_half_window']),
                              'decreasing': config['remove_baseline'].getboolean('decreasing'),
                              'smooth_half_window': None,
                              'filter_order': int(config['remove_baseline']['filter_order']),
                              'sigma': None,
                              'increment': int(config['remove_baseline']['increment']),
                              'max_hits': int(config['remove_baseline']['max_hits']),
                              'window_tol': float(config['remove_baseline']['window_tol']),
                              'lambda_': int(config['remove_baseline']['lambda_']),
                              'porder': int(config['remove_baseline']['porder']),
                              'repetition': None,
                              'degree': int(config['remove_baseline']['degree']),
                              'gradient': float(config['remove_baseline']['gradient'])}
    normalize_intensity_params = {'method': config['normalize_intensity']['method']}
    bin_spectrum_params = {'n_bins': int(config['bin_spectrum']['n_bins']),
                           'lower_mass_range': int(config['bin_spectrum']['lower_mass_range']),
                           'upper_mass_range': int(config['bin_spectrum']['upper_mass_range'])}
    peak_picking_params = {'method': config['peak_picking']['method'],
                           'snr': int(config['peak_picking']['snr']),
                           'widths': None,
                           'deisotope': config['peak_picking'].getboolean('deisotope'),
                           'fragment_tolerance': float(config['peak_picking']['fragment_tolerance']),
                           'fragment_unit_ppm': config['peak_picking'].getboolean('fragment_unit_ppm'),
                           'min_charge': int(config['peak_picking']['min_charge']),
                           'max_charge': int(config['peak_picking']['max_charge']),
                           'keep_only_deisotoped': config['peak_picking'].getboolean('keep_only_deisotoped'),
                           'min_isopeaks': int(config['peak_picking']['min_isopeaks']),
                           'max_isopeaks': int(config['peak_picking']['max_isopeaks']),
                           'make_single_charged': config['peak_picking'].getboolean('make_single_charged'),
                           'annotate_charge': config['peak_picking'].getboolean('annotate_charge'),
                           'annotate_iso_peak_count': config['peak_picking'].getboolean('annotate_iso_peak_count'),
                           'use_decreasing_model': config['peak_picking'].getboolean('use_decreasing_model'),
                           'start_intensity_check': int(config['peak_picking']['start_intensity_check']),
                           'add_up_intensity': config['peak_picking'].getboolean('add_up_intensity')}
    if config['remove_baseline']['smooth_half_window'] != 'None':
        remove_baseline_params['smooth_half_window'] = int(config['remove_baseline']['smooth_half_window'])
    if config['remove_baseline']['sigma'] != 'None':
        remove_baseline_params['sigma'] = float(config['rmeove_baseline']['sigma'])
    if config['remove_baseline']['repetition'] != 'None':
        remove_baseline_params['repetition'] = int(config['remove_baseline']['repetition'])
    if config['peak_picking']['widths'] != 'None':
        peak_picking_params['widths'] = int(config['peak_picking']['widths'])
    return {'TRIM_SPECTRUM': trim_spectrum_params,
            'TRANSFORM_INTENSITY': transform_intensity_params,
            'SMOOTH_BASELINE': smooth_baseline_params,
            'REMOVE_BASELINE': remove_baseline_params,
            'NORMALIZE_INTENSITY': normalize_intensity_params,
            'BIN_SPECTRUM': bin_spectrum_params,
            'PEAK_PICKING': peak_picking_params}


def get_spectrum(spectrum, label_peaks=False):
    """
    Plot the spectrum to a plotly.express.line plot wrapped by plotly_resampler.FigureResampler.

    :param spectrum: Spectrum object whose data is used to generate the figure.
    :type spectrum: pymaldiproc.classes.OpenMALDISpectrum|pymaldiproc.classes.PMPTsfSpectrum|pymaldiproc.classes.PMP2DTdfSpectrum
    :param label_peaks: Whether to label the peak based on peak picking that has been performed.
    :type label_peaks: bool
    :return: Plotly figure containing mass spectrum.
    """
    spectrum_df = pd.DataFrame({'m/z': copy.deepcopy(spectrum.preprocessed_mz_array),
                                'Intensity': copy.deepcopy(spectrum.preprocessed_intensity_array)})

    if label_peaks:
        labels = copy.deepcopy(np.round(copy.deepcopy(spectrum.preprocessed_mz_array), decimals=4).astype(str))
        mask = np.ones(labels.size, dtype=bool)
        mask[spectrum.peak_picking_indices] = False
        labels[mask] = ''
        fig = FigureResampler(px.line(data_frame=spectrum_df,
                                      x='m/z',
                                      y='Intensity',
                                      hover_data={'m/z': ':.4f',
                                                  'Intensity': ':.1f'},
                                      text=labels))
        fig.update_traces(textposition='top center')
    else:
        fig = FigureResampler(px.line(data_frame=spectrum_df,
                                      x='m/z',
                                      y='Intensity',
                                      hover_data={'m/z': ':.4f',
                                                  'Intensity': ':.1f'}))
    fig.update_layout(xaxis_tickformat='d',
                      yaxis_tickformat='~e')

    return fig


def cleanup_file_system_backend():
    """
    Clean up temp files generated by plotly_resampler.
    """
    for filename in os.listdir(FILE_SYSTEM_BACKEND):
        os.remove(os.path.join(FILE_SYSTEM_BACKEND, filename))


def toggle_savitzky_golay_style():
    """
    Toggle the visibility of baseline smoothing parameters if the 'Savitzky Golay' method is selected.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN)]


def toggle_apodization_style():
    """
    Toggle the visibility of baseline smoothing parameters if the 'apodization' method is selected.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN)]


def toggle_rebin_style():
    """
    Toggle the visibility of baseline smoothing parameters if the 'rebin' method is selected.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN)]


def toggle_fast_change_style():
    """
    Toggle the visibility of baseline smoothing parameters if the 'fast change' method is selected.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(SHOWN)]


def toggle_smoothing_median_style():
    """
    Toggle the visibility of baseline smoothing parameters if the 'median' method is selected.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN)]


def toggle_snip_style():
    """
    Toggle the visibility of baseline removal parameters if the 'SNIP' method is selected.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(HIDDEN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN)]


def toggle_tophat_style():
    """
    Toggle the visibility of baseline removal parameters if the 'TopHat' method is selected.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN)]


def toggle_removal_median_style():
    """
    Toggle the visibility of baseline removal parameters if the 'median' method is selected.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(HIDDEN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN)]


def toggle_zhangfit_style():
    """
    Toggle the visibility of baseline removal parameters if the 'ZhangFit' method is selected.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN)]


def toggle_modpoly_style():
    """
    Toggle the visibility of baseline removal parameters if the 'ModPoly' method is selected.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN)]


def toggle_imodpoly_style():
    """
    Toggle the visibility of baseline removal parameters if the 'IModPoly' method is selected.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN)]


def toggle_locmax_style():
    """
    Toggle the visibility of peak picking parameters if the 'locmax' method is selected.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN)]


def toggle_cwt_style():
    """
    Toggle the visibility of peak picking parameters if the 'CWT' method is selected.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN)]


def toggle_deisotope_on_style():
    """
    Toggle the visibility of deisotoping parameters if deisotoping is enabled.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN)]


def toggle_deisotope_off_style():
    """
    Toggle the visibility of deisotoping parameters if deisotoping is disabled.

    :return: List of dictionaries containing style template to show or hide parameters.
    """
    return [copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN)]