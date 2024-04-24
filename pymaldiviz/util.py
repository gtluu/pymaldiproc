import os
import copy
import configparser
import numpy as np
import pandas as pd
import plotly.express as px
from plotly_resampler import FigureResampler


SHOWN = {'margin': '10px',
         'display': 'flex'}
HIDDEN = {'margin': '10px',
          'display': 'none'}


def get_preprocessing_params():
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
    for filename in os.listdir(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend')):
        os.remove(os.path.join(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend'),
                               filename))


def toggle_savitzky_golay_style():
    return [copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN)]


def toggle_apodization_style():
    return [copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN)]


def toggle_rebin_style():
    return [copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN)]


def toggle_fast_change_style():
    return [copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(SHOWN)]


def toggle_smoothing_median_style():
    return [copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN),
            copy.deepcopy(HIDDEN)]


def toggle_snip_style():
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
    return [copy.deepcopy(SHOWN),
            copy.deepcopy(HIDDEN)]


def toggle_cwt_style():
    return [copy.deepcopy(SHOWN),
            copy.deepcopy(SHOWN)]


def toggle_deisotope_on_style():
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