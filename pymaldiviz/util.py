import os
import copy
import configparser
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly_resampler import FigureResampler


SHOWN = {'margin': '10px',
         'display': 'flex'}
HIDDEN = {'margin': '10px',
          'display': 'none'}


def get_preprocessing_params(config_file=''):
    """
    Parse preprocessing parameters from the configuration file provided with pyMALDIproc.

    :param config_file: Configuration file containing default preprocessing parameters.
    :type config_file: str
    :return: Nest dictionaries containing preprocessing parameters for each preprocessing step.
    :rtype: dict
    """
    config = configparser.ConfigParser()
    if config_file == '':
        config.read(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'etc', 'preprocessing.cfg'))
    else:
        config.read(config_file)
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
    peak_picking_3d_params = {'min_distance': int(config['peak_picking_3d']['min_distance']),
                              'noise': None,
                              'snr': int(config['peak_picking_3d']['snr']),
                              'exclude_border': int(config['peak_picking_3d']['exclude_border'])}
    if config['remove_baseline']['smooth_half_window'] != 'None':
        remove_baseline_params['smooth_half_window'] = int(config['remove_baseline']['smooth_half_window'])
    if config['remove_baseline']['sigma'] != 'None':
        remove_baseline_params['sigma'] = float(config['rmeove_baseline']['sigma'])
    if config['remove_baseline']['repetition'] != 'None':
        remove_baseline_params['repetition'] = int(config['remove_baseline']['repetition'])
    if config['peak_picking']['widths'] != 'None':
        peak_picking_params['widths'] = int(config['peak_picking']['widths'])
    if config['peak_picking_3d']['noise'] != 'None':
        peak_picking_3d_params['noise'] = int(config['peak_picking_3d']['noise'])
    return {'TRIM_SPECTRUM': trim_spectrum_params,
            'TRANSFORM_INTENSITY': transform_intensity_params,
            'SMOOTH_BASELINE': smooth_baseline_params,
            'REMOVE_BASELINE': remove_baseline_params,
            'NORMALIZE_INTENSITY': normalize_intensity_params,
            'BIN_SPECTRUM': bin_spectrum_params,
            'PEAK_PICKING': peak_picking_params,
            'PEAK_PICKING_3D': peak_picking_3d_params}


def blank_figure():
    """
    Obtain a blank figure wrapped by plotly_resampler.FigureResampler to be used as a placeholder.

    :return: Blank figure.
    """
    fig = FigureResampler(go.Figure(go.Scatter(x=[], y=[])))
    fig.update_layout(template=None)
    fig.update_xaxes(showgrid=False, showticklabels=False, zeroline=False)
    fig.update_yaxes(showgrid=False, showticklabels=False, zeroline=False)
    return fig


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
        fig.update_traces(marker=dict(color='rgba(0,0,0,0)', size=1))
    else:
        fig = FigureResampler(px.line(data_frame=spectrum_df,
                                      x='m/z',
                                      y='Intensity',
                                      hover_data={'m/z': ':.4f',
                                                  'Intensity': ':.1f'}))
    fig.update_layout(xaxis_tickformat='d',
                      yaxis_tickformat='~e')

    return fig


def get_mirror_spectrum(spectrum_1, spectrum_2, label_peaks=False):
    """
    Plot the two spectra in a mirror plot to a plotly.express.line plot wrapped by plotly_resampler.FigureResampler.

    :param spectrum_1: Spectrum object whose data is used to generate the figure.
    :type spectrum_1: pymaldiproc.classes.OpenMALDISpectrum|pymaldiproc.classes.PMPTsfSpectrum|pymaldiproc.classes.PMP2DTdfSpectrum
    :param spectrum_2: Spectrum object whose data is used to generate the figure.
    :type spectrum_2: pymaldiproc.classes.OpenMALDISpectrum|pymaldiproc.classes.PMPTsfSpectrum|pymaldiproc.classes.PMP2DTdfSpectrum
    :param label_peaks: Whether to label the peak based on peak picking that has been performed.
    :type label_peaks: bool
    :return: Plotly figure containing mass spectrum.
    """
    spectrum_df_1 = pd.DataFrame({'m/z': copy.deepcopy(spectrum_1.preprocessed_mz_array),
                                  'Intensity': copy.deepcopy(spectrum_1.preprocessed_intensity_array)})
    spectrum_df_2 = pd.DataFrame({'m/z': copy.deepcopy(spectrum_2.preprocessed_mz_array),
                                  'Intensity': copy.deepcopy(spectrum_2.preprocessed_intensity_array)})

    def get_labels(spectrum):
        labels = copy.deepcopy(np.round(copy.deepcopy(spectrum.preprocessed_mz_array), decimals=4).astype(str))
        mask = np.ones(labels.size, dtype=bool)
        mask[spectrum.peak_picking_indices] = False
        labels[mask] = ''
        return labels

    if label_peaks:
        labels_1 = get_labels(spectrum_1)
        labels_2 = get_labels(spectrum_2)
        fig_1 = px.line(data_frame=spectrum_df_1,
                        x='m/z',
                        y='Intensity',
                        hover_data={'m/z': ':.4f',
                                    'Intensity': ':.1f'},
                        text=labels_1)
        fig_1.update_traces(textposition='top center')
        fig_1.update_traces(marker=dict(color='rgba(0,0,0,0)', size=1))
        fig_2 = px.line(data_frame=spectrum_df_2,
                        x='m/z',
                        y='Intensity',
                        hover_data={'m/z': ':.4f',
                                    'Intensity': ':.1f'},
                        text=labels_2)
        fig_2.update_traces(textposition='top center')
        fig_2.update_traces(marker=dict(color='rgba(0,0,0,0)', size=1))
    else:
        fig_1 = px.line(data_frame=spectrum_df_1,
                        x='m/z',
                        y='Intensity',
                        hover_data={'m/z': ':.4f',
                                    'Intensity': ':.1f'})
        fig_2 = px.line(data_frame=spectrum_df_2,
                        x='m/z',
                        y='Intensity',
                        hover_data={'m/z': ':.4f',
                                    'Intensity': ':.1f'})

    fig = FigureResampler(go.Figure(data=fig_1.data+fig_2.data))
    fig.update_layout(xaxis_tickformat='d',
                      yaxis_tickformat='~e')

    return fig


def get_peakmap(spectrum, use_log_intensity=True, label_peaks=False):
    spectrum_df = pd.DataFrame({'m/z': copy.deepcopy(spectrum.preprocessed_mz_array),
                                '1/K0': copy.deepcopy(spectrum.preprocessed_mobility_array),
                                'Intensity': copy.deepcopy(spectrum.preprocessed_intensity_array)})

    fig = spectrum_df.plot(x='m/z', y='1/K0', z='Intensity',
                           kind='peakmap',
                           add_marginals=True,
                           z_log_scale=use_log_intensity,
                           xlabel='m/z', ylabel='1/K0',
                           aggregation_method='sum',
                           num_x_bins=int((max(spectrum_df['m/z'].values) - min(spectrum_df['m/z'].values)) / 0.01),
                           num_y_bins=int((max(spectrum_df['1/K0'].values) - min(spectrum_df['1/K0'].values)) / 0.001),
                           width=1000, height=1000,
                           backend='ms_plotly',
                           show_plot=False).fig
    fig.update_layout(dragmode='zoom')
    fig.update_layout(showlegend=False)
    fig['layout']['annotations'] = []
    if label_peaks:
        for index, row in spectrum.get_feature_list().iterrows():
            fig.add_annotation(text=f"m/z {round(row['m/z'], 4)}<br>1/K0 {round(row['1/K0'], 3)}",
                               x=row['m/z'], y=row['1/K0'],
                               xref='x4', yref='y4',
                               showarrow=False,
                               font={'size': 16, 'color': 'Black'})
    fig.data[0].marker.colorscale = 'balance'
    fig['layout']['xaxis']['domain'] = [0.0, 0.25]
    fig['layout']['xaxis2']['domain'] = [0.25, 1.0]
    fig['layout']['xaxis3']['domain'] = [0.0, 0.25]
    fig['layout']['xaxis4']['domain'] = [0.25, 1.0]
    fig['layout']['yaxis']['domain'] = [0.75, 1.0]
    fig['layout']['yaxis2']['domain'] = [0.75, 1.0]
    fig['layout']['yaxis3']['domain'] = [0.75, 1.0]
    fig['layout']['yaxis4']['domain'] = [0.0, 0.75]
    fig['layout']['yaxis5']['domain'] = [0.0, 0.75]
    fig['layout']['width'] = np.inf
    fig['layout']['height'] = np.inf

    return fig


def cleanup_file_system_backend(file_system_backend):
    """
    Clean up temp files generated by plotly_resampler.
    """
    for filename in os.listdir(file_system_backend):
        os.remove(os.path.join(file_system_backend, filename))


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