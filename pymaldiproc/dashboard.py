import os
import gc
import configparser
from pymaldiproc.data_import import import_mzml, import_timstof_raw_data
from pymaldiproc.layout import *
from dash import State, callback_context, no_update
from dash_extensions.enrich import Input, Output, DashProxy, MultiplexerTransform, Serverside, ServersideOutputTransform
import dash_bootstrap_components as dbc
import tkinter
from tkinter.filedialog import askopenfilenames, askdirectory, asksaveasfilename

# will be a dictionary of MALDISpectrum objects used for the spectrum plot
INDEXED_DATA = {}
# default processing parameters from config file
config = configparser.ConfigParser()
config.read(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'etc', 'preprocessing.cfg'))
TRIM_SPECTRUM_PARAMS = {'lower_mass_range': int(config['trim_spectrum']['lower_mass_range']),
                        'upper_mass_range': int(config['trim_spectrum']['upper_mass_range'])}
TRANSFORM_INTENSITY_PARAMS = {'method': config['transform_intensity']['method']}
SMOOTH_BASELINE_PARAMS = {'method': config['smooth_baseline']['method'],
                          'window_length': int(config['smooth_baseline']['window_length']),
                          'polyorder': int(config['smooth_baseline']['polyorder']),
                          'delta_mz': float(config['smooth_baseline']['delta_mz']),
                          'diff_thresh': float(config['smooth_baseline']['diff_thresh'])}
REMOVE_BASELINE_PARAMS = {'method': config['remove_baseline']['method'],
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
NORMALIZE_INTENSITY_PARAMS = {'method': config['normalize_intensity']['method']}
BIN_SPECTRUM_PARAMS = {'n_bins': int(config['bin_spectrum']['n_bins']),
                       'lower_mass_range': int(config['bin_spectrum']['lower_mass_range']),
                       'upper_mass_range': int(config['bin_spectrum']['upper_mass_range'])}
PEAK_PICKING_PARAMS = {'method': config['peak_picking']['method'],
                       'snr': int(config['peak_picking']['snr']),
                       'widths': None}
if config['remove_baseline']['smooth_half_window'] != 'None':
    REMOVE_BASELINE_PARAMS['smooth_half_window'] = int(config['remove_baseline']['smooth_half_window'])
if config['remove_baseline']['sigma'] != 'None':
    REMOVE_BASELINE_PARAMS['sigma'] = float(config['rmeove_baseline']['sigma'])
if config['remove_baseline']['repetition'] != 'None':
    REMOVE_BASELINE_PARAMS['repetition'] = int(config['remove_baseline']['repetition'])
if config['peak_picking']['widths'] != 'None':
    PEAK_PICKING_PARAMS['widths'] = int(config['peak_picking']['widths'])
PREPROCESSING_PARAMS = {'TRIM_SPECTRUM': TRIM_SPECTRUM_PARAMS,
                        'TRANSFORM_INTENSITY': TRANSFORM_INTENSITY_PARAMS,
                        'SMOOTH_BASELINE': SMOOTH_BASELINE_PARAMS,
                        'REMOVE_BASELINE': REMOVE_BASELINE_PARAMS,
                        'NORMALIZE_INTENSITY': NORMALIZE_INTENSITY_PARAMS,
                        'BIN_SPECTRUM': BIN_SPECTRUM_PARAMS,
                        'PEAK_PICKING': PEAK_PICKING_PARAMS}

# Use DashProxy instead of Dash to allow for multiple callbacks to the same plot
app = DashProxy(prevent_initial_callbacks=True,
                transforms=[MultiplexerTransform(), ServersideOutputTransform()],
                external_stylesheets=[dbc.themes.BOOTSTRAP])
app.layout = get_dashboard_layout(PREPROCESSING_PARAMS)


@app.callback(Output('dropdown', 'children'),
              [Input('upload_mzml', 'n_clicks'),
               Input('upload_d', 'n_clicks')])
def upload_data(n_clicks_mzml, n_clicks_d):
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
            data = import_timstof_raw_data(dirname, mode='profile')
            for spectrum in data:
                INDEXED_DATA[spectrum.spectrum_id] = spectrum
    return get_dropdown_layout(INDEXED_DATA)


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('spectrum_id', 'value'))
def plot_spectrum(value):
    global INDEXED_DATA
    fig = get_spectrum(INDEXED_DATA[value])
    for filename in os.listdir(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend')):
        os.remove(os.path.join(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend'),
                               filename))
    return [get_spectrum_plot_layout(fig)], Serverside(fig)


@app.callback(Output('edit_processing_parameters_modal', 'is_open'),
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
               Input('peak_picking_widths_value', 'value')],
              State('edit_processing_parameters_modal', 'is_open'))
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
                                               is_open):
    global PREPROCESSING_PARAMS
    changed_id = [i['prop_id'] for i in callback_context.triggered][0]
    if (changed_id == 'edit_preprocessing_parameters.n_clicks' or
            changed_id == 'edit_processing_parameters_save.n_clicks' or
            changed_id == 'edit_processing_parameters_cancel.n_clicks'):
        if changed_id == 'edit_processing_parameters_save.n_clicks':
            PREPROCESSING_PARAMS['TRIM_SPECTRUM']['lower_mass_range'] = trim_spectrum_lower_mass_range
            PREPROCESSING_PARAMS['TRIM_SPECTRUM']['upper_mass_range'] = trim_spectrum_upper_mass_range
            PREPROCESSING_PARAMS['TRANSFORM_INTENSITY']['method'] = transform_intensity_method
            PREPROCESSING_PARAMS['SMOOTH_BASELINE']['method'] = smooth_baseline_method
            PREPROCESSING_PARAMS['SMOOTH_BASELINE']['window_length'] = smooth_baseline_window_length
            PREPROCESSING_PARAMS['SMOOTH_BASELINE']['polyorder'] = smooth_baseline_polyorder
            PREPROCESSING_PARAMS['SMOOTH_BASELINE']['delta_mz'] = smooth_baseline_delta_mz
            PREPROCESSING_PARAMS['SMOOTH_BASELINE']['diff_thresh'] = smooth_baseline_diff_thresh
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['method'] = remove_baseline_method
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['min_half_window'] = remove_baseline_min_half_window
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['max_half_window'] = remove_baseline_max_half_window
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['decreasing'] = remove_baseline_decreasing
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['smooth_half_window'] = remove_baseline_smooth_half_window
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['filter_order'] = remove_baseline_filter_order
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['sigma'] = remove_baseline_sigma
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['increment'] = remove_baseline_increment
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['max_hits'] = remove_baseline_max_hits
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['window_tol'] = remove_baseline_window_tol
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['lambda_'] = remove_baseline_lambda_
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['porder'] = remove_baseline_porder
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['repetition'] = remove_baseline_repetition
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['degree'] = remove_baseline_degree
            PREPROCESSING_PARAMS['REMOVE_BASELINE']['gradient'] = remove_baseline_gradient
            PREPROCESSING_PARAMS['NORMALIZE_INTENSITY']['method'] = normalize_intensity_method
            PREPROCESSING_PARAMS['BIN_SPECTRUM']['n_bins'] = bin_spectrum_n_bins
            PREPROCESSING_PARAMS['BIN_SPECTRUM']['lower_mass_range'] = bin_spectrum_lower_mass_range
            PREPROCESSING_PARAMS['BIN_SPECTRUM']['upper_mass_range'] = bin_spectrum_upper_mass_range
            PREPROCESSING_PARAMS['PEAK_PICKING']['method'] = peak_picking_method
            PREPROCESSING_PARAMS['PEAK_PICKING']['snr'] = peak_picking_snr
            PREPROCESSING_PARAMS['PEAK_PICKING']['widths'] = peak_picking_widths
        return not is_open
    return is_open


@app.callback(Output('edit_processing_parameters_modal_saved', 'is_open'),
              [Input('edit_processing_parameters_save', 'n_clicks'),
               Input('edit_processing_parameters_modal_saved_close', 'n_clicks')],
              State('edit_processing_parameters_modal_saved', 'is_open'))
def toggle_edit_processing_parameters_saved_modal(n_clicks_save, n_clicks_close, is_open):
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
    global PREPROCESSING_PARAMS
    shown = {'margin': '10px',
             'display': 'flex'}
    hidden = {'margin': '10px',
              'display': 'none'}
    if value == 'SavitzkyGolay':
        return [copy.deepcopy(shown),
                copy.deepcopy(shown),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden)]
    elif value == 'apodization':
        return [copy.deepcopy(shown),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden)]
    elif value == 'rebin':
        return [copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(shown),
                copy.deepcopy(hidden)]
    elif value == 'fast_change':
        return [copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(shown)]
    elif value == 'median':
        return [copy.deepcopy(shown),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden)]


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
    global PREPROCESSING_PARAMS
    shown = {'margin': '10px',
             'display': 'flex'}
    hidden = {'margin': '10px',
              'display': 'none'}
    if value == 'SNIP':
        return [copy.deepcopy(hidden),
                copy.deepcopy(shown),
                copy.deepcopy(shown),
                copy.deepcopy(shown),
                copy.deepcopy(shown),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden)]
    elif value == 'TopHat':
        return [copy.deepcopy(shown),
                copy.deepcopy(shown),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(shown),
                copy.deepcopy(shown),
                copy.deepcopy(shown),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden)]
    elif value == 'Median':
        return [copy.deepcopy(hidden),
                copy.deepcopy(shown),
                copy.deepcopy(hidden),
                copy.deepcopy(shown),
                copy.deepcopy(hidden),
                copy.deepcopy(shown),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden)]
    elif value == 'ZhangFit':
        return [copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(shown),
                copy.deepcopy(shown),
                copy.deepcopy(shown),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden)]
    elif value == 'ModPoly':
        return [copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(shown),
                copy.deepcopy(shown),
                copy.deepcopy(shown)]
    elif value == 'IModPoly':
        return [copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(hidden),
                copy.deepcopy(shown),
                copy.deepcopy(shown),
                copy.deepcopy(shown)]


@app.callback([Output('peak_picking_snr', 'style'),
               Output('peak_picking_widths', 'Style')],
              [Input('edit_preprocessing_parameters', 'n_clicks'),
               Input('peak_picking_method', 'value')])
def toggle_peak_picking_method_parameters(n_clicks, value):
    global PREPROCESSING_PARAMS
    shown = {'margin': '10px',
             'display': 'flex'}
    hidden = {'margin': '10px',
              'display': 'none'}
    if value == 'locmax':
        return [copy.deepcopy(shown),
                copy.deepcopy(hidden)]
    elif value == 'cwt':
        return [copy.deepcopy(shown),
                copy.deepcopy(shown)]


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('trim_spectrum', 'n_clicks'),
              State('spectrum_id', 'value'))
def trim_spectrum_button(n_clicks, value):
    global INDEXED_DATA
    global PREPROCESSING_PARAMS
    INDEXED_DATA[value].trim_spectrum(**PREPROCESSING_PARAMS['TRIM_SPECTRUM'])
    fig = get_spectrum(INDEXED_DATA[value])
    for filename in os.listdir(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend')):
        os.remove(os.path.join(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend'),
                               filename))
    return [get_spectrum_plot_layout(fig)], Serverside(fig)


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('transform_intensity', 'n_clicks'),
              State('spectrum_id', 'value'))
def transform_intensity_button(n_clicks, value):
    global INDEXED_DATA
    global PREPROCESSING_PARAMS
    INDEXED_DATA[value].transform_intensity(**PREPROCESSING_PARAMS['TRANSFORM_INTENSITY'])
    fig = get_spectrum(INDEXED_DATA[value])
    for filename in os.listdir(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend')):
        os.remove(os.path.join(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend'),
                               filename))
    return [get_spectrum_plot_layout(fig)], Serverside(fig)


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('smooth_baseline', 'n_clicks'),
              State('spectrum_id', 'value'))
def smooth_baseline_button(n_clicks, value):
    global INDEXED_DATA
    global PREPROCESSING_PARAMS
    INDEXED_DATA[value].smooth_baseline(**PREPROCESSING_PARAMS['SMOOTH_BASELINE'])
    fig = get_spectrum(INDEXED_DATA[value])
    for filename in os.listdir(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend')):
        os.remove(os.path.join(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend'),
                               filename))
    return [get_spectrum_plot_layout(fig)], Serverside(fig)


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('remove_baseline', 'n_clicks'),
              State('spectrum_id', 'value'))
def remove_baseline_button(n_clicks, value):
    global INDEXED_DATA
    global PREPROCESSING_PARAMS
    INDEXED_DATA[value].remove_baseline(**PREPROCESSING_PARAMS['REMOVE_BASELINE'],
                                        )
    fig = get_spectrum(INDEXED_DATA[value])
    for filename in os.listdir(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend')):
        os.remove(os.path.join(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend'),
                               filename))
    return [get_spectrum_plot_layout(fig)], Serverside(fig)


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('normalize_intensity', 'n_clicks'),
              State('spectrum_id', 'value'))
def normalize_intensity_button(n_clicks, value):
    global INDEXED_DATA
    global PREPROCESSING_PARAMS
    INDEXED_DATA[value].normalize_intensity(**PREPROCESSING_PARAMS['NORMALIZE_INTENSITY'])
    fig = get_spectrum(INDEXED_DATA[value])
    for filename in os.listdir(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend')):
        os.remove(os.path.join(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend'),
                               filename))
    return [get_spectrum_plot_layout(fig)], Serverside(fig)


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('bin_spectrum', 'n_clicks'),
              State('spectrum_id', 'value'))
def bin_spectrum_button(n_clicks, value):
    global INDEXED_DATA
    global PREPROCESSING_PARAMS
    INDEXED_DATA[value].bin_spectrum(**PREPROCESSING_PARAMS['BIN_SPECTRUM'])
    fig = get_spectrum(INDEXED_DATA[value])
    for filename in os.listdir(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend')):
        os.remove(os.path.join(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend'),
                               filename))
    return [get_spectrum_plot_layout(fig)], Serverside(fig)


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('peak_picking', 'n_clicks'),
              State('spectrum_id', 'value'))
def peak_picking_button(n_clicks, value):
    global INDEXED_DATA
    global PREPROCESSING_PARAMS
    INDEXED_DATA[value].peak_picking(**PREPROCESSING_PARAMS['PEAK_PICKING'])
    fig = get_spectrum(INDEXED_DATA[value], label_peaks=True)
    for filename in os.listdir(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend')):
        os.remove(os.path.join(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend'),
                               filename))
    return [get_spectrum_plot_layout(fig)], Serverside(fig)


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('undo_peak_picking', 'n_clicks'),
              State('spectrum_id', 'value'))
def undo_peak_picking(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value].peak_picked_mz_array = None
    INDEXED_DATA[value].peak_picked_intensity_array = None
    INDEXED_DATA[value].peak_picking_indices = None
    del INDEXED_DATA[value].data_processing['peak picking']
    gc.collect()
    fig = get_spectrum(INDEXED_DATA[value])
    for filename in os.listdir(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend')):
        os.remove(os.path.join(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend'),
                               filename))
    return [get_spectrum_plot_layout(fig)], Serverside(fig)


@app.callback(Output('dummy', 'children'),
              Input('export_peak_list', 'n_clicks'),
              State('spectrum_id', 'value'))
def export_peak_list(n_clicks, value):
    global INDEXED_DATA
    if INDEXED_DATA[value].peak_picked_mz_array is None and INDEXED_DATA[value].peak_picked_intensity_array is None:
        INDEXED_DATA[value].peak_picking()
    spectrum_df = pd.DataFrame(data={'m/z': copy.deepcopy(INDEXED_DATA[value].peak_picked_mz_array),
                                     'Intensity': copy.deepcopy(INDEXED_DATA[value].peak_picked_intensity_array)})
    main_tk_window = tkinter.Tk()
    main_tk_window.attributes('-topmost', True, '-alpha', 0)
    csv_filename = asksaveasfilename(confirmoverwrite=True,
                                     filetypes=[('Comma Separated Values', '*.csv')],
                                     defaultextension='csv')
    main_tk_window.destroy()
    spectrum_df.to_csv(csv_filename, index=False)
    return []


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('undo_preprocessing', 'n_clicks'),
              State('spectrum_id', 'value'))
def undo_preprocessing(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value].undo_all_processing()
    fig = get_spectrum(INDEXED_DATA[value])
    for filename in os.listdir(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend')):
        os.remove(os.path.join(os.path.join(os.path.split(os.path.dirname(__file__))[0], 'file_system_backend'),
                               filename))
    return [get_spectrum_plot_layout(fig)], Serverside(fig)


@app.callback(Output('spectrum_plot', 'figure', allow_duplicate=True),
              Input('spectrum_plot', 'relayoutData'),
              State('store_plot', 'data'),
              prevent_initial_call=True,
              memoize=True)
def resample_spectrum(relayoutdata: dict, fig: FigureResampler):
    if fig is None:
        return no_update
    return fig.construct_update_data_patch(relayoutdata)


if __name__ == '__main__':
    app.run_server(debug=False)
