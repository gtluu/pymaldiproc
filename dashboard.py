from pymaldiproc.classes import *
from pymaldiproc.data_import import *
from pymaldiproc.preprocessing import *
from pymaldiproc.layout import *
import plotly.express as px
from dash import Dash, dcc, html, State, callback_context, no_update
from dash_extensions.enrich import Input, Output, DashProxy, MultiplexerTransform, Serverside, ServersideOutputTransform
import dash_bootstrap_components as dbc
import base64
import configparser

# will be a dictionary of MALDISpectrum objects used for the spectrum plot
INDEXED_DATA = {}
# default processing parameters from config file
config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), 'etc', 'preprocessing.cfg'))
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
# relative path for directory where uploaded data is stored
UPLOAD_DIR = 'data'
if not os.path.exists(UPLOAD_DIR):
    os.makedirs(UPLOAD_DIR)

# Use DashProxy instead of Dash to allow for multiple callbacks to the same plot
app = DashProxy(prevent_initial_callbacks=True,
                transforms=[MultiplexerTransform(), ServersideOutputTransform()],
                external_stylesheets=[dbc.themes.BOOTSTRAP])
app.layout = get_dashboard_layout(PREPROCESSING_PARAMS)


@app.callback(Output('dropdown', 'children'),
              Input('upload', 'contents'),
              State('upload', 'filename'))
def upload_data(list_of_contents, list_of_filenames):
    global INDEXED_DATA
    if list_of_contents is not None:
        for contents, filename in zip(list_of_contents, list_of_filenames):
            content_type, content_string = contents.split(',')
            decoded = base64.b64decode(content_string)
            with open(os.path.join(UPLOAD_DIR, filename), 'w') as mzml_file:
                mzml_file.write(decoded.decode('utf-8'))
            data = import_mzml(os.path.join(UPLOAD_DIR, filename))
            for spectrum in data:
                INDEXED_DATA[spectrum.spectrum_id] = spectrum
    return get_dropdown_layout(INDEXED_DATA)


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('spectrum_id', 'value'))
def plot_spectrum(value):
    global INDEXED_DATA
    fig = get_spectrum(INDEXED_DATA[value])
    for filename in os.listdir('file_system_backend'):
        os.remove(os.path.join('file_system_backend', filename))
    return [get_spectrum_plot_layout(fig)], Serverside(fig)


@app.callback(Output('edit_processing_parameters_modal', 'is_open'),
              [Input('edit_preprocessing_parameters', 'n_clicks'),
               Input('edit_processing_parameters_save', 'n_clicks'),
               Input('edit_processing_parameters_cancel', 'n_clicks')],
              State('edit_processing_parameters_modal', 'is_open'))
def toggle_edit_preprocessing_parameters_modal(n_clicks_button, n_clicks_save, n_clicks_cancel, is_open):
    changed_id = [i['prop_id'] for i in callback_context.triggered][0]
    if n_clicks_button or n_clicks_save or n_clicks_cancel:
        if changed_id == 'edit_processing_parameters_save.n_clicks':
            # TODO: add code to parse edited parameters.
            # TODO: need a global dict to store processing parameters; use code from flex_maldi_dda_automation config file
            print('save')
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


@app.callback(Output('smooth_baseline_method_parameters', 'children'),
              [Input('edit_preprocessing_parameters', 'n_clicks'),
               Input('smooth_baseline_method', 'value')])
def toggle_smooth_baseline_method_parameters(n_clicks, value):
    global PREPROCESSING_PARAMS
    if value == 'SavitzkyGolay':
        return get_smooth_baseline_savitzky_golay_parameters(PREPROCESSING_PARAMS)
    elif value == 'apodization':
        return get_smooth_baseline_apodization_parameters(PREPROCESSING_PARAMS)
    elif value == 'rebin':
        return get_smooth_baseline_rebin_parameters(PREPROCESSING_PARAMS)
    elif value == 'fast_change':
        return get_smooth_baseline_fast_change_parameters(PREPROCESSING_PARAMS)
    elif value == 'median':
        return get_smooth_baseline_median_parameters(PREPROCESSING_PARAMS)


@app.callback(Output('remove_baseline_method_parameters', 'children'),
              [Input('edit_preprocessing_parameters', 'n_clicks'),
               Input('remove_baseline_method', 'value')])
def toggle_remove_baseline_method_parameters(n_clicks, value):
    global PREPROCESSING_PARAMS
    if value == 'SNIP':
        return get_remove_baseline_snip_parameters(PREPROCESSING_PARAMS)
    elif value == 'TopHat':
        return get_remove_baseline_tophat_parameters(PREPROCESSING_PARAMS)
    elif value == 'Median':
        return get_remove_baseline_median_parameters(PREPROCESSING_PARAMS)
    elif value == 'ZhangFit':
        return get_remove_baseline_zhangfit_parameters(PREPROCESSING_PARAMS)
    elif value == 'ModPoly':
        return get_remove_baseline_modpoly_parameters(PREPROCESSING_PARAMS)
    elif value == 'IModPoly':
        return get_remove_baseline_imodpoly_parameters(PREPROCESSING_PARAMS)


@app.callback(Output('peak_picking_method_parameters', 'children'),
              [Input('edit_preprocessing_parameters', 'n_clicks'),
               Input('peak_picking_method', 'value')])
def toggle_peak_picking_method_parameters(n_clicks, value):
    global PREPROCESSING_PARAMS
    if value == 'locmax':
        return get_peak_picking_locmax_parameters(PREPROCESSING_PARAMS)
    elif value == 'cwt':
        return get_peak_picking_cwt_parameters(PREPROCESSING_PARAMS)


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('trim_spectrum', 'n_clicks'),
              State('spectrum_id', 'value'))
def trim_spectrum_button(n_clicks, value):
    global INDEXED_DATA
    global PREPROCESSING_PARAMS
    INDEXED_DATA[value].trim_spectrum(**PREPROCESSING_PARAMS['TRIM_SPECTRUM'])
    fig = get_spectrum(INDEXED_DATA[value])
    for filename in os.listdir('file_system_backend'):
        os.remove(os.path.join('file_system_backend', filename))
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
    for filename in os.listdir('file_system_backend'):
        os.remove(os.path.join('file_system_backend', filename))
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
    for filename in os.listdir('file_system_backend'):
        os.remove(os.path.join('file_system_backend', filename))
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
    for filename in os.listdir('file_system_backend'):
        os.remove(os.path.join('file_system_backend', filename))
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
    for filename in os.listdir('file_system_backend'):
        os.remove(os.path.join('file_system_backend', filename))
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
    for filename in os.listdir('file_system_backend'):
        os.remove(os.path.join('file_system_backend', filename))
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
    for filename in os.listdir('file_system_backend'):
        os.remove(os.path.join('file_system_backend', filename))
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
    for filename in os.listdir('file_system_backend'):
        os.remove(os.path.join('file_system_backend', filename))
    return [get_spectrum_plot_layout(fig)], Serverside(fig)


@app.callback(Output('peak_list', 'data'),
              Input('export_peak_list', 'n_clicks'),
              State('spectrum_id', 'value'))
def export_peak_list(n_clicks, value):
    global INDEXED_DATA
    if INDEXED_DATA[value].peak_picked_mz_array is None and INDEXED_DATA[value].peak_picked_intensity_array is None:
        INDEXED_DATA[value].peak_picking()
    spectrum_df = pd.DataFrame(data={'m/z': copy.deepcopy(INDEXED_DATA[value].peak_picked_mz_array),
                                     'Intensity': copy.deepcopy(INDEXED_DATA[value].peak_picked_intensity_array)})
    return dcc.send_data_frame(spectrum_df.to_csv, value + '|peak_list.csv', index=False)


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('undo_preprocessing', 'n_clicks'),
              State('spectrum_id', 'value'))
def undo_preprocessing(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value].undo_all_processing()
    fig = get_spectrum(INDEXED_DATA[value])
    for filename in os.listdir('file_system_backend'):
        os.remove(os.path.join('file_system_backend', filename))
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
