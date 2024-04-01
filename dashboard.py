from pymaldiproc.classes import *
from pymaldiproc.data_import import *
from pymaldiproc.preprocessing import *
from pymaldiproc.layout import *
import plotly.express as px
from dash import Dash, dcc, html, State, callback_context, no_update
from dash_extensions.enrich import Input, Output, DashProxy, MultiplexerTransform, Serverside, ServersideOutputTransform
import dash_bootstrap_components as dbc
import base64

# will be a dictionary of MALDISpectrum objects used for the spectrum plot
INDEXED_DATA = {}
# relative path for directory where uploaded data is stored
UPLOAD_DIR = 'data'
if not os.path.exists(UPLOAD_DIR):
    os.makedirs(UPLOAD_DIR)

# Use DashProxy instead of Dash to allow for multiple callbacks to the same plot
app = DashProxy(prevent_initial_callbacks=True,
                transforms=[MultiplexerTransform(), ServersideOutputTransform()],
                external_stylesheets=[dbc.themes.BOOTSTRAP])
app.layout = get_dashboard_layout()


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
              [State('edit_processing_parameters_modal', 'is_open')])
def toggle_edit_preprocessing_parameters_modal(n_clicks_button, n_clicks_save, n_clicks_cancel, is_open):
    if n_clicks_button or n_clicks_save or n_clicks_cancel:
        if n_clicks_save:
            # TODO: add code to parse edited parameters.
            # TODO: need a global dict to store processing parameters; use code from flex_maldi_dda_automation config file
            pass
        return not is_open
    return is_open


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('trim_spectrum', 'n_clicks'),
              State('spectrum_id', 'value'))
def trim_spectrum_button(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value].trim_spectrum(100, 2000)
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
    INDEXED_DATA[value].transform_intensity()
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
    INDEXED_DATA[value].smooth_baseline()
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
    INDEXED_DATA[value].remove_baseline()
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
    INDEXED_DATA[value].normalize_intensity()
    fig = get_spectrum(INDEXED_DATA[value])
    for filename in os.listdir('file_system_backend'):
        os.remove(os.path.join('file_system_backend', filename))
    return [get_spectrum_plot_layout(fig)], Serverside(fig)


@app.callback([Output('spectrum', 'children'),
               Output('store_plot', 'data')],
              Input('bin_spectrum', 'n_clicks'),
              State('spectrum_id', 'value'))
def normalize_intensity_button(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value].bin_spectrum(10000, 100, 2000)
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
    INDEXED_DATA[value].peak_picking()
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
