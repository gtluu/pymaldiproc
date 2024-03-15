from pymaldiproc.classes import *
from pymaldiproc.data_import import *
from pymaldiproc.preprocessing import *
from pymaldiproc.layout import *
import plotly.express as px
from dash import Dash, dcc, html, State, callback_context
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
app = DashProxy(prevent_initial_callbacks=True, transforms=[MultiplexerTransform(), ServersideOutputTransform()])
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


@app.callback(Output('spectrum', 'children'),
              Input('spectrum_id', 'value'))
def plot_spectrum(value):
    global INDEXED_DATA
    fig = get_spectrum(INDEXED_DATA[value])
    return [get_spectrum_plot_layout(fig)]


@app.callback(Output('spectrum', 'children'),
              Input('trim_spectrum', 'n_clicks'),
              State('spectrum_id', 'value'))
def trim_spectrum_button(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value].trim_spectrum(100, 2000)
    return get_spectrum(INDEXED_DATA[value])


@app.callback(Output('spectrum', 'children'),
              Input('transform_intensity', 'n_clicks'),
              State('spectrum_id', 'value'))
def transform_intensity_button(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value].transform_intensity()
    return get_spectrum(INDEXED_DATA[value])


@app.callback(Output('spectrum', 'children'),
              Input('smooth_baseline', 'n_clicks'),
              State('spectrum_id', 'value'))
def smooth_baseline_button(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value].smooth_baseline()
    return get_spectrum(INDEXED_DATA[value])


@app.callback(Output('spectrum', 'children'),
              Input('remove_baseline', 'n_clicks'),
              State('spectrum_id', 'value'))
def remove_baseline_button(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value].remove_baseline()
    return get_spectrum(INDEXED_DATA[value])


@app.callback(Output('spectrum', 'children'),
              Input('normalize_intensity', 'n_clicks'),
              State('spectrum_id', 'value'))
def normalize_intensity_button(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value].normalize_intensity()
    return get_spectrum(INDEXED_DATA[value])


@app.callback(Output('spectrum', 'children'),
              Input('bin_spectrum', 'n_clicks'),
              State('spectrum_id', 'value'))
def normalize_intensity_button(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value].bin_spectrum(5000, 100, 2000)
    return get_spectrum(INDEXED_DATA[value])


@app.callback(Output('spectrum', 'children'),
              Input('peak_picking', 'n_clicks'),
              State('spectrum_id', 'value'))
def peak_picking_button(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value].peak_picking()
    return get_spectrum(INDEXED_DATA[value])


# TODO: add a button and callback to remove peak picking


@app.callback(Output('peak_list', 'data'),
              Input('export_current_peak_list', 'n_clicks'),
              State('spectrum_id', 'value'))
def export_peak_list(n_clicks, value):
    global INDEXED_DATA
    if INDEXED_DATA[value].peak_picked_mz_array is not None and \
            INDEXED_DATA[value].peak_picked_intensity_array is not None:
        spectrum_df = pd.DataFrame(data={'m/z': copy.deepcopy(INDEXED_DATA[value].peak_picked_mz_array),
                                         'Intensity': copy.deepcopy(INDEXED_DATA[value].peak_picked_intensity_array)})
    else:
        # TODO: add message that says no centroided peak list found
        spectrum_df = pd.DataFrame(data={'m/z': copy.deepcopy(INDEXED_DATA[value].preprocessed_mz_array),
                                         'Intensity': copy.deepcopy(INDEXED_DATA[value].preprocessed_intensity_array)})
    return dcc.send_data_frame(spectrum_df.to_csv, value + '|peak_list.csv', index=False)


@app.callback(Output('spectrum', 'children'),
              Input('undo_preprocessing', 'n_clicks'),
              State('spectrum_id', 'value'))
def undo_preprocessing(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value].undo_all_processing()
    return get_spectrum(INDEXED_DATA[value])


if __name__ == '__main__':
    app.run_server(debug=False)
