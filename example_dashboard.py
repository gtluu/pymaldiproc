from pymaldiproc.classes import *
from pymaldiproc.data_import import *
from pymaldiproc.preprocessing import *
import plotly.express as px
from dash import Dash, dcc, html, State, callback_context
from dash.exceptions import PreventUpdate   
from dash_extensions.enrich import Input, Output, DashProxy, MultiplexerTransform
import dash_bootstrap_components as dbc
import base64

# will be a dictionary of MALDISpectrum objects used for the spectrum plot
INDEXED_DATA = {}

# relative path for directory where uploaded data is stored
UPLOAD_DIR = 'data'
if not os.path.exists(UPLOAD_DIR):
    os.makedirs(UPLOAD_DIR)

# Use DashProxy instead of Dash to allow for multiple callbacks to the same plot
app = DashProxy(prevent_initial_callbacks=True, transforms=[MultiplexerTransform()])

# barebones initial app layout. html "children" elements returned by callback functions and added to this on the fly
app.layout = html.Div([
    html.Div(
        dcc.Upload(
            id='upload',
            children=html.Div(
                [
                    'Drag and Drop or ',
                    html.A('Select mzML Files')
                ]
            ),
            style={
                'width': '97%',
                'height': '100px',
                'lineHeight': '100px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '20px'
            },
            multiple=True
        )
    ),

    html.Div(
        [
            html.Div(id='dropdown', className='one column')
        ]
    ),

    html.Div(
        id='spectrum',
        className='row'
    )
])


# get a spectrum df with m/z and intensity columns and plot it using a px.line plot
# also add UI elements for preprocessing buttons
def get_spectrum_graph(spectrum):
    spectrum_df = pd.DataFrame(data={'m/z': spectrum.get_mz_array(),
                                     'Intensity': spectrum.get_intensity_array()})

    fig = px.line(data_frame=spectrum_df,
                  x='m/z',
                  y='Intensity')
    fig.update_layout(xaxis_tickformat='d',
                      yaxis_tickformat='~e')

    children = [
        html.Div(
            dcc.Graph(
                id='spectrum_plot',
                figure=fig,
                style={
                    'width': '100%',
                    'height': '600px'
                }
            )
        ),
        html.Div(
            html.H1('Preprocessing', className='row')
        ),
        html.Div(
            [
                #html.Button('Trim Spectrum', id='trim_spectrum'),
                html.Button('Transform Intensity', id='transform_intensity'),
                html.Button('Smooth Baseline', id='smooth_baseline'),
                html.Button('Remove Baseline', id='remove_baseline'),
                #html.Button('Normalize Intensity', id='normalize_intensity'),
                #html.Button('Bin Spectrum', id='bin_spectrum'),
                html.Button('Peak Picking', id='peak_picking'),
                html.Button('Export Peak List', id='export_peak_list'),
                html.Button('Undo Preprocessing', id='undo_preprocessing'),
                dcc.Download(id='peak_list')
            ]
        ),

        html.Div(
            html.H1('Modifications (optional)', className='row')
        ),
        html.Div(
            [   
                html.H4("Transform Intensity Method:", 
                    style={'marginTop': '10px', 'marginBottom': '5px'}
                ),
                dcc.RadioItems(
                    id='transform_options',
                    options=[
                        {'label': 'sqrt', 'value': 'sqrt'},
                        {'label': 'ln', 'value': 'log'},
                        {'label': 'log2', 'value': 'log2'},
                        {'label': 'log10', 'value': 'log10'}
                    ],
                    value='sqrt',  # default value
                    labelStyle={'display': 'inline-block', 'margin-right': '20px'}  # display inline
                )
            ], style={'padding': '10px', 'textAlign': 'left'}
        ),
        html.Div(
            [
                html.H4("Smooth Baseline Method:",
                    style={'marginTop': '5px', 'marginBottom': '10px'}
                ),
                dcc.RadioItems(
                    id='smooth_options',
                    options=[
                        {'label': 'Savitzky-Golay', 'value': 'SavitzkyGolay'},
                        {'label': 'Apodization', 'value': 'apodization'},
                        {'label': 'Rebin', 'value': 'rebin'},
                        {'label': 'Fast Change', 'value': 'fast_change'},
                        {'label': 'Median', 'value': 'median'}
                    ],
                    value='SavitzkyGolay',  # default value
                    labelStyle={'display': 'inline-block', 'margin-right': '20px'}
                )
            ], style={'padding': '10px'}
        ),
        html.Div(
            [
                html.H4("Remove Baseline Method:",
                    style={'marginTop': '5px', 'marginBottom': '10px'}
                ),
                dcc.RadioItems(
                    id='remove_options',
                    options=[
                        {'label': 'SNIP', 'value': 'SNIP'},
                        {'label': 'Top Hat', 'value': 'TopHat'},
                        {'label': 'Median', 'value': 'Median'},
                        {'label': 'Zhang Fit', 'value': 'ZhangFit'},
                        {'label': 'ModPoly', 'value': 'ModPoly'},
                        {'label': 'IModPoly', 'value': 'IModPoly'}
                    ],
                    value='SNIP',  # default value
                    labelStyle={'display': 'inline-block', 'margin-right': '20px'}
                )
            ], style={'padding': '10px'}
        ),
    ]

    return children


# after uploading data, add a dropdown UI element to select a spectrum
@app.callback(Output('dropdown', 'children'),
              Input('upload', 'contents'),
              State('upload', 'filename'))
def upload_data(list_of_contents, list_of_filenames):
    # uploaded data is added to the global variable INDEXED_DATA here
    # the following callbacks below here also call "global INDEXED_DATA" to use this same persistent data and modify
    # the spectrum plot
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
                # INDEXED_DATA[spectrum.spectrum_id] = {'original': spectrum} #, 'transformed': None}

        children = [
            html.Div(
                html.H1('Spectrum ID', className='row')
            ),
            html.Div(
                dcc.Dropdown(id='spectrum_id',
                             multi=False,
                             options=[{'label': i, 'value': i} for i in INDEXED_DATA.keys()],
                             value=[i for i in INDEXED_DATA.keys()])
            )
        ]

        return children


@app.callback(Output('spectrum', 'children'),
              Input('spectrum_id', 'value'))
def graph_spectrum(value):
    print("graph started")
    if value is not None:
        #global INDEXED_DATA
        spectrum = INDEXED_DATA[value]
        children = get_spectrum_graph(spectrum)
        print("graph done")
        return children


# @app.callback(Output('spectrum', 'children'),
#               Input('transform_intensity', 'n_clicks'),
#               State('spectrum_id', 'value'))
# def transform_intensity_button(n_clicks, value):
#     print("transform button clicked")
#     #global INDEXED_DATA
#     INDEXED_DATA[value] = transform_intensity([INDEXED_DATA[value]])[0]
#     spectrum = INDEXED_DATA[value]
#     children = get_spectrum_graph(spectrum)
#     print("transform button done")
#     return children

# @app.callback(
#     Output('spectrum', 'children'),
#     Input('transform_intensity', 'n_clicks'),
#     State('spectrum_id', 'value')
# )
# def transform_intensity(n_clicks, value):
#     if n_clicks is not None:
#         global INDEXED_DATA
#         spectrum = INDEXED_DATA.get(value)
#         if spectrum:
#             transformed_spectra = transform_intensity([spectrum], method='sqrt')  # todo: check for method and apply
#             INDEXED_DATA[value] = transformed_spectra[0]  # update spectrum
#             children = get_spectrum_graph(transformed_spectra[0])   #update graph with new data
#             print("Intensity transformation applied")
#             return children
#     raise PreventUpdate #prevent other callbacks from randonly being called

@app.callback(
    Output('spectrum', 'children'),
    [Input('transform_intensity', 'n_clicks'),
     State('spectrum_id', 'value'),
     State('transform_options', 'value')]
)
def transform_intensity_callback(n_clicks, spectrum_id, selected_method):
    if n_clicks is not None:
        global INDEXED_DATA
        spectrum = INDEXED_DATA.get(spectrum_id)
        if spectrum:
            transformed_spectrum = transform_intensity([spectrum], selected_method)[0]  
            INDEXED_DATA[spectrum_id] = transformed_spectrum
            return get_spectrum_graph(transformed_spectrum)
    raise PreventUpdate

# @app.callback(Output('spectrum', 'children'),
#               Input('smooth_baseline', 'n_clicks'),
#               State('spectrum_id', 'value'))
# def smooth_baseline_button(n_clicks, value):
#     print("smooth button clicked")
#     #global INDEXED_DATA
#     INDEXED_DATA[value] = smooth_baseline([INDEXED_DATA[value]])[0]
#     spectrum = INDEXED_DATA[value]
#     children = get_spectrum_graph(spectrum)
#     print("smooth button done")
#     return children

@app.callback(
    Output('spectrum', 'children'),  # Assuming this updates your spectrum graph
    Input('smooth_baseline', 'n_clicks'),  # Button for smoothing baseline
    State('spectrum_id', 'value')  # Dropdown or input that specifies which spectrum to operate on
)
def smooth_baseline_button(n_clicks, spectrum_id):
    if n_clicks is None:
        raise PreventUpdate
    spectrum = INDEXED_DATA.get(spectrum_id)
    if not spectrum:
        raise PreventUpdate
    smoothed_spectrum = smooth_baseline([spectrum], method='SavitzkyGolay')[0] 
    INDEXED_DATA[spectrum_id] = smoothed_spectrum
    return get_spectrum_graph(smoothed_spectrum)


# @app.callback(Output('spectrum', 'children'),
#               Input('remove_baseline', 'n_clicks'),
#               State('spectrum_id', 'value'))
# def remove_baseline_button(n_clicks, value):
#     print("remove button clicked")
#     #global INDEXED_DATA
#     INDEXED_DATA[value] = remove_baseline([INDEXED_DATA[value]])[0]
#     spectrum = INDEXED_DATA[value]
#     children = get_spectrum_graph(spectrum)
#     print("remove button done")
#     return children

@app.callback(
    Output('spectrum', 'children'),
    Input('remove_baseline', 'n_clicks'),
    [State('spectrum_id', 'value'),
     State('remove_options', 'value')] 
)
def on_remove_baseline_click(n_clicks, spectrum_id, baseline_method):
    if n_clicks is None:
        raise PreventUpdate
    original_spectrum = INDEXED_DATA[spectrum_id]['original'] # retrieve original data
    processed_spectrum = remove_baseline([original_spectrum], method=baseline_method)[0]    #use correct method
    INDEXED_DATA[spectrum_id]['processed'] = processed_spectrum  # update w processed data
    return get_spectrum_graph(processed_spectrum)



# '''@app.callback(Output('spectrum', 'children'),
#               Input('normalize_intensity', 'n_clicks'),
#               State('spectrum_id', 'value'))
# def normalize_intensity_button(n_clicks, value):
#     global INDEXED_DATA
#     INDEXED_DATA[value] = normalize_intensity([INDEXED_DATA[value]])[0]
#     spectrum = INDEXED_DATA[value]
#     children = get_spectrum_graph(spectrum)
#     return children'''


# @app.callback(Output('spectrum', 'children'),
#               Input('peak_picking', 'n_clicks'),
#               State('spectrum_id', 'value'))
# def peak_picking_button(n_clicks, value):
#     print("peak button clicked")
#     #global INDEXED_DATA
#     INDEXED_DATA[value] = peak_picking([INDEXED_DATA[value]])[0]
#     spectrum = INDEXED_DATA[value]
#     children = get_spectrum_graph(spectrum)
#     print("peak button done")
#     return children


# @app.callback(Output('peak_list', 'data'),
#               Input('export_peak_list', 'n_clicks'),
#               State('spectrum_id', 'value'))
# def export_peak_list(n_clicks, value):
#     print("export button clicked")
#     spectrum = INDEXED_DATA[value]
#     spectrum_df = pd.DataFrame(data={'m/z': spectrum.get_mz_array(),
#                                      'Intensity': spectrum.get_intensity_array()})
#     print("export button done")
#     return dcc.send_data_frame(spectrum_df.to_csv, value + '|peak_list.csv', index=False)


# @app.callback(Output('spectrum', 'children'),
#               Input('undo_preprocessing', 'n_clicks'),
#               State('spectrum_id', 'value'))
# def undo_preprocessing(n_clicks, value):
#     print("undo button clicked")
#     #global INDEXED_DATA
#     INDEXED_DATA[value].preprocessed_mz_array = None
#     INDEXED_DATA[value].preprocessed_intensity_array = None
#     INDEXED_DATA[value].peak_picked_mz_array = None
#     INDEXED_DATA[value].peak_picked_intensity_array = None
#     INDEXED_DATA[value].data_processing = {}
#     spectrum = INDEXED_DATA[value]
#     children = get_spectrum_graph(spectrum)
#     print("undo button clicked")
#     return children
@app.callback(
    Output('spectrum', 'children'),
    Input('undo_preprocessing', 'n_clicks'),
    State('spectrum_id', 'value')
)
def undo_preprocessing(n_clicks, value):
    if n_clicks is None:
        raise PreventUpdate
    spectrum_info = INDEXED_DATA.get(value)
    if spectrum_info:
        return get_spectrum_graph(spectrum_info['original'])
    raise PreventUpdate  # display original graph


if __name__ == '__main__':
    #app.run_server(debug=False)
    app.run_server(port=8051, debug=False)
