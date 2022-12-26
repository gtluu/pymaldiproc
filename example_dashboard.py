from pymaldiproc.classes import *
from pymaldiproc.data_import import *
from pymaldiproc.preprocessing import *
import plotly.express as px
from dash import Dash, dcc, html, State, callback_context
from dash_extensions.enrich import Input, Output, DashProxy, MultiplexerTransform
import dash_bootstrap_components as dbc
import base64

INDEXED_DATA = {}
UPLOAD_DIR = 'data'
if not os.path.exists(UPLOAD_DIR):
    os.makedirs(UPLOAD_DIR)

#app = Dash(__name__)
app = DashProxy(prevent_initial_callbacks=True, transforms=[MultiplexerTransform()])

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
    if value is not None:
        global INDEXED_DATA
        spectrum = INDEXED_DATA[value]
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
                    html.Button('Normalize Intensity', id='normalize_intensity'),
                    #html.Button('Bin Spectrum', id='bin_spectrum'),
                    html.Button('Peak Picking', id='peak_picking'),
                    html.Button('Export Peak List', id='export_peak_list'),
                    html.Button('Undo Preprocessing', id='undo_preprocessing'),
                    dcc.Download(id='peak_list')
                ]
            ),
        ]

        return children


@app.callback(Output('spectrum', 'children'),
              Input('transform_intensity', 'n_clicks'),
              State('spectrum_id', 'value'))
def transform_intensity_button(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value] = transform_intensity([INDEXED_DATA[value]])[0]
    spectrum = INDEXED_DATA[value]
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
                # html.Button('Trim Spectrum', id='trim_spectrum'),
                html.Button('Transform Intensity', id='transform_intensity'),
                html.Button('Smooth Baseline', id='smooth_baseline'),
                html.Button('Remove Baseline', id='remove_baseline'),
                html.Button('Normalize Intensity', id='normalize_intensity'),
                #html.Button('Bin Spectrum', id='bin_spectrum'),
                html.Button('Peak Picking', id='peak_picking'),
                html.Button('Export Peak List', id='export_peak_list'),
                html.Button('Undo Preprocessing', id='undo_preprocessing'),
                dcc.Download(id='peak_list')
            ]
        ),
    ]

    return children


@app.callback(Output('spectrum', 'children'),
              Input('smooth_baseline', 'n_clicks'),
              State('spectrum_id', 'value'))
def smooth_baseline_button(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value] = smooth_baseline([INDEXED_DATA[value]])[0]
    spectrum = INDEXED_DATA[value]
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
                # html.Button('Trim Spectrum', id='trim_spectrum'),
                html.Button('Transform Intensity', id='transform_intensity'),
                html.Button('Smooth Baseline', id='smooth_baseline'),
                html.Button('Remove Baseline', id='remove_baseline'),
                html.Button('Normalize Intensity', id='normalize_intensity'),
                #html.Button('Bin Spectrum', id='bin_spectrum'),
                html.Button('Peak Picking', id='peak_picking'),
                html.Button('Export Peak List', id='export_peak_list'),
                html.Button('Undo Preprocessing', id='undo_preprocessing'),
                dcc.Download(id='peak_list')
            ]
        ),
    ]

    return children


@app.callback(Output('spectrum', 'children'),
              Input('remove_baseline', 'n_clicks'),
              State('spectrum_id', 'value'))
def remove_baseline_button(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value] = remove_baseline([INDEXED_DATA[value]])[0]
    spectrum = INDEXED_DATA[value]
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
                # html.Button('Trim Spectrum', id='trim_spectrum'),
                html.Button('Transform Intensity', id='transform_intensity'),
                html.Button('Smooth Baseline', id='smooth_baseline'),
                html.Button('Remove Baseline', id='remove_baseline'),
                html.Button('Normalize Intensity', id='normalize_intensity'),
                #html.Button('Bin Spectrum', id='bin_spectrum'),
                html.Button('Peak Picking', id='peak_picking'),
                html.Button('Export Peak List', id='export_peak_list'),
                html.Button('Undo Preprocessing', id='undo_preprocessing'),
                dcc.Download(id='peak_list')
            ]
        ),
    ]

    return children


@app.callback(Output('spectrum', 'children'),
              Input('normalize_intensity', 'n_clicks'),
              State('spectrum_id', 'value'))
def normalize_intensity_button(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value] = normalize_intensity([INDEXED_DATA[value]])[0]
    spectrum = INDEXED_DATA[value]
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
                # html.Button('Trim Spectrum', id='trim_spectrum'),
                html.Button('Transform Intensity', id='transform_intensity'),
                html.Button('Smooth Baseline', id='smooth_baseline'),
                html.Button('Remove Baseline', id='remove_baseline'),
                html.Button('Normalize Intensity', id='normalize_intensity'),
                #html.Button('Bin Spectrum', id='bin_spectrum'),
                html.Button('Peak Picking', id='peak_picking'),
                html.Button('Export Peak List', id='export_peak_list'),
                html.Button('Undo Preprocessing', id='undo_preprocessing'),
                dcc.Download(id='peak_list')
            ]
        ),
    ]

    return children


@app.callback(Output('spectrum', 'children'),
              Input('peak_picking', 'n_clicks'),
              State('spectrum_id', 'value'))
def peak_picking_button(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value] = peak_picking([INDEXED_DATA[value]])[0]
    spectrum = INDEXED_DATA[value]
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
                # html.Button('Trim Spectrum', id='trim_spectrum'),
                html.Button('Transform Intensity', id='transform_intensity'),
                html.Button('Smooth Baseline', id='smooth_baseline'),
                html.Button('Remove Baseline', id='remove_baseline'),
                html.Button('Normalize Intensity', id='normalize_intensity'),
                #html.Button('Bin Spectrum', id='bin_spectrum'),
                html.Button('Peak Picking', id='peak_picking'),
                html.Button('Export Peak List', id='export_peak_list'),
                html.Button('Undo Preprocessing', id='undo_preprocessing'),
                dcc.Download(id='peak_list')
            ]
        ),
    ]

    return children


@app.callback(Output('peak_list', 'data'),
              Input('export_peak_list', 'n_clicks'),
              State('spectrum_id', 'value'))
def export_peak_list(n_clicks, value):
    spectrum = INDEXED_DATA[value]
    spectrum_df = pd.DataFrame(data={'m/z': spectrum.get_mz_array(),
                                     'Intensity': spectrum.get_intensity_array()})

    return dcc.send_data_frame(spectrum_df.to_csv, value + '|peak_list.csv', index=False)


@app.callback(Output('spectrum', 'children'),
              Input('undo_preprocessing', 'n_clicks'),
              State('spectrum_id', 'value'))
def undo_preprocessing(n_clicks, value):
    global INDEXED_DATA
    INDEXED_DATA[value].preprocessed_mz_array = None
    INDEXED_DATA[value].preprocessed_intensity_array = None
    INDEXED_DATA[value].peak_picked_mz_array = None
    INDEXED_DATA[value].peak_picked_intensity_array = None
    INDEXED_DATA[value].data_processing = {}
    spectrum = INDEXED_DATA[value]
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
                # html.Button('Trim Spectrum', id='trim_spectrum'),
                html.Button('Transform Intensity', id='transform_intensity'),
                html.Button('Smooth Baseline', id='smooth_baseline'),
                html.Button('Remove Baseline', id='remove_baseline'),
                html.Button('Normalize Intensity', id='normalize_intensity'),
                # html.Button('Bin Spectrum', id='bin_spectrum'),
                html.Button('Peak Picking', id='peak_picking'),
                html.Button('Export Peak List', id='export_peak_list'),
                html.Button('Undo Preprocessing', id='undo_preprocessing'),
                dcc.Download(id='peak_list')
            ]
        ),
    ]

    return children


if __name__ == '__main__':
    app.run_server(debug=True)
