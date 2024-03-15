from pymaldiproc.classes import *
from pymaldiproc.data_import import *
from pymaldiproc.preprocessing import *
import plotly.express as px
from dash import Dash, dcc, html, State, callback_context
from dash_extensions.enrich import Input, Output, DashProxy, MultiplexerTransform
import dash_bootstrap_components as dbc
import base64


# barebones initial app layout. html "children" elements returned by callback functions and added to this on the fly
def get_dashboard_layout():
    dashboard_layout = html.Div([
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
    return dashboard_layout


def get_spectrum_plot_layout(fig):
    spectrum_plot = html.Div(
        dcc.Graph(
            id='spectrum_plot',
            figure=fig,
            style={
                'width': '100%',
                'height': '600px'
            }
        )
    )
    return spectrum_plot


def get_preprocessing_layout():
    preprocessing_title = html.Div(
        html.H1('Preprocessing', className='row')
    )

    preprocessing_buttons = html.Div(
        [
            html.Button('Trim Spectrum', id='trim_spectrum'),
            html.Button('Transform Intensity', id='transform_intensity'),
            html.Button('Smooth Baseline', id='smooth_baseline'),
            html.Button('Remove Baseline', id='remove_baseline'),
            html.Button('Normalize Intensity', id='normalize_intensity'),
            html.Button('Bin Spectrum', id='bin_spectrum'),
            html.Button('Peak Picking', id='peak_picking'),
            html.Button('Export Current Peak List', id='export_current_peak_list'),
            html.Button('Undo Preprocessing', id='undo_preprocessing'),
            dcc.Download(id='peak_list')
        ]
    )

    return [preprocessing_title, preprocessing_buttons]


def get_dropdown_layout(data):
    dropdown = [
        html.Div(
            html.H1('Spectrum ID', className='row')
        ),
        html.Div(
            dcc.Dropdown(id='spectrum_id',
                         multi=False,
                         options=[{'label': i, 'value': i} for i in data.keys()],
                         #options=[{'label': '|'.join(i.split('|')[:-1]), 'value': i} for i in data.keys()],
                         value=[i for i in data.keys()])
        )
    ]
    return dropdown + get_preprocessing_layout()


def get_spectrum(spectrum):
    spectrum_df = pd.DataFrame({'m/z': copy.deepcopy(spectrum.preprocessed_mz_array),
                                'Intensity': copy.deepcopy(spectrum.preprocessed_intensity_array)})

    fig = px.line(data_frame=spectrum_df,
                  x='m/z',
                  y='Intensity')
    fig.update_layout(xaxis_tickformat='d',
                      yaxis_tickformat='~e')

    return [get_spectrum_plot_layout(fig)]
