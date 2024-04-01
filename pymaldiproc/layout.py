from pymaldiproc.classes import *
from pymaldiproc.data_import import *
from pymaldiproc.preprocessing import *
import plotly.express as px
from plotly_resampler import register_plotly_resampler, FigureResampler
from dash import Dash, dcc, html, State, callback_context
from dash_extensions.enrich import Input, Output, DashProxy, MultiplexerTransform
import dash_bootstrap_components as dbc
import base64


# barebones initial app layout. html "children" elements returned by callback functions and added to this on the fly
def get_dashboard_layout():
    dashboard_layout = html.Div([
        html.Div(
            # TODO: replace this with a button that gets a file path for input.
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
        ),
        dcc.Loading(dcc.Store(id='store_plot'))
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


def get_preprocessing_parameters_layout():
    trim_spectrum_parameters = html.Div([
        html.H5('Spectrum Trimming Parameters')
    ], style={'margin': '20px'})

    transform_intensity_parameters = html.Div([
        html.H5('Intensity Transformation Parameters'),
        html.P('Method'),
        dcc.RadioItems(
            id='transform_options',
            options=[
                {'label': 'Sqrt', 'value': 'sqrt'},
                {'label': 'Natural Log', 'value': 'ln'},
                {'label': 'Log Base 2', 'value': 'log2'},
                {'label': 'Log Base 10', 'value': 'log10'}
            ],
            value='sqrt',
            labelStyle={'display': 'inline-block', 'marginRight': '20px'},
            inputStyle={'margin-right': '10px'}
        )
    ], style={'margin': '20px'})

    smooth_baseline_parameters = html.Div([
        html.H5('Baseline Smoothing Parameters'),
        html.P('Method'),
        dcc.RadioItems(
            id='smooth_options',
            options=[
                {'label': 'Savitzky-Golay', 'value': 'SavitzkyGolay'},
                {'label': 'Apodization', 'value': 'apodization'},
                {'label': 'Rebin', 'value': 'rebin'},
                {'label': 'Fast Change', 'value': 'fast_change'},
                {'label': 'Median', 'value': 'median'}
            ],
            value='SavitzkyGolay',
            labelStyle={'display': 'inline-block', 'marginRight': '20px'},
            inputStyle={'margin-right': '10px'}
        )
    ], style={'margin': '20px'})

    remove_baseline_parameters = html.Div([
        html.H5('Baseline Removal Parameters'),
        html.P('Method'),
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
            value='SNIP',
            labelStyle={'display': 'inline-block', 'marginRight': '20px'},
            inputStyle={'margin-right': '10px'}
        )
    ], style={'margin': '20px'})

    normalize_intensity_parameters = html.Div([
        html.H5('Intensity Normalization Parameters'),
        html.P('Method'),
        dcc.RadioItems(
            id='normalize_options',
            options=[
                {'label': 'TIC', 'value': 'tic'},
                {'label': 'RMS', 'value': 'rms'},
                {'label': 'MAD', 'value': 'mad'},
                {'label': 'Sqrt', 'value': 'sqrt'}
            ],
            value='tic',
            labelStyle={'display': 'inline-block', 'marginRight': '20px'},
            inputStyle={'margin-right': '10px'}
        )
    ], style={'margin': '20px'})

    bin_spectrum_parameters = html.Div([
        html.H5('Spectrum Binning Parameters')
    ], style={'margin': '20px'})

    peak_picking_parameters = html.Div([
        html.H5('Peak Picking Parameters'),
        html.P('Method'),
        dcc.RadioItems(
            id='remove_options',
            options=[
                {'label': 'LocMax', 'value': 'locmax'},
                {'label': 'CWT', 'value': 'cwt'}
            ],
            value='locmax',
            labelStyle={'display': 'inline-block', 'marginRight': '20px'},
            inputStyle={'margin-right': '10px'}
        )
    ], style={'margin': '20px'})

    return [trim_spectrum_parameters,
            transform_intensity_parameters,
            smooth_baseline_parameters,
            remove_baseline_parameters,
            normalize_intensity_parameters,
            bin_spectrum_parameters,
            peak_picking_parameters]


def get_preprocessing_layout():
    preprocessing_title = html.Div(html.H1('Preprocessing', className='row'))

    preprocessing_buttons = html.Div([
        dbc.Button('Trim Spectrum', id='trim_spectrum', style={'margin': '5px'}),
        dbc.Button('Transform Intensity', id='transform_intensity', style={'margin': '5px'}),
        dbc.Button('Smooth Baseline', id='smooth_baseline', style={'margin': '5px'}),
        dbc.Button('Remove Baseline', id='remove_baseline', style={'margin': '5px'}),
        dbc.Button('Normalize Intensity', id='normalize_intensity', style={'margin': '5px'}),
        dbc.Button('Bin Spectrum', id='bin_spectrum', style={'margin': '5px'}),
        dbc.Button('Label Peaks', id='peak_picking', style={'margin': '5px'}),
        dbc.Button('Export Peak List from Labeled Peaks', id='export_peak_list', style={'margin': '5px'}),
        dbc.Button('Undo Preprocessing', id='undo_preprocessing', style={'margin': '5px'}),
        dbc.Button('Undo Peak Labeling', id='undo_peak_picking', style={'margin': '5px'}),
        dbc.Button('Edit Preprocessing Parameters', id='edit_preprocessing_parameters', style={'margin': '5px'}),
        dbc.Modal([
            dbc.ModalHeader(dbc.ModalTitle('Preprocessing Parameters')),
            dbc.ModalBody(get_preprocessing_parameters_layout()),
            dbc.ModalFooter(dbc.ButtonGroup([
                dbc.Button('Cancel', id='edit_processing_parameters_cancel', className='ms-auto'),
                dbc.Button('Save', id='edit_processing_parameters_save', className='ms-auto')
            ]))
        ],
        id='edit_processing_parameters_modal',
        size='lg',
        backdrop='static',
        scrollable=True,
        centered=True,
        is_open=False),
        dcc.Download(id='peak_list')
    ])

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
                         # options=[{'label': '|'.join(i.split('|')[:-1]), 'value': i} for i in data.keys()],
                         value=[i for i in data.keys()])
        )
    ]
    return dropdown + get_preprocessing_layout()


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

    # TODO: add code to add peak labels
    if spectrum.peak_picked_mz_array is not None and spectrum.peak_picked_intensity_array is not None:
        pass

    return fig
