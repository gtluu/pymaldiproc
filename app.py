from pymaldiproc.classes import *
import plotly.express as px
from dash import Dash, Input, Output, dcc, html, State, callback_context
#from dash_extensions.enrich import Input, Output, DashProxy, MultiplexerTransform

dataset = MALDIDataset('prototype')
indexed_dataset = {spectrum.spectrum_id: spectrum for spectrum in dataset.spectra}

app = Dash(__name__)
#app = DashProxy(prevent_initial_callbacks=False, transforms=[MultiplexerTransform()])

app.layout = html.Div([
    html.Div(dcc.Graph(id='spectrum', figure={}), className='row'),
    html.Div(dcc.Graph(id='preprocessed_spectrum', figure={}), className='row'),
    html.Div([
        html.Div(
            html.H1('Spectrum ID', className='row')
        ),
        html.Div(dcc.Dropdown(id='spectrum_id',
                              multi=False,
                              options=[{'label': i.spectrum_id, 'value': i.spectrum_id} for i in dataset.spectra],
                              value=[i.spectrum_id for i in dataset.spectra]),
                 className='one column')
    ]),
    html.Div([
        html.Div(
            html.H1('Preprocessing', className='row')
        ),
        html.Div(
            html.H2('Preprocessing Steps/Order', className='row')
        ),
        html.Div(dcc.Dropdown(id='preprocessing',
                              multi=True,
                              options=[{'label': i, 'value': i} for i in ['trim spectrum',
                                                                          'transform intensity',
                                                                          'smooth baseline',
                                                                          'remove baseline',
                                                                          'normalize intensity',
                                                                          'peak picking']],
                              value=[i for i in ['trim spectrum',
                                                 'transform intensity',
                                                 'smooth baseline',
                                                 'remove baseline',
                                                 'normalize intensity',
                                                 'peak picking']]),
                 className='one column'),
        html.Div(
            html.H2('Preprocessing Parameters', className='row')
        ),
        html.Div(
            html.H3('Trim Spectra', className='row')
        ),
        html.Div([dcc.Input(id='trim_spectra_lower_mass_range',
                            placeholder='Lower Mass Range',
                            type='text',
                            value=''),
                  dcc.Input(id='trim_spectra_upper_mass_range',
                            placeholder='Upper Mass Range',
                            type='text',
                            value='')],
                 className='one column'),
        html.Div(
            html.H3('Transform Intensity', className='row')
        ),
        html.Div(dcc.Dropdown(id='transform_intensity_method',
                              multi=False,
                              options=[{'label': i, 'value': i} for i in ['sqrt', 'log', 'log10', 'log2']],
                              value=[i for i in ['sqrt', 'log', 'log10', 'log2']]),
                 className='one column'),
        html.Div(
            html.H3('Smooth Baseline', className='row')
        ),
        # method, window_length, polyorder, delta_mz, diff_thresh
        html.Div([dcc.Dropdown(id='smooth_baseline_method',
                               multi=False,
                               options=[{'label': i, 'value': i} for i in ['SavitzkyGolay',
                                                                           #'MovingAverage',
                                                                           'apodization',
                                                                           'rebin',
                                                                           'fast_change',
                                                                           'median']],
                               value=[i for i in ['SavitzkyGolay',
                                                  #'MovingAverage',
                                                  'apodization',
                                                  'rebin',
                                                  'fast_change',
                                                  'median']]),
                  dcc.Input(id='smooth_baseline_window_length',
                            placeholder='Window Length',
                            type='text',
                            value=20),
                  dcc.Input(id='smooth_baseline_polyorder',
                            placeholder='Polyorder (SavitzkyGolay)',
                            type='text',
                            value=3),
                  dcc.Input(id='smooth_baseline_delta_mz',
                            placeholder='Delta m/z (rebin)',
                            type='text',
                            value=0.2),
                  dcc.Input(id='smooth_baseline_diff_thresh',
                            placeholder='Difference Threshold (fast change)',
                            type='text',
                            value=0.01)],
                 className='one column'),
        html.Div(
            html.H3('Remove Baseline', className='row')
        ),
        html.Div(dcc.Dropdown(id='remove_baseline_method',
                              multi=False,
                              options=[{'label': i, 'value': i} for i in ['ZhangFit']],
                              value=[i for i in ['ZhangFit']]),
                 className='one column'),
        html.Div(
            html.H3('Normalize Intensity', className='row')
        ),
        html.Div(dcc.Dropdown(id='normalize_intensity_method',
                              multi=False,
                              options=[{'label': i, 'value': i} for i in ['tic', 'rms', 'mad', 'sqrt']],
                              value=[i for i in ['tic', 'rms', 'mad', 'sqrt']]),
                 className='one column'),
        html.Div(
            html.H3('Peak Picking', className='row')
        ),
        html.Div([dcc.Dropdown(id='peak_picking_method',
                               multi=False,
                               options=[{'label': i, 'value': i} for i in ['locmax', 'cwt']],
                               value=[i for i in ['locmax', 'cwt']]),
                  dcc.Input(id='peak_picking_widths_start',
                            placeholder='Peak Picking Widths Start',
                            type='text',
                            value='None'),
                  dcc.Input(id='peak_picking_widths_stop',
                            placeholder='Peak Picking Widths Stop',
                            type='text',
                            value='None'),
                  dcc.Input(id='peak_picking_snr',
                            placeholder='Signal to Noise Ratio',
                            type='text',
                            value=3)],
                 className='one column'),
        html.Div([
            html.H1('', className='row'),
            html.Button('Preprocess', id='preprocess_submit')
        ])
    ])
])


@app.callback(Output(component_id='spectrum', component_property='figure'),
              Input(component_id='spectrum_id', component_property='value'),)
def graph_spectrum(value):
    # get spectrum
    if value is not None:
        spectrum = indexed_dataset[value]
        spectrum_df = pd.DataFrame(data={'mz': spectrum.get_mz_array(),
                                         'intensity': spectrum.get_intensity_array()})

        fig = px.line(data_frame=spectrum_df,
                      x='mz',
                      y='intensity')
        fig.update_xaxes(range=[2000, 20000])
        return fig


@app.callback(Output(component_id='preprocessed_spectrum', component_property='figure'),
              Input(component_id='preprocess_submit', component_property='n_clicks'),
              State(component_id='spectrum_id', component_property='value'),
              State(component_id='preprocessing', component_property='value'),
              State(component_id='trim_spectra_lower_mass_range', component_property='value'),
              State(component_id='trim_spectra_upper_mass_range', component_property='value'),
              State(component_id='transform_intensity_method', component_property='value'),
              State(component_id='smooth_baseline_method', component_property='value'),
              State(component_id='smooth_baseline_window_length', component_property='value'),
              State(component_id='smooth_baseline_polyorder', component_property='value'),
              State(component_id='smooth_baseline_delta_mz', component_property='value'),
              State(component_id='smooth_baseline_diff_thresh', component_property='value'),
              State(component_id='remove_baseline_method', component_property='value'),
              State(component_id='normalize_intensity_method', component_property='value'),
              State(component_id='peak_picking_method', component_property='value'),
              State(component_id='peak_picking_widths_start', component_property='value'),
              State(component_id='peak_picking_widths_stop', component_property='value'),
              State(component_id='peak_picking_snr', component_property='value'),)
def spectrum(n_clicks,
             spectrum_id,
             preprocessing,
             trim_spectra_lower_mass_range,
             trim_spectra_upper_mass_range,
             transform_intensity_method,
             smooth_baseline_method,
             smooth_baseline_window_length,
             smooth_baseline_polyorder,
             smooth_baseline_delta_mz,
             smooth_baseline_diff_thresh,
             remove_baseline_method,
             normalize_intensity_method,
             peak_picking_method,
             peak_picking_widths_start,
             peak_picking_widths_stop,
             peak_picking_snr):
    changed_id = [i['prop_id'] for i in callback_context.triggered][0]

    if 'preprocess_submit' in changed_id:
        spectrum = indexed_dataset[spectrum_id]
        for step in preprocessing:
            '''if step == 'trim spectra':
                params = {'lower mass range': trim_spectra_lower_mass_range,
                          'upper mass range': trim_spectra_upper_mass_range}
            elif step == 'transform intensity':
                params = {'method': transform_intensity_method}
            elif step == 'smooth baseline':
                params = {'method': smooth_baseline_method,
                          'window length': smooth_baseline_window_length,
                          'polyorder': smooth_baseline_polyorder,
                          'delta mz': smooth_baseline_delta_mz,
                          'diff thresh': smooth_baseline_diff_thresh}
            elif step == 'remove baseline':
                params = {'method': remove_baseline_method}
            elif step == 'normalize intensity':
                params = {'method': normalize_intensity_method}
            elif step == 'peak picking':
                params = {'method': peak_picking_method,
                          'widths': np.arange(peak_picking_widths_start,
                                              peak_picking_widths_stop,
                                              (peak_picking_widths_stop - peak_picking_widths_start) / 10),
                          'snr': peak_picking_snr}'''
            if step == 'trim spectrum' and 'trim spectrum' not in spectrum.data_processing.keys():
                spectrum.trim_spectrum(lower_mass_range=int(trim_spectra_lower_mass_range),
                                       upper_mass_range=int(trim_spectra_upper_mass_range))
            elif step == 'transform intensity' and 'intensity transformation' not in spectrum.data_processing.keys():
                spectrum.transform_intensity(method=transform_intensity_method)
            elif step == 'smooth baseline' and 'baseline smoothing' not in spectrum.data_processing.keys():
                spectrum.smooth_baseline(method=smooth_baseline_method,
                                         window_length=int(smooth_baseline_window_length),
                                         polyorder=int(smooth_baseline_polyorder),
                                         delta_mz=float(smooth_baseline_delta_mz),
                                         diff_thresh=float(smooth_baseline_diff_thresh))
            elif step == 'remove baseline' and 'baseline removal' not in spectrum.data_processing.keys():
                spectrum.remove_baseline(method=remove_baseline_method)
            elif step == 'normalize intensity' and 'intensity normalization' not in spectrum.data_processing.keys():
                spectrum.normalize_intensity(method=normalize_intensity_method)
            elif step == 'peak picking' and 'peak picking' not in spectrum.data_processing.keys():
                if peak_picking_widths_start != 'None' and peak_picking_widths_stop != 'None':
                    peak_picking_widths = np.arange(float(peak_picking_widths_start),
                                                    float(peak_picking_widths_stop),
                                                    float((peak_picking_widths_stop - peak_picking_widths_start) / 10))
                else:
                    peak_picking_widths = None
                spectrum.peak_picking(method=peak_picking_method,
                                      widths=peak_picking_widths,
                                      snr=int(peak_picking_snr))

    spectrum_df = pd.DataFrame(data={'mz': spectrum.get_mz_array(),
                                     'intensity': spectrum.get_intensity_array()})

    fig = px.line(data_frame=spectrum_df,
                  x='mz',
                  y='intensity')
    # TODO: dynamic m/z axis range
    fig.update_xaxes(range=[2000, 20000])

    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
