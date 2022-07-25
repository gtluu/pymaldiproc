from pymaldiproc.classes import *
from pymaldiproc.data_import import *
from pymaldiproc.preprocessing import *
import plotly.express as px
from dash import Dash, dcc, html, Input, Output
from pymaldiproc.classes import *

dataset = MALDIDataset('prototype')

app = Dash(__name__)

app.layout = html.Div([
    html.Div(dcc.Graph(id='spectrum', figure={}), className='row'),
    html.Div([
        html.Div(dcc.Dropdown(id='spectrum_id',
                              multi=True,
                              options=[{'label': i, 'value': i} for i in dataset.spectra_ids],
                              value=[i for i in dataset.spectra_ids]),
                 className='one column')
    ])
])


@app.callback(Output(component_id='spectrum', component_property='figure'),
              [Input(component_id='spectrum_id', component_property='value')],)
def graph_spectrum(value):
    # get spectrum
    spectrum_id = value[0].split('|')[0]
    spectrum_replicate = value[0].split('|')[1]
    spectrum = dataset.spectra[spectrum_id][spectrum_replicate]
    spectrum_df = pd.DataFrame(data={'mz': spectrum.raw_mz_array, 'intensity': spectrum.raw_intensity_array})

    fig = px.line(data_frame=spectrum_df,
                  x='mz',
                  y='intensity')
    fig.update_xaxes(range=[2000, 20000])
    return fig


if __name__ == '__main__':
    app.run_server()
