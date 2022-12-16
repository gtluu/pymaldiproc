from pymaldiproc.classes import *
from pymaldiproc.data_import import *
from pymaldiproc.preprocessing import *
import plotly.express as px
from dash import Dash, Input, Output, dcc, html, State, callback_context
import base64

INDEXED_DATA = {}
UPLOAD_DIR = 'data'
if not os.path.exists(UPLOAD_DIR):
    os.makedirs(UPLOAD_DIR)

app = Dash(__name__)

app.layout = html.Div([
    html.Div(
        dcc.Upload(
            id='upload',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select mzML Files')
            ]),
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
        id='spectrum',
        className='row'
    ),

    html.Div([
        html.Div(id='dropdown', className='one column')
    ])
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
        spectrum = INDEXED_DATA[value]
        spectrum_df = pd.DataFrame(data={'m/z': spectrum.get_mz_array(),
                                         'Intensity': spectrum.get_intensity_array()})

        fig = px.line(data_frame=spectrum_df,
                      x='m/z',
                      y='Intensity')
        fig.update_layout(xaxis_tickformat='d',
                          yaxis_tickformat='~e')

        children = [
            dcc.Graph(
                id='spectrum_plot',
                figure=fig,
                style={
                    'width': '100%',
                    'height': '800px'
                }
            )
        ]

        return children


if __name__ == '__main__':
    app.run_server(debug=True)
