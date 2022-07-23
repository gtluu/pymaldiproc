from pymaldiproc.classes import *
from pymaldiproc.data_import import *
from pymaldiproc.preprocessing import *
import plotly.express as px
from dash import Dash, dcc, html, Input, Output


app = Dash(__name__)

app.layout = 