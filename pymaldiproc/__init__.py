import os
import copy
import gc
import tempfile
import configparser
import numpy as np
import pandas as pd
from uuid import uuid4
from functools import reduce
from contextlib import redirect_stdout
from io import StringIO

from pyteomics import mzml as pyt_mzml
from pyteomics import mzxml as pyt_mzxml
from pyTDFSDK.classes import TsfData, TdfData, TsfSpectrum, TdfSpectrum
from pyTDFSDK.init_tdf_sdk import init_tdf_sdk_api

from scipy.signal import savgol_filter, find_peaks, find_peaks_cwt, peak_widths
from pyMSpec.smoothing import apodization, rebin, fast_change, median
from pybaselines.smooth import snip, noise_median
from pybaselines.morphological import tophat
from BaselineRemoval import BaselineRemoval
from pyMSpec.normalisation import tic, rms, mad, sqrt
from icoshift import icoshift
from pyopenms import *

import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
from plotly_resampler import FigureResampler

from dash import State, callback_context, no_update, dcc, html
from dash_extensions.enrich import (Input, Output, DashProxy, MultiplexerTransform, Serverside,
                                    ServersideOutputTransform, FileSystemBackend)
import dash_bootstrap_components as dbc
import webview

import tkinter
from tkinter.filedialog import askopenfilenames, askdirectory, asksaveasfilename

from pymaldiproc.classes import *
from pymaldiproc.data_import import *
from pymaldiproc.preprocessing import *

VERSION = '0.3.0a13'
