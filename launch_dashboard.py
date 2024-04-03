import os
import configparser
from contextlib import redirect_stdout
from io import StringIO
import webview
import webview.menu as wm
from pymaldiproc import VERSION
from pymaldiproc.data_import import import_timstof_raw_data, import_mzml, import_mzxml
from pymaldiproc.layout import get_dropdown
from pymaldiproc.dashboard import *
import tkinter
from tkinter.filedialog import askopenfilenames, askdirectory


# will be a dictionary of MALDISpectrum objects used for the spectrum plot
INDEXED_DATA = {}
# default processing parameters from config file
config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), 'etc', 'preprocessing.cfg'))
TRIM_SPECTRUM_PARAMS = {'lower_mass_range': int(config['trim_spectrum']['lower_mass_range']),
                        'upper_mass_range': int(config['trim_spectrum']['upper_mass_range'])}
TRANSFORM_INTENSITY_PARAMS = {'method': config['transform_intensity']['method']}
SMOOTH_BASELINE_PARAMS = {'method': config['smooth_baseline']['method'],
                          'window_length': int(config['smooth_baseline']['window_length']),
                          'polyorder': int(config['smooth_baseline']['polyorder']),
                          'delta_mz': float(config['smooth_baseline']['delta_mz']),
                          'diff_thresh': float(config['smooth_baseline']['diff_thresh'])}
REMOVE_BASELINE_PARAMS = {'method': config['remove_baseline']['method'],
                          'min_half_window': int(config['remove_baseline']['min_half_window']),
                          'max_half_window': int(config['remove_baseline']['max_half_window']),
                          'decreasing': config['remove_baseline'].getboolean('decreasing'),
                          'smooth_half_window': None,
                          'filter_order': int(config['remove_baseline']['filter_order']),
                          'sigma': None,
                          'increment': int(config['remove_baseline']['increment']),
                          'max_hits': int(config['remove_baseline']['max_hits']),
                          'window_tol': float(config['remove_baseline']['window_tol']),
                          'lambda_': int(config['remove_baseline']['lambda_']),
                          'porder': int(config['remove_baseline']['porder']),
                          'repetition': None,
                          'degree': int(config['remove_baseline']['degree']),
                          'gradient': float(config['remove_baseline']['gradient'])}
NORMALIZE_INTENSITY_PARAMS = {'method': config['normalize_intensity']['method']}
BIN_SPECTRUM_PARAMS = {'n_bins': int(config['bin_spectrum']['n_bins']),
                       'lower_mass_range': int(config['bin_spectrum']['lower_mass_range']),
                       'upper_mass_range': int(config['bin_spectrum']['upper_mass_range'])}
PEAK_PICKING_PARAMS = {'method': config['peak_picking']['method'],
                       'snr': int(config['peak_picking']['snr']),
                       'widths': None}
if config['remove_baseline']['smooth_half_window'] != 'None':
    REMOVE_BASELINE_PARAMS['smooth_half_window'] = int(config['remove_baseline']['smooth_half_window'])
if config['remove_baseline']['sigma'] != 'None':
    REMOVE_BASELINE_PARAMS['sigma'] = float(config['rmeove_baseline']['sigma'])
if config['remove_baseline']['repetition'] != 'None':
    REMOVE_BASELINE_PARAMS['repetition'] = int(config['remove_baseline']['repetition'])
if config['peak_picking']['widths'] != 'None':
    PEAK_PICKING_PARAMS['widths'] = int(config['peak_picking']['widths'])
PREPROCESSING_PARAMS = {'TRIM_SPECTRUM': TRIM_SPECTRUM_PARAMS,
                        'TRANSFORM_INTENSITY': TRANSFORM_INTENSITY_PARAMS,
                        'SMOOTH_BASELINE': SMOOTH_BASELINE_PARAMS,
                        'REMOVE_BASELINE': REMOVE_BASELINE_PARAMS,
                        'NORMALIZE_INTENSITY': NORMALIZE_INTENSITY_PARAMS,
                        'BIN_SPECTRUM': BIN_SPECTRUM_PARAMS,
                        'PEAK_PICKING': PEAK_PICKING_PARAMS}


SPECTRA_LIST = app.layout.children[0]


def get_dash_element(element, id_query):
    print(element)
    try:
        if not isinstance(element, list):
            if element is not None:
                if element.id == id_query:
                    return element
        elif element != [] and element is not None:
            if isinstance(element, list):
                for child in element:
                    if child.id == id_query:
                        return child
                    elif not isinstance(child.children, str):
                        get_dash_element(child.children, id_query)
            elif isinstance(element, html.Div):
                get_dash_element(element, id_query)
    except AttributeError:
        print(f'Skipping {element}')


def load_bruker_d_files():
    global INDEXED_DATA
    main_tk_window = tkinter.Tk()
    main_tk_window.attributes('-topmost', True, '-alpha', 0)
    dirname = askdirectory(mustexist=True)
    main_tk_window.destroy()
    print(dirname)
    if dirname.endswith('.d'):
        data = import_timstof_raw_data(dirname, mode='profile')
        for spectrum in data:
            INDEXED_DATA[spectrum.spectrum_id] = spectrum
    print(INDEXED_DATA)
    print(SPECTRA_LIST.children)
    print(get_dropdown(INDEXED_DATA))
    app.layout.children[0].children = get_dropdown(INDEXED_DATA)
    webview.active_window().load_url(app.server)


def load_mzml_files():
    global INDEXED_DATA
    main_tk_window = tkinter.Tk()
    main_tk_window.attributes('-topmost', True, '-alpha', 0)
    filenames = askopenfilenames(filetypes=[('mzML Files', '*.mzML')])
    main_tk_window.destroy()
    print(filenames)
    for filename in filenames:
        data = import_mzml(filename)
        for spectrum in data:
            INDEXED_DATA[spectrum.spectrum_id] = spectrum
    print(INDEXED_DATA)
    print(SPECTRA_LIST.children)
    print(get_dropdown(INDEXED_DATA))
    SPECTRA_LIST.children = get_dropdown(INDEXED_DATA)


def load_mzxml_files():
    global INDEXED_DATA
    main_tk_window = tkinter.Tk()
    main_tk_window.attributes('-topmost', True, '-alpha', 0)
    filenames = askopenfilenames(filetypes=[('mzXML Files', '*.mzXML')])
    main_tk_window.destroy()
    print(filenames)
    for filename in filenames:
        data = import_mzxml(filename)
        for spectrum in data:
            INDEXED_DATA[spectrum.spectrum_id] = spectrum
    print(INDEXED_DATA)
    print(SPECTRA_LIST.children)
    print(get_dropdown(INDEXED_DATA))
    SPECTRA_LIST.children = get_dropdown(INDEXED_DATA)


def export_spectrum():
    pass


def export_peak_list():
    pass


def edit_preprocessing_parameters():
    pass


def trim_spectrum():
    pass


def transform_intensity():
    pass


def smooth_baseline():
    pass


def remove_baseline():
    pass


def normalize_intensity():
    pass


def bin_spectrum():
    pass


def align_spectra():
    pass


def label_peaks():
    pass


def undo_preprocessing():
    pass


def undo_peak_labeling():
    pass


def get_pymaldiproc_menu():
    return [
        wm.Menu(
            'File',
            [
                wm.MenuAction('Load Bruker *.d File', load_bruker_d_files),
                wm.MenuAction('Load mzML File(s)', load_mzml_files),
                wm.MenuAction('Load mzXML File(s)', load_mzxml_files),
                wm.MenuSeparator(),
                wm.Menu(
                    'Menu',
                    [
                        wm.MenuAction('Export Current Spectrum to mzML File', export_spectrum),
                        wm.MenuAction('Export Current Peak List to CSV File', export_peak_list)
                    ]
                )
            ],
        ),
        wm.Menu(
            'Preprocessing',
            [
                wm.MenuAction('Edit Preprocessing Parameters', edit_preprocessing_parameters),
                wm.MenuSeparator(),
                wm.MenuAction('Trim Spectrum', trim_spectrum),
                wm.MenuAction('Transform Intensity', transform_intensity),
                wm.MenuAction('Smooth Baseline', smooth_baseline),
                wm.MenuAction('Remove Baseline', remove_baseline),
                wm.MenuAction('Normalize Intensity', normalize_intensity),
                wm.MenuAction('Bin Spectrum', bin_spectrum),
                wm.MenuAction('Align Spectra', align_spectra),
                wm.MenuAction('Label Peaks', label_peaks),
                wm.MenuSeparator(),
                wm.MenuAction('Undo Preprocessing', undo_preprocessing),
                wm.MenuAction('Undo Peak Labeling', undo_peak_labeling)
            ]
        )
    ]


if __name__ == '__main__':
    stream = StringIO()
    with redirect_stdout(stream):
        window = webview.create_window(f'pyMALDIproc Dashboard {VERSION}', app.server)
        webview.start(menu=get_pymaldiproc_menu())
