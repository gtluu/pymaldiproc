from contextlib import redirect_stdout
from io import StringIO
import webview
import webview.menu as wm
from pymaldiproc.dashboard import app
from pymaldiproc import VERSION


def load_bruker_d_files():
    pass


def load_mzml_files():
    pass


def load_mzxml_files():
    pass


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
