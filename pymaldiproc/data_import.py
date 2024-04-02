import os
from pyteomics import mzml as pyt_mzml
from pyteomics import mzxml as pyt_mzxml
from pymaldiproc.classes import OpenMALDISpectrum, PMPTsfSpectrum, PMPTdfSpectrum
from pyTDFSDK.classes import TsfData, TdfData
from pyTDFSDK.init_tdf_sdk import init_tdf_sdk_api


def schema_detection(bruker_dot_d_file):
    """
    Detect the schema used by the raw data in the Bruker .d directory.

    :param bruker_dot_d_file: Path to the .d directory of interest.
    :type: str
    :return: Capitalized schema extension (TDF, TSF, or BAF).
    :rtype: str
    """
    exts = [os.path.splitext(fname)[1] for dirpath, dirnames, filenames in os.walk(bruker_dot_d_file)
            for fname in filenames]
    if '.tdf' in exts and '.tsf' not in exts and '.baf' not in exts:
        return 'TDF'
    elif '.tsf' in exts and '.tdf' not in exts and '.baf' not in exts:
        return 'TSF'
    elif '.baf' in exts and '.tdf' not in exts and '.tsf' not in exts:
        return 'BAF'


def import_timstof_raw_data(input_path, mode, profile_bins=0, encoding=64, exclude_mobility=True):
    # find Bruker .d directories
    if input_path.endswith('.d'):
        input_files = [input_path]
    elif not input_path.endswith('.d') and os.path.isdir(input_path):
        input_files = [os.path.join(dirpath, directory) for dirpath, dirnames, filenames in os.walk(input_path)
                       for directory in dirnames if directory.endswith('.d')]
    # read in data with pyTDFSDK
    list_of_spectra = []
    for dot_d_directory in input_files:
        print('Importing ' + dot_d_directory)
        if schema_detection(dot_d_directory) == 'TSF':
            data = TsfData(dot_d_directory, init_tdf_sdk_api())
            for frame in range(1, data.analysis['Frames'].shape[0] + 1):
                list_of_spectra.append(PMPTsfSpectrum(data, frame, mode, profile_bins, encoding))
        elif schema_detection(dot_d_directory) == 'TDF':
            data = TdfData(dot_d_directory, init_tdf_sdk_api())
            for frame in range(1, data.analysis['Frames'].shape[0] + 1):
                list_of_spectra.append(PMPTdfSpectrum(data, frame, mode, profile_bins=profile_bins, encoding=encoding,
                                                      exclude_mobility=exclude_mobility))
    return list_of_spectra


def import_mzml(input_path):
    # find mzML files
    if input_path.endswith('.mzML'):
        input_files = [input_path]
    elif not input_path.endswith('.mzML') and os.path.isdir(input_path):
        input_files = [os.path.join(dirpath, filename) for dirpath, dirnames, filenames in os.walk(input_path)
                       for filename in filenames if filename.endswith('.mzML')]
    # read in data with pyteomics
    list_of_spectra = []
    for mzml_filename in input_files:
        print('Importing ' + mzml_filename)
        mzml_data = list(pyt_mzml.read(mzml_filename))
        for scan_dict in mzml_data:
            list_of_spectra.append(OpenMALDISpectrum(scan_dict, mzml_filename))
    return list_of_spectra


def import_mzxml(input_path):
    # find mzML files
    if input_path.endswith('.mzXML'):
        input_files = [input_path]
    elif not input_path.endswith('.mzxML') and os.path.isdir(input_path):
        input_files = [os.path.join(dirpath, filename) for dirpath, dirnames, filenames in os.walk(input_path)
                       for filename in filenames if filename.endswith('.mzXML')]
    # read in data with pyteomics
    list_of_spectra = []
    for mzxml_filename in input_files:
        print('Importing ' + mzxml_filename)
        mzxml_data = list(pyt_mzxml.read(mzxml_filename))
        for scan_dict in mzxml_data:
            list_of_spectra.append(OpenMALDISpectrum(scan_dict, mzxml_filename))
    return list_of_spectra
