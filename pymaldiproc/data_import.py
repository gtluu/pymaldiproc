import os
from pyteomics import mzml as pyt_mzml
from pyteomics import mzxml as pyt_mzxml
from pymaldiproc.classes import OpenMALDISpectrum, PMPTsfSpectrum, PMP2DTdfSpectrum, PMP3DTdfSpectrum
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


def import_timstof_raw_data(input_path, mode, profile_bins=0, encoding=64, exclude_mobility=False):
    """
    Import spectra from Bruker TSF and TDF files. Due to incompatibility with current preprocessing methods, ion
    mobility data is not parsed by default.

    :param input_path: Path to the directory to be searched. If the path is a Bruker .d directory, spectra from that
        dataset will be imported. If the path is a directory containing multiple Bruker TSF/TDF files, all data will be
        loaded.
    :type input_path: str
    :param mode: Data array mode, either "profile", "centroid", or "raw".
    :type mode: str
    :param profile_bins: Number of bins to bin spectrum to.
    :type profile_bins: int
    :param encoding: Encoding bit mode, either "64" or "32"
    :type encoding: int
    :param exclude_mobility: Whether to include mobility data in the output files, defaults to True.
    :type exclude_mobility: bool | None
    :return: List of spectra.
    :rtype: list[pymaldiproc.classes.PMPTsfSpectrum|pymaldiproc.classes.PMP2DTdfSpectrum]
    """
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
            if exclude_mobility:
                for frame in range(1, data.analysis['Frames'].shape[0] + 1):
                    list_of_spectra.append(PMP2DTdfSpectrum(data, frame, mode, profile_bins=profile_bins,
                                                            encoding=encoding))
            elif not exclude_mobility:
                for frame in range(1, data.analysis['Frames'].shape[0] + 1):
                    list_of_spectra.append(PMP3DTdfSpectrum(data, frame, mode='centroid', profile_bins=profile_bins,
                                                            encoding=encoding))
    return list_of_spectra


def import_mzml(input_path):
    """
    Import spectra from mzML files.

    :param input_path: Path to the directory to be searched. If the path is an mzML file, spectra from that dataset
        will be imported. If the path is a directory containing multiple mzML files, all data will be loaded.
    :type input_path: str
    :return: List of spectra.
    :rtype: list[pymaldiproc.classes.OpenMALDISpectrum]
    """
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
    """
    Import spectra from mzXML files.

    :param input_path: Path to the directory to be searched. If the path is an mzXML file, spectra from that dataset
        will be imported. If the path is a directory containing multiple mzXML files, all data will be loaded.
    :type input_path: str
    :return: List of spectra.
    :rtype: list[pymaldiproc.classes.OpenMALDISpectrum]
    """
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
