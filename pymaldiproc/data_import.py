import os
from pyteomics import mzml as pyt
from lxml.etree import parse, XMLParser
from pymaldiproc.classes import MALDISpectrum


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
        mzml_data = list(pyt.read(mzml_filename))
        for scan_dict in mzml_data:
            list_of_spectra.append(MALDISpectrum(scan_dict, mzml_filename))

    return list_of_spectra
