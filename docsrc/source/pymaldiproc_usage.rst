pyMALDIproc Usage
=================
pyMALDIproc provides functions and classes to easily read in MALDI datasets. Furthermore, all data arrays are stored as
``numpy`` arrays, making them easy to work with and cross compatible with other Python libraries.

Importing Data from mzML Files
------------------------------

.. code-block::

    data = import_mzml('brevi_brachy.mzML')

Importing Data from mzXML Files
-------------------------------

.. code-block::

    data = import_mzxml('brevi_brachy.mzXML')

Importing Data from Raw Bruker *.d Files (TSF Format Files)
-----------------------------------------------------------
For more information on TSF files, see
`here <https://gtluu.github.io/timsconvert/introduction.html#timstof-file-formats>`_. TIMS on datasets are not
currently supported but may be in the future.

.. code-block::

    data = import_timstof_raw_data('brevi_brachy.d', mode='profile')

All data is imported as an ``OpenMALDISpectrum``, ``PMPTsfSpectrum``, or ``PMPTdfSpectrum``.

.. code-block::

    for key, value in data[0].info().items():
        print(f'{key}: {value}')

.. code-block::

    Name: brevi_brachy_0
    Coordinate: None
    MS Level: 1
    Mode: Profile
    Source File: [path to]\brevi_brachy.mzML
    UUID: a15ae270-004c-4593-a048-56921a682b6c

Data Preprocessing and Peak Picking
-----------------------------------
All data preprocessing can be performed using the appropriate methods from the ``OpenMALDISpectrum``,
``PMPTsfSpectrum``, or ``PMPTdfSpectrum`` classes. Refer to the methods in the ``PMPMethods`` class for more details.

.. code-block::

    for spectrum in data:
        spectrum.trim_spectrum(lower_mass_range=100, upper_mass_range=3000)
        spectrum.transform_intensity(method='sqrt')
        spectrum.smooth_baseline(method='SavitzkyGolay')
        spectrum.remove_baseline(method='SNIP')
        spectrum.normalize_intensity(method='tic')
        spectrum.bin_spectrum(n_bins=8000, lower_mass_range=100, upper_mass_range=3000)
        spectrum.peak_picking(method='locmax', snr=4, deisotope=True)

Feature Matrix
--------------
To get a feature matrix containing all peak lists for each spectrum in the dataset, simply run the
``get_feature_matrix`` function. This will return a ``pandas.DataFrame`` with an consensus m/z column and intensity
columns for each spectrum.

.. code-block::

    feature_matrix = get_feature_matrix(data, tolerance=0.05, decimals=4, missing_value_imputation=True)

To export the feature matrix, a simple wrapper around ``pandas.DataFrame.to_csv()`` method is included.

.. code-block::

    export_feature_matrix(feature_matrix, 'feature_matrix.csv')
