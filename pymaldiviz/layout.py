from dash import dcc, html
import dash_bootstrap_components as dbc
from pymaldiviz import VERSION
from pymaldiviz.util import blank_figure, get_preprocessing_params


def get_preprocessing_parameters_layout(param_dict):
    """
    Obtain the layout for the preprocessing parameters modal window body.

    :param param_dict: Dictionary of parameters used to populate default values.
    :return: List of divs containing the layout for the preprocessing parameters modal window.
    """
    trim_spectrum_parameters = html.Div(
        [
            html.H5('Spectrum Trimming Parameters'),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Lower Mass Range'),
                    dbc.Input(id='trim_spectrum_lower_mass_range_value',
                              placeholder=param_dict['TRIM_SPECTRUM']['lower_mass_range'],
                              value=param_dict['TRIM_SPECTRUM']['lower_mass_range'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='trim_spectrum_lower_mass_range',
                style={'margin': '10px'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Upper Mass Range'),
                    dbc.Input(id='trim_spectrum_upper_mass_range_value',
                              placeholder=param_dict['TRIM_SPECTRUM']['upper_mass_range'],
                              value=param_dict['TRIM_SPECTRUM']['upper_mass_range'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='trim_spectrum_upper_mass_range',
                style={'margin': '10px'}
            )
        ],
        id='trim_spectrum_parameters',
        style={'margin': '20px'}
    )

    transform_intensity_parameters = html.Div(
        [
            html.H5('Intensity Transformation Parameters'),
            html.P('Method'),
            dbc.RadioItems(
                id='transform_intensity_method',
                options=[
                    {'label': 'Square Root', 'value': 'sqrt'},
                    {'label': 'Natural Log', 'value': 'log'},
                    {'label': 'Log Base 2', 'value': 'log2'},
                    {'label': 'Log Base 10', 'value': 'log10'}
                ],
                value=param_dict['TRANSFORM_INTENSITY']['method'],
                labelStyle={'display': 'inline-block', 'marginRight': '20px'},
                inputStyle={'margin-right': '10px'},
                className='btn-group',
                inputClassName='btn-check',
                labelClassName='btn btn-outline-primary',
                labelCheckedClassName='active',
            )
        ],
        id='transform_intensity_parameters',
        style={'margin': '20px'}
    )

    smooth_baseline_parameters = html.Div(
        [
            html.H5('Baseline Smoothing Parameters'),
            html.P('Method'),
            dbc.RadioItems(
                id='smooth_baseline_method',
                options=[
                    {'label': 'Savitzky-Golay', 'value': 'SavitzkyGolay'},
                    {'label': 'Apodization', 'value': 'apodization'},
                    {'label': 'Rebin', 'value': 'rebin'},
                    {'label': 'Fast Change', 'value': 'fast_change'},
                    {'label': 'Median', 'value': 'median'}
                ],
                value=param_dict['SMOOTH_BASELINE']['method'],
                labelStyle={'display': 'inline-block', 'marginRight': '20px'},
                inputStyle={'margin-right': '10px'},
                className='btn-group',
                inputClassName='btn-check',
                labelClassName='btn btn-outline-primary',
                labelCheckedClassName='active',
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Window Length'),
                    dbc.Input(id='smooth_baseline_window_length_value',
                              placeholder=param_dict['SMOOTH_BASELINE']['window_length'],
                              value=param_dict['SMOOTH_BASELINE']['window_length'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='smooth_baseline_window_length',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Polyorder'),
                    dbc.Input(id='smooth_baseline_polyorder_value',
                              placeholder=param_dict['SMOOTH_BASELINE']['polyorder'],
                              value=param_dict['SMOOTH_BASELINE']['polyorder'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='smooth_baseline_polyorder',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Delta m/z'),
                    dbc.Input(id='smooth_baseline_delta_mz_value',
                              placeholder=param_dict['SMOOTH_BASELINE']['delta_mz'],
                              value=param_dict['SMOOTH_BASELINE']['delta_mz'],
                              type='number',
                              min=0,
                              step=0.001)
                ],
                id='smooth_baseline_delta_mz',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Difference Threshold'),
                    dbc.Input(id='smooth_baseline_diff_thresh_value',
                              placeholder=param_dict['SMOOTH_BASELINE']['diff_thresh'],
                              value=param_dict['SMOOTH_BASELINE']['diff_thresh'],
                              type='number',
                              min=0,
                              step=0.001)
                ],
                id='smooth_baseline_diff_thresh',
                style={'margin': '10px',
                       'display': 'flex'}
            )
        ],
        id='smooth_baseline_parameters',
        style={'margin': '20px'}
    )

    remove_baseline_parameters = html.Div(
        [
            html.H5('Baseline Removal Parameters'),
            html.P('Method'),
            dbc.RadioItems(
                id='remove_baseline_method',
                options=[
                    {'label': 'SNIP', 'value': 'SNIP'},
                    {'label': 'Top Hat', 'value': 'TopHat'},
                    {'label': 'Median', 'value': 'Median'},
                    {'label': 'Zhang Fit', 'value': 'ZhangFit'},
                    {'label': 'ModPoly', 'value': 'ModPoly'},
                    {'label': 'IModPoly', 'value': 'IModPoly'}
                ],
                value=param_dict['REMOVE_BASELINE']['method'],
                labelStyle={'display': 'inline-block', 'marginRight': '20px'},
                inputStyle={'margin-right': '10px'},
                className='btn-group',
                inputClassName='btn-check',
                labelClassName='btn btn-outline-primary',
                labelCheckedClassName='active',
            ),
            dbc.RadioItems(
                id='remove_baseline_decreasing',
                options=[
                    {'label': 'Use Decreasing Iterative Window Sizes', 'value': True},
                    {'label': 'Do Not Use Decreasing Iterative Window Sizes', 'value': False}
                ],
                value=param_dict['REMOVE_BASELINE']['decreasing'],
                labelStyle={'display': 'inline-block', 'marginRight': '20px'},
                inputStyle={'margin-right': '10px'},
                className='btn-group',
                inputClassName='btn-check',
                labelClassName='btn btn-outline-primary',
                labelCheckedClassName='active',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Minimum Half Window'),
                    dbc.Input(id='remove_baseline_min_half_window_value',
                              placeholder=param_dict['REMOVE_BASELINE']['max_half_window'],
                              value=param_dict['REMOVE_BASELINE']['max_half_window'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='remove_baseline_min_half_window',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Maximum Half Window'),
                    dbc.Input(id='remove_baseline_max_half_window_value',
                              placeholder=param_dict['REMOVE_BASELINE']['max_half_window'],
                              value=param_dict['REMOVE_BASELINE']['max_half_window'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='remove_baseline_max_half_window',

                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Smoothing Half Window'),
                    dbc.Input(id='remove_baseline_smooth_half_window_value',
                              placeholder=param_dict['REMOVE_BASELINE']['smooth_half_window'],
                              value=param_dict['REMOVE_BASELINE']['smooth_half_window'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='remove_baseline_smooth_half_window',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Filter Order'),
                    dbc.Input(id='remove_baseline_filter_order_value',
                              placeholder=param_dict['REMOVE_BASELINE']['filter_order'],
                              value=param_dict['REMOVE_BASELINE']['filter_order'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='remove_baseline_filter_order',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Sigma'),
                    dbc.Input(id='remove_baseline_sigma_value',
                              placeholder=param_dict['REMOVE_BASELINE']['sigma'],
                              value=param_dict['REMOVE_BASELINE']['sigma'],
                              type='number',
                              min=0,
                              step=0.001)
                ],
                id='remove_baseline_sigma',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Increment'),
                    dbc.Input(id='remove_baseline_increment_value',
                              placeholder=param_dict['REMOVE_BASELINE']['increment'],
                              value=param_dict['REMOVE_BASELINE']['increment'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='remove_baseline_increment',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Max Hits'),
                    dbc.Input(id='remove_baseline_max_hits_value',
                              placeholder=param_dict['REMOVE_BASELINE']['max_hits'],
                              value=param_dict['REMOVE_BASELINE']['max_hits'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='remove_baseline_max_hits',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Window Tolerance'),
                    dbc.Input(id='remove_baseline_window_tol_value',
                              placeholder=param_dict['REMOVE_BASELINE']['window_tol'],
                              value=param_dict['REMOVE_BASELINE']['window_tol'],
                              type='number',
                              min=0,
                              step=0.000001)
                ],
                id='remove_baseline_window_tol',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Lambda'),
                    dbc.Input(id='remove_baseline_lambda__value',
                              placeholder=param_dict['REMOVE_BASELINE']['lambda_'],
                              value=param_dict['REMOVE_BASELINE']['lambda_'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='remove_baseline_lambda_',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('p-order'),
                    dbc.Input(id='remove_baseline_porder_value',
                              placeholder=param_dict['REMOVE_BASELINE']['porder'],
                              value=param_dict['REMOVE_BASELINE']['porder'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='remove_baseline_porder',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Repetition'),
                    dbc.Input(id='remove_baseline_repetition_value',
                              placeholder=param_dict['REMOVE_BASELINE']['repetition'],
                              value=param_dict['REMOVE_BASELINE']['repetition'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='remove_baseline_repetition',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Degree'),
                    dbc.Input(id='remove_baseline_degree_value',
                              placeholder=param_dict['REMOVE_BASELINE']['degree'],
                              value=param_dict['REMOVE_BASELINE']['degree'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='remove_baseline_degree',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Gradient'),
                    dbc.Input(id='remove_baseline_gradient_value',
                              placeholder=param_dict['REMOVE_BASELINE']['gradient'],
                              value=param_dict['REMOVE_BASELINE']['gradient'],
                              type='number',
                              min=0,
                              step=0.001)
                ],
                id='remove_baseline_gradient',
                style={'margin': '10px',
                       'display': 'flex'}
            )
        ],
        id='remove_baseline_parameters',
        style={'margin': '20px'}
    )

    normalize_intensity_parameters = html.Div(
        [
            html.H5('Intensity Normalization Parameters'),
            html.P('Method'),
            dbc.RadioItems(
                id='normalize_intensity_method',
                options=[
                    {'label': 'Total Ion Count', 'value': 'tic'},
                    {'label': 'Root Mean Square', 'value': 'rms'},
                    {'label': 'Mean Absolute Deviation', 'value': 'mad'},
                    {'label': 'Square Root', 'value': 'sqrt'}
                ],
                value=param_dict['NORMALIZE_INTENSITY']['method'],
                labelStyle={'display': 'inline-block', 'marginRight': '20px'},
                inputStyle={'margin-right': '10px'},
                className='btn-group',
                inputClassName='btn-check',
                labelClassName='btn btn-outline-primary',
                labelCheckedClassName='active',
            )
        ],
        id='normalize_intensity_parameters',
        style={'margin': '20px'}
    )

    bin_spectrum_parameters = html.Div(
        [
            html.H5('Spectrum Binning Parameters'),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Number of Bins'),
                    dbc.Input(id='bin_spectrum_n_bins_value',
                              placeholder=param_dict['BIN_SPECTRUM']['n_bins'],
                              value=param_dict['BIN_SPECTRUM']['n_bins'],
                              type='number',
                              min=0,
                              step=100)
                ],
                id='bin_spectrum_n_bins',
                style={'margin': '10px'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Lower Mass Range'),
                    dbc.Input(id='bin_spectrum_lower_mass_range_value',
                              placeholder=param_dict['BIN_SPECTRUM']['lower_mass_range'],
                              value=param_dict['BIN_SPECTRUM']['lower_mass_range'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='bin_spectrum_lower_mass_range',
                style={'margin': '10px'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Upper Mass Range'),
                    dbc.Input(id='bin_spectrum_upper_mass_range_value',
                              placeholder=param_dict['BIN_SPECTRUM']['upper_mass_range'],
                              value=param_dict['BIN_SPECTRUM']['upper_mass_range'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='bin_spectrum_upper_mass_range',
                style={'margin': '10px'}
            )
        ],
        id='bin_spectrum_parameters',
        style={'margin': '20px'}
    )

    peak_picking_parameters = html.Div(
        [
            html.H5('2D Peak Picking Parameters'),
            html.P('Method'),
            dbc.RadioItems(
                id='peak_picking_method',
                options=[
                    {'label': 'Local Maxima', 'value': 'locmax'},
                    {'label': 'Continuous Wavelet Transform', 'value': 'cwt'}
                ],
                value=param_dict['PEAK_PICKING']['method'],
                labelStyle={'display': 'inline-block', 'marginRight': '20px'},
                inputStyle={'margin-right': '10px'},
                className='btn-group',
                inputClassName='btn-check',
                labelClassName='btn btn-outline-primary',
                labelCheckedClassName='active',
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Signal-to-Noise Ratio'),
                    dbc.Input(id='peak_picking_snr_value',
                              placeholder=param_dict['PEAK_PICKING']['snr'],
                              value=param_dict['PEAK_PICKING']['snr'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='peak_picking_snr',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Widths (Optional)'),
                    dbc.Input(id='peak_picking_widths_value',
                              placeholder=param_dict['PEAK_PICKING']['widths'],
                              value=param_dict['PEAK_PICKING']['widths'],
                              type='number',
                              min=0,
                              step=0.001)
                ],
                id='peak_picking_widths',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.RadioItems(
                id='peak_picking_deisotope',
                options=[
                    {'label': 'Deisotope Peak List', 'value': True},
                    {'label': 'Do Not Deisotope', 'value': False}
                ],
                value=param_dict['PEAK_PICKING']['deisotope'],
                labelStyle={'display': 'inline-block', 'marginRight': '20px'},
                inputStyle={'margin-right': '10px'},
                className='btn-group',
                inputClassName='btn-check',
                labelClassName='btn btn-outline-primary',
                labelCheckedClassName='active',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Deisotoping Fragment Tolerance'),
                    dbc.Input(id='peak_picking_deisotope_fragment_tolerance_value',
                              placeholder=param_dict['PEAK_PICKING']['fragment_tolerance'],
                              value=param_dict['PEAK_PICKING']['fragment_tolerance'],
                              type='number',
                              min=0,
                              step=0.0001)
                ],
                id='peak_picking_deisotope_fragment_tolerance',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            html.P('Deisotoping Fragment_Tolerance Unit',
                   id='peak_picking_deisotope_fragment_unit_ppm_label',
                   style={'margin': '10px',
                          'display': 'flex'}),
            dbc.RadioItems(
                id='peak_picking_deisotope_fragment_unit_ppm',
                options=[
                    {'label': 'PPM', 'value': True},
                    {'label': 'Da', 'value': False}
                ],
                value=param_dict['PEAK_PICKING']['fragment_unit_ppm'],
                labelStyle={'display': 'inline-block', 'marginRight': '20px'},
                inputStyle={'margin-right': '10px'},
                className='btn-group',
                inputClassName='btn-check',
                labelClassName='btn btn-outline-primary',
                labelCheckedClassName='active',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Deisotoping Minimum Charge'),
                    dbc.Input(id='peak_picking_deisotope_min_charge_value',
                              placeholder=param_dict['PEAK_PICKING']['min_charge'],
                              value=param_dict['PEAK_PICKING']['min_charge'],
                              type='number',
                              min=1,
                              step=1)
                ],
                id='peak_picking_deisotope_min_charge',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Deisotoping Maximum Charge'),
                    dbc.Input(id='peak_picking_deisotope_max_charge_value',
                              placeholder=param_dict['PEAK_PICKING']['max_charge'],
                              value=param_dict['PEAK_PICKING']['max_charge'],
                              type='number',
                              min=1,
                              step=1)
                ],
                id='peak_picking_deisotope_max_charge',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.RadioItems(
                id='peak_picking_deisotope_keep_only_deisotoped',
                options=[
                    {'label': 'Retain Only Deisotoped Peaks', 'value': True},
                    {'label': 'Retain Monoisotopic and All Other Peaks', 'value': False}
                ],
                value=param_dict['PEAK_PICKING']['keep_only_deisotoped'],
                labelStyle={'display': 'inline-block', 'marginRight': '20px'},
                inputStyle={'margin-right': '10px'},
                className='btn-group',
                inputClassName='btn-check',
                labelClassName='btn btn-outline-primary',
                labelCheckedClassName='active',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Deisotoping Minimum Number of Isotopic Peaks'),
                    dbc.Input(id='peak_picking_deisotope_min_isopeaks_value',
                              placeholder=param_dict['PEAK_PICKING']['min_isopeaks'],
                              value=param_dict['PEAK_PICKING']['min_isopeaks'],
                              type='number',
                              min=2,
                              step=1)
                ],
                id='peak_picking_deisotope_min_isopeaks',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Deisotoping Maximum Number of Isotopic Peaks'),
                    dbc.Input(id='peak_picking_deisotope_max_isopeaks_value',
                              placeholder=param_dict['PEAK_PICKING']['max_isopeaks'],
                              value=param_dict['PEAK_PICKING']['max_isopeaks'],
                              type='number',
                              min=2,
                              step=1)
                ],
                id='peak_picking_deisotope_max_isopeaks',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.RadioItems(
                id='peak_picking_deisotope_make_single_charged',
                options=[
                    {'label': 'Convert Deisotoped Monoisotopic Peak to Single Charge', 'value': True},
                    {'label': 'Retain Original Charge', 'value': False}
                ],
                value=param_dict['PEAK_PICKING']['make_single_charged'],
                labelStyle={'display': 'inline-block', 'marginRight': '20px'},
                inputStyle={'margin-right': '10px'},
                className='btn-group',
                inputClassName='btn-check',
                labelClassName='btn btn-outline-primary',
                labelCheckedClassName='active',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.RadioItems(
                id='peak_picking_deisotope_annotate_charge',
                options=[
                    {'label': 'Annotate Charge', 'value': True},
                    {'label': 'Do Not Annotate Charge', 'value': False}
                ],
                value=param_dict['PEAK_PICKING']['annotate_charge'],
                labelStyle={'display': 'inline-block', 'marginRight': '20px'},
                inputStyle={'margin-right': '10px'},
                className='btn-group',
                inputClassName='btn-check',
                labelClassName='btn btn-outline-primary',
                labelCheckedClassName='active',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.RadioItems(
                id='peak_picking_deisotope_annotate_iso_peak_count',
                options=[
                    {'label': 'Annotate Number of Isotopic Peaks', 'value': True},
                    {'label': 'Do Not Annotate Number of Isotopic Peaks', 'value': False}
                ],
                value=param_dict['PEAK_PICKING']['annotate_iso_peak_count'],
                labelStyle={'display': 'inline-block', 'marginRight': '20px'},
                inputStyle={'margin-right': '10px'},
                className='btn-group',
                inputClassName='btn-check',
                labelClassName='btn btn-outline-primary',
                labelCheckedClassName='active',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.RadioItems(
                id='peak_picking_deisotope_use_decreasing_model',
                options=[
                    {'label': 'Use Decreasing/Averagine Model', 'value': True},
                    {'label': 'Do Not Perform Peak Intensity Check', 'value': False}
                ],
                value=param_dict['PEAK_PICKING']['use_decreasing_model'],
                labelStyle={'display': 'inline-block', 'marginRight': '20px'},
                inputStyle={'margin-right': '10px'},
                className='btn-group',
                inputClassName='btn-check',
                labelClassName='btn btn-outline-primary',
                labelCheckedClassName='active',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Deisotoping Intensity Check Starting Peak'),
                    dbc.Input(id='peak_picking_deisotope_start_intensity_check_value',
                              placeholder=param_dict['PEAK_PICKING']['start_intensity_check'],
                              value=param_dict['PEAK_PICKING']['start_intensity_check'],
                              type='number',
                              min=1,
                              step=1)
                ],
                id='peak_picking_deisotope_start_intensity_check',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.RadioItems(
                id='peak_picking_deisotope_add_up_intensity',
                options=[
                    {'label': 'Sum Isotopic Pattern Intensities Into Monoisotopic Peak', 'value': True},
                    {'label': 'Do Not Sum Intensities', 'value': False}
                ],
                value=param_dict['PEAK_PICKING']['add_up_intensity'],
                labelStyle={'display': 'inline-block', 'marginRight': '20px'},
                inputStyle={'margin-right': '10px'},
                className='btn-group',
                inputClassName='btn-check',
                labelClassName='btn btn-outline-primary',
                labelCheckedClassName='active',
                style={'margin': '10px',
                       'display': 'flex'}
            )
        ],
        id='peak_picking_parameters',
        style={'margin': '20px'}
    )

    peak_picking_3d_parameters = html.Div(
        [
            html.H5('3D Peak Picking Parameters'),
            html.P('Method'),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Minimum Distance'),
                    dbc.Input(id='peak_picking_3d_min_distance_value',
                              placeholder=param_dict['PEAK_PICKING_3D']['min_distance'],
                              value=param_dict['PEAK_PICKING_3D']['min_distance'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='peak_picking_3d_min_distance',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Noise (Optional)'),
                    dbc.Input(id='peak_picking_3d_noise_value',
                              placeholder=param_dict['PEAK_PICKING_3D']['noise'],
                              value=param_dict['PEAK_PICKING_3D']['noise'],
                              type='number',
                              min=0,
                              step=0.1)
                ],
                id='peak_picking_3d_noise',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Signal-to-Noise Ratio'),
                    dbc.Input(id='peak_picking_3d_snr_value',
                              placeholder=param_dict['PEAK_PICKING_3D']['snr'],
                              value=param_dict['PEAK_PICKING_3D']['snr'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='peak_picking_3d_snr',
                style={'margin': '10px',
                       'display': 'flex'}
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupText('Exclude Border'),
                    dbc.Input(id='peak_picking_3d_exclude_border_value',
                              placeholder=param_dict['PEAK_PICKING_3D']['exclude_border'],
                              value=param_dict['PEAK_PICKING_3D']['exclude_border'],
                              type='number',
                              min=0,
                              step=1)
                ],
                id='peak_picking_3d_exclude_border',
                style={'margin': '10px',
                       'display': 'flex'}
            )
        ],
        id='peak_picking_3d_parameters',
        style={'margin': '20px'}
    )

    return [trim_spectrum_parameters,
            transform_intensity_parameters,
            smooth_baseline_parameters,
            remove_baseline_parameters,
            normalize_intensity_parameters,
            bin_spectrum_parameters,
            peak_picking_parameters,
            peak_picking_3d_parameters]


def get_dashboard_layout(param_dict):
    """
    Obtain the main dashboard layout. html.Div children elements are returned by various callback functions and added
    to the layout on the fly.

    :param param_dict: Dictionary of parameters used to populate default values.
    :return: Div containing the main dashboard layout.
    """
    dashboard_layout = html.Div(
        [
            dcc.Loading(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Div(
                                    [
                                        dbc.Button('Load mzML File',
                                                   id='upload_mzml',
                                                   style={'margin': '20px'},
                                                   disabled=False),
                                        dbc.Button('Load Bruker *.d Directory',
                                                   id='upload_d',
                                                   style={'margin': '20px'},
                                                   disabled=False)
                                    ],
                                    style={'justify-content': 'center',
                                           'display': 'flex'}
                                ),
                                width={'size': 4, 'offset': 4}
                            )
                        ]
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Label('Spectrum ID'),
                                width=12
                            )
                        ]
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                dcc.Dropdown(id='spectrum_id',
                                             multi=False,
                                             options=[],
                                             value=[],
                                             disabled=False,
                                             style={'margin': '5px'}),
                                width=12
                            )
                        ]
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                dcc.Dropdown(id='spectrum_id_2',
                                             multi=False,
                                             options=[],
                                             value=[],
                                             disabled=False,
                                             style={'margin': '5px',
                                                    'display': 'none'}),
                                width=12
                            )
                        ]
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Div(
                                    [
                                        dbc.Button('Trim Spectrum',
                                                   id='trim_spectrum',
                                                   style={'margin': '5px'},
                                                   disabled=False),
                                        dbc.Button('Transform Intensity',
                                                   id='transform_intensity',
                                                   style={'margin': '5px'},
                                                   disabled=False),
                                        dbc.Button('Smooth Baseline',
                                                   id='smooth_baseline',
                                                   style={'margin': '5px'},
                                                   disabled=False),
                                        dbc.Button('Remove Baseline',
                                                   id='remove_baseline',
                                                   style={'margin': '5px'},
                                                   disabled=False),
                                        dbc.Button('Normalize Intensity',
                                                   id='normalize_intensity',
                                                   style={'margin': '5px'},
                                                   disabled=False),
                                        dbc.Button('Bin Spectrum',
                                                   id='bin_spectrum',
                                                   style={'margin': '5px'},
                                                   disabled=False),
                                        dbc.Button('Label Peaks',
                                                   id='peak_picking',
                                                   style={'margin': '5px'},
                                                   disabled=False)
                                    ],
                                    style={'justify-content': 'center',
                                           'display': 'flex'}
                                ),
                                width=12
                            )
                        ]
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Div(
                                    [
                                        dbc.Button('Export Peak List from Labeled Peaks',
                                                   id='export_peak_list',
                                                   style={'margin': '5px'},
                                                   disabled=False),
                                        dbc.Button('Undo Preprocessing',
                                                   id='undo_preprocessing',
                                                   style={'margin': '5px'},
                                                   disabled=False),
                                        dbc.Button('Undo Peak Labeling',
                                                   id='undo_peak_picking',
                                                   style={'margin': '5px'},
                                                   disabled=False),
                                        dbc.Button('Toggle log Intensity Scale',
                                                   id='toggle_log_intensity',
                                                   style={'margin': '5px'},
                                                   disabled=True),
                                        dbc.Button('Toggle Mirror Plot',
                                                   id='toggle_mirror_plot',
                                                   style={'margin': '5px'},
                                                   disabled=False),
                                        dbc.Button('Edit Preprocessing Parameters',
                                                   id='edit_preprocessing_parameters',
                                                   style={'margin': '5px'},
                                                   disabled=False),
                                        dbc.Button('About',
                                                   id='about',
                                                   style={'margin': '5px'},
                                                   disabled=False)
                                    ],
                                    style={'justify-content': 'center',
                                           'display': 'flex'}
                                ),
                                width=12
                            )
                        ]
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                dcc.Graph(
                                    id='spectrum_plot',
                                    figure=blank_figure(),
                                    style={'width': '100%',
                                           'height': '600px'}
                                ),
                                width=12
                            )
                        ]
                    ),

                    dbc.Modal(
                        [
                            dbc.ModalHeader(dbc.ModalTitle('Preprocessing Parameters')),
                            dbc.ModalBody(children=get_preprocessing_parameters_layout(param_dict),
                                          id='edit_processing_parameters_modal_body'),
                            dbc.ModalFooter(
                                [
                                    dbc.ButtonGroup(
                                        [
                                            dbc.Button('Load Parameters',
                                                       id='edit_processing_parameters_load',
                                                       className='ms-auto'),
                                            dbc.Button('Export Parameters',
                                                       id='edit_processing_parameters_export',
                                                       className='ms-auto')
                                        ]
                                    ),
                                    dbc.ButtonGroup(
                                        [
                                            dbc.Button('Cancel',
                                                       id='edit_processing_parameters_cancel',
                                                       className='ms-auto'),
                                            dbc.Button('Save',
                                                       id='edit_processing_parameters_save',
                                                       className='ms-auto')
                                        ]
                                    )
                                ]
                            )
                        ],
                        id='edit_processing_parameters_modal',
                        fullscreen=True,
                        backdrop='static',
                        scrollable=True,
                        centered=True,
                        is_open=False
                    ),
                    dbc.Modal(
                        [
                            dbc.ModalHeader(dbc.ModalTitle('Preprocessing parameters have been saved.')),
                            dbc.ModalFooter(dbc.Button('Close',
                                                       id='edit_processing_parameters_modal_saved_close',
                                                       className='ms-auto'))
                        ],
                        id='edit_processing_parameters_modal_saved',
                        centered=True,
                        is_open=False
                    ),
                    dbc.Modal(
                        [
                            dbc.ModalHeader(dbc.ModalTitle('Error')),
                            dbc.ModalBody('Mirror spectrum cannot be visualized with a three dimensional spectrum '
                                          'containing ion mobility data.'),
                            dbc.ModalFooter(dbc.Button('Close',
                                                       id='plot_mirror_spectrum_error_modal_close',
                                                       className='ms-auto'))
                        ],
                        id='plot_mirror_spectrum_error_modal',
                        centered=True,
                        is_open=False
                    ),
                    dbc.Modal(
                        [
                            dbc.ModalHeader(dbc.ModalTitle(f'pyMALDIviz {VERSION}')),
                            dbc.ModalBody('Included software components: Copyright © 2022 by Bruker Daltonics GmbH & '
                                          'Co. KG. All rights reserved.'),
                            dbc.ModalFooter(dbc.Button('Close',
                                                       id='about_close',
                                                       className='ms-auto'))
                        ],
                        id='about_modal',
                        centered=True,
                        is_open=False
                    ),

                    dcc.Store(id='store_preprocessing_params',
                              data=get_preprocessing_params()),
                    dcc.Store(id='store_use_log_intensity',
                              data=True),
                    dcc.Store(id='store_plot_type',
                              data='single'),
                    dcc.Store(id='store_plot'),

                    html.Div(
                        id='dummy',
                        children=[],
                        style={'display': 'none'}
                    )
                ],
                overlay_style={'visibility': 'visible', 'opacity': 0.5}
            )
        ],
        style={'margin': '20px'}
    )

    return dashboard_layout
