from Metallicity_Stack_Commons.analysis import error_prop

deep2_dir = 'tests_data/DEEP2_Ly2015/'


def test_fluxes_derived_prop():

    # Even though these are individual galaxies, we can trick it to think
    # it's binned data.
    params_list = [
        {'raw': True,  # No MC, no dust, use default valid_table
         'binned_data': True,
         'apply_dust': False,
         'revised': False},
        {'raw': True,  # No MC, no dust, use revised valid_table
         'binned_data': True,
         'apply_dust': False,
         'revised': True},
        {'raw': True,  # No MC, apply dust, use default valid_table
         'binned_data': True,
         'apply_dust': True,
         'revised': False},
        {'raw': True,  # No MC, apply dust, use revised valid_table
         'binned_data': True,
         'apply_dust': True,
         'revised': True},
        {'raw': False,  # MC, no dust, use default valid_table
         'binned_data': True,
         'apply_dust': False,
         'revised': False},
        {'raw': False,  # MC, no dust, use revised valid_table
         'binned_data': True,
         'apply_dust': False,
         'revised': True},
        {'raw': False,  # MC, apply dust, use default valid_table
         'binned_data': True,
         'apply_dust': True,
         'revised': False},
        {'raw': False,  # MC, apply dust, use revised valid_table
         'binned_data': True,
         'apply_dust': True,
         'revised': True},
    ]
    for params in params_list:
        error_prop.fluxes_derived_prop(deep2_dir, **params, verbose=True)
