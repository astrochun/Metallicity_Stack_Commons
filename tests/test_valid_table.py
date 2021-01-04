from Metallicity_Stack_Commons import valid_table

deep2_dir = 'tests_data/DEEP2_Ly2015/'


def test_valid_table():

    valid_table.make_validation_table(deep2_dir)
