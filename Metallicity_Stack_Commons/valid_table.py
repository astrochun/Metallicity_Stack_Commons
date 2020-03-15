import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import Table

from Metallicity_Stack_Commons.column_names import filename_dict, bin_names0

def make_validation_table(fitspath): #bin_type_str
    '''
    Purpose: 
        This function creates a validation table for a given binning set. The validation table
        contains a OIII4363 detection column where 1.0 means detection, 0.5 means non-detection with
        reliable OIII5007, and 0.0 means unreliable non-detection.
        This function will be run every time the analysis is completed and will create a validation 
        table for every analysis. 
        
    Usage:
        valid_table.make_validation_table(fitspath, bin_type_str)
        
    Params:
        fitspath --> a string of the file path where the input file is and where the output file
            will be placed.
        not in use: 
        bin_type_str --> a string describing the binning type. (e.g. 'massLHbetabin' or 'massbin') 
                     --> This is dataset for Zcalbase_gal analysis
        
    Returns: 
        None
        
    Outputs:
        fitspath + 'validation.tbl' --> a validation table containing bin IDs;
                                        number of galaxies in each bin; and 
                                        column indicating  OIII4363 detection/non-detection
                                        OIII4363_Flux_Observed
                                        OIII4363_S/N
    '''

    em_table = asc.read(fitspath + filename_dict['bin_fit'])  #this is combine_flux_ascii

    bin_ID = em_table['bin_ID'].data
    raw_OIII4363 = em_table['OIII_4363_Flux_Observed'].data
    O_4363_SN = em_table['OIII_4363_S/N'].data
    O_4363_sigma = em_table['OIII_4363_Sigma'].data
    O_5007_SN = em_table['OIII_5007_S/N'].data
    
    N_stack = em_table['N_stack'].data
    Hgamma_SN    = em_table['HGAMMA_S/N'].data
    Hgamma = em_table['HGAMMA_Flux_Observed'].data

    detection  = np.zeros(len(bin_ID))
    OIII4363 = np.zeros(len(bin_ID))
    up_limit = (Hgamma/Hgamma_SN) *3


    valid_stacks_idx = np.where((O_4363_SN >= 3) & (O_5007_SN > 100) & (O_4363_sigma < 2))[0] 
    reliable_5007_stacks = np.where((O_4363_SN < 3) & (O_5007_SN > 100))[0]
    wide_lines_valid = np.where((O_4363_SN >= 3) & (O_5007_SN > 100) & (O_4363_sigma >= 2))[0]
    detection[valid_stacks_idx] = 1
    detection[reliable_5007_stacks] = 0.5
    detection[wide_lines_valid] = 0.5
    print(detection)

    for ii in range(len(OIII4363)):
            
            if detection[ii] == 1:
                OIII4363[ii]= raw_OIII4363[ii]
            if detection[ii] == 0.5:
                OIII4363[ii]= up_limit[ii] 
            if detection[ii] ==0: 
                OIII4363[ii]= up_limit[ii]

    ver_tab = fitspath+ filename_dict['bin_valid']
    n =('bin_ID','N_stack','Detection', 'OIII_4363_Flux_Observed', 'OIII_4363_S/N')
    tab1 = Table([bin_ID, N_stack, detection, OIII4363, O_4363_SN], names=n)
    asc.write(tab1, ver_tab, format='fixed_width_two_line')

def quality_assurance(bin_type_str, QA_flag):
    '''
    Purpose:
        This function allows for manual flagging of sources for quality assurance of OIII4363 detections 
        for the mass luminosity analysis. Based on the bin_type_str keyword, the user can override the 
        detection flag by setting a specificbin index equal to the desired flag (1.0, 0.5, or 0.0).
        
    Usage:
        valid_table.quality_assurance(bin_type_str, QA_flag)
        
    Params:
        bin_type_str --> a string describing the binning type. (e.g. 'massLHbetabin' or 'massbin')
        QA_flag --> a numpy zeros array the size of the number of bins. This is used to flag sources by
            changing the value at a specific index to the desired flag.
        
    Returns: 
        QA_flag --> the updated flag array.
        
    Outputs:
        None
    '''
    if bin_type_str == 'mass_LHbeta_bin':
        QA_flag[10] = 1.0    #has large line width on OIII4363
        QA_flag[11] = 1.0    #has large line width on OIII4363
    elif bin_type_str == 'massbin':   
        QA_flag[5] = 1.0     #has large line width on OIII4363
        
    return QA_flag
