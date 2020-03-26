from os.path import exists
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import Table, Column

from Metallicity_Stack_Commons.column_names import filename_dict, bin_names0, remove_from_list

def make_validation_table(fitspath):
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
        fitspath + 'bin_validation.tbl' --> a validation table containing bin IDs;
                                        number of galaxies in each bin; and 
                                        column indicating  OIII4363 detection/non-detection
                                        OIII4363_Flux_Observed
                                        OIII4363_S/N
    '''

    bin_table = asc.read(fitspath + filename_dict['bin_info'])
    em_table = asc.read(fitspath + filename_dict['bin_fit']) 

    bin_ID = em_table['bin_ID'].data
    raw_OIII4363 = em_table['OIII_4363_Flux_Observed'].data
    O_4363_SN = em_table['OIII_4363_S/N'].data
    O_4363_sigma = em_table['OIII_4363_Sigma'].data
    O_5007_SN = em_table['OIII_5007_S/N'].data
    
    N_stack = bin_table['N_stack'].data
    Hgamma_SN    = em_table['HGAMMA_S/N'].data
    Hgamma = em_table['HGAMMA_Flux_Observed'].data

    detection  = np.zeros(len(bin_ID))
    OIII4363 = np.zeros(len(bin_ID))
    up_limit = (Hgamma/Hgamma_SN) *3


    valid_stacks_idx = np.where((O_4363_SN >= 3) & (O_5007_SN > 100) & (O_4363_sigma < 1.6))[0] 
    reliable_5007_stacks = np.where((O_4363_SN < 3) & (O_5007_SN > 100))[0]
    wide_lines_valid = np.where((O_4363_SN >= 3) & (O_5007_SN > 100) & (O_4363_sigma >= 1.6))[0]
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
    tab1 = Table([bin_ID, N_stack, detection, OIII4363, O_4363_SN], names = valid_table_names0)   #names=n)
    asc.write(tab1, ver_tab, format='fixed_width_two_line')

    # Write revised file for human editing
    ver_tab_revised = fitspath + filename_dict['bin_valid_rev']
    if not exists(ver_tab_revised):
        asc.write(tab1, ver_tab_revised, format='fixed_width_two_line')
        print("   ")
        print("URGENT!!! HUMAN EDITING OF FILE NEEDED : "+ver_tab_revised)
        print("   ")
    else:
        print("   ")
        print("ERROR!!! FILE EXISTS!!!  WILL NOT OVERWRITE !!!")
        print("ERROR!!! PLEASE RENAME/DELETE FILE TO REGENERATE !!!")
        print("   ")


def compare_to_by_eye(fitspath,dataset):
    '''
    Purpose -> This function takes the automated validation table and checks it against 
               inputted measurement that are determined by eye. These inputted measurements 
               are in the np.where statements. It outputs a revised validation table based 
               on the inputted measurements. 
    
    Usage -> valid_table.make_validation_table(fitspath, dataset)

    Params:
        fitspath --> a string of the file path where the input file is and where the output file
            will be placed
        dataset  --> a string that is used to determine which eye measurements to use
    
    Returns: 
        None
        
    Outputs:
        fitspath + 'bin_validation_revised.tbl' and '.csv' --> 
                                        a validation table containing bin IDs;
                                        number of galaxies in each bin; and 
                                        column indicating  OIII4363 detection/non-detection
                                        OIII4363_Flux_Observed
                                        OIII4363_S/N
                                        Notes

    '''
    ver_table = fitspath+ filename_dict['bin_valid']
    ver_tab = asc.read(ver_table)
    indicate = ver_tab['Detection']
    ID = ver_tab['bin_ID']

    #Detections By Eye
    if dataset== 'Voronoi20': det_4363 = np.where((ID ==0)| (ID ==2)| (ID ==3) | (ID ==5) | (ID==6))[0]
    if dataset== 'Voronoi14': det_4363 = np.where((ID ==0)| (ID ==7)| (ID ==10) | (ID ==11) | (ID==12))[0]
    if dataset== 'Voronoi10': det_4363 = np.where((ID ==1)| (ID ==9)| (ID ==18) | (ID ==21))[0]
    if dataset== 'Grid': det_4363 = np.where((ID ==11)| (ID ==13)| (ID ==19) | (ID ==20) | (ID == 21))[0]
    if dataset== 'R23_Grid': det_4363 = np.where((ID ==0)| (ID ==4)| (ID ==5) | (ID ==6))[0]
    if dataset== 'O32_Grid': det_4363 = np.where((ID ==6))[0]
    if dataset== 'Double_Bin': det_4363 = np.where((ID ==0)| (ID ==1)| (ID ==2) | (ID ==7) | (ID ==9) | (ID ==10) | (ID ==11) | (ID ==13))[0]
    if dataset== 'n_Bins':
        det_4363= np.where( (ID ==10) | (ID ==11) |(ID==14) | (ID ==15) | (ID ==20) | (ID==23) | (ID ==26))[0]
        rlimit = np.where( (ID ==5) | (ID ==7) |(ID==8) | (ID ==13) | (ID ==16) | (ID==17) | (ID ==19)| (ID ==22))[0]

    ######Caroline: Add you conditions here#######

    check_ID = np.zeros(len(ID))
  
    check_ID[det_4363]=1
    check_ID[rlimit] =0.5

    for ii in range(len(ID)):
        if check_ID[ii] == indicate[ii]:
            print(ID[ii], 'matches with by eye validation')
        else:
            print('*****', ID[ii], 'does not match calculated values. Please check!')

    ###This is where I need to add the column for notes
    if dataset == 'n_Bins':
        notes = ['N/A','N/A','N/A','N/A','N/A',
                 'N/A','N/A','Broad features, but reliable OIII5007 and HGAMMA',
                 'Bad fit, but good OIII5007','N/A',
                 'N/A','N/A','N/A','N/A','N/A',
                 'High Temperature','not fit well, but reliable OIII5007 and HGAMMA',
                 'N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A']
        note_add = Column(name='Notes', data = notes)
        ver_tab.add_column(note_add,5)

    #######Caroline: Add your notes column here and copy the note_add and ver_tab.add_column lines to your if statement######

    
    ver_tab.remove_column('Detection')
    
    detect_add = Column(name='Detections', data =check_ID)
    ver_tab.add_column(detect_add, 2)

    asc.write(ver_tab, fitspath+ 'bin_validation_revised.tbl',  format = 'fixed_width_two_line') 
    asc.write(ver_tab, fitspath+ 'bin_validation_revised.csv', format = 'csv')    
    
        
