# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 11:52:54 2025

@author: Angel.BAUDON
"""

import pandas as pd, numpy as np, matplotlib.pyplot as plt, glob, scipy.stats as stat, os, glob
from scipy.signal import savgol_filter, find_peaks


Motherfolder = r'C:\Angel.BAUDON\Exp\Data\RImp Wounding\Data\RImp Wound Mutants'

# def func(x, Y0 = 65.03, Plateau = -2088301, K = 4.148e-009): # For Glut puff
def func(x, Y0 = 98.22, Plateau = 5.745, K = 0.0006842): # For wounding
    #https://www.graphpad.com/guides/prism/latest/curve-fitting/reg_exponential_decay_1phase.htm
    return (Y0 - Plateau)*np.exp(-K*x) + Plateau

for folder in glob.glob(f'{Motherfolder}\*'):
    folder_name = folder.split('\\')[-1]
    print(folder_name)

    distance_file = glob.glob(rf'{folder}\raw\*.xlsx')[0]
    distances = pd.read_excel(distance_file, header=None)
    predicted = pd.Series([func(distance) for distance in distances.iloc[:,1]])
    
    amp_file = rf'{folder}\analysis\{folder_name}.xlsx'
    amps = pd.read_excel(amp_file, sheet_name='Amp').iloc[:,1]
    
    pd.Series(amps)
    
    dMPP = pd.Series([M/P for P, M in zip(predicted, amps)])
    
    out = pd.concat([distances, predicted, amps, dMPP], axis=1)
    out.columns = ['ID', 'Distance', 'Predicted', 'Measured', 'M-P/P']
    
    writer = pd.ExcelWriter(rf'{folder}/analysis/Prediction {folder_name}.xlsx')
    out.to_excel(writer)
    writer.save()  

    

