# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 15:28:27 2025

@author: Angel.BAUDON
"""
import pandas as pd, numpy as np, matplotlib.pyplot as plt, glob, scipy.stats as stat, os
from scipy.signal import savgol_filter, find_peaks


folder = r"C:\Angel.BAUDON\Exp\Data\Root CaIm wounding\Root MCA-GCaMP"
if not os.path.exists(rf'{folder}\analysis'): os.makedirs(rf'{folder}\analysis')

file = glob.glob(f'{folder}\*.xlsx')[0]
file_name = file.split('\\')[-1]
sampling_Hz, rec_len = .5, 289
data, Amps = [], []
xl = pd.ExcelFile(file)


for sheet_name in xl.sheet_names:
    
    raw = pd.read_excel(file, sheet_name=sheet_name).to_numpy()
    _, n_rec = raw.shape
    
    dFF0, amps = [], []
    
    for i in range(int(n_rec/2)):
        
        camera_background = raw[:,i*2]
        F = raw[:,i*2+1] - camera_background
        
        baseline = np.nanmean(F[:30])
        dff0 = (F[:rec_len]-baseline)/baseline
        
        fltr = savgol_filter(dff0, 5, 2)
       
        # plt.figure(), plt.title(f'{sheet_name} Rec n°{i}')
        # plt.plot(dff0), plt.plot(fltr)

        dFF0.append(fltr), amps.append(max(dff0[30:]))
    data.append(np.asarray(dFF0)), Amps.append(amps)

x_ax = np.linspace(0, rec_len/sampling_Hz, rec_len)

plt.figure()
for i, d in enumerate(data):
    m, s = np.nanmean(d, axis=0), stat.sem(d, axis=0, nan_policy='omit')

    # m, s = savgol_filter(m, 3, 1), savgol_filter(s, 3, 1)
    plt.plot(x_ax, m, label=xl.sheet_names[i]), plt.fill_between(x_ax, m-s, m+s, alpha=0.5)
plt.xlabel('Time(s)'), plt.ylabel('dF/F0'), plt.legend()
plt.savefig(rf'{folder}/analysis/{file_name[:-5]}.pdf')

      
writer = pd.ExcelWriter(rf'{folder}/analysis/{file_name[:-5]} analysis.xlsx')
for d, name in zip(data, xl.sheet_names): pd.DataFrame(d).to_excel(writer, sheet_name = f'{name} dFF0')
for a, name in zip(Amps, xl.sheet_names): pd.DataFrame(a).to_excel(writer, sheet_name = f'{name} Amp')
writer.save()  

    