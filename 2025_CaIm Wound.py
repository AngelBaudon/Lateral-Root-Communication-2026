# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 16:08:51 2025

@author: Angel.BAUDON
"""
import pandas as pd, numpy as np, matplotlib.pyplot as plt, glob, scipy.stats as stat, os
from scipy.signal import savgol_filter, find_peaks
from scipy.optimize import curve_fit



folder = r"C:\Angel.BAUDON\Exp\Data\Root CaIm wounding\Root GCaMP wounding 12µm pipette 10x"
if not os.path.exists(rf'{folder}\analysis'): os.makedirs(rf'{folder}\analysis')

file = glob.glob(f'{folder}\*.xlsx')[0]
file_name = file.split('\\')[-1]
wound, sampling_Hz, rec_len = 28, .5, 330
data, Amps, Max_indexes, Rise_taus, Decay_taus, HMFW = [], [], [], [], [], []
xl = pd.ExcelFile(file)

for s, sheet_name in enumerate(xl.sheet_names):
    
    raw = pd.read_excel(file, sheet_name=sheet_name, header=None).to_numpy()
    n_roi = raw.shape[1]-1

    camera_background = raw[:,0]

    dFF0, amps, max_indexes, rise_taus, decay_taus, hmfw = [], [], [], [], [], []

    # plt.figure()

    for i in range(n_roi):
        F = raw[:,i+1] - camera_background
        
        baseline = np.nanmean(F[:wound])
        dff0 = (F[:rec_len]-baseline)/baseline
        fltr = savgol_filter(dff0, 5, 2)

        # plt.subplot(n_roi, 3, i*3+1), plt.plot(dff0), plt.plot(fltr)

        
        amp = max(fltr[wound:])
        max_index = np.where(fltr[wound:] == amp)[0][0] + wound
    
    
        #Find return to baseline
        try: return_to_baseline = np.where(fltr[max_index:] < amp/10)[0][0]+max_index
        except IndexError: return_to_baseline = len(fltr)
        
        #Define and calculate rise and decay
        rise = fltr[wound:max_index]
        decay = fltr[max_index:return_to_baseline]
        
        plt.figure()
        plt.subplot(211), plt.plot(fltr)
        plt.subplot(223), plt.plot(rise)
        plt.subplot(224), plt.plot(decay)
        
        
#         try:
#             start, stop = np.where(fltr[:max_index] > amp/2)[0][0], np.where(fltr[max_index:] < amp/2)[0][0]+max_index
#             hmfw.append((stop-start)/sampling_Hz)
#         except IndexError: hmfw.append(np.nan)

#         taus = []
#         def expo(x, a, b, c): return a * np.exp(b * x) + c
        
#         for k, kin in enumerate((rise, decay)):
#             if len(kin)>2:
#                 kin = np.array([x for x in kin if str(x) != 'nan'])
#                 kin= kin - min(kin)
        
#                 p0 = [-100, 0, 100] if k else [100, 0, -100]
#                 popt, _ = curve_fit(expo, np.arange(len(kin)), kin, p0 = p0, maxfev=5000)
#                 fit = expo(np.arange(len(kin)), *popt)
                
#                 # plt.subplot(n_roi, 3, i*3+k+2)
#                 # plt.plot(kin)
#                 # plt.plot(np.arange(len(kin)), fit)
                    
#                 fit = fit - min(fit)
#                 tau = np.where(fit < 0.368*max(fit)) if k else np.where(fit > 0.632*max(fit))
#                 taus.append(tau[0][0]/sampling_Hz)
            
#             else: taus.append(np.nan)
            
#         dFF0.append(fltr), amps.append(amp), max_indexes.append(max_index/sampling_Hz)
#         rise_taus.append(taus[0]), decay_taus.append(taus[1])

    
#     plt.figure(), plt.title(sheet_name)
#     for k, kaboom in enumerate(dFF0): plt.plot(kaboom, label=f"{k*100}µm")
#     plt.legend()
#     plt.savefig(rf'{folder}/analysis/{sheet_name}.pdf')
#     plt.close()

    
#     data.append(np.asarray(dFF0)), Amps.append(amps), Max_indexes.append(max_indexes)
#     Rise_taus.append(rise_taus), Decay_taus.append(decay_taus), HMFW.append(hmfw)

# x_ax = np.linspace(0, rec_len/sampling_Hz, rec_len)
# max_roi, max_img = max([d.shape[0] for d in data]), max([d.shape[1] for d in data])

# ar = np.empty((len(data), max_roi, max_img))
# ar[:] = np.nan

# for d, dat in enumerate(data):
#     roi, img = dat.shape
#     ar[d, :roi, :img] = dat


# plt.figure()
# means, sems = np.nanmean(ar, axis=0), stat.sem(ar, axis=0, nan_policy='omit')
# for i, (m, s) in enumerate(zip(means, sems)): 
#     plt.plot(x_ax, m, label=f'{i*100} µm'), plt.fill_between(x_ax, m-s, m+s, alpha=0.5)
# plt.xlabel('Time(s)'), plt.ylabel('dF/F0'), plt.legend()
# plt.savefig(rf'{folder}/analysis/{file_name[:-5]}.pdf')

      
# writer = pd.ExcelWriter(rf'{folder}/analysis/{file_name[:-5]} analysis.xlsx')
# pd.DataFrame(Amps).to_excel(writer, sheet_name = 'Amp')
# pd.DataFrame(Max_indexes).to_excel(writer, sheet_name = 'Index')
# pd.DataFrame(Rise_taus).to_excel(writer, sheet_name = 'Rise')
# pd.DataFrame(Decay_taus).to_excel(writer, sheet_name = 'Decay')
# pd.DataFrame(HMFW).to_excel(writer, sheet_name = 'HMFW')
# writer.save()  









