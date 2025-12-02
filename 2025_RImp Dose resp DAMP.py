# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 11:59:07 2024

@author: Angel.BAUDON
"""
import pyabf, matplotlib.pyplot as plt, numpy as np, glob, pandas as pd, os, scipy.stats as stat
from pyabf.filter import gaussian
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter



folder = r'C:\Users\Angel.BAUDON\Desktop\NaCl vs KCl'
sub_folders = [x for x in glob.glob(rf'{folder}\*') if x.split('\\')[-1] not in ('analysis', 'raw')]
if not os.path.exists(rf'{folder}\analysis'): os.makedirs(rf'{folder}\analysis')
All_time_course = {}
show_fig = False
index_puff, cut_off_early_response = 20, 30
separate_early_and_late_response = True
writer = pd.ExcelWriter(f'{folder}/analysis/DR analysis.xlsx')


for sub_folder in sub_folders:
    subfo = sub_folder.split('\\')[-1]

    print('\n'*2, ':'*30, '\n', subfo, '\n', ':'*30, '\n'*2)
    cells = [x for x in glob.glob(rf'{sub_folder}\*.abf') if x.split('\\')[-1] not in ('analysis')]
    if not os.path.exists(rf'{folder}\analysis\{subfo}'): os.makedirs(rf'{folder}\analysis\{subfo}')

    Amp, time_course = {}, []
    for cell in cells:
        cell_id = cell.split('\\')[-1]
        print('\n'*2, '='*30, '\n', cell_id, '\n', '='*30, '\n'*2)
    
        abf = pyabf.ABF(cell)
        raw = list(abf.sweepY)[10000:300000:100]
        while len(raw) < 2900: raw = np.append(raw, np.nan)
        
        # gaussian(abf, 1000)
        gaussian(abf, 100)
        fltr, time, sampling = abf.sweepY[10000:300000:100], np.linspace(0, 290, 2900), 10
        idx_puff, cut_off_early = index_puff*sampling, cut_off_early_response*sampling
        while len(fltr) < 2900: fltr = np.append(fltr, np.nan)
        baseline = np.nanmean(fltr[(index_puff-10)*sampling:idx_puff])
        data = fltr - baseline
        
        index = np.where(data == np.nanmax(data[index_puff*sampling:]))[0][0]
        
        #Early response
        early = data[idx_puff:idx_puff + cut_off_early]
        index_early = np.where(data == np.nanmax(early))[0][0]
        
        #Late response
        late = data[idx_puff + cut_off_early:]
        index_late = np.where(data == np.nanmax(late))[0][0]
        
        
        
        
        rise, decay, taus = data[idx_puff:index], data[index:], []
        def expo(x, a, b, c): return a * np.exp(b * x) + c
        # plt.figure()
        for k, kin in enumerate((rise, decay)):
            try:

                p0 = [100, 0, -100] if k else [-100, 0, 100]
                kin = np.array([x for x in kin if str(x) != 'nan'])
                kin= kin - min(kin)
    
                popt, _ = curve_fit(expo, np.arange(len(kin)), kin, p0 = p0, maxfev=5000)
                fit = expo(np.arange(len(kin)), *popt)
                
                fit = fit - min(fit)
                tau = np.where(fit < 0.368*max(fit)) if k else np.where(fit > 0.632*max(fit))
                taus.append(tau[0][0]/sampling)

                # plt.subplot(1, 2, k+1), plt.plot(kin)
                # plt.plot(np.arange(len(kin)), fit, label=tau), plt.legend()
                
            except TypeError: taus.append(np.nan)
            except ValueError: taus.append(np.nan)
        
        
        #Stock section
        Amp[cell_id] = [np.nanmax(data), np.nanmax(early), np.nanmax(late),
                        *[x/sampling - 20 for x in (index, index_early, index_late)], *taus]
        time_course.append(data)
        
        #Plot section
        plt.figure(), plt.title(cell_id)
        plt.plot(time, raw), plt.plot(time, fltr), plt.ylim(-220, -20)
        Vinit = np.nanmean(fltr[:index_puff*sampling])    
        plt.axvline(index_puff, color='purple')
        plt.axvspan(index_puff, index_puff + cut_off_early_response, color='gold', alpha=.5)
        plt.text(index_puff+1, -190, 'Early', size=15, weight='bold')
        plt.text(0, -10, f'Vinit = {Vinit}', size=15, weight='bold')
        for index in (index_early, index_late): plt.plot(index/sampling, fltr[index], 'xr')
        plt.savefig(rf'{folder}\analysis\{subfo}\{cell_id}.pdf')
        if not show_fig: plt.close()
    
    #Stock section
    All_time_course[subfo] = time_course
    out = pd.DataFrame.from_dict(Amp, orient='index',
                                  columns = ('Max', 'Early', 'Late', 'Idx Mx', 'Idx E', 'Idx L', 'Rise', 'Deacy'))
    out.to_excel(writer, sheet_name = subfo)
writer.save()    


plt.figure()
for key in All_time_course.keys():
    Time_course = All_time_course[key]
    if not len(Time_course): continue
    mean, sem = np.nanmean(Time_course, axis=0), stat.sem(Time_course, axis=0, nan_policy='omit')
    m, s = savgol_filter(mean, 51, 3), savgol_filter(sem, 51, 3)
    plt.plot(time, m), plt.fill_between(time, m-s, m+s, alpha=0.5, zorder=1, label=key)
plt.legend()



