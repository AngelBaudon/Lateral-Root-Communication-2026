# -*- coding: utf-8 -*-
"""
Created on Tue May 13 12:02:38 2025

@author: Angel.BAUDON
"""
import pyabf, matplotlib.pyplot as plt, numpy as np, glob, pandas as pd, os, scipy.signal, scipy.stats as stat
from pyabf.filter import gaussian
from scipy.optimize import curve_fit


folder = r"C:\Users\Angel.BAUDON\Desktop\To analyse\OptoPuff Glut"

for sub_folder in glob.glob(rf'{folder}\*/'):
    subfo = sub_folder.split('\\')[-2]
    print('\n'*2, ':'*30, '\n', subfo, '\n', ':'*30, '\n'*2)
    cells = glob.glob(rf'{sub_folder}\*.abf')
    if not os.path.exists(rf'{sub_folder}\analysis'): os.makedirs(rf'{sub_folder}\analysis')

    Raw_Vm, Data, Amp, Indexes = [], [], [], []
    for cell in cells:
        cell_id = cell.split('\\')[-1]
        print('\n'*2, '='*30, '\n', cell_id, '\n', '='*30, '\n'*2)
    
        abf = pyabf.ABF(cell)
        gaussian(abf, 10)

        fltr = abf.sweepY[:1000000:100]
        while len(fltr) < 10000: fltr = np.append(fltr, np.nan)

        sampling = 10
        time = np.linspace(0, len(fltr)/sampling, len(fltr))
      
        #Find peak
        baseline = np.nanmean(fltr[:29*sampling])
        data = fltr - baseline
        max_amp = np.nanmax(data[620*sampling:])
        max_index = np.where(data[620*sampling:] == max_amp)[0][0]+620*sampling
        
        #Calculate the local depolarization
        local_baseline = np.nanmin(data[600*sampling:max_index])
        local_amp = max_amp-local_baseline
        
        #Find return to baseline
        wash = np.nanmin(data[max_index:])
            
        #PLot section
        plt.figure(), plt.title(cell_id)
        plt.plot(time, fltr), plt.ylim(-200, 0)
        plt.axvline(max_index/sampling, c='r', label='Glut puff')
        plt.plot(max_index/sampling, fltr[max_index], 'og', label='Amplitude max')
        plt.legend()
        plt.savefig(rf'{sub_folder}\analysis\{cell_id}.pdf')
        # plt.close()
        
        Raw_Vm.append([x+baseline for x in (0, max_amp, wash)])
        Data.append(data), Amp.append((max_amp, local_amp))
        Indexes.append(max_index/sampling)

    
    writer = pd.ExcelWriter(f'{sub_folder}/analysis/{subfo}.xlsx')
    for x, y in zip((Raw_Vm, Amp, Indexes), ('Raw Vm', 'Amp', 'Indexes')):
        df = pd.DataFrame(x)
        df.rename(index=pd.Series([x.split('\\')[-1] for x in cells])).to_excel(writer, sheet_name=y)
    writer.save()
    
        
    plt.figure(), plt.title(subfo)
    mean, sem = np.nanmean(Data, axis=0), stat.sem(Data, nan_policy='omit', axis=0)
    plt.plot(time, mean), plt.fill_between(time, mean-sem, mean+sem, alpha=.5, zorder=1)
    plt.legend(), plt.savefig(rf'{sub_folder}\analysis\Mean trace.pdf')
    plt.close()



