# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 17:56:11 2025

@author: Angel.BAUDON
"""

import pyabf, matplotlib.pyplot as plt, numpy as np, glob, pandas as pd, os, scipy.stats as stat
from pyabf.filter import gaussian


folder = r"C:\Users\Angel.BAUDON\Desktop\To analyse\OptoPuff Sorb"
index_puff = (30, 1290)

for sub_folder in glob.glob(rf'{folder}\*/'):
    subfo = sub_folder.split('\\')[-2]
    print('\n'*2, ':'*30, '\n', subfo, '\n', ':'*30, '\n'*2)
    cells = glob.glob(rf'{sub_folder}\*.abf')
    if not os.path.exists(rf'{sub_folder}\analysis'): os.makedirs(rf'{sub_folder}\analysis')

    Raw_Vm, Data, Amp, Amp_local, Indexes = [], [], [], [], []
    for cell in cells:
        cell_id = cell.split('\\')[-1]
        print('\n'*2, '='*30, '\n', cell_id, '\n', '='*30, '\n'*2)

        if cell[-13:] == 'to merge1.abf':
            abf = pyabf.ABF(cell[:-5]+'1.abf')
            gaussian(abf, 10)
            raw1 = abf.sweepY
            
            abf = pyabf.ABF(cell[:-5]+'2.abf')
            gaussian(abf, 10)
            raw2 = abf.sweepY

            raw = np.concatenate((raw1, raw2))
            
        elif cell[-13:] == 'to merge2.abf':
            continue
        
        else:
            abf = pyabf.ABF(cell)
            gaussian(abf, 10)
            raw = abf.sweepY
            
    
        fltr = raw[:1900000:100]
            
        while len(fltr) < 19000: fltr = np.append(fltr, np.nan)
        
        sampling = 10
        time = np.linspace(0, len(fltr)/sampling, len(fltr))
     
        #Find peak
        baseline = np.nanmean(fltr[:index_puff[0]*sampling])
        data = fltr - baseline
        
        baselines, amps, amps_loc, amp_indexes, washs = [], [], [], [], []
        for indx_puff in index_puff:
            start, stop = indx_puff*sampling, (indx_puff+60)*sampling
            amp = np.nanmax(data[start:stop])
            max_index = np.where(data[start:stop] == amp)[0][0]+start

            local_baseline = np.nanmin(data[(indx_puff-20)*sampling:indx_puff*sampling])
            local_amp = amp-local_baseline
            wash = np.nanmin(data[max_index:max_index+100*sampling])
            
            baselines.append(local_baseline), amps.append(amp), amps_loc.append(local_amp)
            amp_indexes.append(max_index), washs.append(wash)
            
        #PLot section
        plt.figure(), plt.title(cell_id)
        plt.plot(time, fltr), plt.ylim(-200, 0)
        for indx in index_puff: plt.axvline(indx, c='r', label='Sorb puff')
        for amp_indx in amp_indexes: plt.plot(amp_indx/sampling, fltr[amp_indx], 'og', label='Amplitude max')
        plt.legend()
        plt.savefig(rf'{sub_folder}\analysis\{cell_id}.pdf')
        # plt.close()        
        
        Raw_Vm.append([x+baseline for x in (baselines[0], amps[0], washs[0], baselines[1], amps[1], washs[1])])
        Data.append(data), Amp.append(amps), Amp_local.append(amps_loc)
        Indexes.append(max_index/sampling)

    
    writer = pd.ExcelWriter(f'{sub_folder}/analysis/{subfo}.xlsx')
    for x, y in zip((Raw_Vm, Amp, Amp_local, Indexes),
                    ('Raw Vm', 'Amp', 'Amp loc', 'Indexes')):
        df = pd.DataFrame(x)
        df.rename(index=pd.Series([x.split('\\')[-1] for x in cells])).to_excel(writer, sheet_name=y)
    writer.save()
    
        
    plt.figure(), plt.title(subfo)
    mean, sem = np.nanmean(Data, axis=0), stat.sem(Data, nan_policy='omit', axis=0)
    plt.plot(time, mean), plt.fill_between(time, mean-sem, mean+sem, alpha=.5, zorder=1)
    plt.legend(), plt.savefig(rf'{sub_folder}\analysis\Mean trace.pdf')
    plt.close()



