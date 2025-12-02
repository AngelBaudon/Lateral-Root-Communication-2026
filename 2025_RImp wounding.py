# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 16:52:09 2024

@author: Angel.BAUDON
"""
import pyabf, matplotlib.pyplot as plt, numpy as np, glob, pandas as pd, os, scipy.signal, scipy.stats as stat
from pyabf.filter import gaussian
from scipy.optimize import curve_fit


folder = r'C:\Angel.BAUDON\Exp\Data\RImp Wounding\Data\RImp Wound Mutants'
sub_folders = [x for x in glob.glob(rf'{folder}\*') if x.split('\\')[-1] != 'analysis']
show_fig = False

for sub_folder in sub_folders:
    subfo = sub_folder.split('\\')[-1]
    print('\n'*2, ':'*30, '\n', subfo, '\n', ':'*30, '\n'*2)
    cells = glob.glob(rf'{sub_folder}\*.abf')
    if not os.path.exists(rf'{sub_folder}\analysis'): os.makedirs(rf'{sub_folder}\analysis')

    Raw_Vm, Data, Amp, Indexes, Taus, HMFW = [], [], [], [], [], []
    for cell in cells:
        cell_id = cell.split('\\')[-1]
        print('\n'*2, '='*30, '\n', cell_id, '\n', '='*30, '\n'*2)
    
        abf = pyabf.ABF(cell)
        gaussian(abf, 10)
        
        if 'SCW' in sub_folder:
            fltr, time = abf.sweepY[:200000:100], abf.sweepX[:200000:100]
            while len(fltr) < 2000: fltr = np.append(fltr, np.nan)

        else:
            fltr = abf.sweepY[:500000:100]
            while len(fltr) < 5000: fltr = np.append(fltr, np.nan)
            
        sampling = 10
        time = np.linspace(0, len(fltr)/sampling, len(fltr))
            
        #Find the wounding
        deriv2 = np.diff(np.diff(fltr)).clip(0)
        _, ppts = scipy.signal.find_peaks(deriv2[270:310], height=3, prominence=3)
        deriv_peak = ppts['left_bases']
        wound = deriv_peak[0]+271 if len(deriv_peak) else 300
        
        #Find peak
        baseline = np.nanmean(fltr[wound-10*sampling:wound])
        data = fltr - baseline
        initial_amp = np.nanmax(data[wound:wound+60])
        max_amp = np.nanmax(data[wound:wound+3000])
        index, max_index = [np.where(data[wound:] == x)[0][0]+wound for x in (initial_amp, max_amp)]
        
        #Find return to baseline
        try: return_to_baseline = np.where(data[index:] < max_amp/10)[0][0]+index
        except IndexError: return_to_baseline = 3000
        
        #Define and calculate rise and decay
        rise = data[wound:index]
        decay = data[index:return_to_baseline]
        depol_time = (return_to_baseline - wound)/sampling

        taus = []
        def expo(x, a, b, c): return a * np.exp(b * x) + c
        # plt.figure()
        for k, kin in enumerate((rise, decay)):
            try:

                p0 = [100, 0, -100] if k else [-100, 0, 100]
                kin = np.array([x for x in kin if str(x) != 'nan'])
                kin= kin - min(kin)
    
                # plt.subplot(1, 2, k+1), plt.plot(kin)
                popt, _ = curve_fit(expo, np.arange(len(kin)), kin, p0 = p0, maxfev=5000)
                fit = expo(np.arange(len(kin)), *popt)
                # plt.plot(np.arange(len(kin)), fit)
                
                fit = fit - min(fit)
                tau = np.where(fit < 0.368*max(fit)) if k else np.where(fit > 0.632*max(fit))
                taus.append(tau[0][0]/sampling)
                
            except TypeError: taus.append(np.nan)
            except ValueError: taus.append(np.nan)
            
        
        #Half Max Full Width        
        try: hmfw = len(np.where(data[wound:] >= max_amp/2)[0])/sampling
        except: hmfw = np.nan

            
            
        #PLot section
        plt.figure(), plt.title(cell_id)
        plt.plot(time, fltr), plt.ylim(-200, 0)
        plt.plot(time[2:], deriv2+baseline, c='gold', label='2nd derivative')
        plt.axvline(wound/sampling, c='r', label='Wounding')
        plt.plot(max_index/sampling, fltr[max_index], 'og', label='Amplitude max')
        plt.plot(index/sampling, fltr[index], 'xr', label='Initial max')
        for start, stop, color, label in ((wound-100, wound, 'plum', 'baseline'),
                                          (wound, index, 'moccasin', 'rise'),
                                          (index, return_to_baseline, 'mediumpurple', 'decay')):
            plt.axvspan(start/sampling, stop/sampling, color=color, alpha=.5, label=label)
        plt.legend()
        plt.savefig(rf'{sub_folder}\analysis\{cell_id}.pdf')
        if not show_fig: plt.close()        
        
        Raw_Vm.append([x+baseline for x in (0, initial_amp, max_amp, np.nanmean(data[-10*sampling:]))])
        Data.append(data), Amp.append([initial_amp, max_amp])
        Indexes.append([(x/sampling)-30 for x in [index, max_index]])
        Taus.append(taus), HMFW.append(hmfw)

    
    writer = pd.ExcelWriter(f'{sub_folder}/analysis/{subfo}.xlsx')
    for x, y in zip((Raw_Vm, Amp, Indexes, Taus, HMFW),
                    ('Raw Vm', 'Amp', 'Indexes', 'Taus', 'HMFW')):
        df = pd.DataFrame(x)
        df.rename(index=pd.Series([x.split('\\')[-1] for x in cells])).to_excel(writer, sheet_name=y)
    writer.save()
    
        
    plt.figure(), plt.title(subfo)
    mean, sem = np.nanmean(Data, axis=0), stat.sem(Data, nan_policy='omit', axis=0)
    plt.plot(time, mean), plt.fill_between(time, mean-sem, mean+sem, alpha=.5, zorder=1)
    plt.legend(), plt.savefig(rf'{sub_folder}\analysis\Mean trace.pdf')
    plt.close()



