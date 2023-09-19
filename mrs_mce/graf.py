# graf.py - plotting sim.py
#
# Author: Konstantin E Bosbach <konstantin.bosbach@mars.uni-freiburg.de>

import glob
import re
import os
from pathlib import Path

import pandas as pd
import seaborn as sns
import numpy as np

def mcGraf(dex, explicit={}): 
    # load all data-files from path dex, matching ptrn, 
    # explicit={} constraints: list of acceptable values for respective match_obj[:]
	# cave: explicit{} can only have list values! 
	# e.g. expli = {"parameter":[para], "noise-SD":[n], …}
    dfs = pd.DataFrame()
    ptrn = re.compile(r'(\w{0,5})-(\w{3,})~([-\d\.]+)_noiseSD-([\d\.]+)_REP-(\d+)\.csv')
    
    #print(os.getcwd())  # use relative file path
    file_paths = [os.path.basename(x) for x in dex.glob('*.csv')]
    
    paras = []
    modes = []
    delts = []
    noiSD = []
    runss = []
    datas = []
    
    for fip in file_paths:
        match_obj = ptrn.match(fip)
        
        readIn = {
            "parameter" : match_obj[1],
            "mode" : match_obj[2],
            "change" : float(match_obj[3]),
            "noise-SD" : float(match_obj[4]),
            "runs" : int(match_obj[5]),
        }
        
        RtA = True # RtA: Return to Array
        # check if all explicit{…} conditions are met
        for key in explicit:
            if readIn[key] not in explicit[key]:
                        RtA = False
                    
        if RtA:
            paras.append(readIn["parameter"])
            modes.append(readIn["mode"])
            delts.append(readIn["change"])
            noiSD.append(readIn["noise-SD"])
            runss.append(readIn["runs"])
        
            # read the csv files
            datas.append(pd.read_csv(Path(dex/fip), index_col=0).stack())
    
    mindex = pd.MultiIndex.from_arrays([paras, modes, delts, noiSD, runss], 
                                 names=['transmitter', 'abs/rel', 'change', 'noise-SD', 'runs'])
    return datas, mindex
	
def mcLoad(subdir="", metabolite="GABA", return_format_data=False, return_SDs=False, explicit={}):
    ### uses mcGraf to load data, sorting for *metabolite*
	# cave: explicit{} can only have list values! 
	# e.g. expli = {"parameter":[para], "noise-SD":[n], …}
    loc = Path(subdir)
    datas, mindex = mcGraf(loc, explicit)
    data = pd.concat(datas, axis=1)
    data.columns = mindex
    format_data = data.unstack(level=1).stack(level="runs").droplevel(axis=0, level=1).reset_index(drop=True)
    format_data.columns = format_data.columns.set_names('metabolites', level=4)

    # list metabolite data
    index = format_data.melt().metabolites==metabolite
    to_plot = format_data.melt().loc[index,:]
    SDs=to_plot['noise-SD'].unique()
    
    return_tuple = index, to_plot            
    if return_format_data:
        return_tuple += (format_data,)
    if return_SDs:
        return_tuple += (SDs,)
    return return_tuple 
    
def mcLinFit(index, to_plot, metabolite, SDs, abszissa='change', ordinate='value', verbose=False):
    ### Linear fit on input data, per SD    
    # trendline  numpy.polyfit
    trend_a = []
    trend_b = []
    # individual fit per sd
    for s in SDs:
        # work copy
        ind = index.copy()
        pot = to_plot.copy()
        if verbose:
            print(pot)
        # only use data with specific noise
        ind = pot['noise-SD']==s
        pot = pot.loc[ind,:]
        # filter out nan
        ind = pot['value'].notna()
        pot = pot.loc[ind,:]
        # fitting (linear)
        x = np.array(pot['change'])
        y = np.array(pot['value'])
        t = np.polyfit(x, y, 1) # polyfit deprecated? np.polynomial?
        if verbose:
            print(s,t)
        trend_a.append(t[0])
        trend_b.append(t[1])
    fit  = pd.DataFrame(data={'noise-SD': SDs, 'slope': trend_a, 'offset':trend_b})
    return fit
    
def mcTexture(transmitter, transmitters=[], fits=[], to_plot_s=[], explicit={}, subdir=""):
    ### Appends three different data structs to given lists
    ### for quick drawing
    index, to_plot, format_data, SDs = mcLoad(metabolite=transmitter, 
                                              return_format_data=True, return_SDs=True, 
                                              explicit=explicit,
					      subdir=subdir)
    fit = mcLinFit(index, to_plot, transmitter, SDs)
    transmitters.append(transmitter)
    to_plot_s.append(to_plot)
    fits.append(fit)
    return transmitters, fits, to_plot_s
    
def mcMultiTexture(transmitter, transmitters=[], fits=[], to_plot_s=[], subdir=""):
    ### try to load two transmitters, add values, and return as mcTexture does.
    plot2, data2 = [],[]
    trans_name = str(transmitter[0])
    for trans in transmitter:
        index, plot, data, SDs = mcLoad(metabolite=trans, 
					return_format_data=True, return_SDs=True,
					subdir=subdir)
        plot2.append(plot)
        data2.append(data)
        if trans != transmitter[0]: #combine name
            trans_name += str("+"+trans)
    to_plot = plot2[0].copy()
    to_plot["value"]=plot2[0]["value"].values+plot2[1]["value"].values
    to_plot["metabolites"]=trans_name
    
    fit = mcLinFit(index, to_plot, transmitter, SDs)
    
    transmitters.append(trans_name)
    to_plot_s.append(to_plot)
    fits.append(fit)
        
    return transmitters, fits, to_plot_s

def mcPaint(index, to_plot, xlim=False, ylim=False):
    sns.lineplot(data=to_plot, x='change', y='value', hue='noise-SD', 
                  legend="full")

    if xlim:
        mtplot.xlim(xlim)
    if ylim:
        mtplot.ylim(ylim)
	 
