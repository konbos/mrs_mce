# graf.py - plotting sim.py
#
# Author: Konstantin E Bosbach <konstantin.bosbach@mars.uni-freiburg.de>


import seaborn as sb
import glob
import re
import os
import pandas as pd

def mcGraf(dex):
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
        
        # parameter from file-name for labeling data
        paras.append(match_obj[1])
        modes.append(match_obj[2])
        delts.append(float(match_obj[3]))
        noiSD.append(float(match_obj[4]))
        runss.append(int(match_obj[5]))
        
        # read the csv files
        datas.append(pd.read_csv(Path(dex/fip), index_col=0).stack())
    
    mindex = pd.MultiIndex.from_arrays([paras, modes, delts, noiSD, runss], 
                                 names=['transmitter', 'abs/rel', 'change', 'noise-SD', 'runs'])
    return datas, mindex

def mcLoad(subdir, metab="NAA"):
    datas, mindex = mcGraf(subdir)
    data = pd.concat(datas, axis=1)
    data.columns = mindex
    format_data = data.unstack(level=1).stack(level="runs").droplevel(axis=0, level=1).reset_index(drop=True)
    format_data.columns = format_data.columns.set_names('metabolites', level=4)

    # list data for metabolite
    index = format_data.melt().metabolites==metab
    to_plot = format_data.melt().loc[index,:]

    return index, to_plot

def mcPaint(index, to_plot, xlim, ylim):
    sns.lineplot(data=to_plot, x='change', y='value', hue='noise-SD', 
                  legend="full")

    if xlim:
        mtplot.xlim(xlim)
    if ylim:
        mtplot.ylim(ylim)
	 
