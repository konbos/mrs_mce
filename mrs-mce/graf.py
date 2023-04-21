# graf.py - plotting sim.py
#
# Author: Konstantin E Bosbach <konstantin.bosbach@mars.uni-freiburg.de>


import seaborn as sb
import glob
import re
import os

def mcGraf(dex):
    dfs = pd.DataFrame()
    ptrn = re.compile(r'(\w{0,5})-(\w{3,})~([-\d\.]+)_noiseSD-([\d\.]+)_REP-(\d+)\.csv')
    
    #print(os.getcwd())
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
