# sim.py - Monte-Carlo simulation on fsl_mrs
#
# Author: Konstantin E Bosbach <konstantin.bosbach@mars.uni-freiburg.de> 

import numpy as np
import pandas as pd
import time

from fsl_mrs.core import MRS
from fsl_mrs.utils import mrs_io
from fsl_mrs.utils.synthetic import synthetic_from_basis as synth
from pathlib import Path

def synth_and_ana(noise_sd, fit_parameter, fit_snr, fit_sdnoise,
                  syn_parameter_dict, basis_path, printIm = None):
    """Synthetise spectra with given noise-covariance, analyse the data and
    return table of fitting parameters with a fit-plot."""
    
    # #Generate synthetic data
    fidS, mrsA, concentrationsS = synth.syntheticFromBasisFile(
        basis_path,
        concentrations=syn_parameter_dict,
        # correct for complex noise
        noisecovariance=[[np.divide(np.square(noise_sd), 2)]],
        bandwidth=6000,
        ind_scaling=['Mac'],
        broadening=(syn_parameter_dict['gamma_0'], 0.0)
    )
    metab_groups = mrsA.parse_metab_groups("combine_all")
    mrsA.FID = fidS

    ## FIT: Voigt line broadening, between .2 and 4.2,
    #  with a 2nd order polynomial baseline
    FitArgs = {
        'model': 'voigt',
        'metab_groups': metab_groups,
        'ppmlim': (.2, 4.2),
        'baseline_order': 2}
    res = mrsA.fit(**FitArgs)
    
    # Combine highly correlated metabolites
    combinationList = [['Glu', 'Gln'],
                       ['GPC', 'PCh'],    #changed "PCho"->"PCh"
                       ['Cr', 'PCr'],
                       ['Glc', 'Tau'],
                       ["NAA", "NAAG"]]
    res.combine(combinationList)
    
    ## Write results
    noise_sd = np.std(mrsA.FID[1000:1600])   
    fit_sdnoise.append(noise_sd)
    fit_snr.append(res.SNR.spectrum) 
    fit_parameter.append(res.fitResults)
    
    # generate result-df
    df_params = pd.concat(fit_parameter, ignore_index=True)
    df_params['SNR'] = fit_snr
    df_params['noise_sd'] = fit_sdnoise
    
    if printIm is not None:
        res.plot(mrsA, out=printIm)
        
    return df_params

def mc(noise_sd, syn_parameter_dict, basis_path, 
       file_out_path, n, shortage=True):
    """Function for calling synth_and_ana repeatedly,
     as in Monte-Carlo approach"""
    runtime = time.time()       # timer for inspection

    if shortage:
        try:
            noise_fit = pd.read_csv(file_out_path)
            print("Found existing file/s, latest: ", file_out_path, end='\r')
            return noise_fit, file_out_path
        except:
            None

    print(
        "Starting  noise_sd", round(noise_sd, 2), " with ",
        round(n, 2), "repetitions", end='\r'
    )
    
    # structs used to gather input back
    fit_parameter=[]
    fit_snr=[]
    fit_sdnoise=[]
    fit_varnoise=[]

    # Call function generation the desired amount of times
    ## first run save image
    noise_fit = synth_and_ana(
            noise_sd, 
            fit_parameter=fit_parameter,
            fit_snr=fit_snr, fit_sdnoise=fit_sdnoise,
            syn_parameter_dict=syn_parameter_dict, 
            basis_path=basis_path,
            printIm = Path(str(file_out_path)+"-precusor.png")
        )
    
    for k in range(0, n-1):
        noise_fit = synth_and_ana(
            noise_sd, 
            fit_parameter=fit_parameter,
            fit_snr=fit_snr, fit_sdnoise=fit_sdnoise,
            syn_parameter_dict=syn_parameter_dict, 
            basis_path=basis_path
        )

    # Print results to file
    noise_fit.to_csv(file_out_path)

    print("Finishing noise_sd", round(noise_sd, 2), " with ",
          round(n, 2), "repetitions, Runtime took ",
          round(time.time()-runtime, 2), "[s]")

    return noise_fit, file_out_path

def mcCall(n, noise_sd, para="NAA", step=[-0.1, 1, 10], absolute=False, output="default"):
    '''
    Pipeline mc call : get simulation running in folder
	assumes : output in folder of script 
	and: workspace/[fsl_mrs_mce]/[experiment X]/[script.ipynb.py]
	and: workspace/fit_conc_result.csv     # input data
	and: workspace/basis-steam/  	       # spectrum basis
    Helper code
    '''
    # path
    workspace_path =  Path('../../')
    basis_path = Path(workspace_path / 'basis-steam')
    output_path = Path('./')
    
    #read in data
    csv_path = Path(workspace_path / 'fit_conc_result.csv')
    df_parameter_synth = pd.read_csv(csv_path, index_col=0)
    syn_parameter_dic = df_parameter_synth.mean().to_dict()
    syn_parameter_dic['gamma_0'] = 24.49517127663636       # in vivo
    syn_parameter_dic['PCh'] = syn_parameter_dic['PCho']
    syn_parameter_dic['Mac'] = syn_parameter_dic['mm']
        
    # init mc
    for s in step:
        spd = dict(syn_parameter_dic)
        if absolute:
            spd[para] = s
            context = "abs"
        else:
            spd[para] *= s
            context = "rel"

        file_name = str(para+"-"+context+"~"+str(round(s, 3))
                        +"_"+"noiseSD-"+str(round(noise_sd, 3))
                        +"_"+"REP-"+str(n)+".csv")
        file_out_path = Path(output_path / file_name)
        
        # call mc
        noise_fit, out_path = mc(noise_sd, spd, basis_path, 
           file_out_path, n, shortage=True)     
    return
