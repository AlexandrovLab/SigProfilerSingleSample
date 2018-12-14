# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 11:22:45 2018

@author: compactmatter
"""
import numpy as np
import time
def remove_all_single_signatures(exposures_init, allSignatures, genome, saOptions, sigNames):

    from remove_one_sig import remove_one_sig
    
    from eval_single_sample import eval_single_sample
    # Check initial sample
    if ( sum(exposures_init) > 0):
        #print(exposures_init)
        [exposures, accr_first, kl_div_first, frob_rel_div_first, norm_one_dif_first] = eval_single_sample(exposures_init, allSignatures, genome, saOptions)
        #time.sleep(50) 
        #print(accr_first)
        #print(frob_rel_div_first)
        exposuresOutput = exposures     
        accr = accr_first 
        kl_div = kl_div_first 
        frob_rel_div = frob_rel_div_first
        norm_one_dif =  norm_one_dif_first

        # Removing singature one by one
        numSig = int(sum((np.array([exposures])>0).ravel()))
        for j in range(numSig - 1): #BEWARE code modification
        #for j in range(numSig):
            #print(exposuresOutput)
            exposures_temp = [i for i in exposuresOutput] #This step is done to ensure remove_one_sig does not modify exposuresOutput
            [accr_temp, kl_div_temp, frob_rel_div_temp, norm_one_dif_temp, exposuresSample_temp, fID] = remove_one_sig(exposuresOutput, allSignatures, genome, saOptions)
            exposuresOutput = np.asarray(exposures_temp)
            #print(accr_temp)
            if ( (accr_first - accr_temp < 0.01) ):
                exposuresOutput = exposuresSample_temp
                accr = accr_temp
                kl_div = kl_div_temp 
                frob_rel_div = frob_rel_div_temp
                norm_one_dif =  norm_one_dif_temp
            else: 
                break
    else:
      exposuresOutput = exposures_init
      accr = 0
      kl_div = 0
      frob_rel_div = 0
      norm_one_dif = 0
    
    return [exposuresOutput, accr, kl_div, frob_rel_div, norm_one_dif]
