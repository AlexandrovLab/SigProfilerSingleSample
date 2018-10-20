# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 17:29:41 2018

@author: compactmatter
"""

import numpy as np

def add_all_single_signatures(exposures_init, allSignatures, genome, saOptions, sigNames):
    
    from eval_single_sample import eval_single_sample
    from add_one_sig import add_one_sig
    #Check initial sample
    [exposuresOutput, accr, kl_div, frob_rel_div, norm_one_dif] = eval_single_sample(exposures_init, allSignatures, genome, saOptions)
    
    #Removing singature one by one
    numSig = len(exposuresOutput) - sum(exposuresOutput>0)
    for j in range (numSig):
        [accr_temp, kl_div_temp, frob_rel_div_temp, norm_one_dif_temp, exposuresSample_temp, fID] = add_one_sig(exposuresOutput, allSignatures, genome, saOptions)
               
        if ( (accr_temp - accr) > 0.05 ):
            exposuresOutput = exposuresSample_temp   
            accr = accr_temp; kl_div = kl_div_temp
            frob_rel_div = frob_rel_div_temp
            norm_one_dif =  norm_one_dif_temp
        else: 
            break
    
    #return [exposuresOutput, accr, kl_div, frob_rel_div, norm_one_dif] #BEWARE modification code
    return exposuresOutput