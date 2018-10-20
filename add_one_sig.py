# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 17:42:40 2018

@author: compactmatter
"""
import numpy as np
import scipy

def add_one_sig(exposures, allSignatures, genome, saOptions):
    
    from eval_single_sample import eval_single_sample
    
    numberOfSignatures = len(exposures) - sum(exposures>0)
    sig_IDs = [i for i, x in enumerate(exposures) if x <= 0] #BEWARE code modification originally x == 0; #MATLAB CODE    sig_IDs = find(exposures==0)
    totalSignatures = allSignatures.shape[1] #MATLAB CODE size(allSignatures, 2)
    #print(exposures)
    #print(len(exposures))
    #print(sum(exposures > 0))
    #print(numberOfSignatures)
    #print(sig_IDs)
    #print(len(sig_IDs))
    if ( numberOfSignatures > 1 ):
        exposuresSample_with = np.zeros((totalSignatures, numberOfSignatures))
        accr_with = np.zeros((numberOfSignatures, 1), dtype = int) #MATLAB CODE accr_with = zeros(numberOfSignatures, 1);
        kl_div_with = np.zeros((numberOfSignatures, 1), dtype = int) #MATLAB CODE kl_div_with = zeros(numberOfSignatures, 1);
        frob_rel_div_with = np.zeros((numberOfSignatures, 1), dtype = int) #MATLAB CODE frob_rel_div_with = zeros(numberOfSignatures, 1);
        norm_one_dif_with = np.zeros((numberOfSignatures, 1), dtype = int) #MATLAB CODE norm_one_dif_with = zeros(numberOfSignatures, 1);
        
        for j in range(numberOfSignatures):
            exposuresAddOne = exposures
            exposuresAddOne[sig_IDs[j]] = 1
            #print(j)
            #print(exposuresAddOne.shape)
            #print(sig_IDs[j])
            [exposuresSample_with[:,j], accr_with[j], kl_div_with[j], frob_rel_div_with[j], norm_one_dif_with[j]] = eval_single_sample(exposuresAddOne, allSignatures, genome, saOptions)

        fVal =  min(scipy.stats.hmean([(1-accr_with), frob_rel_div_with])) #BEWARE
        fID = np.argmin(fVal)
        #MATLAB CODE [fVal, fID] = min(harmmean([1-accr_with frob_rel_div_with],2));
        accr = accr_with[fID]
        kl_div = kl_div_with[fID]
        frob_rel_div = frob_rel_div_with[fID]
        norm_one_dif = norm_one_dif_with[fID]
        exposuresSample = exposuresSample_with[:, fID]
        fID = sig_IDs[fID]
    else:
        [exposuresSample, accr, kl_div, frob_rel_div, norm_one_dif] = eval_single_sample(exposures, allSignatures, genome, saOptions)
        fID = 0
        
    return [accr, kl_div, frob_rel_div, norm_one_dif, exposuresSample, fID]