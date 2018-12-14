# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 17:42:40 2018

@author: compactmatter
"""
import numpy as np
import scipy.stats

def add_one_sig(exposures_onesig, allSignatures, genome, saOptions):
    
    from eval_single_sample import eval_single_sample
    #print(exposures)
    numberOfSignatures = len(exposures_onesig) - sum(exposures_onesig>0)
    #print(numberOfSignatures)
    sig_IDs = [i for i, x in enumerate(exposures_onesig) if x == 0] #BEWARE code modification originally x == 0; #MATLAB CODE    sig_IDs = find(exposures==0)
    #print(sig_IDs)
    totalSignatures = allSignatures.shape[1]

    if ( numberOfSignatures > 1 ):
        exposuresSample_with = np.zeros((totalSignatures, numberOfSignatures))
        accr_with = np.zeros((numberOfSignatures, 1))
        kl_div_with = np.zeros((numberOfSignatures, 1))
        frob_rel_div_with = np.zeros((numberOfSignatures, 1))
        norm_one_dif_with = np.zeros((numberOfSignatures, 1))
        #print(exposures_onesig)
        for j in range(numberOfSignatures):
            exposuresAddOne = exposures_onesig
            exposuresAddOne[sig_IDs[j]] = 1
            #print(j)
            #print(exposuresAddOne)
            #print(exposures_onesig)
            [exposuresSample_with[:,j], accr_with[j], kl_div_with[j], frob_rel_div_with[j], norm_one_dif_with[j]] = eval_single_sample(exposuresAddOne, allSignatures, genome, saOptions)
            #print(exposuresSample_with[:,j])
        #print(exposuresSample_with.shape)
        #print(frob_rel_div_with)
        #print(1-accr_with)
        #print(1 - accr_with)
        fVal =  min(scipy.stats.hmean([(1-accr_with), frob_rel_div_with])) #BEWARE
        #print(scipy.stats.hmean([(1-accr_with), frob_rel_div_with]))
        fID = np.argmin(scipy.stats.hmean([(1-accr_with), frob_rel_div_with]))
        #fID = 0
        #MATLAB CODE [fVal, fID] = min(harmmean([1-accr_with frob_rel_div_with],2));
        accr = accr_with[fID]
        kl_div = kl_div_with[fID]
        frob_rel_div = frob_rel_div_with[fID]
        norm_one_dif = norm_one_dif_with[fID]
        exposuresSample = exposuresSample_with[:, fID]
        #print(exposuresSample_with.shape)
        #print(exposuresSample.shape)
        fID = sig_IDs[fID]
        #print('exposuresSample')
        #print(exposuresSample)
    else:
        [exposuresSample, accr, kl_div, frob_rel_div, norm_one_dif] = eval_single_sample(exposures_onesig, allSignatures, genome, saOptions)
        fID = 0
        
    return [accr, kl_div, frob_rel_div, norm_one_dif, exposuresSample, fID]
