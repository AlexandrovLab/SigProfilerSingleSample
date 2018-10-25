# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 11:26:10 2018

@author: compactmatter
"""
import numpy as np
import scipy

def remove_one_sig(exposures, allSignatures, genome, saOptions):
    
    from eval_single_sample import eval_single_sample
    
    numberOfSignatures = int(sum((np.array([exposures])>0).ravel()))
    print(exposures.shape)
    sig_IDs = [i for i, x in enumerate(exposures) if x > 0]
    totalSignatures = allSignatures.shape[1]
    
    if ( numberOfSignatures > 1 ):
        exposuresSample_wo = np.zeros((totalSignatures, numberOfSignatures))
        accr_wo = np.zeros((numberOfSignatures, 1))
        kl_div_wo = np.zeros((numberOfSignatures, 1))
        frob_rel_div_wo = np.zeros((numberOfSignatures, 1))
        norm_one_dif_wo = np.zeros((numberOfSignatures, 1))
        
        for j in range(numberOfSignatures):
            exposuresRemoveOne = exposures

            # exposuresRemoveOne[sig_IDs[j]] = 0 #BEWARE code modification this line has been sent 2 lines below
            
            [exposuresSample_wo[:,j], accr_wo[j], kl_div_wo[j], frob_rel_div_wo[j], norm_one_dif_wo[j]] = eval_single_sample(exposuresRemoveOne, allSignatures, genome, saOptions) #BEWARE
            #[exposuresSample_wo[:,j], accr_wo[j], kl_div_wo[j], frob_rel_div_wo[j], norm_one_dif_wo[j]] = eval_single_sample(exposuresRemoveOne, allSignatures, genome, saOptions)
            
            #exposuresRemoveOne[sig_IDs[j]] = 0 
            #exposuresSample_wo[:,j] = exposuresSample_wo_j
        fVal = min(scipy.stats.hmean([1-accr_wo, frob_rel_div_wo])) #BEWARE
        fID = np.argmin(fVal)
        accr = accr_wo[fID]
        kl_div = kl_div_wo[fID]
        frob_rel_div = frob_rel_div_wo[fID]
        norm_one_dif = norm_one_dif_wo[fID]
        exposuresSample = exposuresSample_wo[:,fID]
        fID = sig_IDs[fID]
    else:
        [exposuresSample, accr, kl_div, frob_rel_div, norm_one_dif] = eval_single_sample(exposures, allSignatures, genome, saOptions)
        fID = 0

    return [accr, kl_div, frob_rel_div, norm_one_dif, exposuresSample, fID]