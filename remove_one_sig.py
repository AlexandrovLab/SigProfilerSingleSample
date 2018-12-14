# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 11:26:10 2018

@author: compactmatter
"""
import numpy as np
import scipy.stats

def remove_one_sig(exposures, allSignatures, genome, saOptions):
    
    from eval_single_sample import eval_single_sample
    
    numberOfSignatures = int(sum((np.array([exposures])>0).ravel()))
    sig_IDs = [i for i, x in enumerate(exposures) if x > 0]
    totalSignatures = allSignatures.shape[1]
    #print(sig_IDs )
    if ( numberOfSignatures > 1 ):
        exposuresSample_wo = np.zeros((totalSignatures, numberOfSignatures))
        accr_wo = np.zeros((numberOfSignatures, 1))
        kl_div_wo = np.zeros((numberOfSignatures, 1))
        frob_rel_div_wo = np.zeros((numberOfSignatures, 1))
        norm_one_dif_wo = np.zeros((numberOfSignatures, 1))
        #print(exposures.shape)
        #exposuresRemoveOne = np.zeros((exposures.shape[0],)) #BEWARE code addition
        #print(exposures)
        for j in range(numberOfSignatures):
            
            exposuresRemoveOne = np.asarray([i for i in exposures]).ravel()
            #print(exposuresRemoveOne)
            #print(exposures)
            exposuresRemoveOne[sig_IDs[j]] = 0 #BEWARE code modification this line has been sent 2 lines below
            #print(j)
            #print(exposuresRemoveOne)
            [exposuresSample_wo[:,j], accr_wo[j], kl_div_wo[j], frob_rel_div_wo[j], norm_one_dif_wo[j]] = eval_single_sample(exposuresRemoveOne, allSignatures, genome, saOptions) #BEWARE
            #[exposuresSample_wo[:,j], accr_wo[j], kl_div_wo[j], frob_rel_div_wo[j], norm_one_dif_wo[j]] = eval_single_sample(exposuresRemoveOne, allSignatures, genome, saOptions)
            
            #exposuresRemoveOne[sig_IDs[j]] = 0 
            #exposuresSample_wo[:,j] = exposuresSample_wo_j
        #print(accr_wo)
        #print(frob_rel_div_wo)
        #if frob_rel_div_wo[numberOfSignatures] == 0:
        #    frob_rel_div_wo[numberOfSignatures] = 1E-10
        
        fVal = min(scipy.stats.hmean([1-accr_wo, frob_rel_div_wo])) #BEWARE
        fID = np.argmin(scipy.stats.hmean([1-accr_wo, frob_rel_div_wo]))
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
