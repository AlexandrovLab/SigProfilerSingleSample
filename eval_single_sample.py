# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 18:49:40 2018

@author: compactmatter
"""

import numpy as np
import scipy
from scipy import optimize
from scipy.spatial.distance import pdist


def eval_single_sample(exposures, allSignatures, genome, saOptions):

    from parameterized_objective2_custom import parameterized_objective2_custom
    from KLDiv import KLDiv
    
    if  ( int(sum(exposures.ravel())) > 0 ):
        #print(exposures.ravel())
        exposuresSample = np.zeros((allSignatures.shape[1], 1))
        numberOfSignatures = int(sum((exposures>0).ravel()))
        maxMutations = sum(genome)
        
        x0 = maxMutations * np.random.rand(numberOfSignatures, 1)
        x0 = x0 / sum(x0) * maxMutations
        A = np.ones((1, numberOfSignatures ))
        lb = np.zeros(( numberOfSignatures , 1))
        ub = maxMutations * np.ones((numberOfSignatures, 1))

        #print(allSignatures.shape) #DEBUG
        #print((exposures > 0).shape) #DEBUG
        
        subSignatures = allSignatures[: , (exposures > 0).ravel()]
        
        #ObjectiveFunction = @(x) parameterized_objective2_custom(x, subSignatures, genome);
        x = (scipy.optimize.minimize(parameterized_objective2_custom, x0, args = (subSignatures, genome))).x
        #[x, minFunction] = fmincon(ObjectiveFunction, x0, [], [], A, maxMutations, lb, ub, [], saOptions) # simulannealbnd(ObjectiveFunction, X0, lb, ub, saOptions);
        #x = [round(i) for i in x] #BEWARE
        if ( sum(x) != maxMutations):
             #[A B] = max(x);
             B = np.argmax(x) 
             #B = 0
             x[B] = x[B] + maxMutations - sum(x)

        #print(exposuresSample.shape) DEBUG
        #print((exposures > 0).shape) DEBUG
        #print(np.transpose(np.array([x])).shape) DEBUG
        
        exposuresSample[(exposures>0).ravel(),:] = np.transpose(np.array([x]))
        
        #recon = allSignatures * exposuresSample 
        recon = np.matmul(allSignatures, exposuresSample)
        accr = 1 - pdist(np.transpose(np.concatenate((recon, np.transpose(np.array([genome]))), axis = 1)), 'cosine') #BEWARE
        
        #accr = 0
        #print(np.array([genome]).shape)
        #print(recon.shape)
        kl_div = KLDiv(genome, recon)
        #frob_rel_div = np.divide(np.linalg.norm(recon - genome,'fro'), np.linalg.norm(genome,'fro')) #Use np.divide #BEWARE
        frob_rel_div = np.divide(np.linalg.norm(recon - genome), np.linalg.norm(genome)) #Use np.divide
        norm_one_dif = np.linalg.norm(recon - genome, 1)
        
    else:
        exposuresSample = exposures
        accr = 0
        kl_div = 0
        frob_rel_div = 0
        norm_one_dif = 0 
        
    #print(exposuresSample.shape)
    #print(accr.shape)
    #print(kl_div.shape)
    #print(frob_rel_div.shape)
    #print(norm_one_dif.shape)
    return [exposuresSample.ravel(), accr, kl_div, frob_rel_div, norm_one_dif]