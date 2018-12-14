# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 18:49:40 2018

@author: compactmatter
"""

import numpy as np
import scipy
from scipy import optimize
from scipy.spatial.distance import pdist
from create_bounds import create_bounds
from constraints1 import constraints1

def eval_single_sample(exposures, allSignatures, genome, saOptions):

    from parameterized_objective2_custom import parameterized_objective2_custom
    from KLDiv import KLDiv
    
    if  ( int(sum(exposures.ravel())) > 0 ):
        exposuresSample = np.zeros((allSignatures.shape[1], 1))
        numberOfSignatures = int(sum((exposures>0).ravel()))
        maxMutations = sum(genome)
        
        x0 = maxMutations * np.random.rand(numberOfSignatures, 1)
        x0 = x0 / sum(x0) * maxMutations
        A = np.ones((1, numberOfSignatures ))
        lb = np.zeros(( numberOfSignatures , 1))
        ub = maxMutations * np.ones((numberOfSignatures, 1))

        
        subSignatures = allSignatures[: , (exposures > 0).ravel()]
        #print(x0.shape)
        #print(subSignatures.shape)
        #print(genome.shape)
        #set the bounds and constraints
        #bnds = create_bounds([], genome, 1) 
        cons1 ={'type': 'eq', 'fun': constraints1, 'args':[genome]} 
        #sprint(allSignatures)
        #ObjectiveFunction = @(x) parameterized_objective2_custom(x, subSignatures, genome);
        # print(x0)
        #print(sum(x0))
        x = (scipy.optimize.minimize(parameterized_objective2_custom, x0, args = (subSignatures, genome), bounds = [(0,maxMutations)]*len(x0), constraints = cons1, tol = 1e-30)).x
        
        #[x, minFunction] = fmincon(ObjectiveFunction, x0, [], [], A, maxMutations, lb, ub, [], saOptions) # simulannealbnd(ObjectiveFunction, X0, lb, ub, saOptions);
        x = [round(i) for i in x]
        if ( sum(x) != maxMutations):
             #[A B] = max(x);
             B = np.argmax(x) 
             #B = 0
             x[B] = x[B] + maxMutations - sum(x)
        #print(x)
        #print(len(x))
        #print(exposuresSample)
        exposuresSample[(exposures>0).ravel(),:] = np.transpose(np.array([x]))
        #print(exposuresSample)
        
        recon = np.matmul(allSignatures, exposuresSample)
        accr = 1 - pdist(np.transpose(np.concatenate((recon, np.transpose(np.array([genome]))), axis = 1)), 'cosine') #BEWARE
        
        kl_div = KLDiv(genome, recon)
        frob_rel_div = np.divide(np.linalg.norm(recon.ravel() - genome), np.linalg.norm(genome)) #Use np.divide
        norm_one_dif = np.linalg.norm(recon - genome, 1)
        #print(np.linalg.norm(recon - genome))
        #print(np.linalg.norm(genome))
        #print(genome)
        
        #print(accr)
        #print(frob_rel_div)
    else:
        #print('y')
        exposuresSample = exposures
        accr = 0
        kl_div = 0
        frob_rel_div = 1E-10 #BEWARE
        norm_one_dif = 0 

    return [exposuresSample.ravel(), accr, kl_div, frob_rel_div, norm_one_dif]
