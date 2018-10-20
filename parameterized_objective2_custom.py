# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 11:12:50 2018

@author: compactmatter
"""

#%% The parameterization of objective function
import numpy as np

def parameterized_objective2_custom(x, signatures, samples):

    rec = np.matmul(signatures, x) # reconstructed genome using signatures and their activities
      
    #y = (samples - rec).max() #BEWARE        # L^2-Norm 
    y = np.linalg.norm(samples - rec,2)
    # MATLAB CODE y = norm(samples - rec, 2);        # L^2-Norm 
#     y = KLDiv(samples', rec');         # Kullback?Leibler divergence
#     y = norm(samples - rec, 'fro');    # Frobenius Norm

    return y