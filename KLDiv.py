# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 10:51:36 2018

@author: compactmatter
"""

import numpy as np

def KLDiv(P,Q):
# dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
# distributions
# P and Q  are automatically normalised to have the sum of one on rows
# have the length of one at each 
# P =  n x nbins
# Q =  1 x nbins or n x nbins(one to one)
# dist = n x 1

    P = np.transpose(np.array([P])) #BEWARE new line of code added

    if P.shape[1] != Q.shape[1]:
        print('the number of columns in P and Q should be the same')

    if (np.all(~np.isfinite(P[:]))):
        print('the inputs contain non-finite values!') 
        
    if (np.all(~np.isfinite(Q[:]))):
        print('the inputs contain non-finite values!') 
        
# normalizing the P and Q
    if Q.shape[0] == 1:
        Q = np.array(Q) /sum(Q)
        P = np.array(P) /np.matlibrepmat(sum(P), 1, P.shape[1]) #BEWARE sum(P,2) 2 has been removed
        temp =  P* np.log(P /np.matlibrepmat(Q, P.shape[0], 1))
        temp[np.isnan(temp)]=0    # resolving the case when P(i)==0
        dist = sum(temp)

    else:
        Q = Q /np.matlib.repmat(sum(Q),1, Q.shape[1])
        P = P /np.matlib.repmat(sum(P),1, P.shape[1])
        for p_i in range(len(P)):
            if P[p_i] == 0:
                P[p_i] = Q[p_i]
        temp =  P * np.log(np.divide(P,Q))
        temp[np.isnan(temp)]=0 # resolving the case when P(i)==0
        dist = sum(temp)
        
    return dist
