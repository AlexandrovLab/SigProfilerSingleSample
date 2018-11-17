#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 15:38:49 2018

@author: mishugeb
"""

from numpy import linalg as LA
import numpy as np
from scipy.optimize import minimize
import math
from scipy import spatial


#function to calculate the cosine similarity
def cos_sim(a, b):
	"""Takes 2 vectors a, b and returns the cosine similarity according 
	to the definition of the dot product
    
    Dependencies: 
        *Requires numpy library. 
        *Does not require any custom function (constructed by me)
        
    Required by:
        * pairwise_cluster_raw
	"""
	dot_product = np.dot(a, b)
	norm_a = np.linalg.norm(a)
	norm_b = np.linalg.norm(b)
	return dot_product / (norm_a * norm_b)
    
def parameterized_objective2_custom(x, signatures, samples):
    #print(signatures[:1])
    #print(samples[:1])
    rec = np.dot(signatures, x)
    y = LA.norm(samples-rec)
    #print(y)
    #print(x)
    #print("\n")
    return y


def constraints1(x, samples):
    sumOfSamples=np.sum(samples, axis=0)
    #print(sumOfSamples)
    #print(x)
    result = sumOfSamples-(np.sum(x))
    #print (result)
    return result


def create_bounds(idxOfZeros, samples, numOfSignatures):
    total = np.sum(samples)
    b = (0.0, float(total))
    lst =[b]*numOfSignatures
    
    for i in idxOfZeros:
        lst[i] = (0.0, 0.0) 
    
    return lst 



def add_signatures(W, genome):
    
    # This function takes an array of signature and a single genome as input, returns a dictionray of cosine similarity, exposures and presence 
    # of signatures according to the indices of the original signature array
    
    originalSimilarity = -1 # it can be also written as oldsimilarity
    maxmutation = round(np.sum(genome))
    init_listed_idx = []
    init_nonlisted_idx = list(range(W.shape[1]))
    finalRecord = [["similarity place-holder" ], ["newExposure place-holder"], ["signatures place-holder"]] #for recording the cosine difference, similarity, the new exposure and the index of the best signauture
    
    
    while True:
        bestDifference = -1 
        bestSimilarity = -1
        loopRecord = [["newExposure place-holder"], ["signatures place-holder"], ["best loop signature place-holder"]]
        for sig in init_nonlisted_idx:
            
            
            if len(init_listed_idx)!=0:
                loop_liststed_idx=init_listed_idx+[sig]
                loop_liststed_idx.sort()
                #print(loop_liststed_idx)
                W1 = W[:,loop_liststed_idx]
                #print (W1.shape)
                #initialize the guess
                x0 = np.random.rand(W1.shape[1], 1)*maxmutation
                x0= x0/np.sum(x0)*maxmutation
                
                #set the bounds and constraints
                bnds = create_bounds([], genome, W1.shape[1]) 
                cons1 ={'type': 'eq', 'fun': constraints1, 'args':[genome]} 
            # for the first time addintion  
            else:
                W1 = W[:,sig][:,np.newaxis]
                #print (W1.shape)        
                #initialize the guess
                x0 = np.ones((1,1))*maxmutation    
            
                #set the bounds and constraints
                bnds = create_bounds([], genome, 1) 
                cons1 ={'type': 'eq', 'fun': constraints1, 'args':[genome]} 
            
            #the optimization step
            sol = minimize(parameterized_objective2_custom, x0, args=(W1, genome),  bounds=bnds, constraints =cons1, tol=1e-30)
            
            #print(W1)
            #convert the newExposure vector into list type structure
            newExposure = list(sol.x)
            
            # get the maximum value of the new Exposure
            maxcoef = max(newExposure)
            idxmaxcoef = newExposure.index(maxcoef)
            
            newExposure = np.round(newExposure)
            
            # We may need to tweak the maximum value of the new exposure to keep the total number of mutation equal to the original mutations in a genome
            if np.sum(newExposure)!=maxmutation:
                newExposure[idxmaxcoef] = round(newExposure[idxmaxcoef])+maxmutation-sum(newExposure)
             
            # compute the estimated genome
            est_genome = np.dot(W1, newExposure)
            newSimilarity = cos_sim(genome, est_genome)
            
            difference = newSimilarity - originalSimilarity 
            
            # record the best values so far
            if difference>bestDifference:
                bestDifference = difference
                bestSimilarity = newSimilarity
                loopRecord = [newExposure, W1, sig]  #recording the cosine difference, the new exposure and the index of the best signauture
                #print(newSimilarity)
        
        # 0.01 is the thresh-hold for now 
        if bestSimilarity-originalSimilarity>0.011:
            originalSimilarity = bestSimilarity
            init_listed_idx.append(loopRecord[2])
            init_nonlisted_idx.remove(loopRecord[2])
            init_listed_idx.sort()
            #print(originalSimilarity)
            finalRecord = [originalSimilarity, loopRecord[0], init_listed_idx, loopRecord[1], genome]
            #print (finalRecord)
            
            if len(init_nonlisted_idx)!= 0:
                
                continue
            else:
                break
        else:
            break
        
    #print(finalRecord)
    return {"similarity":finalRecord[0], "exposures":finalRecord[1], "signatures": finalRecord[2]}  
     
if __name__ == "__main__":  
    
    # test the function
    mat = scipy.io.loadmat('21_breast_WGS_substitutions.mat')
    mat = extract_input(mat)
    genomes = mat[1]
    
    
    W = pd.read_csv('processes.txt', sep="\t")
    W = np.array(W)
    W = W[:,1:]
    
    
    for i in range(genomes.shape[1]):
        results = add_signatures(W, genomes[:, i])
        print (results["exposures"])
        print(results["signatures"])
        print('\n')
        
        
    
     
    
    
    
        
        
    for i in range(0,4):
        #print (i)
        #print ("Sample", i)
        Hnew[:,i], a, b, c, d = remove_all_single_signatures(W, H[:,i], genomes[:,i])
        #print (cos_sim(genomes[:,i], np.dot(W, H[:,i])))
        #print (cos_sim(genomes[:,i], np.dot(W, Hnew[:,i])))
        #print ("\n\n\n\n")
        
    
    results = add_signatures(W, c)
    
    print(Hnew[:,i])
    print(results[1])    
        
        