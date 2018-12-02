"""
Created on Thu Nov  8 15:38:49 2018

@author: mishugeb
"""

import xml.etree.ElementTree as ET  
from lxml import etree
from numpy import linalg as LA
import numpy as np
from scipy.optimize import minimize
import math
from scipy import spatial
import multiprocessing as mp


# Functions to calculate the cosine similarity, parameters, constrains and bounds
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




#################################################################### Function to check signatures rules ########################################################


def check_signature_rules(signaturesInSample,
                          signatures,
                          sampleID,
                          seqType,
                          totalMutations,
                          strandBias,
                          C_to_A_p_file,
                          C_to_A_d_file,
                          C_to_G_p_file,
                          C_to_G_d_file,
                          C_to_T_p_file,
                          C_to_T_d_file,
                          T_to_A_p_file,
                          T_to_A_d_file,
                          T_to_C_p_file,
                          T_to_C_d_file,
                          T_to_C_ATN_p_file,
                          T_to_C_ATN_d_file,
                          T_to_G_p_file,
                          T_to_G_d_file,
                          signatureRulesXML,
                          signatureRulesXMLSchema):
    
    totalSignatures = len(signaturesInSample)  
    
    C_to_A_p = [float(i) for i in open(C_to_A_p_file,'r').read().split('\n')]
    C_to_A_d = [float(i) for i in open(C_to_A_d_file,'r').read().split('\n')]
    C_to_G_p = [float(i) for i in open(C_to_G_p_file,'r').read().split('\n')]
    C_to_G_d = [float(i) for i in open(C_to_G_d_file,'r').read().split('\n')]
    C_to_T_p = [float(i) for i in open(C_to_T_p_file,'r').read().split('\n')]
    C_to_T_d = [float(i) for i in open(C_to_T_d_file,'r').read().split('\n')]
    T_to_A_p = [float(i) for i in open(T_to_A_p_file,'r').read().split('\n')]
    T_to_A_d = [float(i) for i in open(T_to_A_d_file,'r').read().split('\n')]
    T_to_C_p = [float(i) for i in open(T_to_C_p_file,'r').read().split('\n')]
    T_to_C_d = [float(i) for i in open(T_to_C_d_file,'r').read().split('\n')]
    T_to_C_ATN_p = [float(i) for i in open(T_to_C_ATN_p_file,'r').read().split('\n')]
    T_to_C_ATN_d = [float(i) for i in open(T_to_C_ATN_d_file,'r').read().split('\n')]
    T_to_G_p = [float(i) for i in open(T_to_G_p_file,'r').read().split('\n')]
    T_to_G_d = [float(i) for i in open(T_to_G_d_file,'r').read().split('\n')]

    tree = ET.parse(signatureRulesXML)  
    
    #Validate XML file against schema
    xml = etree.parse(signatureRulesXML)

    schema_doc = etree.parse(signatureRulesXMLSchema)
    schema = etree.XMLSchema(schema_doc)

    if (schema.validate(xml) == False):
        print(schema.assertValid(xml))
    
    root = tree.getroot()

    signaturesList = []    
    for i in range(totalSignatures):
        for signature in root.findall('signatureName'):
            signaturesList.append(signature.get("signatureSBS"))
            
        if(signatures[i] in signaturesList):
            if (root.find(".//strandBias/..[@signatureSBS='" + signatures[i] + "']") is not None):
                strandbias = root.find(".//strandBias/..[@signatureSBS='" + signatures[i] + "']").find('strandBias').findall('mutationType')
                for sb in range(len(strandbias)):
                    if( eval(strandbias[sb].get('mutType').replace('>','_to_') + "_p[sampleID]") > float(strandbias[sb][1].text) or
                       eval(strandbias[sb].get('mutType').replace('>','_to_') + "_d[sampleID]") != float(strandbias[sb][0].text)):
                       signaturesInSample[i] = 0
                       
            if (root.find(".//totalMutations/..[@signatureSBS='" + signatures[i] + "']") is not None):
                totalmutations = root.find(".//totalMutations/..[@signatureSBS='" + signatures[i] + "']").find('totalMutations').findall('seqType')
                for st in range(len(totalmutations)):
                    if( totalmutations[st].get('type') == seqType[sampleID]):
                        if( float(totalMutations[sampleID]) < float(totalmutations[st][0].text) ):
                            signaturesInSample[i] = 0
                
    return signaturesInSample

"""EXAMPLE

# Test the "check_signature_rules" function with an example:

Take the necessary input data

result = check_signature_rules(signaturesInSample,
                          signatures,
                          sampleID,
                          seqType,
                          totalMutations,
                          strandBias,
                          C_to_A_p_file,
                          C_to_A_d_file,
                          C_to_G_p_file,
                          C_to_G_d_file,
                          C_to_T_p_file,
                          C_to_T_d_file,
                          T_to_A_p_file,
                          T_to_A_d_file,
                          T_to_C_p_file,
                          T_to_C_d_file,
                          T_to_C_ATN_p_file,
                          T_to_C_ATN_d_file,
                          T_to_G_p_file,
                          T_to_G_d_file,
                          signatureRulesXML,
                          signatureRulesXMLSchema)
print(results)
"""


#################################################################### Function to get the initial Exposures with all the signatures #############################
def initial_optimization(processes, genomes):
    exposure = np.zeros([processes.shape[1], genomes.shape[1]] )
    # get the total mutations for the given sample
    maxmutation = np.round(np.sum(genomes, axis=0))
    
    for i in range(genomes.shape[1]):
        Gi = genomes[:,i]
    
        #initialize the guess
        x0 = np.random.rand(processes.shape[1], 1)*maxmutation[i]
        x0= x0/np.sum(x0)*maxmutation[i]
        
        #set the bounds and constraints
        bnds = create_bounds([], Gi, processes.shape[1]) 
        cons1 ={'type': 'eq', 'fun': constraints1, 'args':[Gi]} 
        
        #the optimization step
        sol = minimize(parameterized_objective2_custom, x0, args=(processes, Gi),  bounds=bnds, constraints =cons1, tol=1e-15)
        exposure[:,i]=sol.x
        
    return exposure

"""
# EXAMPLE 1
#Test the "initial_optimization" function with an example:
processes = np.random.rand(96, 4)
sample = np.random.rand(96, 1)
# Here, we get the initial optimized exposure for the processes on the sample
exposure = initial_optimization(processes, sample)
print(exposure)
"""

#################################################################### Function to remomove signatures from samples  #############################
def remove_signatures(indices, W, exposures, totoalgenomes):
    i = indices
    H = exposures[:,i]
    genomes= totoalgenomes[:,i]
    # make the empty list of the successfull combinations
    successList = [0,[],0] 
    # get the cos_similarity with sample for the oringinal W and H[:,i]
    originalSimilarity= cos_sim(genomes, np.dot(W, H))
    # make the original exposures of specific sample round
    oldExposures = np.round(H)
    
    # set the flag for the while loop
    if len(oldExposures[np.nonzero(oldExposures)])>1:
        Flag = True
    else: 
        Flag = False
        return oldExposures
    # The while loop starts here
    while Flag: 
        
        # get the list of the indices those already have zero values
        if len(successList[1]) == 0:
            initialZerosIdx = list(np.where(oldExposures==0)[0]) 
            #get the indices to be selected
            selectableIdx = list(np.where(oldExposures>0)[0]) 
        elif len(successList[1]) > 1: 
            initialZerosIdx = list(np.where(successList[1]==0)[0]) 
            #get the indices to be selected
            selectableIdx = list(np.where(successList[1]>0)[0])
        else:
            print("iteration is completed")
            #break
        
        
        # get the total mutations for the given sample
        maxmutation = round(np.sum(genomes))
        
        # new signature matrix omiting the column for zero
        #Winit = np.delete(W, initialZerosIdx, 1)
        Winit = W[:,selectableIdx]
        
        # set the initial cos_similarity
        record  = [0.11, []] 
        # get the number of current nonzeros
        l= Winit.shape[1]
        
        for i in range(l):
            #print(i)
            loopSelection = list(range(l))
            del loopSelection[i]
            #print (loopSelection)
            W1 = Winit[:,loopSelection]
           
            
            #initialize the guess
            x0 = np.random.rand(l-1, 1)*maxmutation
            x0= x0/np.sum(x0)*maxmutation
            
            #set the bounds and constraints
            bnds = create_bounds([], genomes, W1.shape[1]) 
            cons1 ={'type': 'eq', 'fun': constraints1, 'args':[genomes]} 
            
            #the optimization step
            sol = minimize(parameterized_objective2_custom, x0, args=(W1, genomes),  bounds=bnds, constraints =cons1, tol=1e-15)
            
            #print (sol.success)
            #print (sol.x)
            
            #convert the newExposure vector into list type structure
            newExposure = list(sol.x)
            
            #insert the loopZeros in its actual position 
            newExposure.insert(i, 0)
            
            #insert zeros in the required position the newExposure matrix
            initialZerosIdx.sort()
            for zeros in initialZerosIdx:
                newExposure.insert(zeros, 0)
            
            # get the maximum value the new Exposure
            maxcoef = max(newExposure)
            idxmaxcoef = newExposure.index(maxcoef)
            
            newExposure = np.round(newExposure)
            
            
            if np.sum(newExposure)!=maxmutation:
                newExposure[idxmaxcoef] = round(newExposure[idxmaxcoef])+maxmutation-sum(newExposure)
                
            
            newSample = np.dot(W, newExposure)
            newSimilarity = cos_sim(genomes, newSample) 
             
            difference = originalSimilarity - newSimilarity
            #print(originalSimilarity)
            #print(newSample)
            #print(newExposure)
            #print(newSimilarity)
            
            #print(difference)
            #print (newExposure)
            #print (np.round(H))
            #print ("\n\n")
             
            if difference<record[0]:
                record = [difference, newExposure, newSimilarity]
            
            
        #print ("This loop's selection is {}".format(record))
        
        if record[0]>0.01:   
            Flag=False
        elif len(record[1][np.nonzero(record[1])])==1:
            successList = record 
            Flag=False
        else:
            successList = record
        #print("The loop selection is {}".format(successList))
        
        #print (Flag)
        #print ("\n\n")
    
    #print ("The final selection is {}".format(successList))
    
    if len(successList[1])==0:
        successList = [0.0, oldExposures, originalSimilarity]
    
    #print ("one sample completed")
    return successList[1]

"""
# EXAMPLE 2
#Test the "remove_signatures" function with an example:
# Continue the processes, sample and exposure varialbes from the "initial_optimization" function example (EXAMPLE 1)

pool = mp.Pool(processes=1)
results = [pool.apply_async(remove_signatures, args=(x,processes,exposure,sample,)) for x in range(sample.shape[1])]
output = [p.get() for p in results]
#print(results)

exposureAfterRemovingSignatures = np.zeros(exposure.shape)
for i in range(len(output)):
    #print(results[i])
    exposureAfterRemovingSignatures[:,i]=output[i] 
print(exposureAfterRemovingSignatures)
"""

#################################################################### Function to add signatures to samples from database #############################
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
            newSimilarity = cos_sim(genome[:,0], est_genome)
            
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
    dictExposure= {"similarity":finalRecord[0], "exposures":finalRecord[1], "signatures": finalRecord[2]}  
    addExposure = np.zeros([W.shape[1]])
    addExposure[dictExposure["signatures"]]=dictExposure["exposures"]
    
    return  addExposure

"""
# EXAMPLE 3
#Test the "add_signatures" function with an example:
# Continue the processes, sample varialbes from the "initial_optimization" function example (EXAMPLE 1)

exposureAfterAddingSignatures = add_signatures(processes, sample)
print(exposureAfterAddingSignatures)

"""




if __name__ == "__main__":

    ################################################################################################################################################################
    ######################################################## The SigProfilerSingleSample codes starts here #########################################################
    ################################################################################################################################################################
    
    # INPUTS:
        # signaturesInSample = database of Signatures
        # sample  : Data for a single sample: The list should contain single point mutations vector, double point mutations vector, indel vector, strand_bias, etc.
        # Rules: An input file containing all the rules         
    
    ##################### The check_signature_rules function will go here ##################################################################################################### 
    """
    result = check_signature_rules(signaturesInSample, sample, rules,..... )
    
    result will be used in the downstream analysis
    
    """
    
    
    
    
    #################################################################### get the initial exposures from a given set signatures and genome ##########################
    exposure =   initial_optimization(processes, sample)      
    ################################################################################################################################################################       
    
    ######################################################### remove signatures ####################################################################################
    pool = mp.Pool(processes=1)
    results = [pool.apply_async(remove_signatures, args=(x,processes,exposure,sample,)) for x in range(sample.shape[1])]
    output = [p.get() for p in results]
    #print(results)
    
    exposureAfterRemovingSignatures = np.zeros(exposure.shape)
    for i in range(len(output)):
        #print(results[i])
        exposureAfterRemovingSignatures[:,i]=output[i] 
    print(exposureAfterRemovingSignatures)
    #################################################################################################################################################################
    
    
    #################################################################### add signatures #############################################################################
    exposureAfterAddingSignatures = add_signatures(processes, sample)
    print(exposureAfterAddingSignatures)
    #################################################################################################################################################################    
