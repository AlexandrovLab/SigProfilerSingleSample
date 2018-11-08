# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 13:33:06 2018

@author: compactmatter
"""
import numpy as np

from check_signature_rules import check_signature_rules
from remove_all_single_signatures import remove_all_single_signatures
from add_all_single_signatures import add_all_single_signatures
from eval_single_sample import eval_single_sample

def parallel_for_loop(iSample, inputSamples, allGenomeSignatures, allExomeSignatures, cancerType_file,
                      seqType_file, totalMutations_file, sampleNames, signaturesInSamples_file,
                      longSampleNames, totalSignatures, useRules, sigNames, idsToAdd, connected,
                      saOptions, exposuresNew, accr_org, kl_div_org, frob_rel_div_org, norm_one_dif_org):
    
    genome = inputSamples['originalGenomes'][:,iSample]
    if (inputSamples['seqType'][iSample] == 'WGS'):
        allSignatures =  allGenomeSignatures
    else:
        allSignatures =  allExomeSignatures
        
        #Identify signatures in cancer type and apply signature rules
    cancerType = open(cancerType_file,'r').read().split('\n')
    seqType = open(seqType_file,'r').read().split('\n')
    totalMutations = [float(i) for i in open(totalMutations_file,'r').read().split('\n')]               

    longSampleName = cancerType[iSample] + '::' + sampleNames[iSample]

    signaturesInSamples = np.loadtxt(signaturesInSamples_file)
        
    exposures = signaturesInSamples[:, longSampleNames.index(longSampleName)]
    if ( exposures.size == 0 ):
        exposures = np.zeros((totalSignatures, 1))
    
    strandBias = inputSamples['strandBias']
    
    if ( useRules == 1 ):
        exposures = check_signature_rules(exposures,
                                          sigNames,
                                          iSample,
                                          seqType,
                                          totalMutations,
                                          strandBias,
                                          'input/C_to_A_p.txt',
                                          'input/C_to_A_d.txt',
                                          'input/C_to_G_p.txt',
                                          'input/C_to_G_d.txt',
                                          'input/C_to_T_p.txt',
                                          'input/C_to_T_d.txt',
                                          'input/T_to_A_p.txt',
                                          'input/T_to_A_d.txt',
                                          'input/T_to_C_p.txt',
                                          'input/T_to_C_d.txt',
                                          'input/T_to_C_ATN_p.txt',
                                          'input/T_to_C_ATN_d.txt',
                                          'input/T_to_G_p.txt',
                                          'input/T_to_G_d.txt')
            
    #Remove all signatures in the sample: one by one
    #Add signatures that are allowed everywhere
    #if ( len(idsToAdd) > 0 ): #BEWARE modification to code in a few more ifs down as well
    if (idsToAdd != False):
        exposures[idsToAdd] = 1 
    
    #Add connected signatures
    if (len(connected) > 0):
        for iConnect in range(len(connected)):
            if ( sum(exposures[connected[iConnect]]) > 0 ):
                exposures[connected[iConnect]] = 1
    
    [exposuresOutput, accr, kl_div, frob_rel_div, norm_one_dif] = remove_all_single_signatures(exposures, allSignatures, genome, saOptions, sigNames)
       
    #Add all remaining signatures to the sample: one by one
    #Add signatures that are allowed everywhere
    if (idsToAdd != False):
        for i in idsToAdd:
            exposures[i] = 1
            
    #Add connected signatures
    if (len(connected) > 0):
        for iConnect in range(len(connected)):
            if ( sum(exposures[connected[iConnect]]) > 0 ):
                exposures[connected[iConnect]] = 1
                
    exposures = add_all_single_signatures(exposures, allSignatures, genome, saOptions, sigNames)
        
    #Add signatures that are allowed everywhere
    #Add signatures that are allowed everywhere
    if (idsToAdd != False):
        exposures[idsToAdd] = 1

    #Add connected signatures
    if (len(connected) > 0):
        for iConnect in range(len(connected)):
            if ( sum(exposures[connected[iConnect]]) > 0 ):
                exposures[connected[iConnect]] = 1
        
    [exposures, accr_org[iSample], kl_div_org[iSample], frob_rel_div_org[iSample], norm_one_dif_org[iSample]] = eval_single_sample(exposures, allSignatures, genome, saOptions)
    exposuresNew[:, iSample] = exposures

    #Dispay summary output data
    out_str = ['Sample #' + str(iSample) + ': ' + str(inputSamples['cancerType'][iSample]) + ' ' + str(inputSamples['sampleNames'][iSample]) + ' with ' + str(sum(exposuresNew[:, iSample]>0)) + ' signatures and an accuracy of ' + str(accr_org[iSample])] #BEWARE ,'%.2f' removed
    return(out_str)               
