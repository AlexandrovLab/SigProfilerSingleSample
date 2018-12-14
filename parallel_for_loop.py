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
    
    #print(sum(exposures))
    if ( useRules == 1 ):
        exposures = check_signature_rules(exposures,
                                          sigNames,
                                          iSample,
                                          seqType,
                                          totalMutations,
                                          #strandBias,
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
                                          'input/T_to_G_d.txt',
                                          'input/20181108_Signature_Rules.xml',
                                          'input/Signature_Rules_Schema.xsd')
    
    #print('signaturesInsample')
    #print(exposures)
    #print('signatures')
    #print(sigNames)
    #print('sampleID')
    #print(iSample)
    #print('seqType')
    #print(seqType)
    #print('totalMutations')
    #print(totalMutations)
    #print('strandBias')
    #print(strandBias)
       
    #Remove all signatures in the sample: one by one
    #Add signatures that are allowed everywhere
    #if ( len(idsToAdd) > 0 ): #BEWARE modification to code in a few more ifs down as well
    
    if (idsToAdd != False):
        exposures[idsToAdd] = 1 
    #print(exposures)
    
    #Add connected signatures
    if (len(connected) > 0):
        for iConnect in range(len(connected)):
            if ( sum(exposures[connected[iConnect]]) > 0 ):
                exposures[connected[iConnect]] = 1
    
    #[exposuresOutput, accr, kl_div, frob_rel_div, norm_one_dif] = remove_all_single_signatures(exposures, allSignatures, genome, saOptions, sigNames)
    #print(exposures)
    [exposures, accr, kl_div, frob_rel_div, norm_one_dif] = remove_all_single_signatures(exposures, allSignatures, genome, saOptions, sigNames)
    #print(exposures)
    #print(exposures)   
    #Add all remaining signatures to the sample: one by one
    #Add signatures that are allowed everywhere
    if (idsToAdd != False):
        for i in idsToAdd:
            exposures[i] = 1
    #print(exposures)        
    #Add connected signatures
    if (len(connected) > 0):
        for iConnect in range(len(connected)):
            if ( sum(exposures[connected[iConnect]]) > 0 ):
                exposures[connected[iConnect]] = 1
    #print(exposures)      
    #print(sum(exposures)) #BEWARE DIFFERENT OUTPUT AS COMPARED TO MATLAB      
    exposures = add_all_single_signatures(exposures, allSignatures, genome, saOptions, sigNames)
    #print(sum(exposures))    
    #Add signatures that are allowed everywhere
    #Add signatures that are allowed everywhere
    #print(exposures)
    if (idsToAdd != False):
        exposures[idsToAdd] = 1

    #Add connected signatures
    if (len(connected) > 0):
        for iConnect in range(len(connected)):
            if ( sum(exposures[connected[iConnect]]) > 0 ):
                exposures[connected[iConnect]] = 1
    
    #print("1")
    #print(exposures)
    #print(sum(exposures))     
    [exposures, accr_org[iSample], kl_div_org[iSample], frob_rel_div_org[iSample], norm_one_dif_org[iSample]] = eval_single_sample(exposures, allSignatures, genome, saOptions)
    #print(sum(exposures))
    #print(exposures)
    exposuresNew[:, iSample] = exposures

    #Dispay summary output data
    out_str = ['Sample #' + str(iSample+1) + ': ' + str(inputSamples['cancerType'][iSample][0]) + ' ' + str(inputSamples['sampleNames'][iSample][0]) + ' with ' + str(sum(exposuresNew[:, iSample]>0)) + ' signatures and an accuracy of ' + str(round(accr_org[iSample][0],2))] #BEWARE ,'%.2f' removed
    return(out_str)               
