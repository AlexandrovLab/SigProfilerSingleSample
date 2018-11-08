# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 18:17:39 2018

@author: compactmatter
"""

def analysis_individual_samples(signaturesSet,
                                signaturesInSamples,
                                samplesForAnalysis,
                                outputFolder,
                                outputFile,
                                useRules,
                                allowSigsEverywhere,
                                connected,
                                newSignatures_file,
                                genomeSignatures_file,
                                exomeSignatures_file,
                                sampleNames_file,
                                sampleCancerTypes_file,
                                cancerType_file,
                                seqType_file,
                                totalMutations_file,
                                signaturesInSamples_file):
    import numpy as np
    import scipy.io as sio
    from joblib import Parallel, delayed
    import multiprocessing
    #Analysis of data
    #Loading samples
    
    from check_signature_rules import check_signature_rules
    from remove_all_single_signatures import remove_all_single_signatures
    from eval_single_sample import eval_single_sample
    from parameterized_objective2_custom import parameterized_objective2_custom
    from add_all_single_signatures import add_all_single_signatures
    from parallel_for_loop import parallel_for_loop
    
    inputSamples = sio.loadmat(samplesForAnalysis)
    totalSamples = inputSamples['originalGenomes'].shape[1] #MATLAB CODE totalSamples = size(inputSamples.originalGenomes, 2)
    accr_org = np.zeros((totalSamples, 1))
    kl_div_org = np.zeros((totalSamples, 1))  
    frob_rel_div_org = np.zeros((totalSamples, 1))
    norm_one_dif_org = np.zeros((totalSamples, 1)) 

    #Loading signatures
    newSignatures = open(newSignatures_file,'r').read().split('\n')
    
    if(dir().count('allowSigsEverywhere') > 0 and len(allowSigsEverywhere) > 0):
        print('Allow these signatures in all samples:')
        for i in range(len(allowSigsEverywhere)):
            print(newSignatures[allowSigsEverywhere[i] - 1])
            idsToAdd = allowSigsEverywhere
    else:
        idsToAdd = False #BEWARE modification to code
    
    if (dir().count('connected') > 0 and len(connected) > 0):
        print('These mutational signatures are connected:')
        for i in range(len(connected)):
            dispStr = 'Set ' + str(i+1) + ':' + newSignatures[connected[i][0] - 1]
            for j in range(1, len(connected[i])):
                dispStr = dispStr + '; ' + newSignatures[connected[i][j] - 1]
            print(dispStr)
    
    allGenomeSignatures = np.loadtxt(genomeSignatures_file)  
    allExomeSignatures  = np.loadtxt(exomeSignatures_file)
    sigNames = newSignatures
    totalSignatures = allGenomeSignatures.shape[1]
    exposuresNew = np.zeros((totalSignatures, totalSamples))

    #Loading signatures in samples
    #sigsInCanType = load(signaturesInSamples)
    
    sampleNames = open(sampleNames_file,'r').read().split('\n')
    sampleCancerTypes = open(sampleCancerTypes_file,'r').read().split('\n')
    
    longSampleNames = ["" for x in range(len(sampleNames))]
    
    for i in range(len(sampleNames)): #MATLAB CODE for i = 1 : size(sigsInCanType.sampleNames, 1)
       longSampleNames[i] = sampleCancerTypes[i] + '::' + sampleNames[i]

    #Minization options
    
    saOptions = [ 'Display', 'off', 'TolFun', 1e-100,
                          'MaxFunEvals', float("inf"), 'MaxIter', 100000,
                          'Algorith', 'interior-point', 'FinDiffType', 'central',
                          'TolCon', 1e-100, 'TolX', 1e-100 ]

    #parfor iSample = 1 : totalSamples
    
    num_cores = multiprocessing.cpu_count()
    #print(Parallel(n_jobs = num_cores)(delayed(parallel_for_loop)(iSample, inputSamples, allGenomeSignatures, allExomeSignatures, cancerType_file, seqType_file, totalMutations_file, sampleNames, signaturesInSamples_file, longSampleNames, totalSignatures, useRules, sigNames, idsToAdd, connected, saOptions, exposuresNew, accr_org, kl_div_org, frob_rel_div_org, norm_one_dif_org) for iSample in range(totalSamples)))
    
    for iSample in range(totalSamples):
        print(parallel_for_loop(iSample, inputSamples, allGenomeSignatures, allExomeSignatures, cancerType_file,
                                                                seqType_file, totalMutations_file, sampleNames, signaturesInSamples_file,
                                                                longSampleNames, totalSignatures, useRules, sigNames, idsToAdd, connected,
                                                                saOptions, exposuresNew, accr_org, kl_div_org, frob_rel_div_org, norm_one_dif_org))
        #Select cancer sample and set of signatures      
        #COME BACK LATER TO CREATE OUTPUTFORDER AUTOMATICALLY
        #if ( exist(outputFolder,'dir') == 0 ):
        #   mkdir(outputFolder)
        #save([outputFolder filesep outputFile])