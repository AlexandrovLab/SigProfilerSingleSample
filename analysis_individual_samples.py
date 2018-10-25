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
    #Analysis of data
    #Loading samples
    
    from check_signature_rules import check_signature_rules
    from remove_all_single_signatures import remove_all_single_signatures
    from eval_single_sample import eval_single_sample
    from parameterized_objective2_custom import parameterized_objective2_custom
    from add_all_single_signatures import add_all_single_signatures
    
    input = sio.loadmat(samplesForAnalysis)
    totalSamples = input['originalGenomes'].shape[1] #MATLAB CODE totalSamples = size(input.originalGenomes, 2)
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
    for iSample in range(totalSamples):
        #Select cancer sample and set of signatures
        
        genome = input['originalGenomes'][:,iSample]
        if (input['seqType'][iSample] == 'WGS'):
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
        
        strandBias = input['strandBias']
        
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
        print(['Sample #' + str(iSample) + ': ' + str(input['cancerType'][iSample]) + ' ' + str(input['sampleNames'][iSample]) + ' with ' + str(sum(exposuresNew[:, iSample]>0)) + ' signatures and an accuracy of ' + str(accr_org[iSample])]) #BEWARE ,'%.2f' removed
        #COME BACK LATER TO CREATE OUTPUTFORDER AUTOMATICALLY
        #if ( exist(outputFolder,'dir') == 0 ):
        #   mkdir(outputFolder)
        #save([outputFolder filesep outputFile])






