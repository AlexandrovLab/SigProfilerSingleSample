# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 13:33:06 2018

@author: compactmatter
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from sigproSS import signatures_optimizer as sigopt
import os
import shutil
import multiprocessing as mp
import xml.etree.ElementTree as ET  
from lxml import etree
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from sigproextractor import subroutines as sub
import sigproSS as cosmic
import sigProfilerPlotting as plot


def get_connected_alwaysinclude_signatures(signatureRulesXML, signatures):
    tree = ET.parse(signatureRulesXML)      
    root = tree.getroot()
    
    #Connected
    csList = [] 
    csIDs = []      
    setList = []
    csLength = []
    connectedSignatures = []
    connectedSignaturesIDs = []

    for conn in root.find('connectedSignatures').findall('connected'):
        csLength.append(len(conn))
        setList.append(conn.get("Set"))
        for i in conn.findall('signature'):
            csList.append(i.text)

    csLengthCumSum = np.cumsum(csLength)
    for j in range(len(csList)):
        csIDs.append(signatures.index(csList[j]))
    for i in range(len(csLength)):
        if i == 0:
            connectedSignatures.append(csList[0:csLengthCumSum[i]])
            connectedSignaturesIDs.append(csIDs[0:csLengthCumSum[i]])
        else:
            connectedSignatures.append(csList[csLengthCumSum[i-1]:csLengthCumSum[i]])
            connectedSignaturesIDs.append(csIDs[csLengthCumSum[i-1]:csLengthCumSum[i]])
    
    #Always Include
    allowEverywhere = []
    for ai in root.find('alwaysIncludeSignatures').findall('signature'):
        allowEverywhere.append(ai.text)
    
    allowEverywhereIDs = [i for i, x in enumerate(signatures) if x in allowEverywhere]
    
    return[connectedSignaturesIDs, allowEverywhereIDs]
    


def check_signature_rules(signaturesInSample,
                          signatures,
                          sampleID,
                          seqType,
                          totalMutations,
                          C_to_A_p,
                          C_to_A_d,
                          C_to_G_p,
                          C_to_G_d,
                          C_to_T_p,
                          C_to_T_d,
                          T_to_A_p,
                          T_to_A_d,
                          T_to_C_p,
                          T_to_C_d,
                          T_to_C_ATN_p,
                          T_to_C_ATN_d,
                          T_to_G_p,
                          T_to_G_d,
                          signatureRulesXML,
                          signatureRulesXMLSchema):
    
    totalSignatures = len(signaturesInSample)  
    
    
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
                    if( eval(strandbias[sb].get('mutType').replace('>','_to_') + "_p") > float(strandbias[sb][1].text) or
                       eval(strandbias[sb].get('mutType').replace('>','_to_') + "_d") != float(strandbias[sb][0].text)):
                       signaturesInSample[i] = 0
                       
            if (root.find(".//totalMutations/..[@signatureSBS='" + signatures[i] + "']") is not None):
                totalmutations = root.find(".//totalMutations/..[@signatureSBS='" + signatures[i] + "']").find('totalMutations').findall('seqType')
                for st in range(len(totalmutations)):
                    if( totalmutations[st].get('type') == seqType):
                        if( float(totalMutations) < float(totalmutations[st][0].text) ):
                            signaturesInSample[i] = 0
                
    return signaturesInSample



def parallel_for_loop(iSample, inputSamples, allGenomeSignatures, allExomeSignatures, cancerType,
                      seqType, totalMutations, sampleNames, samplePvalue,
                      longSampleNames, totalSignatures, useRules, sigNames, idsToAdd, connected,
                      exposuresNew, accr_org, kl_div_org, frob_rel_div_org, norm_one_dif_org):
    
    
    
    #get the path for files
    paths = cosmic.__path__[0]
    
    genome = inputSamples[:,iSample]
    genome = genome.reshape(96,1) #Uncomment for other functions
    if (seqType[iSample] == 'WGS'):
        allSignatures =  allGenomeSignatures
    else:
        allSignatures =  allExomeSignatures
        
    #Identify signatures in cancer type and apply signature rules
    
    
                  

    longSampleName = cancerType[iSample] + '::' + sampleNames[iSample]

    signaturesInSamples = np.zeros([allSignatures.shape[1], inputSamples.shape[1]])
        
    exposures = signaturesInSamples[:, longSampleNames.index(longSampleName)]
    #print(exposures.shape)
    if ( exposures.size == 0 ):
        exposures = np.zeros((totalSignatures, 1))
    
      
    
    C_to_A_p = samplePvalue.iloc[0][0]
    C_to_A_d = samplePvalue.iloc[0][1]
    C_to_G_p = samplePvalue.iloc[1][0]
    C_to_G_d = samplePvalue.iloc[1][1]
    C_to_T_p = samplePvalue.iloc[2][0]
    C_to_T_d = samplePvalue.iloc[2][1]
    T_to_A_p = samplePvalue.iloc[3][0]
    T_to_A_d = samplePvalue.iloc[3][1]
    T_to_C_p = samplePvalue.iloc[4][0]
    T_to_C_d = samplePvalue.iloc[4][1]
    T_to_C_ATN_p = samplePvalue.iloc[5][0]
    T_to_C_ATN_d = samplePvalue.iloc[5][1]
    T_to_G_p = samplePvalue.iloc[6][0]
    T_to_G_d = samplePvalue.iloc[6][1]
    
    if ( useRules == True ):
        #print("Using rules")
        exposures = check_signature_rules(exposures,
                                          sigNames,
                                          iSample,
                                          seqType,
                                          totalMutations,
                                          C_to_A_p,
                                          C_to_A_d,
                                          C_to_G_p,
                                          C_to_G_d,
                                          C_to_T_p,
                                          C_to_T_d,
                                          T_to_A_p,
                                          T_to_A_d,
                                          T_to_C_p,
                                          T_to_C_d,
                                          T_to_C_ATN_p,
                                          T_to_C_ATN_d,
                                          T_to_G_p,
                                          T_to_G_d,
                                          paths+'/input/20181108_Signature_Rules.xml',
                                          paths+'/input/Signature_Rules_Schema.xsd')
    
    
    #print("Total Mutations", totalMutations[0])
    
    #Remove all signatures in the sample: one by one
    #Add signatures that are allowed everywhere    
    if (idsToAdd != False):
        exposures[idsToAdd] = 1 
    
    #Add connected signatures
    if (len(connected) > 0):
        for iConnect in range(len(connected)):
            if ( sum(exposures[connected[iConnect]]) > 0 ):
                exposures[connected[iConnect]] = 1
    
    #print(exposures)
    ##########################################################[exposures, accr, kl_div, frob_rel_div, norm_one_dif] = remove_all_single_signatures(exposures, allSignatures, genome, sigNames)
    
    
    IDs = [i for i,x in enumerate(exposures) if x > 0]
    #!print("Signatures Selected to Remove after checking the rules:")
    #!print(IDs)
    processes = allSignatures[:,IDs]

    exposure =   sigopt.initial_optimization(processes, genome) #processes, sample 
    #exposure, similarity =   sigopt.fit_signatures(processes, genome) #processes, sample
    #exposure = exposure[:,np.newaxis]

    pool = mp.Pool(processes=1)
    results = [pool.apply_async(sigopt.remove_signatures, args=(x,processes,exposure,genome,)) for x in range(genome.shape[1])]
    
    output = [p.get() for p in results]
    pool.close()
    
    exposureAfterRemovingSignatures = np.zeros(exposure.shape)
    for i in range(len(output)):
        exposureAfterRemovingSignatures[:,i]=output[i] 

    np.put(exposures, IDs, exposureAfterRemovingSignatures.ravel())
    #print(np.nonzero(exposures)[0])
    
    #Add all remaining signatures to the sample: one by one
    #Add signatures that are allowed everywhere
    #print(exposures)
    
    if (idsToAdd != False):
        for i in idsToAdd:
            exposures[i] = 1
       
    #Add connected signatures
    if (len(connected) > 0):
        for iConnect in range(len(connected)):
            if ( sum(exposures[connected[iConnect]]) > 0 ):
                exposures[connected[iConnect]] = 1

    #print(exposures)
    IDs = [i for i,x in enumerate(exposures) if x > 0] #Changed i to i+1 
    #print(IDs)
    
    #print(exposureAfterRemovingSignatures) ###############
    
    
    ######################################################exposures = add_all_single_signatures(exposures, allSignatures, genome, sigNames) #No eval single sample now
    ## ## ## ##processes = allSignatures[:,IDs]
    
    ####################print(IDs) 
    #pd.DataFrame(processes).to_csv("processes_file.txt", sep = "\t")
    #pd.DataFrame(genome).to_csv("genome_file.txt", sep = "\t")
    #print("done")
    
    exposureAfterAddingSignatures = sigopt.add_signatures(allSignatures, genome, 0.025, IDs)[0]
    #print(len(exposureAfterAddingSignatures))
    #Add signatures that are allowed everywhere
    exposures = exposureAfterAddingSignatures
    ## ## ## ## np.put(exposures, IDs, exposureAfterAddingSignatures.ravel())
   
    
    
    #print(exposures)#####################
    
    
    if (idsToAdd != False):
        exposures[idsToAdd] = 1

    #Add connected signatures
    if (len(connected) > 0):
        for iConnect in range(len(connected)):
            if ( sum(exposures[connected[iConnect]]) > 0 ):
                exposures[connected[iConnect]] = 1
    
    #print(exposures)
    IDs = [i for i,x in enumerate(exposures) if x > 0]
   
    
    ################################################################[exposures, accr_org[iSample], kl_div_org[iSample], frob_rel_div_org[iSample], norm_one_dif_org[iSample]] = eval_single_sample(exposures, allSignatures, genome)

    #IDs = [i for i,x in enumerate(exposures) if x > 0] #BEWARE
    #processes = allSignatures[:,IDs]#BEWARE

    #exposure =   initial_optimization(processes, genome)#BEWARE
    #exposures = np.zeros((totalSignatures, 1)) #BEWARE
    #np.put(exposures, IDs, exposure.ravel()) #BEWARE
    #print(exposures)
    exposuresNew[:, iSample] = exposures.ravel() #exposures
    #print(exposuresNew)
    indices = np.nonzero(exposuresNew)[0] 
    #print(exposuresNew[list(indices)])
    
    
    #fit the signatures 
    print("\n")
    sigNames=np.array(sigNames)
    exposures, similarity = sigopt.fit_signatures(allSignatures[:,list(indices)], genome)
    
    
    
    
    
    #Dispay summary output data
    #out_str = ['Sample #' + str(iSample+1) + ': ' + str(inputSamples['cancerType'][iSample][0]) + ' ' + str(inputSamples['sampleNames'][iSample][0]) + ' with ' + str(sum(exposuresNew[:, iSample]>0)) + ' signatures and an accuracy of ' + str(round(accr_org[iSample][0],2))] #BEWARE ,'%.2f' removed
    #out_str = ['Sample #' + str(iSample+1) + ': [' + str(cancerType[iSample]) + '] [' + str(sampleNames[iSample]) + '] with ' + str(sum(exposuresNew[:, iSample]>0)) + ' signatures']
    
    print(indices)
    print(exposures)
    print(sigNames)
    print(allSignatures.shape)
    print(similarity)
    return([indices,exposures,sigNames, allSignatures, similarity])               


def analysis_individual_samples(#signaturesSet,
                                #signaturesInSamples,
                                samples,
                                #outputFolder,
                                #outputFile,
                                useRules,
                                #allowSigsEverywhere,
                                #connected,
                                newSignatures_file,
                                genomeSignatures_file,
                                exomeSignatures_file,
                                sampleNames,
                                sampleCancerTypes,
                                cancerType,
                                seqType,
                                totalMutations,
                                p_value,
                                signatureRulesXML):
    
   
    #Analysis of data
    #Loading samples
    
    inputSamples = samples
    
    
    
    #print(inputSamples)
    #print("\n")
    try:
        totalSamples = inputSamples.shape[1] 
    except:
        inputSamples = inputSamples[:,np.newaxis]
        totalSamples = inputSamples.shape[1] 
        
    accr_org = np.zeros((totalSamples, 1))
    kl_div_org = np.zeros((totalSamples, 1))  
    frob_rel_div_org = np.zeros((totalSamples, 1))
    norm_one_dif_org = np.zeros((totalSamples, 1)) 
    #print(inputSamples)
    #Loading signatures
    newSignatures = open(newSignatures_file,'r').read().split('\n')
    

    
    allGenomeSignatures = np.loadtxt(genomeSignatures_file)  
    allExomeSignatures  = np.loadtxt(exomeSignatures_file)
    sigNames = newSignatures
    totalSignatures = allGenomeSignatures.shape[1]
    exposuresNew = np.zeros((totalSignatures, totalSamples))
    #Loading signatures in samples
    
    
    
    
    
    longSampleNames = ["" for x in range(len(sampleNames))]
    
    for i in range(len(sampleNames)): #MATLAB CODE for i = 1 : size(sigsInCanType.sampleNames, 1)
       longSampleNames[i] = sampleCancerTypes[i] + '::' + sampleNames[i]

    

    #print(Parallel(n_jobs = num_cores)(delayed(parallel_for_loop)(iSample, inputSamples, allGenomeSignatures, allExomeSignatures, cancerType_file, seqType_file, totalMutations_file, sampleNames, signaturesInSamples_file, longSampleNames, totalSignatures, useRules, sigNames, idsToAdd, connected, saOptions, exposuresNew, accr_org, kl_div_org, frob_rel_div_org, norm_one_dif_org) for iSample in range(totalSamples)))
    #print(sigNames)
    [connected, allowSigsEverywhere] = get_connected_alwaysinclude_signatures(signatureRulesXML, sigNames)
    #print(connected)
    #print(allowSigsEverywhere)
    
    if(dir().count('allowSigsEverywhere') > 0 and len(allowSigsEverywhere) > 0):
        ###!print('Allow these signatures in all samples:')
        for i in range(len(allowSigsEverywhere)):
            ###!print(newSignatures[allowSigsEverywhere[i]])
            idsToAdd = allowSigsEverywhere
    else:
        idsToAdd = False 
    
    if (dir().count('connected') > 0 and len(connected) > 0):
        ###!print('These mutational signatures are connected:')
        for i in range(len(connected)):
            dispStr = 'Set ' + str(i+1) + ':' + newSignatures[connected[i][0]]
            for j in range(1, len(connected[i])):
                dispStr = dispStr + '; ' + newSignatures[connected[i][j]]
            ###!print(dispStr)

    
    for iSample in range(totalSamples):
        samplePvalue = p_value.iloc[:,iSample]
        results = parallel_for_loop(iSample, inputSamples, allGenomeSignatures, allExomeSignatures, cancerType,
                                seqType, totalMutations, sampleNames, samplePvalue,
                                longSampleNames, totalSignatures, useRules, sigNames, idsToAdd, connected
                                , exposuresNew, accr_org, kl_div_org, frob_rel_div_org, norm_one_dif_org)
        
    
        return results


def decomposition_profile(exposures, similarity, signatureList , sampleName):

  
    
    signames = signatureList 
    
    
    exposure_percentages = exposures[np.nonzero(exposures)]/np.sum(exposures[np.nonzero(exposures)])*100
    listofinformation = list("0"*len(np.nonzero(exposures)[0])*3)
    
    count =0
    for j in np.nonzero(exposures)[0]:
        listofinformation[count*3] = signames[j]
        listofinformation[count*3+1] = round(exposure_percentages[count],2)
        listofinformation[count*3+2]="%"
        count+=1
    ListToTumple = tuple([sampleName] +listofinformation+ [similarity])
    strings ="Sample %s,"+" %s (%0.2f%s) &"*(len(np.nonzero(exposures)[0])-1)+" %s (%0.2f%s), %0.2f\n" 
    return(strings%(ListToTumple))
       



def single_sample(vcf, outputdir, ref="GRCh37", exome=False):
    
    
    """
    Decompose the query samples into the global signatures.
    
    parameters
    ----------
    vcf: A string. The name of the folder containing the vcf files. The folder should be present in the current working directory. 
    outputdir: A string. The name of the output folder. The output folder will be generated in the current working directory according to name provided in the current working directory. 
    ref: A string. The name of the reference genome file. The file should be installed previously through "SigProfilerMatrixGenerator". 
    Please see the "INSTALLATION" part of the README.md file. The default reference genome is "GRCh37".
    exome: A boolean. If the agrument is True, that will genearate the mutational profile only for the exomes. If False, the profile 
    for the whole genome sequence will be generated. 
          
    
    Returns:
    -------
    After the single_sample function is successfully executed, an output directory will be generated in the current working directory. 
    
    The output folder will contain the following files:
        -exposure.txt 
        -signature.txt 
        -probabilities.txt 
        -signature plot pdf 
        -dendrogram plot
        -decomposition profile.csv

    Example: 
    -------
    >>> from sigproSS import spss 
    >>> data = spss.importdata()
    >>> spss.single_sample(data, "results", ref="GRCh37", exome=False)
    
    
        
    """
    
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
        
        
    
    #get the path for files
    paths = cosmic.__path__[0]
    
    if vcf[-1] != "/":
        vcf_name = vcf.split("/")[-1]
    else:
        vcf_name = vcf.split("/")[-2]
    

    data = matGen.SigProfilerMatrixGeneratorFunc(vcf_name, ref, vcf, exome=exome)
    
    # make the totalExposure dataframe which have dimention of totalsignatures and totalsamples
    
    totalExposures = np.zeros([65,data["96"].shape[1]])
    listOfSamples = list(data["96"].head(0))
    
    # open a file to profile the signatures
    fh = open(outputdir+"/decomposition profile.csv", "w")
    fh.write("Sample Names, Global NMF Signatures, Similarity\n")
    fh.close()
    
    for i in range(data["96"].shape[1]):
        print("##########################################################")
        print("Exacting Profile for "+"Sample " +str(i+1))
        index = i
        samples = data["96"].iloc[:,index:index+1]
        p_value = data["7_pvalue"]
        #print(p_value)
        samples = np.array(samples)
        sampleNames = list(data["96"].head(0))[index:index+1]
        cancerType = ['Breast Cancer']*samples.shape[1]
        seqType = ['WGS']*samples.shape[1]
        totalMutations= np.sum(samples, axis=0)
        
      
        
        
        #results variable contains [indices,exposures, signatureNames, allSignatures, similarity]
        results = analysis_individual_samples(samples, 
                                    True,
                                    paths+'/input/signaturesSet.txt',
                                    paths+'/input/genomeSignatures.txt',
                                    paths+'/input/exomeSignatures.txt',
                                    sampleNames,
                                    cancerType,
                                    cancerType,
                                    seqType,
                                    totalMutations,
                                    p_value,
                                    paths+'/input/20181108_Signature_Rules.xml') 
        totalExposures[results[0],i]=results[1]
        listOfSignatures = results[2]
        signatures = pd.DataFrame(results[3])
        profile = decomposition_profile(totalExposures[:,i],  results[4], results[2], sampleNames[0])
        
        #write the profiles into file
        fh = open(outputdir+"/decomposition profile.csv", "a")
        fh.write(profile)
        fh.close()
        
        
    
    
    
    
    #prepare the exposures dataframe
    
    totalExposures = pd.DataFrame(totalExposures)
    totalExposures = totalExposures.set_index(listOfSignatures)
    totalExposures.columns = listOfSamples
    totalExposures = totalExposures.rename_axis("Samples", axis="columns")
    #Convert the floats to integers
    totalExposures[listOfSamples] = totalExposures[listOfSamples].applymap(np.int64)
    
    
    
    #remove the rows with all zeros to create the final exposure dataframe
    exposures = totalExposures.loc[~(totalExposures==0).all(axis=1)]
    
    #presure the signatures dataframe
    signatures = pd.DataFrame(results[3])
    signatures.columns = listOfSignatures
    signatures = signatures.set_index(data["96"].index)
    signatures = signatures.rename_axis("Signatures", axis="columns")
    
    #Filter the signatures by the exposures rows to get the final signature dataframe
    signatures = signatures.loc[:,list(exposures.index)] 
    
    #create the probalities 
    probability = sub.probabilities(signatures, exposures)
    probability = probability.set_index("Sample")
    probability = probability.rename_axis("", axis="columns")
    
    try:
        #create the dedrogrames
        Y, dn = sub.dendrogram(exposures, 0.05, outputdir)
    except:
        pass
        
    #export results
    signatures.to_csv(outputdir+"/signatures.txt", "\t", index_label=[signatures.columns.name]) 
    exposures.to_csv(outputdir+"/exposures.txt", "\t", index_label=[exposures.columns.name]) 
    probability.to_csv(outputdir+"/probabilities.txt", "\t")
    plot.plotSBS(outputdir+"/signatures.txt", outputdir+"/Signature_plot", "", "96", True)

def importdata(inpute_type="vcf"):
    
    """
    Imports the path of example data.
    
    parameters
    ----------
            
    
    Returns:
    -------

    The path of the example data.

    Example: 
    -------
    >>> from sigproSS import spss 
    >>> data = spss.importdata()
    
    This "data" variable can be used as a parameter of the first argument of the sigpross function.
        
    """
    
    
    
    paths = cosmic.__path__[0]
    directory = os.getcwd()
    if inpute_type=="pcwag96":
        data=paths+"/input/csv_example96.csv"
    elif inpute_type=="pcwag192":
        data=paths+"/input/csv_example192.csv"
    
    elif inpute_type=="pcwag192":
        dataold = paths+"/input/vcf"
        datanew = directory+"/vcf"
        if not os.path.exists(datanew):
            shutil.copytree(dataold , datanew) 
        data="vcf"
    return data


