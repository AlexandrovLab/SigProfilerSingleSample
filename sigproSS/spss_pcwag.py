#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 15:22:11 2019

@author: mishugeb
"""


import os
import warnings
import pandas as pd
import numpy as np
import copy
import time
from scipy.stats import binom_test
import multiprocessing
from functools import partial
from sigproSS import signatures_optimizer as sigopt
from sigproSS import spss
from sigproextractor import subroutines as sub
import sigproSS as cosmic

mutation_categories = {
    'SBS':['C>A','C>G','C>T','T>A','T>C','T>G'],
    'DBS':['AC>','AT>','CC>','CG>','CT>','GC>','TA>','TC>','TG>','TT>'],
    'ID':['DEL_C_1', 'DEL_T_1', 'INS_C_1', 'INS_T_1', 'DEL_repeats_2', 'DEL_repeats_3', 'DEL_repeats_4', 'DEL_repeats_5+',
    'INS_repeats_2', 'INS_repeats_3', 'INS_repeats_4', 'INS_repeats_5+', 'DEL_MH_2', 'DEL_MH_3', 'DEL_MH_4', 'DEL_MH_5+'],
}

# each signature corresponds to a dictionary with categories and relevant TSB signs:
# -1 means a negative TSB, i.e. more mutations on untranscribed strand,
# 1 means a positive TBS, i.e. more mutations on transcribed strand

def get_indeces(a, b):
    
    """ 
    Extracts the indices multiple items in a list.
    
    Parameters:
        a: list. where we want to get the index of the items.
        b; list. the items we want to get index of. 
    """

    indeces = []
    for i in b:
        try: 
            idx = a.index(i)
            indeces.append(idx)
        except: 
            next

    return indeces 

    """
    #example: 
    x = [1,3,5,8]
    y = [3, 8]
    get_indeces(x, y)
    #result
    >>> [1,3]
    """
    

def get_items_from_index(x,y):
    """ decipher the values of items in a list from their indices.
    """
    z = []
    for i in y:
        z.append(x[i])
    return z

    """
    #example: 
    x = [1,3,5,8]
    y = [1, 3]
    get_items_from_index(x, y)
    #result
    >>> [3,8]
    """
    
#############################################################################################################
##################################### Get Input From CSV Files ##############################################
#############################################################################################################
def read_csv(filename, folder = False):
    
    # if filename is not a string then get the genomes from the object
    if type(filename) != str:
        genomes = filename
        
    else:
        
        if folder==False:
            genomes = pd.read_csv(filename, sep=",").iloc[:, :]    
        else:
        
            count = 1
            for file in os.listdir(filename):
                #print(count) 
                df = pd.read_csv(filename+"/"+file)
                if count==1:
                        genomes = df
                else:
                        genomes = pd.merge(genomes, df, on=["Mutation type", "Trinucleotide"])
                count += 1


    

    
    
    mtypes = [str(genomes.shape[0])]
    
    if mtypes == ['96']:
        # create the list to sort the mutation types
        orderlist = ['A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'A[C>G]A', 'A[C>G]C', 'A[C>G]G', 'A[C>G]T', 'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T',
                     'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T', 'A[T>C]A', 'A[T>C]C', 'A[T>C]G', 'A[T>C]T', 'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T',
                     'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T', 'C[C>G]A', 'C[C>G]C', 'C[C>G]G', 'C[C>G]T', 'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T',
                     'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T', 'C[T>C]A', 'C[T>C]C', 'C[T>C]G', 'C[T>C]T', 'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T',
                     'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T', 'G[C>G]A', 'G[C>G]C', 'G[C>G]G', 'G[C>G]T', 'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T',
                     'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T', 'G[T>C]A', 'G[T>C]C', 'G[T>C]G', 'G[T>C]T', 'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T',
                     'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T', 'T[C>G]A', 'T[C>G]C', 'T[C>G]G', 'T[C>G]T', 'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T',
                     'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T', 'T[T>C]A', 'T[T>C]C', 'T[T>C]G', 'T[T>C]T', 'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T']
        #Contruct the indeces of the matrix
        #setting index and columns names of processAvg and exposureAvg
        index1 = genomes.iloc[:,1]
        index2 = genomes.iloc[:,0]
        index = []
        for i, j in zip(index1, index2):
            index.append(i[0]+"["+j+"]"+i[2])
    
    
        index = np.array(pd.Series(index))
        genomes["index"] = index
    
    
        #sort the mutationa types
        genomes["index"]= pd.Categorical(genomes["index"], orderlist)
        genomes = genomes.sort_values("index")
        
    
        genomes = genomes.iloc[:,2:genomes.shape[1]]
    
    
    
        # set the index 
        genomes = genomes.set_index("index")
    
        # prepare the index and colnames variables
        index = np.array(orderlist)
        colnames = np.array(pd.Series(genomes.columns.tolist()))
    
    else:
        
        index = np.array(genomes.iloc[:,0])
        genomes = genomes.iloc[:,1:]
        genomes = genomes.loc[:, (genomes != 0).any(axis=0)]
        colnames = genomes.columns
       
        
        
   
    
    
    return(genomes, index, colnames, mtypes)


def make_folder_if_not_exists(folder):
    if not os.path.exists(folder):
        try:
            os.makedirs(folder)
        except:
            warnings.warn("Could not create a folder ", folder)

def calculate_strand_bias(list_of_values):
    """
    Perform a test to compare two Poisson means (mutation counts on transcribed
    and untranscribed strands) usng binom_test function from scipy.stats.
    Here, an exact binomial test is used which is preferred for low mutation counts.

    Parameters:
    ----------
    list_of_values: array_like
        list of two values, e.g. [23, 32]

    Returns:
    -------
    The p-value of the hypothesis test
    -------
    """
    if len(list_of_values)!=2:
        raise ValueError("Please provide a list of two values (transcribed and untranscribed strands counts)")

    return binom_test(list_of_values)

def annotate_strand_bias(input_data, p_value=0.05):
    #print("###################  INPUT DATA ########################")
    #print(input_data)
    """
    Perform a test of strand bias on each row of an input dataframe.

    Parameters:
    ----------
    input_data: pandas dataframe
        input dataframe that has exactly two columns (transcribed and untranscribed strand mutation counts)
    p_value: float, optional
        threshold p-value to deem the strand bias significant. Default is 0.05

    Returns:
    -------
    dataframe, strand_bias_flag, affected_categories: tuple

    dataframe: pandas dataframe
    Dataframe with star-marked rows with significant TSB

    strand_bias_flag: boolean
    Boolean flag showing if significant TSB is observed in any row

    affected_categories: dict
        Dictionary of row indexes with significant TSB with absolute differences (transcribed minus untransribed)
    -------
    """
    dataframe = copy.deepcopy(input_data)
    strand_bias_flag = False
    if not isinstance(dataframe, pd.DataFrame):
        raise ValueError("Input not a dataframe, can not annotate for strand bias. Given type:",type(dataframe))
    affected_categories = {}
    for category in list(dataframe.index):
        transcribed_strand_count, untranscribed_strand_count = dataframe.loc[category].tolist()
        strand_bias_p_value = calculate_strand_bias([transcribed_strand_count, untranscribed_strand_count])
        if strand_bias_p_value < p_value:
            # check which strand has more mutations
            difference = transcribed_strand_count-untranscribed_strand_count
            dataframe.rename(index={category: '$\star$'+category}, inplace=True)
            strand_bias_flag = True
            affected_categories[category] = difference
    #print("################### DATA FRAME ########################")
    #print(dataframe)
    #print("################### strand_bias_flag ########################")
    #print(strand_bias_flag)
    #print("################### affected_categories ########################")
    #print(affected_categories)
    return dataframe, strand_bias_flag, affected_categories

def rearrange_stranded_data(input_data):
    
    """
    Rearrange stranded data of a single sample into two different columns of the dataframe

    Parameters:
    ----------
    input_data: pandas dataframe
        Input single-column dataframe with a 3-level MultiIndex of 192-context PCAWG format (e.g. WGS_PCAWG.192.csv on Synapse)
        Essential to have a 'Strand' level of the MultiIndex, with 'T' and 'U' values

    Returns:
    -------
    rearranged_data: pandas dataframe
        Dataframe with removed 'Strand' level, but with two columns (transcribed and untranscribed strands)
    -------
    """
    data = copy.deepcopy(input_data)
    
    T_strand, U_strand = pd.DataFrame(), pd.DataFrame()
    T_strand = data.iloc[data.index.get_level_values('Strand') == 'T']
    U_strand = data.iloc[data.index.get_level_values('Strand') == 'U']
    T_strand = T_strand.droplevel(0)
    U_strand = U_strand.droplevel(0)
    rearranged_data = pd.concat([T_strand, U_strand], names=['transcribed strand','untranscribed strand'], axis=1, sort=True)
    rearranged_data.columns = ['transcribed strand','untranscribed strand']
    return rearranged_data

def categorise_mutations(input_data, categories, condensed=False, strand_bias=False):
    """
    Categorise the dataframe by a given a list of categories

    Parameters:
    ----------
    input_data: pandas dataframe
        Input dataframe. If strand_bias flag is True, a single-column dataframe is expected with
        a 3-level MultiIndex of 192-context PCAWG format (see rearrange_stranded_data function)

    categories: list of strings
        Mutation categories to condense mutations to, e.g.: ['C>A','C>G','C>T','T>A','T>C','T>G']

    condensed: boolean
        If True, the function returns a dataframe with mutations summed up for given categories,
        otherwise, a dictionary of dataframes is returned with keys for each category

    strand_bias: boolean
        If True, rearrange_stranded_data function is called on a dataframe to obtain a two-column format.

    Returns:
    -------
    categorised_mutations: pandas dataframe
        Dataframe with categorised mutations if condensed flag is True, otherwise
        a dictionary of dataframes with keys for each category listed in input categories
    -------
    """
  
    
    data = copy.deepcopy(input_data)
    if strand_bias:
        data = rearrange_stranded_data(data)

    # dataframe in case of condensed format, dictionary of dataframes otherwise
    if condensed:
        if strand_bias:
            categorised_mutations = pd.DataFrame(columns=['transcribed strand','untranscribed strand'])
        else:
            if isinstance(data, pd.Series):
                categorised_mutations = pd.DataFrame(columns=['N'])
            else:
                categorised_mutations = pd.DataFrame(columns=data.columns)
    else:
        categorised_mutations = {}

    for category in categories:
        mutations_per_category = data.filter(like=category, axis=0)

        # sum up categories for condensed format, otherwise leave it in dictionary
        if condensed:
            categorised_mutations.loc[category] = mutations_per_category.sum()
        else:
            categorised_mutations[category] = mutations_per_category
    
    
    return categorised_mutations

def apply_TSB_rules(input_mutations, input_signatures, mutation_type, verbose=False):
    
    #print("###################  input_mutations ########################")
    #print(input_mutations)
    #print("###################  input_signatures ########################")
    #print(input_signatures)
    #print("###################  mutation_type ########################")
    #print(mutation_type)
    """
    Apply transcriptional strand bias rules for a single-column dataframe (single sample)
    If a strand bias is observed in a sample but not a signature, it is excluded from the input list
    A 'sign' of TSB is respected (has to be the same in both samples and signatures)

    Parameters:
    ----------
    input_mutations: pandas dataframe
        Input dataframe. A single-column dataframe is expected with a 3-level MultiIndex
        of 192-context PCAWG format (see rearrange_stranded_data function)

    input_signatures: list of strings
        List of input signatures to apply the rules for

    mutation_type: string
        Mutation type (e.g. SBS, DBS, ID), used to pull the list of relevant categories

    verbose: boolean
        Verbosity flag: lots of output for debugging if set to True

    Returns:
    -------
    signatures_after_TSB_rules: list of strings
        The list of signatures passing the strand bias rules
    -------
    """
    signatures_with_TSB = {
    'SBS4':{'C>A':1,'T>A':1},
    'SBS5':{'T>C':1},
    'SBS7a':{'C>T':-1},
    'SBS7b':{'C>T':-1},
    'SBS7c':{'T>A':-1},
    'SBS7d':{'T>C':-1},
    'SBS8':{'C>A':1},
    'SBS16':{'T>C':1},
    'SBS19':{'C>T':1},
    'SBS21':{'T>C':-1},
    'SBS22':{'T>A':1},
    'SBS23':{'C>T':1},
    'SBS24':{'C>A':1},
    'SBS25':{'T>A':1},
    'SBS26':{'T>C':-1},
    'SBS27':{'T>A':1},
    'SBS29':{'C>A':1},
    'SBS31':{'C>T':1},
    'SBS32':{'C>T':1},
    'SBS33':{'T>C':-1},
    'SBS35':{'C>A':1,'C>T':1},
    'SBS42':{'C>A':1,'C>T':1},
    'SBS45':{'C>A':1,'T>C':1},
    'SBS47':{'T>A':-1,'T>G':1},
    'SBS50':{'C>A':-1},
    'SBS51':{'C>A':1,'C>T':-1,'T>A':1},
    'SBS57':{'T>C':-1,'T>G':-1},
    'SBS60':{'T>G':-1}
     }
    
    categorised_mutations = categorise_mutations(input_mutations, mutation_categories[mutation_type],
                            condensed=True, strand_bias=True)

    _, strand_bias_present, categories_with_TSB = annotate_strand_bias(categorised_mutations)

    if verbose:
        print('TSB found: ', strand_bias_present, categories_with_TSB)
        print('Initial signatures:', input_signatures.columns.tolist())

    signatures_after_TSB_rules = input_signatures
    for signature in input_signatures.columns.tolist():
        if signature in signatures_with_TSB.keys():
            if set(signatures_with_TSB[signature].keys())<=set(categories_with_TSB.keys()):
                # all necessary categories with significant strand bias are present for a given signature
                # now need to check the TSB sign equivalence
                for category in signatures_with_TSB[signature].keys():
                    if signatures_with_TSB[signature][category] != np.sign(categories_with_TSB[category]):
                        if verbose:
                            print('Excluding signature',signature,'as for category',category,'the sign',signatures_with_TSB[signature][category],'differs from',np.sign(categories_with_TSB[category]))
                        signatures_after_TSB_rules.drop(signature, axis=1, inplace=True)
                        # break out of the loop as the signature only needs to be excluded once
                        break
            else:
                signatures_after_TSB_rules.drop(signature, axis=1, inplace=True)
                if verbose:
                    print('Excluding signature',signature,'as',signatures_with_TSB[signature],' is not equal or a subset of ',categories_with_TSB)

    if verbose:
        print('Signatures passing TSB rules:',signatures_after_TSB_rules.columns.tolist())
    
    #print(signatures_after_TSB_rules)
    #Just return the column names to the remaining signatures 
    return signatures_after_TSB_rules.columns



def assign_signatures96(s, genome=0, genome192=0, sig_pcwag=0, par_one=0.01, par_two=0.025):
        
        """
        samples = a single sample
        genome = all samples in 96 context
        genome192 = all samples in 192 context
        sig_pcwaf = the global signatures in 96 context
        
        """
        
        #get the sample's name
        samples = genome.columns[s]
        # get the Transcription bias
        initial_signatures = sig_pcwag.copy()
        #print updates
        print("\n#################################################\n")
        print("Processing sample "+samples+"\n" ) 
        
        if len(genome192) == 0:
            selected_signatures = initial_signatures
        else:
            sample192 = genome192.loc[:,samples]    
            
            # Apply the strand bias rule
            selected_signatures = apply_TSB_rules(sample192, initial_signatures, "SBS", verbose=False)
            
           
        #get the selected signatures indices
        selected_indeces  = get_indeces(list(sig_pcwag.columns), list(selected_signatures))
        sample = np.array(genome.loc[:,samples])[:,np.newaxis]
        #print(len(selected_indeces))
        
        # sigs means the whole signatures
        sigs = sig_pcwag.reset_index()
        sigs, idx, col, _ = sub.read_csv(sigs)
        sigs = sigs.set_index(idx)
        sigs.columns = col
        sigs = np.array(sigs) 
        
        # adding the sigtures in the first round       
        try:    #if the signature database is default
            #get the initial exposure
            preselected = ["SBS1", "SBS5"]
            preselected_idx = get_indeces(list(col), preselected)
            #exposures, similarity = fit_signatures(sigs[:,0:5], sample)
        except:  # if the signature database is custom
            preselected_idx = []
            
        # Add the signatures in the first round. a = exposures, b= similarity
        a, b = sigopt.add_signatures(sigs, sample, cutoff=par_one, presentSignatures=preselected_idx, toBeAdded=selected_indeces)
        present = np.nonzero(a)[0]
        
        #added_first = get_items_from_index(list(col),list(present))
        
        
        
        
        
        
        # adding the sigtures in the second round 
        # Add the signatures in the second round. a = exposures, b= similarity
        nonselected_indeces = list(set(list(range(sigs.shape[1])))-set(selected_indeces))
        a,b = sigopt.add_signatures(sigs, sample, cutoff=par_two, presentSignatures=list(present), toBeAdded= nonselected_indeces)
        present = np.nonzero(a)[0]
        added_second = get_items_from_index(list(col),list(present))
       
       
        
        try:   #if the signature database is default
            #get the connected signatures
            connected_singnatures = [["SBS17a", "SBS17b"], ["SBS7a","SBS7b", "SBS7c","SBS7d"], ["SBS2", "SBS13"], ["SBS10a", "SBS10b"]]
            connected = []     
            for i in connected_singnatures:
                if any(j in added_second for j in i):
                    connected = connected + i
        except:   # if the signature database is custom
            connected = [] 
        
        #get the final signatures
        finalSignatures = list(set(added_second).union(set(connected)))
        finalSinaturesIdx = get_indeces(list(col), finalSignatures)
        finalSinaturesIdx .sort()
        exposures, similarity = sigopt.fit_signatures(sigs[:,finalSinaturesIdx], sample)
        
       
        
        return[finalSinaturesIdx, exposures, list(col), sigs, similarity, genome.columns[s]]
    



def single_sample_pcwag(samples, samples192="", output="results", sigbase="default", par_one=0.01, par_two=0.025, n_cpu=-1):
    
    """ 
    This function takes the csv file of the mutation context SBS96 and SBS192 the same sample/samples. The csv files should be in the PCWAG format.
    The output is the activity of the global signatures of the sample/samples.
    
    Parameters:
        samples:  string or dataframe. if string, this should be a csv file of 96 context PCAWG format. if dataframe, this should be a mutational catalogue where the row 
        index will be the names of mutations and the column names will be the sample names. 
        samples192: string. only valid if the "samples" argument is a 96 contex samples. this should be a csv file in 192 context PCAWG format where the samples should be indentical to the input of the "samples" argument.
        output: string. name of the output folder.
        sigbase: string or dataframe. if string, this should be the name of the csv files in PCAWG format that contains the signature database. The signature database matrix should have the same number
        of rows the "samples" input has. 
        par_one = float. the cut off cosine similarity difference in the first layer of adding signature. defualt value is 0.01.
        par_two = float. the cut off consine similarity difference in the second layer of adding sinature. defualt value is 0.025.
        cpu = integer. the number of cpus to used. default value is -1 which uses all cpus.  
        
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
    >>> from sigproSS import spss, spss_pcwag 
    >>> csv96 = spss.importdata("pcwag96")  #'you can put the path of your own csv96 file here'
    >>> csv192 = spss.importdata("pcwag192") #''you can put the path of your own csv192 file here'
    >>> spss_pcwag.single_sample_pcwag(csv96, csv192, output="example_output")
        
    """
    
    
    # First, get the path for data files
    paths = cosmic.__path__[0]
    
    if not os.path.exists(output):
        os.makedirs(output)
        
           
    # get the genome
    if type(samples)==str:
        genome, idx, col, mtype = sub.read_csv(samples)
        genome = genome.set_index(idx)
        genome.columns = col
    else: 
        genome = samples
        
    if samples192 != "":
        # get the 192 context of the genome     
        genome192 = pd.read_csv(samples192, sep=",", index_col = [0,1,2])
        #Check if the 96 context agrees with the 192 context
        if genome.shape[1] != genome192.shape[1]:
            raise Exception('Please the correct format of the 96 and 192 context files')
        elif genome.shape[0] != 96:
            print("samples192 parameter will not be used for the mutational context of your input sample")
            genome192 = ""
    else:
        genome192 = ""
    
    #get the signature files
    if type(sigbase)==str:
        if sigbase == "default":    
            sig_pcwag = pd.read_csv(paths+"/input/sigProfiler_SBS_signatures.csv", index_col = [0,1])
        else:
            sig_pcwag = pd.read_csv(sigbase, index_col = [0,1])
        
         
    else:
        sig_pcwag = sigbase
        
    #check if the signature matrix is compatible with genome matrix
    if genome.shape[0] != sig_pcwag.shape[0]:
        raise Exception('The size mutation context of samples is not compatible to the signature database')
        
    
    totalExposures = np.zeros([sig_pcwag.shape[1],genome.shape[1]]) 
    listOfSamples = list(col)
    
    # open a file to profile the signatures
    fh = open(output+"/decomposition profile.csv", "w")
    fh.write("Sample Names, Global NMF Signatures, Similarity\n")
    fh.close()
    
    
    
    # Map the signautures using paraller processing 
    if n_cpu==-1:
        pool = multiprocessing.Pool()
    else:
        pool = multiprocessing.Pool(processes=n_cpu)
    
    pool_assign_signatures96 = partial(assign_signatures96, genome=genome, genome192=genome192, sig_pcwag=sig_pcwag, par_one = par_one, par_two = par_two)
    results = pool.map(pool_assign_signatures96, range(len(genome.columns)))
    
    pool.close()
    pool.join()
    
    # Extract the results generated from pool 
    n=0
    similarities = []
    for result in results:
        
        totalExposures[result[0],n]=result[1]
        listOfSignatures = result[2]
        signatures = pd.DataFrame(result[3])
        
        profile = spss.decomposition_profile(totalExposures[:,n],  result[4], result[2], result[5])
        similarities.append(result[4])
        #write the profiles into file
        fh = open(output+"/decomposition profile.csv", "a")
        fh.write(profile)
        fh.close()
        n+=1
    
    
    
    #prepare the exposures dataframe 
    totalExposures = pd.DataFrame(totalExposures)
    totalExposures = totalExposures.set_index(np.array(listOfSignatures))
    totalExposures.columns = listOfSamples
    totalExposures = totalExposures.rename_axis("Cancer Types", axis="columns")
    #Convert the floats to integers
    totalExposures[listOfSamples] = totalExposures[listOfSamples].applymap(np.int64)
    
    
    
    #remove the rows with all zeros to create the final exposure dataframe
    #exposures = totalExposures.loc[~(totalExposures==0).all(axis=1)]
    
    #presure the signatures dataframe
    signatures = pd.DataFrame(result[3])
    signatures.columns = listOfSignatures
    signatures = signatures.set_index(genome.index)
    signatures = signatures.rename_axis("Signatures", axis="columns")
    
    #Filter the signatures by the exposures rows to get the final signature dataframe
    #signatures = signatures.loc[:,list(exposures.index)] 
    
    #create the probalities 
    probability = sub.probabilities(signatures, totalExposures, genome.index, signatures.columns, totalExposures.columns)
    probability = probability.set_index("Sample Names")
    probability = probability.rename_axis("", axis="columns")
    
    
    
    # Transform the exposures matrix and insert the cosine similarity column
    exposures = totalExposures.T    
    exposures.insert(loc=0, column='Similarity', value=similarities )
    
    #export results
    signatures.to_csv(output+"/signatures.txt", "\t", index_label=[signatures.columns.name]) 
    exposures.to_csv(output+"/sig_activities.txt", "\t", index_label=["Cancer Types"]) 
    probability.to_csv(output+"/mutation_probabilities.txt", "\t")
    
    
    try:
        #create the dedrogrames
        Y, dn = sub.dendrogram(totalExposures, 0.05, output)
    except:
        pass
    #plot.plotSBS(outputdir+"/signatures.txt", outputdir+"/Signature_plot", "", "96", True)  
    
    # Finish the process
    print ("\n\nYour Job is Completed! Thank you for using Sigprofiler Single Sample.")
    
        
#single_sample("Biliary-AdenoCA.96.csv", csv192="Biliary-AdenoCA.192.csv", output="results", par_one=0.01, par_two=0.025)        
    
    
    
    
    

            

