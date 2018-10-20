# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 13:06:38 2018

@author: compactmatter
"""
#%%

signaturesSet = 'input/Consensus_subs_mutational_signatures.mat' # set of signatures
signaturesInSamples = 'input/signatures_in_samples_and_cancer_types.mat' # set of signatures in samples
samplesForAnalysis = 'input/Biliary-AdenoCA_example_cancer_samples.mat' # set of individual samples for examination
outputFolder = 'output/' # output folder
outputFile = 'signatures_in_Biliary-AdenoCA_example_cancer_samples.mat' # output file
useRules = 1 # boolean variable indicating whether to use rules or not (1==use rules; 0==do not use rules)
allowSigsEverywhere = [1, 5] # IDs of signatures to be included in all samples regardless of rules or sparsity  
connected = [[2, 17], [7, 8, 9, 10], [13, 14], [21, 22]]