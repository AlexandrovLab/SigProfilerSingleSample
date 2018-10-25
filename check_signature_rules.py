# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 18:36:41 2018

@author: compactmatter
"""

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
                          T_to_G_d_file):
    
    totalSignatures = len(signaturesInSample)
    strand_bias_cutoff = 10**-2
    
    if (seqType[sampleID] == 'WGS'):
        pole_subs_cutoff = 10**5
        msi_subs_cutoff = 10**4
    elif (seqType[sampleID] == 'WES'):
        pole_subs_cutoff = 2 * 10**3
        msi_subs_cutoff = 2 * 10**2
    else:
        print('Invalid type of sequencing data!')    
    
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
    
    for i in range(totalSignatures):

        #Transcriptional Strand Bias
        #[C>A transcribed] Check C>A transcriptional strand bias for signatures 4, 8, 24, 35, 42, and 51
        if ( 'Signature Subs-04' == signatures[i] or
             'Signature Subs-08' == signatures[i] or
             'Signature Subs-24' == signatures[i] or
             'Signature Subs-35' == signatures[i] or
             'Signature Subs-42' == signatures[i] or
             'Signature Subs-51' == signatures[i]):

           if ( C_to_A_p[sampleID] > strand_bias_cutoff or C_to_A_d[sampleID] != -1 ): # damage on G
               signaturesInSample[i] = 0
        
        #[C>A untranscribed] Check C>A transcriptional strand bias for signatures 50
        if ( 'Signature Subs-50' == signatures[i] ):

           if ( C_to_A_p[sampleID] > strand_bias_cutoff or C_to_A_d[sampleID] != 1 ): # damage on C
               signaturesInSample[i] = 0
                
        #[C>T transcribed] Check C>T transcriptional strand bias for signatures 19, 23, 31, 32, 42, and 51
        if ( 'Signature Subs-19' == signatures[i] or
             'Signature Subs-23' == signatures[i] or
             'Signature Subs-31' == signatures[i] or
             'Signature Subs-32' == signatures[i] or
             'Signature Subs-42' == signatures[i] ):
         
            if ( C_to_T_p[sampleID] > strand_bias_cutoff or C_to_T_d[sampleID] != -1): # damage on G
               signaturesInSample[i] = 0
        
        #[C>T untranscribed] Check C>T transcriptional strand bias for signatures 7a and 7b
        if ( 'Signature Subs-07a' == signatures[i] or
             'Signature Subs-07b' == signatures[i] or
             'Signature Subs-51' == signatures[i] ):
         
            if (C_to_T_p[sampleID] > strand_bias_cutoff or C_to_T_d[sampleID] != 1): # damage on C
               signaturesInSample[i] = 0
        
        #[T>A transcribed] Check T>A transcriptional strand bias for signature 22, 25, and 27
        if ( 'Signature Subs-22' == signatures[i] or
             'Signature Subs-25' == signatures[i] or
             'Signature Subs-27' == signatures[i] or
             'Signature Subs-51' == signatures[i] ):
         
            if ( T_to_A_p[sampleID] > strand_bias_cutoff or T_to_A_d[sampleID] != -1): # damage on A
               signaturesInSample[i] = 0
        
        #[T>A untranscribed] Check T>A transcriptional strand bias for signature 7c
        if ( 'Signature Subs-07c' == signatures[i] or 'Signature Subs-47' == signatures[i]):
            
            if ( T_to_A_p[sampleID] > strand_bias_cutoff or T_to_A_d[sampleID] != 1): # damage on T
               signaturesInSample[i] = 0
        
        #[T>C transcribed] Check T>C transcriptional strand bias for signatures 5, 12, and 16
        if ( 'Signature Subs-05' == signatures[i] or 
             'Signature Subs-12' == signatures[i] or 
             'Signature Subs-16' == signatures[i] or 
             'Signature Subs-46' == signatures[i] ):
         
            if ( T_to_C_p[sampleID] > strand_bias_cutoff or T_to_C_d[sampleID] != -1 ): # damage on A
               signaturesInSample[i] = 0
        
        #[T>C untranscribed] Check T>C transcriptional strand bias for signatures 7c, 7d, 21, 26, 33, and 57
        if ( 'Signature Subs-07c' == signatures[i] or
             'Signature Subs-07d' == signatures[i] or 
             'Signature Subs-21' == signatures[i] or
             'Signature Subs-26' == signatures[i] or 
             'Signature Subs-33' == signatures[i] or 
             'Signature Subs-57' == signatures[i]  ):
         
            if ( T_to_C_p[sampleID] > strand_bias_cutoff or T_to_C_d[sampleID] != 1): # damage on T
               signaturesInSample[i] = 0
        
        #[T>G transcribed] Check T>G transcriptional strand bias for signatures 5, 12, and 16
        if ( 'Signature Subs-47' == signatures[i] ): 
         
            if ( T_to_G_p[sampleID] > strand_bias_cutoff or T_to_G_d[sampleID] != -1): # damage on A
               signaturesInSample[i] = 0
        
        #[T>G untranscribed] Check T>G transcriptional strand bias for signatures 7d, 21, 26, and 33
        if ( 'Signature Subs-44' == signatures[i] or 'Signature Subs-57' == signatures[i]): 
         
            if ( T_to_G_p[sampleID] > strand_bias_cutoff or T_to_G_d[sampleID] != 1): # damage on T
               signaturesInSample[i] = 0     
                
        # Check large numbers of mutations for POLE signatures
        if ('Signature Subs-10a' == signatures[i] or 'Signature Subs-10b' == signatures[i]):
            
           #Checking numbers of single base mutations
           if ( float(totalMutations[sampleID]) < pole_subs_cutoff ):
               signaturesInSample[i] = 0
        
        # Check short mutations for MSI signatures
        if ( 'Signature Subs-06' == signatures[i] or 
             'Signature Subs-14' == signatures[i] or 
             'Signature Subs-15' == signatures[i] or 
             'Signature Subs-20' == signatures[i] or 
             'Signature Subs-21' == signatures[i] or 
             'Signature Subs-26' == signatures[i] or 
             'Signature Subs-60' == signatures[i] ):

            #Checking numbers of subs
            if ( float(totalMutations[sampleID]) < msi_subs_cutoff ):
                signaturesInSample[i] = 0
                
    return signaturesInSample