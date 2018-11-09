# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 18:36:41 2018

@author: compactmatter
"""
import xml.etree.ElementTree as ET  

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
                          signatureRulesXML):
    
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

    tree = ET.parse(signatureRulesXML)  
    root = tree.getroot()

    signaturesList = []    
    for i in range(totalSignatures):
        for signature in root.findall('signatureName'):
            signaturesList.append(signature.get("sigName"))
            
        if(signatures[i] in signaturesList):
            if (root.find(".//strandBias/..[@sigName='" + signatures[i] + "']") is not None):
                strandbias = root.find(".//strandBias/..[@sigName='Signature Subs-04']").find('strandBias').findall('mutationType')
                for sb in range(len(strandbias)):
                    if( eval(strandbias[sb].get('type').replace('>','_to_') + "_p[sampleID]") > float(strandbias[sb][1].text) or
                       eval(strandbias[sb].get('type').replace('>','_to_') + "_d[sampleID]") != float(strandbias[sb][0].text)):
                        signaturesInSample[i] = 0
            if (root.find(".//totalMutations/..[@sigName='Signature Subs-04']") is not None):
                totalmutations = root.find(".//totalMutations/..[@sigName='" + signatures[i] + "']").find('totalMutations').findall('seqType')
                for st in range(len(totalmutations)):
                    if( totalmutations[st].get('type') == seqType[sampleID]):
                        if( float(totalMutations[sampleID]) < float(totalmutations[st][0].text) ):
                            signaturesInSample[i] = 0
                
    return signaturesInSample