# This source code file is a part of SigProfilerSingleSample

# SigProfilerSingleSample is a tool included as part of the SigProfiler

# computational framework for comprehensive analysis of mutational

# signatures from next-generation sequencing of cancer genomes.

# SigProfilerSingleSample 

# Copyright (C) 2018 [Gudur Ashrith Reddy]

#

# SigProfilerSingleSample is free software: you can redistribute it and/or modify

# it under the terms of the GNU General Public License as published by

# the Free Software Foundation, either version 3 of the License, or

# (at your option) any later version.

#

# SigProfilerSingleSample is distributed in the hope that it will be useful,

# but WITHOUT ANY WARRANTY; without even the implied warranty of

# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

# GNU General Public License for more details.

#

# You should have received a copy of the GNU General Public License

# along with SigProfilerSingleSample.  If not, see http://www.gnu.org/licenses/

import xml.etree.ElementTree as ET  
from lxml import etree

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
                strandbias = root.find(".//strandBias/..[@signatureSBS='Signature Subs-04']").find('strandBias').findall('mutationType')
                for sb in range(len(strandbias)):
                    if( eval(strandbias[sb].get('mutType').replace('>','_to_') + "_p[sampleID]") > float(strandbias[sb][1].text) or
                       eval(strandbias[sb].get('mutType').replace('>','_to_') + "_d[sampleID]") != float(strandbias[sb][0].text)):
                        signaturesInSample[i] = 0
            if (root.find(".//totalMutations/..[@signatureSBS='" + signatures[i] + "']") is not None):
                totalmutations = root.find(".//totalMutations/..[@sigName='" + signatures[i] + "']").find('totalMutations').findall('seqType')
                for st in range(len(totalmutations)):
                    if( totalmutations[st].get('type') == seqType[sampleID]):
                        if( float(totalMutations[sampleID]) < float(totalmutations[st][0].text) ):
                            signaturesInSample[i] = 0
                
    return signaturesInSample
