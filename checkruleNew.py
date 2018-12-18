import xml.etree.ElementTree as ET  
from lxml import etree

def check_signature_rules(signaturesInSample,
                          signatures,
                          sampleID,
                          seqType,
                          totalMutations,
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
    print(C_to_A_p)
    C_to_A_d = [float(i) for i in open(C_to_A_d_file,'r').read().split('\n')]
    print(C_to_A_d)
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
                strandbias = root.find(".//strandBias/..[@signatureSBS='" + signatures[i] + "']").find('strandBias').findall('mutationType')
                for sb in range(len(strandbias)):
                    if( eval(strandbias[sb].get('mutType').replace('>','_to_') + "_p[sampleID]") > float(strandbias[sb][1].text) or
                       eval(strandbias[sb].get('mutType').replace('>','_to_') + "_d[sampleID]") != float(strandbias[sb][0].text)):
                       signaturesInSample[i] = 0
                       
            if (root.find(".//totalMutations/..[@signatureSBS='" + signatures[i] + "']") is not None):
                totalmutations = root.find(".//totalMutations/..[@signatureSBS='" + signatures[i] + "']").find('totalMutations').findall('seqType')
                for st in range(len(totalmutations)):
                    if( totalmutations[st].get('type') == seqType[sampleID]):
                        if( float(totalMutations[sampleID]) < float(totalmutations[st][0].text) ):
                            signaturesInSample[i] = 0
                        
              
    return signaturesInSample


################################################################################################## Check Rule New Version ##########################################################################

import xml.etree.ElementTree as ET  
from lxml import etree
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen


context = context
project = project 
refgen = refgen 
exome = exome
bed_file = bed_file    


    
# Run the SigProfiler MatrixGenerator Function
data = matGen.SigProfilerMatrixGeneratorFunc(project, refgen, project, exome=False, indel_extended=False, bed_file=None, chrom_based=False, plot=False, gs=False)
#data = a
pValues = data["7_pvalue"]
samples = data["96"]
samples = np.array(samples)
totalMutations = np.sum(samples, axis=0)
signaturesInSample = signaturesInSample
signatures = signatures
if (exome==False and bed_file==None):
    seqType=['WGS']
elif (exome==True and bed_file==None):
    seqType=['Exom']
elif (exome==False and bed_file!=None):
    seqType==['Region']
    
def check_signature_rules(signaturesInSample,
                          signatures,
                          sampleID,
                          seqType,
                          totalMutations,
                          pValues,
                          signatureRulesXML,
                          signatureRulesXMLSchema):
    
    
    

    totalMutations = totalMutations[sampleID]
    samplePvalues = pValues.iloc[:, sampleID]
    
    
    
    totalSignatures = len(signaturesInSample)  
    
    C_to_A_p = samplePvalues.iloc[0][0]
    C_to_A_d = samplePvalues.iloc[0][1]
    C_to_G_p = samplePvalues.iloc[1][0]
    C_to_G_d = samplePvalues.iloc[1][1]
    C_to_T_p = samplePvalues.iloc[2][0]
    C_to_T_d = samplePvalues.iloc[2][1]
    T_to_A_p = samplePvalues.iloc[3][0]
    T_to_A_d = samplePvalues.iloc[3][1]
    T_to_C_p = samplePvalues.iloc[4][0]
    T_to_C_d = samplePvalues.iloc[4][1]
    T_to_C_ATN_p = samplePvalues.iloc[5][0]
    T_to_C_ATN_d = samplePvalues.iloc[5][1]
    T_to_G_p = samplePvalues.iloc[6][0]
    T_to_G_d = samplePvalues.iloc[6][1]

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
                    if( eval(strandbias[sb].get('mutType').replace('>','_to_') + "_p[sampleID]") > float(strandbias[sb][1].text) or
                       eval(strandbias[sb].get('mutType').replace('>','_to_') + "_d[sampleID]") != float(strandbias[sb][0].text)):
                       signaturesInSample[i] = 0
                       
            if (root.find(".//totalMutations/..[@signatureSBS='" + signatures[i] + "']") is not None):
                totalmutations = root.find(".//totalMutations/..[@signatureSBS='" + signatures[i] + "']").find('totalMutations').findall('seqType')
                for st in range(len(totalmutations)):
                    if( totalmutations[st].get('type') == seqType[sampleID]):
                        if( float(totalMutations[sampleID]) < float(totalmutations[st][0].text) ):
                            signaturesInSample[i] = 0
                        
              
    return signaturesInSample







if __name__ == "__main__":
    
    data = a
    strandBias = data["7_pvalue"]
    samples = data["96"]
    samples = np.array(samples)
    totalMutations = np.sum(samples, axis=0)
    if (exome==False and bed_file==None):
        seqType=['WGS']
    elif (exome==True and bed_file==None):
        seqType=['WES']
    
    
    exposures = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    sigNames = ['Signature Subs-01', 'Signature Subs-02', 'Signature Subs-03', 'Signature Subs-04',
                'Signature Subs-05', 'Signature Subs-06', 'Signature Subs-07a', 'Signature Subs-07b',
                'Signature Subs-07c', 'Signature Subs-07d', 'Signature Subs-08', 'Signature Subs-09',
                'Signature Subs-10a', 'Signature Subs-10b', 'Signature Subs-11', 'Signature Subs-12',
                'Signature Subs-13', 'Signature Subs-14', 'Signature Subs-15', 'Signature Subs-16',
                'Signature Subs-17a', 'Signature Subs-17b', 'Signature Subs-18', 'Signature Subs-19',
                'Signature Subs-20', 'Signature Subs-21', 'Signature Subs-22', 'Signature Subs-23',
                'Signature Subs-24', 'Signature Subs-25', 'Signature Subs-26', 'Signature Subs-27',
                'Signature Subs-28', 'Signature Subs-29', 'Signature Subs-30', 'Signature Subs-31',
                'Signature Subs-32', 'Signature Subs-33', 'Signature Subs-34', 'Signature Subs-35',
                'Signature Subs-36', 'Signature Subs-37', 'Signature Subs-38', 'Signature Subs-39',
                'Signature Subs-40', 'Signature Subs-41', 'Signature Subs-42', 'Signature Subs-43',
                'Signature Subs-44', 'Signature Subs-45', 'Signature Subs-46', 'Signature Subs-47',
                'Signature Subs-48', 'Signature Subs-49', 'Signature Subs-50', 'Signature Subs-51',
                'Signature Subs-52', 'Signature Subs-53', 'Signature Subs-54', 'Signature Subs-55',
                'Signature Subs-56', 'Signature Subs-57', 'Signature Subs-58', 'Signature Subs-59',
                'Signature Subs-60']
    iSample = 0
    
    exposures = check_signature_rules(exposures,
                                          sigNames,
                                          iSample,
                                          seqType,
                                          totalMutations,
                                          strandBias,
                                          'input/20181108_Signature_Rules.xml',
                                          'input/Signature_Rules_Schema.xsd')
    print(exposures)
