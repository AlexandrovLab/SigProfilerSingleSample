# SigProfilerSingleSample
SigProfilerSingleSample allows attributing a known set of mutational signatures to an individual sample. The tool identifies the activity of each signature in the sample and assigns the probability for each signature to cause a specific mutation type in the sample. The tool makes use of SigProfilerMatrixGenerator, SigProfilerExtractor and SigProfilerPlotting. 

## INSTALLATION
In the commandline, please type the following line:
```
$pip install sigproSS
```
Install your desired reference genome from the command line/terminal as follows (available reference genomes are: GRCh37, GRCh38, mm9, and mm10):
```
$ python
>> from SigProfilerMatrixGenerator import install as genInstall
>> genInstall.install('GRCh37')
```
This will install the human 37 assembly as a reference genome. You may install as many genomes as you wish.

open a python interpreter and import the SigProfilerExtractor module. Please see the examples of the functions. 

## FUNCTIONS

### importdata 
    
    
    Imports the path of example data.
    
    importdata()

    Example: 
    -------
    >>> from sigproSS import spss 
    >>> data = spss.importdata()
    This "data" variable can be used as a parameter of the first argument of the sigproSS function.
    
    To get help on the parameters and outputs of the "importdata" function, please write down the following line:
    
    >>> help(spss.importdata)
    
### single_sample

    Decompose the query samples into the global signatures.
    
    single_sample(vcf, outputdir, exome=False)
    
    Example: 
    -------
    >>> from sigproSS import spss 
    >>> data = spss.importdata()
    >>> spss.single_sample(data, "results", ref="GRCh37", exome=False)
    
    To get help on the parameters and outputs of the "importdata" function, please write down the following line:
    
    >>> help(spss.single_sample)
## COPYRIGHT
This software and its documentation are copyright 2018 as a part of the sigProfiler project. The SigProfilerExtractor framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## CONTACT INFORMATION
Please address any queries or bug reports to S M Ashiqul Islam (Mishu) at m0islam.ucsd.edu
