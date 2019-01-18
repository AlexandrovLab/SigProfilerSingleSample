# SigProfilerSingleSample
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
