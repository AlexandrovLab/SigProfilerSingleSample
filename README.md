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
>>> from SigProfilerMatrixGenerator import install as genInstall
>>> genInstall.install('GRCh37')
```
This will install the human 37 assembly as a reference genome. You may install as many genomes as you wish.

open a python interpreter and import the SigProfilerExtractor module. Please see the examples of the functions. 

## FUNCTIONS

### importdata 
    
    
    Imports the path of example data.
    
    importdata(inpute_type="vcf")

    Example 1: 
    ----------
    To import an example vcf project provided with the package:
    
    >>> from sigproSS import spss 
    >>> data = spss.importdata("vcf")
    
    This "data" variable can be used as a parameter of the first argument of the single_sample function.
        
    
    Example 2: 
    ----------
    To import an example csv96 file (for more description please the singple_sample_pcwag function) provided with the           package:
        
    >>> from sigproSS import spss 
    >>> data = spss.importdata("pcwag96")
    
    This "data" variable can be used as a parameter of the first argument of the single_sample_pcwag function.
    
    Example 3: 
    ----------
    To import an example csv192 file (for more description please the singple_sample_pcwag function) provided with the           package:
        
    >>> from sigproSS import spss 
    >>> data = spss.importdata("pcwag192")
    
    This "data" variable can be used as a parameter of the second argument of the single_sample_pcwag function.
    
   **To get help on the parameters and outputs of the "importdata" function, please write down the following line:**
    
   _>>> help(spss.importdata)_
    
### single_sample

    Decompose the query samples into the global signatures.
    
    single_sample(vcf, outputdir, exome=False)
    
    Example: 
    -------
    >>> from sigproSS import spss 
    >>> data = spss.importdata() ##'you can put the path of your own vcf project here'
    >>> spss.single_sample(data, "results", ref="GRCh37", exome=False)
    
   **To get help on the parameters and outputs of the "single_sample" function, please write down the following line:**
    
    _>>> help(spss.single_sample)_
    
### single_sample_pcwag
    Decompose the query samples those are in pcwag format into the global signatures.
    
    single_sample_pcwag(csv96, csv192="", output="results", par_one=0.01, par_two=0.025, n_cpu=-1)
    
    Example: 
    -------
    >>> from sigproSS import spss, spss_pcwag 
    >>> csv96 = spss.importdata("pcwag96")  #'you can put the path of your own csv96 file here'
    >>> csv192 = spss.importdata("pcwag192") #''you can put the path of your own csv192 file here'
    >>> spss_pcwag.single_sample_pcwag(csv96, csv192, output="example_output")
    
  **To get help on the parameters and outputs of the "single_sample_pcwag" function, please write down the following line:**
    
  _>>> help(spss_pcwag.single_sample_pcwag)_    
    
## COPYRIGHT
This software and its documentation are copyright 2018 as a part of the sigProfiler project. The SigProfilerExtractor framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## CONTACT INFORMATION
Please address any queries or bug reports to S M Ashiqul Islam (Mishu) at m0islam.ucsd.edu
