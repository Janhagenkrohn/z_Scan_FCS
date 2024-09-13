# z_Scan_FCS
## Scripts for z-Scan FCS: Acquisition at Zeiss LSM980 and analysis

z-Scan FCS is a technique for measuring Fluorescence Correlation Spectroscopy in 2D systems, namely biological or biomimetic membranes.While z-Scan FCS is more time-consuming than standard FCS, it yields more robust results as it circumvents the problems of membrane FCS that come with limited accuracy positioning of the confocal observation volume on the membrane. Inaccurate focus positioning is not only is a source of noise. Any defocus will lead to biased results with reduced apparent diffusion coefficient, increased apparent particle number, and decreased apparent molecular brightness.  

For further details, see _e.g._:

Benda, ..., Hof 2003, doi.org/10.1021/la0270136

Macháň and Hof 2010, doi.org/10.1016/j.bbamem.2010.02.014

Betaneli, Mücksch, and Schwille 2019, doi.org/10.1007/978-1-4939-9512-7_18

The scripts here help perform z-Scan FCS conveniently, with emphasis on automation of multiple subsequent z-scan acquisitions within one sample, and batch-mode analysis. 


#### 01_Acquisition_ZENmacro
The first subfolder contains a macro for (single- or multi-position) z-scan FCS acquisition at a ZEISS LSM980 system with ZEISS FCS module. Note that we know that there are restrictions in terms of which ZEN Versions are compatible with this macro, but I do know exactly what those requirements are.

Instructions on how to use the macro are written directly into the header of the file, you can open and read it in a normal text editor. The macro also creates a .csv file with all the xyz measurement positions included in a single acquisition run, which is needed for the analysis in step 4.


#### 02_Correlation_MATLAB
These MATLAB-based scripts are for batch-mode calculation of auto- and cross-correlation functions from the .raw photon arrival times written by ZEN. 

Batch_correlation_LSM980.m is the script to be edited by the user. Specify a batch of .raw files to correlate, some basic correlator settings, and choose whether to perform calibrated subtraction of afterpulsing (this requires appropriate calibration, see Zhao, ..., Chen 2003, doi.org/10.1364/AO.42.004031) and whether to perform bleaching correction on the data. The output will be automatically written to .cvs files in the 4-column "Kristine" format used by Claus Seidel's lab. 

readConfoCor3.m is the data reader, and has its name from the fact that we wrote it for our old ZEISS ConfoCor3 system, but the data format has not changed. cross_corr.m and cross_corr_weights.m are functions for the correlation itself. cross_corr_weights.m is technically a FLCS/filterdFCS correlator, but used here for bleaching correction. Note that it is a rather naive implementation that is computationally expensive. get_blcorr_weights.m fits a polynomial model to the data to determine the bleaching kinetics, is correction of the latter is desired. detectorsD118_980.mat is a MATLAB data file containing the afterpulsing calibration information for our LSM980 system. 

cross_corr.m, cross_corr_weights.m, and readConfoCor3.m are modified versions of code taken from other sources, please see the copyright and license information specified in the headers of these functions.


#### 03_Individual_ACF_Fitting_Python
acf_fitting_script.py is the Python script meant for user interaction, with rather extensive description of parameters in the code. It expects to find correlation functions in the "Kristine" format produced by the software in 02_Correlation_MATLAB, but also by the software in my other repo "FCS_Fixer" for time-tagged, time-resolved data.

acf_fitting_functions.py contains the functions and classes used called by acf_fitting_script.py. 

Huge Thanks goes to Yusuf Qutbuddin, who wrote the original, more complex version of this code. In this repo I provide a version of the code that is simplified but reduced in functionality.
 
 
#### 04_zScan_Refit_Python
This single, relatively compact, script takes the output table in which 03_Individual_ACF_Fitting_Python compiled the individual correlation function fit parameters, associates them with the coordinate data written by the acquisition macro, and performs the re-fit of the z-Scan profile of count rate, apparent particle number, and apparent diffusion time to get accurate "in-focus" values for all three. 

Note that as described in the publications listed above, this is a self-calibrating analysis that determines the beam waist radius needed for physical interpretation of particle number (-> concentration) and diffusion time (-> diffusion coefficient) directly from the data. Depending on the accuracy of the automated xyz coordinate targeting, this may not be reliable though. If you see that the beam waist diameter is not robust/reproducible, you may be better off interpreting the obtained in-focus diffusion time and particle number (which are more robustly estimated than the beam waist diameter) based on standard FCS calibration, and still profit from automated acquisition and focus-finding.


#### Python environment requirements
The functions used are not especially fancy, the version requirements should be quite relaxed. W
e ran the code on various machines with Anaconda python environments that were inconsistent in module versions.
lmfit, numdifftools and uncertainties, is most likely to cause trouble: Most problems we encounter with the analysis pipeline relate to uncertainty calculations.

python 3.7.11

numpy 1.20.3

scipy 1.7.3

pandas 1.3.5

lmfit 1.2.2

matplotlilb 3.5.1

uncertainties 3.1.6

numdifftools 0.9.41 

