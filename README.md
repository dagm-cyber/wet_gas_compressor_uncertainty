# IMECE2020-23711 Wet Gas Compressor Testing
Documenting code for paper IMECE2020-23711 Wet Gas Compressor Testing - Performance Uncertainty

Cite the code: 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3988097.svg)](https://doi.org/10.5281/zenodo.3988097)


## How to use the code
Set up the case to run in the %% Config section of the main.m script:
* Chose EOS
* Chose wet or dry input array

  .
  
  .
  
  .
  
*  chose number of pertubations e.g. n = 10000
  
  Run the main.m script.

## Results
All compressor performance results will be located in the Compressor structure

All properties from the flash points will be located in the FP structure

Sensitivity tables are saved as .csv files

MCM uncertainty results are saved as .csv files

