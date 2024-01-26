# SpecKS
genome-level speciation simulation to recreate Ks curves

<br>

## Simulation overview
SpecKS is a pipeline application, with the following steps:
1) Make species trees for polyploids (custom code, allo and auto for a range of time frames)
2) Make gene trees (SaGePhy)
3) Relax gene trees (SaGePhy)
4) Evolve sequences through gene trees (PAML: evolver)
5) Prune trees to model post-WGD gene shedding
6) Calculate Ks for polyploid (PAML: codeml)
7) Plot histograms

<br>

## Prerequisites
* Run on an environment with PAML installed
* Have SaGePhy installed (get sagephy.zip from https://compbio.engr.uconn.edu/software/sagephy/)
* Specify the path to the sagephy jar file (ie, /home/Apps/sagephy/sagephy-1.0.0.jar) in the config file

<br>

## OS
* Tested on a mac and linux (ubuntu)
