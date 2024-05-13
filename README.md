# SpecKS
Genome-level polyploid simulation to recreate Ks curves

<br>

## Example command
>python3 SpecKS.py config.xml

<br>

## Simulation overview
SpecKS is a pipeline application, with the following steps:
1) Make species trees for polyploids (custom code, allo and auto for a range of time frames)
2) Make gene trees 
3) Run gene-birth-and-death model
4) Run ohnolog-shedding model
5) Evolve gene sequences along gene trees (PAML: evolver)
6) Calculate Ks for polyploid (PAML: codeml)
7) Plot histograms
8) Prepare the final results

<br>

## Prerequisites
* Run on an environment with PAML 4.10.7 installed

<br>

## OS
* Tested on a mac and linux (ubuntu)
