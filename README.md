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
* Run on an environment with PAML installed
* We suggest v4.10.7 but older versions also work. 

<br>
## Configuration
To automatically configure a conda environment suitable for SpecKS, use "specks_environment.yml" 
included in the sample_configs folder. Example commands to install and do a short run are below. 
"my-test-run.xml" should be a copy of "short-run.xml" with the paths modifed as needed.

> conda env create -f /PathToSpecKS/sample_configs/SpecKs_environment.yml 
> conda activate specks_env
> python3 /PathToSpecKS/SpecKS.py my-test-run.xml 

<br>

## OS
* Tested on a mac and linux (ubuntu)
