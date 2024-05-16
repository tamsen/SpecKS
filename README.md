# SpecKS
Genome-level forward simulation of polyploid evolution. Produces Ks data suitable for generating simulated Ks histograms.

<br>

## Example command
>python3 SpecKS.py config.xml

<br>

## Simulation overview
SpecKS is a pipeline application, with the following steps:
1) Make species trees for polyploids 
2) Make gene trees 
3) Run gene-birth-and-death model
4) Run ohnolog-shedding model
5) Evolve gene sequences along gene trees (PAML: evolver)
6) Calculate Ks for polyploid (PAML: codeml)
7) Plot simple histograms*
8) Prepare the final results (.csv files with Ks data)

* Note, SpecKS is not intended to produce beautiful histograms. It is left to the user to upload the Ks datafiles generated by SpecKS to the visualization tools of choice.
<br>

## Prerequisites
* Run on an environment with PAML installed. We suggest v4.10.7 but have found that versions also work. 
* We give directions for configuring an appropriate conda environment below.

<br>

## Configuration
To automatically configure a conda environment suitable for SpecKS, use the "specks_environment.yml" 
included in the "sample_configs" folder. Example commands to install and do a short run are below.
"my-test-run.xml" should be a copy of "short-run.xml" with the paths modifed as needed.

> conda env create -f /PathToSpecKS/sample_configs/SpecKs_environment.yml

> conda activate specks_env

> python3 /PathToSpecKS/SpecKS.py my-test-run.xml 

<br>

## Help
* There are example config files and an overview of the parameters in the "sample_configs" folder.

## OS
* Tested on a mac and linux (ubuntu)

## Referencing SpecKS
* Publication forthcoming. Please reference "Accurately Inferring Ancient Auto and Allopolyploidization Events using Forward-time Simulations" by T. Dunn and A. Sethuraman
