#!/bin/bash -l
#
#PBS -V
#PBS -l nodes=1:ppn=3
#PBS -N JOBNAME
#PBS -joe
#PBS -q batch
#PBS -m abe
#PBS -M tdunn@sdsu.edu

cd $PBS_O_WORKDIR
NCORES=`wc -w < $PBS_NODEFILE`
DATE=`date`
HOST=`hostname`

echo " "
echo "running on host: $HOSTNAME"
echo "$NCORES cores requested"
echo "job submitted: $DATE"
echo "job STDOUT follows:"
echo " "


source /home/tdunn/bashrc-ana3
export PATH=/usr/bin/conda:$PATH

conda env list
echo "activating my env"
conda activate specks_test
cd $PBS_O_WORKDIR

specks_py=/home/tdunn/git/SpecKS/SpecKS/main.py
config=CONFIGPATH

echo python3 $specks_py $config
python3 $specks_py $config

