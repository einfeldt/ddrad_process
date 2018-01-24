#!/bin/sh

## Submit this script using:

## qsub -pe openmp 16 parallel_stacks.sh

#$ -cwd
#$ -j y
## Runtime
#$ -l h_rt=200:00:00
## Send email when job complete
#$ -M tony.einfeldt@gmail.com
#$ -m e
#$ -pe openmp 16
#$ -R y
#$ -l h_vmem=4G

## Echo details to job output
echo Grid Engine Job began at `date` on `hostname`

## Requires /home/einfeldt/scratch/stacks/05_allreads

source ~/.bashrc
cd /home/einfeldt/scratch/stacks/06_populations
touch /home/einfeldt/scratch/stacks/06_populations/denovo_map.log

SEQLIST=`ls -1 /home/einfeldt/scratch/stacks/05_allreads/* | sed -e 's/^/-s /g'`
denovo_map.pl -S -T 16 -m 5 -M 2 -n 3 -b 1 -D "Corophium_radtags" -o /home/einfeldt/scratch/stacks/06_populations $SEQLIST

## Echo details to job output
echo Grid Engine Job completed at `date` on `hostname`


