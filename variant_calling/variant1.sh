#!/bin/bash

###################################################### DOWNLOAD STEP ####################################################################

# make directories
folder="/Users/andreapassetti/src/CFTR"

mkdir ${folder}
mkdir ${folder}/Aligned ${folder}/Data ${folder}/Reads ${folder}/Results 

prefetch SRR2136533
fasterq-dump -x SRR2136533 -O ${folder}/Reads 
#-t /Users/andreapassetti/tmp/scratch

