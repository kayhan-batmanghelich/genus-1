#!/bin/bash

#SBATCH -t 2-00:00:00
#SBATCH --mem=20G
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --qos=gablab

python logreg.py
