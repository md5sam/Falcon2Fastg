#!/bin/bash

#SBATCH -C new
#SBATCH --nodes=1
#SBATCH -t 14-00:00:00

source /galaxy/home/szr165/.bash_profile


fc_run.py fc_run_filtsub.cfg
