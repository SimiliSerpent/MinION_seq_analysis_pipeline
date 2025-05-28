#!/bin/bash

# Run in sshell to avoid frontal node memory exhaust
cd /env/cnrgh/proj/math_stats/scratch/vacus/Projet_CAPASVIR/scratch/Analysis/micromamba_envs
module load micromamba
micromamba create -p /env/cnrgh/proj/math_stats/scratch/vacus/Projet_CAPASVIR/scratch/Analysis/micromamba_envs/UMI-tools
micromamba activate ./UMI-tools
micromambe install -c bioconda -c conda-forge umi_tools
