#!/bin/bash

# Name of your job
#SBATCH --job-name="su36calibration"
#SBATCH --partition=amarsden

# Specify the name of the output file. The %j specifies the job ID
#SBATCH --output=output.txt

# Specify the name of the error file. The %j specifies the job ID
#SBATCH --error=error.txt

# The walltime you require for your job
#SBATCH --time=24:00:00

# Job priority. Leave as normal for now
#SBATCH --qos=normal

# Number of nodes are you requesting for your job. You can have 24 processors per node
#SBATCH --nodes=1

# Amount of memory you require per node. The default is 4000 MB per node
#SBATCH --mem=24GB

# Number of processors per node
#SBATCH --ntasks-per-node=24

# Send an email to this address when your job starts and finishes
#SBATCH --mail-user=oguzziya@stanford.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# Run normal batch commands

# Name of the executable you want to run on the cluster

$HOME/build-perigee/su36_calibrate/preprocess3d -cpu_size=24
$HOME/build-perigee/su36_calibrate/prepostproc3d -cpu_size=24
python3 calibrate_for_ecg.py

