
#!/bin/bash
#RESCALE_NAME="rubber_bearing_r3"
#RESCALE_ANALYSIS=abaqus
#RESCALE_ANALYSIS_VERSION=2024-golden
#RESCALE_CORE_TYPE=Ruby
#RESCALE_CORES=2
#RESCALE_WALLTIME=96:00
#RESCALE_ENV_LM_LICENSE_FILE=27000@fidelis-abaqus
cd ~/work/shared
abaqus job=rubber_bearing_r3 cpus=2 interactive