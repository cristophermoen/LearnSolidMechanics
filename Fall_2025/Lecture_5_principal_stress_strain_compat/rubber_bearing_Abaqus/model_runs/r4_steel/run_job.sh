
#!/bin/bash
#RESCALE_NAME="steel_bearing_r4"
#RESCALE_ANALYSIS=abaqus
#RESCALE_ANALYSIS_VERSION=2024-golden
#RESCALE_CORE_TYPE=Ruby
#RESCALE_CORES=2
#RESCALE_WALLTIME=96:00
#RESCALE_ENV_LM_LICENSE_FILE=27000@fidelis-abaqus
cd ~/work/shared
abaqus job=steel_bearing_r4 cpus=2 interactive