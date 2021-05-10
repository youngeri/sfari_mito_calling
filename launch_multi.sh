#!/bin/bash
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -J cromwell_wrap
#SBATCH --mem=8000

java -Dconfig.file=slurm.conf -jar \
/home/ml2529/shared/tools/jars/cromwell-56.jar run \
MitochondriaPipeline_multisample.wdl \
-i mito_multisample_pipeline_inputs.json  \
-o workflow.options



