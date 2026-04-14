#!/bin/bash

# Zizhuang Miao
# This script runs fmriprep on the high-performance computing cluster at Dartmouth College
# The raw data following the BIDS format is located at ${BIDSDIR}
# See below or the paper for the details of the preprocessing pipeline

# Slurm parameters ___________________________________________________
#SBATCH --job-name=
#SBATCH --ntasks-per-node=1

#SBATCH --output=
#SBATCH --error=

#SBATCH --mem-per-cpu=

# Array jobs
#SBATCH --array=

# fmriprep parameters ________________________________________________
CONTAINER_IMAGE=""
BIDSDIR=""
WORKDIR=""
OUTPUTDIR=""
FREESURFER_LICENSE=""

subjects=()
PARTICIPANT_LABEL=${subjects[$((SLURM_ARRAY_TASK_ID - 1 ))]}
echo "array id: " ${SLURM_ARRAY_TASK_ID}, "subject id: " ${PARTICIPANT_LABEL}

singularity run \
  --cleanenv \
  -B ${BIDSDIR}:/bids \
  -B ${WORKDIR}:/work \
  -B ${OUTPUTDIR}:/out \
  -B ${FREESURFER_LICENSE}:${FREESURFER_LICENSE} \
  ${CONTAINER_IMAGE} \
  /bids /out participant --participant-label ${PARTICIPANT_LABEL} \
  -w /work/"${SLURM_ARRAY_TASK_ID}_" \
  -vv --notrack \
  --task-id narratives  \
  --ignore slicetiming \
  --fd-spike-threshold 0.9 \
  --dummy-scans 6 \
  --bold2t1w-dof 9 \
  --mem_mb 55000 --nprocs 12 \
  --fs-license-file ${FREESURFER_LICENSE} --skip-bids-validation \
  --random-seed 2022 --skull-strip-fixed-seed \
  --bids-database-dir /bids/bids_database \
  --fs-no-reconall 