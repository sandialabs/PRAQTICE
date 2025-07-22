#!/bin/bash

DATA_DIR=/scratch/gst_resource_data/production_runs_1Q/data_XYI_L1024_rev1

# Venv containing pygsti develop feature-perL-fpr branch
if [[ "$USER" == "mgrace" ]]; then
  source /home/mgrace/Python/pyGSTi_develop/bin/activate
elif [[ "$USER" == "sserita" ]]; then
  source ../../venv/bin/activate
elif [[ "$USER" == "ciostro" ]]; then
  source  /home/ciostro/2_qubit_GST_resource/pyGSTi_fpr_bugfix/.venv/bin/activate
else
  echo "Must set the venv for this user"
  exit 1
fi

# Stage 4: Data Analysis
mprof run -C -M --python --output out_of_memory_testing_one_q_analysis.dat python ../../combined_gst.py --data_dir $DATA_DIR \
  --use_mp_pool --num_mp_workers 27 --maxL_power 10 \
  --run_analysis --df_file XYI_L1024_rev1.partial.h5 \
      1> ./analysis_logs_rev1/analysis.out 2> ./analysis_logs_rev1/analysis.err
