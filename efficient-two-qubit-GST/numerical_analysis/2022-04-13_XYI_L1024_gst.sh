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

# Stage 3: Run a datagen sample for GST for all models
for i in {8..9}; do
  for j in {0..9}; do
    ## XYI overrotation
    python ../../combined_gst.py --data_dir $DATA_DIR \
      --use_mp_pool --num_mp_workers 27 --run_gst \
      --model_file overrot_XYI_1e-2_${i}.pkl --gst_verbosity 0 \
      --num_clicks 1000 --datagen_seed $(( 2021 + ${j} )) \
      --gst_modes "full TP,CPTP,Target,Truth" \
      1> ./gst_logs_rev1/2_overrot_XYI.${i}.${j}.out 2> ./gst_logs_rev1/2_overrot_XYI.${i}.${j}.err

    ## XYI overrot + global depol
    python ../../combined_gst.py --data_dir $DATA_DIR \
      --use_mp_pool --num_mp_workers 27 --run_gst \
      --model_file overrot_XYI_1e-2_global_depol_${i}.pkl --gst_verbosity 0 \
      --num_clicks 1000 --datagen_seed $(( 2021 + ${j} )) \
      --gst_modes "full TP,CPTP,Target,Truth" \
      1> ./gst_logs_rev1/2_overrot_depol_XY.${i}.${j}.out 2> ./gst_logs_rev1/2_overrot_depol_XY.${i}.${j}.err
  done
done

