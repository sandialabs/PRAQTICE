#!/bin/bash

DATA_DIR=/scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L64_Full_Full

# Venv containing pygsti develop feature-perL-fpr branch
if [[ "$USER" == "mgrace" ]]; then
  source ../../venv_new/bin/activate
elif [[ "$USER" == "sserita" ]]; then
  source ../../venv/bin/activate
elif [[ "$USER" == "ciostro" ]]; then
  source  /home/ciostro/2_qubit_GST_resource/pyGSTi_2Q_legacy_1/.venv/bin/activate
else
  echo "Must set the venv for this user"
  exit 1
fi

# Stage 2: Noisy Model Generation
# Will generate 2 noisy models

echo "Generating H + SPAM models."
for i in {0..0}; do
  ## H errors + SPAM
  python -u ../../combined_gst.py --data_dir $DATA_DIR \
    --run_model_gen --modelpack smq2Q_XYCPHASE \
    --model_file H_XYCPHASE_1e-2_normal_${i}.pkl \
    --model_json H_XYCPHASE_1e-2_normal.json \
    tensored_depol_prep_1e-3_povm_1e-3.json \
    --model_gen_seed $((2021 + ${i})) \
    --model_gen_verbose\
    1> ./model_gen_logs/model.out 2> ./model_gen_logs/model.err
done

echo "Generating H + Depol + SPAM models."

for i in {0..0}; do
  ## H errors + depolarization + SPAM
  python -u ../../combined_gst.py --data_dir $DATA_DIR \
    --run_model_gen --modelpack smq2Q_XYCPHASE \
    --model_file H_XYCPHASE_1e-2_normal_depol_${i}.pkl \
    --model_json H_XYCPHASE_1e-2_normal.json \
    depol_XYCPHASE_1e-4_1e-3.json \
    tensored_depol_prep_1e-3_povm_1e-3.json \
    --model_gen_seed $((2031 + ${i})) \
    --model_gen_verbose\
     1>> ./model_gen_logs/model.out 2>> ./model_gen_logs/model.err
done