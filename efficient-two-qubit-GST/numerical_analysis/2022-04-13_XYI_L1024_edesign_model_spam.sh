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

# Stage 1: Edesign Generation
python ../../combined_gst.py --data_dir $DATA_DIR \
   --run_edesign_gen --modelpack smq1Q_XYI \
   --maxL_power 10 --fpr_samples 1 --fpr_types Full PerGerm PerGermPowerRand \
   --gr_types Full Lite Bare  --fpr_fractions .03 .06 .085 .125 .25 .5 \
   --num_per_germ_power_randomizations 500\
   --per_germ_power_min_iterations 50 --per_germ_power_eigval_tol .033 .25 .5\
   --per_germ_power_cond_tol 30 4 2 --use_mp_pool --num_mp_workers 15\
   --num_soln_returned 1 --type_soln_returned best --retry_for_smaller\
   1> ./edesign_gen_logs_rev1/edesign.out 2> ./edesign_gen_logs_rev1/edesign.err

#Stage 2: Noisy Model Generation
# Overotation on XYI only
for i in {0..9}; do
  python ../../combined_gst.py --data_dir $DATA_DIR \
    --run_model_gen --modelpack smq1Q_XYI \
    --model_file overrot_XYI_1e-2_${i}.pkl \
    --model_json spam.json overrot_XYI.json \
    --model_gen_seed $((2021 + ${i})) \
    --model_gen_verbose 1>> ./model_gen_logs_rev1/model1.out 2>> ./model_gen_logs_rev1/model1.err 
done


# Overotation on XYI + global depolarizing on ALL gates
for i in {0..9}; do
  python ../../combined_gst.py --data_dir $DATA_DIR \
    --run_model_gen --modelpack smq1Q_XYI \
    --model_file overrot_XYI_1e-2_global_depol_${i}.pkl \
    --model_json spam.json overrot_XYI.json global_depol_XYI.json \
    --model_gen_seed $((2051 + ${i})) \
    --model_gen_verbose 1>> ./model_gen_logs_rev1/model2.out 2>> ./model_gen_logs_rev1/model2.err
done
