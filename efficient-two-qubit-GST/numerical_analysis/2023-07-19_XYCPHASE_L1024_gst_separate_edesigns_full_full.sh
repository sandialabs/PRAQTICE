#!/bin/bash

# Venv containing pygsti develop feature-perL-fpr branch
if [[ "$USER" == "mgrace" ]]; then
  source /home/mgrace/Python/pyGSTi_develop/bin/activate
elif [[ "$USER" == "sserita" ]]; then
  source ../../venv/bin/activate
elif [[ "$USER" == "ciostro" ]]; then
  source  /home/ciostro/2_qubit_GST_resource/pyGSTi_2Q_legacy_1/.venv/bin/activate
else
  echo "Must set the venv for this user"
  exit 1
fi

# Stage 3: Run a datagen sample for GST for all models

start=$SECONDS

#Full-Full

echo "Running Full-Full Coherent-Only"

mpirun -n 16 python ../../combined_gst.py --data_dir /data/projects/two_qubit_resource/data_XYCPHASE_L1024_Full_Full/ \
    --modelpack smq2Q_XYCPHASE \
    --mem_per_core 43 --use_mpi\
    --run_gst \
    --model_file H_XYCPHASE_1e-2_normal_0.pkl --gst_verbosity 4 \
    --num_clicks 1000 --datagen_seed $(( 2021 )) \
    --gst_modes "full TP,CPTP,Target,Truth" \
    1> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_0_full_full.out 2> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_0_full_full.err 
   
echo "Completed"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"    
    
echo "Running Full-Full Coherent+Depol"    
    
mpirun -n 16 python ../../combined_gst.py --data_dir /data/projects/two_qubit_resource/data_XYCPHASE_L1024_Full_Full/ \
    --modelpack smq2Q_XYCPHASE \
    --mem_per_core 43 --use_mpi\
    --run_gst \
    --model_file H_XYCPHASE_1e-2_normal_depol_0.pkl --gst_verbosity 4 \
    --num_clicks 1000 --datagen_seed $(( 2021 )) \
    --gst_modes "full TP,CPTP,Target,Truth" \
    1> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_depol_0_full_full.out 2> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_depol_0_full_full.err 

echo "Completed"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"    

