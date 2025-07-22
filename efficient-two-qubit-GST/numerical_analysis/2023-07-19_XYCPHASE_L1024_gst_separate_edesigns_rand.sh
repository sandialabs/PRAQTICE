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

#I need this particular run to wait for another round of GST fits to finish so add a delay (I think 2 1/2 hours should work?)
echo "Sleeping while waiting for another job to finish"

sleep 2h 30m

echo "Waking up and running job now"
# Stage 3: Run a datagen sample for GST for all models

start=$SECONDS

#Per Germ Rand .03

echo "Running Per-Germ Random .03-Lite Coherent-Only"

mpirun -n 28 python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_03_Lite/ \
    --modelpack smq2Q_XYCPHASE \
    --mem_per_core 8 --use_mpi\
    --run_gst \
    --model_file H_XYCPHASE_1e-2_normal_0.pkl --gst_verbosity 4 \
    --num_clicks 1000 --datagen_seed $(( 2021 )) \
    --gst_modes "full TP,CPTP,Target,Truth" \
    1> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_0_pergermrand03_lite.out 2> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_0_pergermrand03_lite.err 
    
echo "Completed"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"        
    
echo "Running Per-Germ Random .03-Lite Coherent+Depol"    
    
mpirun -n 28 python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_03_Lite/ \
    --modelpack smq2Q_XYCPHASE \
    --mem_per_core 8 --use_mpi\
    --run_gst \
    --model_file H_XYCPHASE_1e-2_normal_depol_0.pkl --gst_verbosity 4 \
    --num_clicks 1000 --datagen_seed $(( 2021 )) \
    --gst_modes "full TP,CPTP,Target,Truth" \
    1> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_depol_0_pergermrand03_lite.out 2> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_depol_0_pergermrand03_lite.err 

echo "Completed"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"    

#Per Germ Rand .06

echo "Running Per-Germ Random .06-Lite Coherent-Only"

mpirun -n 28 python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_06_Lite/ \
    --modelpack smq2Q_XYCPHASE \
    --mem_per_core 8 --use_mpi\
    --run_gst \
    --model_file H_XYCPHASE_1e-2_normal_0.pkl --gst_verbosity 4 \
    --num_clicks 1000 --datagen_seed $(( 2021 )) \
    --gst_modes "full TP,CPTP,Target,Truth" \
    1> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_0_pergermrand06_lite.out 2> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_0_pergermrand06_lite.err 

echo "Completed"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"    

echo "Running Per-Germ Random .06-Lite Coherent+Depol"     
    
mpirun -n 28 python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_06_Lite/ \
    --modelpack smq2Q_XYCPHASE \
    --mem_per_core 8 --use_mpi\
    --run_gst \
    --model_file H_XYCPHASE_1e-2_normal_depol_0.pkl --gst_verbosity 4 \
    --num_clicks 1000 --datagen_seed $(( 2021 )) \
    --gst_modes "full TP,CPTP,Target,Truth" \
    1> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_depol_0_pergermrand06_lite.out 2> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_depol_0_pergermrand06_lite.err

echo "Completed"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"    

#Per Germ Rand .09

echo "Running Per-Germ Random .09-Lite Coherent-Only"

mpirun -n 28 python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_09_Lite/ \
    --modelpack smq2Q_XYCPHASE \
    --mem_per_core 8 --use_mpi\
    --run_gst \
    --model_file H_XYCPHASE_1e-2_normal_0.pkl --gst_verbosity 4 \
    --num_clicks 1000 --datagen_seed $(( 2021 )) \
    --gst_modes "full TP,CPTP,Target,Truth" \
    1> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_0_pergermrand09_lite.out 2> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_0_pergermrand09_lite.err 

echo "Completed"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"    

echo "Running Per-Germ Random .09-Lite Coherent+Depol"     
    
mpirun -n 28 python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_09_Lite/ \
    --modelpack smq2Q_XYCPHASE \
    --mem_per_core 8 --use_mpi\
    --run_gst \
    --model_file H_XYCPHASE_1e-2_normal_depol_0.pkl --gst_verbosity 4 \
    --num_clicks 1000 --datagen_seed $(( 2021 )) \
    --gst_modes "full TP,CPTP,Target,Truth" \
    1> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_depol_0_pergermrand09_lite.out 2> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_depol_0_pergermrand09_lite.err 

echo "Completed"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"    

#Per Germ Rand .125

echo "Running Per-Germ Random .125-Lite Coherent-Only"

mpirun -n 28 python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_12_Lite/ \
    --modelpack smq2Q_XYCPHASE \
    --mem_per_core 8 --use_mpi\
    --run_gst \
    --model_file H_XYCPHASE_1e-2_normal_0.pkl --gst_verbosity 4 \
    --num_clicks 1000 --datagen_seed $(( 2021 )) \
    --gst_modes "full TP,CPTP,Target,Truth" \
    1> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_0_pergermrand12_lite.out 2> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_0_pergermrand12_lite.err 

echo "Completed"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"    

echo "Running Per-Germ Random .125-Lite Coherent+Depol"     
    
mpirun -n 28 python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_12_Lite/ \
    --modelpack smq2Q_XYCPHASE \
    --mem_per_core 8 --use_mpi\
    --run_gst \
    --model_file H_XYCPHASE_1e-2_normal_depol_0.pkl --gst_verbosity 4 \
    --num_clicks 1000 --datagen_seed $(( 2021 )) \
    --gst_modes "full TP,CPTP,Target,Truth" \
    1> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_depol_0_pergermrand12_lite.out 2> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_depol_0_pergermrand12_lite.err 

echo "Completed"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"    

#Per Germ Rand .25

echo "Running Per-Germ Random .25-Lite Coherent-Only"

mpirun -n 28 python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_25_Lite/ \
    --modelpack smq2Q_XYCPHASE \
    --mem_per_core 8 --use_mpi\
    --run_gst \
    --model_file H_XYCPHASE_1e-2_normal_0.pkl --gst_verbosity 4 \
    --num_clicks 1000 --datagen_seed $(( 2021 )) \
    --gst_modes "full TP,CPTP,Target,Truth" \
    1> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_0_pergermrand25_lite.out 2> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_0_pergermrand25_lite.err 

echo "Completed"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"    

echo "Running Per-Germ Random .25-Lite Coherent+Depol"     
    
mpirun -n 28 python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_25_Lite/ \
    --modelpack smq2Q_XYCPHASE \
    --mem_per_core 8 --use_mpi\
    --run_gst \
    --model_file H_XYCPHASE_1e-2_normal_depol_0.pkl --gst_verbosity 4 \
    --num_clicks 1000 --datagen_seed $(( 2021 )) \
    --gst_modes "full TP,CPTP,Target,Truth" \
    1> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_depol_0_pergermrand25_lite.out 2> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_depol_0_pergermrand25_lite.err 
    
echo "Completed"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"        
    
#Per Germ Rand .50

echo "Running Per-Germ Random .5-Lite Coherent-Only"

mpirun -n 28 python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_50_Lite/ \
    --modelpack smq2Q_XYCPHASE \
    --mem_per_core 8 --use_mpi\
    --run_gst \
    --model_file H_XYCPHASE_1e-2_normal_0.pkl --gst_verbosity 4 \
    --num_clicks 1000 --datagen_seed $(( 2021 )) \
    --gst_modes "full TP,CPTP,Target,Truth" \
    1> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_0_pergermrand50_lite.out 2> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_0_pergermrand50_lite.err 

echo "Completed"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"    

echo "Running Per-Germ Random .5-Lite Coherent+Depol"     
    
mpirun -n 28 python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_50_Lite/ \
    --modelpack smq2Q_XYCPHASE \
    --mem_per_core 8 --use_mpi\
    --run_gst \
    --model_file H_XYCPHASE_1e-2_normal_depol_0.pkl --gst_verbosity 4 \
    --num_clicks 1000 --datagen_seed $(( 2021 )) \
    --gst_modes "full TP,CPTP,Target,Truth" \
    1> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_depol_0_pergermrand50_lite.out 2> ./gst_logs_L1024/H_XYCPHASE_1e-2_normal_depol_0_pergermrand50_lite.err 
    
echo "Completed"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"        












