#!/bin/bash

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

# Stage 4: Data Analysis  
  
start=$SECONDS

#Full-Full

echo "Running Full-Full"

python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Full_Full/ \
    --df_file H_XYCPHASE_L1024_full_full.partial.h5  \
    --maxL_power 10 --run_analysis --sep_edesign_mode \
    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_full_full.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_full_full.err & 
    
#echo "Running Full-Full Coherent+Depol"    
#    
#python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Full_Full/ \
#    --df_file H_XYCPHASE_L10241e-2_normal_depol_0_full_full.partial.h5  \
#    --maxL_power 10 --run_analysis --sep_edesign_mode \
#    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_full_full.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_full_full.err & 

#Full-Lite

echo "Running Full-Lite"

python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Full_Lite/ \
    --df_file H_XYCPHASE_L1024_full_lite.partial.h5  \
    --maxL_power 10 --run_analysis --sep_edesign_mode \
    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_full_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_full_lite.err & 
    
#echo "Running Full-Lite Coherent+Depol"    
#    
#python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Full_Lite/ \
#    --df_file H_XYCPHASE_L1024_1e-2_normal_depol_0_full_lite.partial.h5  \
#    --maxL_power 10 --run_analysis --sep_edesign_mode \
#    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_full_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_full_lite.err & 

#Per Germ E 0.033 C 30 N 1000    

echo "Running Per-Germ E033C30N1000-Lite"
    
python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_E033C30N1000_Lite/ \
    --df_file H_XYCPHASE_L1024_pergermE033C30N1000_lite.partial.h5  \
    --maxL_power 10 --run_analysis --sep_edesign_mode \
    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermE033C30N1000_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermE033C30N1000_lite.err & 

#echo "Running Per-Germ E033C30N1000-Lite Coherent+Depol"    
#    
#python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_E033C30N1000_Lite/ \
#    --df_file H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE033C30N1000_lite.partial.h5  \
#    --maxL_power 10 --run_analysis --sep_edesign_mode \
#    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE033C30N1000_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE033C30N1000_lite.err & 
#
#Per Germ E 0.033 C 30 N 20000    

echo "Running Per-Germ E033C30N20000-Lite"
    
python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_E033C30N20000_Lite/ \
    --df_file H_XYCPHASE_L1024_pergermE033C30N20000_lite.partial.h5  \
    --maxL_power 10 --run_analysis --sep_edesign_mode \
    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermE033C30N20000_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermE033C30N20000_lite.err & 

#echo "Running Per-Germ E033C30N20000-Lite Coherent+Depol"    
#    
#python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_E033C30N20000_Lite/ \
#    --df_file H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE033C30N20000_lite.partial.h5  \
#    --maxL_power 10 --run_analysis --sep_edesign_mode \
#    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE033C30N20000_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE033C30N20000_lite.err & 

#Per Germ E 0.25 C 4 N 1000    
    
echo "Running Per-Germ E25C4N1000-Lite"    
    
python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_E25C4N1000_Lite/ \
    --df_file H_XYCPHASE_L1024_pergermE25C4N1000_lite.partial.h5  \
    --maxL_power 10 --run_analysis --sep_edesign_mode \
    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermE25C4N1000_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermE25C4N1000_lite.err & 

#echo "Running Per-Germ E25C4N1000-Lite Coherent+Depol" 
#
#python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_E25C4N1000_Lite/ \
#    --df_file H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE25C4N1000_lite.partial.h5  \
#    --maxL_power 10 --run_analysis --sep_edesign_mode \
#    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE25C4N1000_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE25C4N1000_lite.err &

#Per Germ E 0.25 C 4 N 20000    

echo "Running Per-Germ E25C4N20000-Lite"   

python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_E25C4N20000_Lite/ \
    --df_file H_XYCPHASE_L1024_pergermE25C4N20000_lite.partial.h5  \
    --maxL_power 10 --run_analysis --sep_edesign_mode \
    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermE25C4N20000_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermE25C4N20000_lite.err & 

#echo "Running Per-Germ E25C4N20000-Lite Coherent+Depol" 
#
#python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_E25C4N20000_Lite/ \
#    --df_file H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE25C4N20000_lite.partial.h5  \
#    --maxL_power 10 --run_analysis --sep_edesign_mode \
#    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE25C4N20000_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE25C4N20000_lite.err &    

#Per Germ E 0.5 C 2 N 1000

echo "Running Per-Germ E5C2N1000-Lite Coherent-Only"   
    
python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_E5C2N1000_Lite/ \
    --df_file H_XYCPHASE_L1024_pergermE5C2N1000_lite.partial.h5  \
    --maxL_power 10 --run_analysis --sep_edesign_mode \
    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermE5C2N1000_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermE5C2N1000_lite.err & 

#echo "Running Per-Germ E5C2N1000-Lite Coherent+Depol"     
#    
#python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_E5C2N1000_Lite/ \
#    --df_file H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE5C2N1000_lite.partial.h5  \
#    --maxL_power 10 --run_analysis --sep_edesign_mode \
#    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE5C2N1000_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE5C2N1000_lite.err &

#Per Germ E 0.5 C 2 N 20000    

echo "Running Per-Germ E5C2N1000-Lite"   

python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_E5C2N20000_Lite/ \
    --df_file H_XYCPHASE_L1024_pergermE5C2N20000_lite.partial.h5  \
    --maxL_power 10 --run_analysis --sep_edesign_mode \
    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermE5C2N20000_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermE5C2N20000_lite.err & 

#echo "Running Per-Germ E5C2N20000-Lite Coherent+Depol"       
#    
#python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_E5C2N20000_Lite/ \
#    --df_file H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE5C2N20000_lite.partial.h5  \
#    --maxL_power 10 --run_analysis --sep_edesign_mode \
#    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE5C2N20000_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermE5C2N20000_lite.err &
 
#Per Germ Rand .03

echo "Running Per-Germ Random .03-Lite"

python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_03_Lite/ \
    --df_file H_XYCPHASE_L1024_pergermrand03_lite.partial.h5  \
    --maxL_power 10 --run_analysis --sep_edesign_mode \
    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermrand03_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermrand03_lite.err & 
     
#echo "Running Per-Germ Random .03-Lite Coherent+Depol"    
#    
#python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_03_Lite/ \
#    --df_file H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand03_lite.partial.h5  \
#    --maxL_power 10 --run_analysis --sep_edesign_mode \
#    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand03_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand03_lite.err & 

#Per Germ Rand .010

echo "Running Per-Germ Random .06-Lite"

python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_06_Lite/ \
    --df_file H_XYCPHASE_L1024_pergermrand06_lite.partial.h5  \
    --maxL_power 10 --run_analysis --sep_edesign_mode \
    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermrand06_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermrand06_lite.err & 

#echo "Running Per-Germ Random .010-Lite Coherent+Depol"     
#    
#python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_06_Lite/ \
#    --df_file H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand06_lite.partial.h5  \
#    --maxL_power 10 --run_analysis --sep_edesign_mode \
#    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand06_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand06_lite.err &

#Per Germ Rand .09

echo "Running Per-Germ Random .09-Lite"

python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_09_Lite/ \
    --df_file H_XYCPHASE_L1024_pergermrand09_lite.partial.h5  \
    --maxL_power 10 --run_analysis --sep_edesign_mode \
    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermrand09_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermrand09_lite.err & 

#echo "Running Per-Germ Random .09-Lite Coherent+Depol"     
#    
#python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_09_Lite/ \
#    --df_file H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand09_lite.partial.h5  \
#    --maxL_power 10 --run_analysis --sep_edesign_mode \
#    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand09_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand09_lite.err & 

#Per Germ Rand .125

echo "Running Per-Germ Random .125-Lite"

python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_12_Lite/ \
    --df_file H_XYCPHASE_L1024_pergermrand12_lite.partial.h5  \
    --maxL_power 10 --run_analysis --sep_edesign_mode \
    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermrand12_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermrand12_lite.err & 

#echo "Running Per-Germ Random .125-Lite Coherent+Depol"     
#    
#python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_12_Lite/ \
#    --df_file H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand12_lite.partial.h5  \
#    --maxL_power 10 --run_analysis --sep_edesign_mode \
#    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand12_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand12_lite.err & 

#Per Germ Rand .25

echo "Running Per-Germ Random .25-Lite"

python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_25_Lite/ \
    --df_file H_XYCPHASE_L1024_pergermrand25_lite.partial.h5  \
    --maxL_power 10 --run_analysis --sep_edesign_mode \
    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermrand25_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermrand25_lite.err & 

#echo "Running Per-Germ Random .25-Lite Coherent+Depol"     
#    
#python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_25_Lite/ \
#    --df_file H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand25_lite.partial.h5  \
#    --maxL_power 10 --run_analysis --sep_edesign_mode \
#    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand25_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand25_lite.err & 
    
#Per Germ Rand .50

echo "Running Per-Germ Random .5-Lite"

python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_50_Lite/ \
    --df_file H_XYCPHASE_L1024_pergermrand50_lite.partial.h5  \
    --maxL_power 10 --run_analysis --sep_edesign_mode \
    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermrand50_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_0_pergermrand50_lite.err & 

#echo "Running Per-Germ Random .5-Lite Coherent+Depol"     
#    
#python ../../combined_gst.py --data_dir /scratch/gst_resource_data/production_runs_2Q/data_XYCPHASE_L1024_Per_Germ_Rand_50_Lite/ \
#    --df_file H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand50_lite.partial.h5  \
#    --maxL_power 10 --run_analysis --sep_edesign_mode \
#    1> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand50_lite.out 2> ./analysis_logs_L1024/H_XYCPHASE_L1024_1e-2_normal_depol_0_pergermrand50_lite.err & 
    
wait    
    
echo "Completed Analysis"
end=$SECONDS
duration=$(( end - start ))
echo "Cumulative Time Elapsed : $duration seconds"
