#!/bin/bash

# conda activate dimet #  <- make sure you activated the env

dir_out="../output"


n_features=1000 #
n_overlapping=500
n_samples_per_group=(20 15 10 7 6 5 4 3)  # (20 10 7 5 4 3)
min_value=10
max_value=999



for i in ${n_samples_per_group[*]}; do

    python3 data_simulator.py --m $n_features --m_with_overlap $n_overlapping --card_a $i --card_b $i --value_min_whole_df $min_value --value_max_whole_df $max_value $dir_out
    
    file_simulated="../output/simulated_data/data_dsta-l_m${n_features}-overlap${n_overlapping}-a${i}-b${i}_${min_value}-${max_value}.tsv"
    
    echo $file_simulated
    
    python3 run_analysis.py $file_simulated $dir_out
    
    file_result_differential="../output/result_dsta-l_m${n_features}-overlap${n_overlapping}-a${i}-b${i}_${min_value}-${max_value}.tsv"

    	
done

echo "now please use now 'run_validation.py' (verify that the names of the files match)"




