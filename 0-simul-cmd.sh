#!/bin/bash
path_out=$HOME/technote_joha

 
# generate pairs of data (same m, overlap, cards) 
# 1.  rationals  [50, ~2000], simulate  abundances
# 2.  fractions [0,1] , simulate fractional contributions

# pair 1
#python3 data_simulator.py --m 50 --m_with_overlap 45 --card_a 11 --card_b 12 --value_min_whole_df 67 --value_max_whole_df 1335 $path_out

#python3 data_simulator.py --m 50 --m_with_overlap 45 --card_a 11 --card_b 12 --value_min_whole_df 0 --value_max_whole_df 1 $path_out


# all the other pairs (9):
mis=(20 22 30 35 40 45 50 55 60)
ovs=(17 29 27 31 36 40 42 50 50)
cardAs=(2 3 4 5 7 10 11 7 5)
cardBs=(2 3 3 5 7 10 10 7 5)
minims=(80 30 60 50 50 50 50 50 50)
maxims=(2200 777 1700 1968 1968 1968 1968 1968 1968)
items=(0 1 3 4 5 6 7 8)


for i in ${items[*]}; do
	echo ${mis[i]}
	echo ${ovs[i]}
	echo ${cardAs[i]}
	echo ${cardBs[i]}
	echo ${minims[i]}
	echo ${maxims[i]}
	echo ""
	
done


#	python3 data_simulator.py --m 19 --m_with_overlap 10 --card_a 6 --card_b 7 --value_min_whole_df 90 --value_max_whole_df 777 $path_out

#	python3 data_simulator.py --m 19 --m_with_overlap 10 --card_a 6 --card_b 7 --value_min_whole_df 0 --value_max_whole_df 1 $path_out
