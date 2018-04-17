#!/bin/bash

evalParamFile="./tsurff.param"
propagateParamFile="./propagate.param"
initialParamFile="./initial.param"

num_phi_surff=$(grep -oP 'num-phi-surff long \K(.+)' $evalParamFile)
num_theta_surff=$(grep -oP 'num-theta-surff long \K(.+)' $evalParamFile)
expansion_scheme=$(grep -oP 'expansion-scheme long \K(.+)' $evalParamFile)
num_theta_partial=$(grep -oP 'num-theta-partial long \K(.+)' $evalParamFile)
num_phi_partial=$(grep -oP 'num-phi-partial long \K(.+)' $evalParamFile)
num_k_surff=$(grep -oP 'num-k-surff long \K(.+)' $evalParamFile)
ell_grid_size=$(grep -oP 'ell-grid-size long \K(.+)' $propagateParamFile)
cache_size_t=$(grep -oP 'cache-size-t long \K(.+)' $evalParamFile)
qprop_dim=$(grep -oP 'qprop-dim long \K(.+)' $initialParamFile)

if [ "$qprop_dim" -eq 34 ]
then ell_m_grid_size=$(echo "$ell_grid_size" | bc -l)
else ell_m_grid_size=$(echo "$ell_grid_size * $ell_grid_size" | bc -l)
fi

#ell_m_grid_size=(qprop_dim==34)?ell_grid_size:ell_grid_size*ell_grid_size;
#const long num_psi_size_surff=(expansion_scheme==1)?(ell_m_grid_size*num_theta_partial*num_phi_partial*num_k_proc):ell_m_grid_size*num_k_proc*ell_m_grid_size;
# the num_k_proc should be replaced by num_k_surff

if [ "$expansion_scheme" -eq 1 ]
then num_psi_size_surff=$(echo "$ell_m_grid_size * $num_theta_surff * $num_phi_surff * $num_k_surff" | bc -l)
else num_psi_size_surff=$(echo "$ell_m_grid_size * $num_k_surff * $ell_m_grid_size" | bc -l)
fi

let totalSize=0
let sizeofComplexDouble=8
let sizeofDouble=4
totalSize=$(echo "$sizeofComplexDouble * ($ell_m_grid_size * $num_k_surff + $num_psi_size_surff * 3 + $ell_m_grid_size * $cache_size_t * 2) + $sizeofDouble * $ell_m_grid_size * $cache_size_t * 4" | bc -l)

if [ "$qprop_dim" -eq 44 ]
then totalSize=$(echo "$totalSize + $sizeofComplexDouble * $num_psi_size_surff" | bc -l)
fi

totalSize_GB=$(echo "$totalSize / (1024.0 *1024.0 * 1024.0)" | bc -l)

#echo totalSize = $totalSize
#echo totalSize_GB = $totalSize_GB

echo $totalSize
#printf "%.5f GB\n" $totalSize_GB

#echo num_psi_size_surff = $num_psi_size_surff
#
#echo $expansion_scheme
#echo $num_theta_partial
#echo $num_phi_partial
#echo $num_k_surff
#echo $ell_grid_size
#echo $cache_size_t
#echo ell_m_grid_size $ell_m_grid_size
