#!/bin/bash
# Grid convergence study script
# Using uniform grid (stretch=1) for proper convergence analysis

echo "N,L2_error_mV,Relative_L2_percent,Linf_error_mV" > results/convergence_data.csv

# Use uniform grid (stretch=1) and lower concentration for larger Debye length
for N in 51 101 201 401 801 1601; do
    output=$(./build/pnp_solver --N $N --c0 0.01 --stretch 1.0 --output results/conv_${N}.dat 2>&1)
    L2=$(echo "$output" | grep "L2 error:" | head -1 | awk '{print $3}')
    L2_rel=$(echo "$output" | grep "Relative L2" | awk '{print $4}')
    Linf=$(echo "$output" | grep "L-inf" | awk '{print $4}')
    echo "$N,$L2,$L2_rel,$Linf" >> results/convergence_data.csv
    echo "N=$N: L2=$L2 mV, Rel_L2=$L2_rel %"
done

echo ""
echo "=== Convergence Data ==="
cat results/convergence_data.csv
