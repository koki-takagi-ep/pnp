#!/bin/bash
# Grid convergence study script
# Tests 2nd-order accuracy of the discretization scheme

cd "$(dirname "$0")/.."

echo "========================================"
echo "  Grid Convergence Study"
echo "========================================"
echo ""
echo "Test conditions:"
echo "  - Bulk concentration: c0 = 0.001 mol/L"
echo "  - Surface potential: phi0 = 50 mV"
echo "  - Domain length: L = 50 nm"
echo "  - Debye length: ~11.9 nm"
echo "  - Grid: uniform (stretch=1.0)"
echo ""

# Create output directory
mkdir -p results

# Output CSV file
echo "N,dx_nm,L2_error_mV,Relative_L2_percent,Linf_error_mV,convergence_order" > results/convergence_data.csv

prev_L2=""
prev_N=""

echo "Running convergence tests..."
echo ""
printf "%-8s %-12s %-14s %-10s\n" "N" "dx [nm]" "L2 [mV]" "Order"
echo "--------------------------------------------"

for N in 51 101 201 401 801 1601; do
    # Calculate grid spacing
    dx=$(python3 -c "print(f'{50.0/($N-1):.4f}')")

    # Run solver
    output=$(./build/pnp_solver --N $N --c0 0.001 --phi0 50 --stretch 1.0 --L 50 --output results/conv_${N}.dat 2>&1)

    # Extract errors
    L2=$(echo "$output" | grep "L2 error:" | head -1 | awk '{print $3}')
    L2_rel=$(echo "$output" | grep "Relative L2" | awk '{print $4}')
    Linf=$(echo "$output" | grep "L-inf" | awk '{print $4}')

    # Calculate convergence order
    if [ -n "$prev_L2" ] && [ -n "$prev_N" ]; then
        order=$(python3 -c "import math; print(f'{math.log(float(\"$prev_L2\")/float(\"$L2\"))/math.log(float(\"$N\")/float(\"$prev_N\")):.2f}')")
    else
        order="--"
    fi

    # Print results
    printf "%-8s %-12s %-14s %-10s\n" "$N" "$dx" "$L2" "$order"

    # Save to CSV
    echo "$N,$dx,$L2,$L2_rel,$Linf,$order" >> results/convergence_data.csv

    prev_L2=$L2
    prev_N=$N
done

echo "--------------------------------------------"
echo ""
echo "Results saved to results/convergence_data.csv"
echo "Run 'python3 scripts/plot_convergence.py' to generate plots."
