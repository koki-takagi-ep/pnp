#!/bin/bash
# Comprehensive grid convergence study script
# Tests 2nd-order accuracy under various conditions

cd "$(dirname "$0")/.."

echo "========================================"
echo "  Comprehensive Grid Convergence Study"
echo "========================================"
echo ""

# Create output directory
mkdir -p results/convergence

# Master CSV file for all results
MASTER_CSV="results/convergence/all_convergence_data.csv"
echo "test_name,c0_M,phi0_mV,L_nm,eps_r,stretch,N,dx_nm,L2_mV,order" > "$MASTER_CSV"

# Function to run convergence test for a given condition
run_convergence_test() {
    local test_name=$1
    local c0=$2
    local phi0=$3
    local L=$4
    local eps=$5
    local stretch=$6

    echo ""
    echo "========================================"
    echo "Test: $test_name"
    echo "  c0 = $c0 M, phi0 = $phi0 mV, L = $L nm, eps_r = $eps, stretch = $stretch"
    echo "========================================"

    # Calculate theoretical Debye length
    debye=$(python3 -c "import math; eps0=8.854e-12; kB=1.381e-23; T=298; e=1.602e-19; NA=6.022e23; c0_SI=$c0*1000; eps=$eps*eps0; lam=math.sqrt(eps*kB*T/(2*e*e*c0_SI*NA))*1e9; print(f'{lam:.3f}')")
    echo "  Debye length: $debye nm"

    # Check if L is sufficient (should be >= 10*lambda_D)
    ratio=$(python3 -c "print(f'{$L/$debye:.1f}')")
    echo "  L/lambda_D = $ratio"

    local prev_L2=""
    local prev_N=""

    printf "\n%-8s %-12s %-14s %-10s\n" "N" "dx [nm]" "L2 [mV]" "Order"
    echo "--------------------------------------------"

    for N in 51 101 201 401 801 1601; do
        # Calculate grid spacing
        dx=$(python3 -c "print(f'{$L/($N-1):.4f}')")

        # Run solver (standard PB model, open system)
        output=$(./build/pnp_solver --N $N --c0 $c0 --phi0 $phi0 --L $L --eps $eps --stretch $stretch --output results/convergence/conv_${test_name}_${N}.dat 2>&1)

        # Extract L2 error
        L2=$(echo "$output" | grep "L2 error:" | head -1 | awk '{print $3}')

        if [ -z "$L2" ]; then
            echo "ERROR: Failed to get L2 error for N=$N"
            echo "Solver output:"
            echo "$output"
            continue
        fi

        # Calculate convergence order
        if [ -n "$prev_L2" ] && [ -n "$prev_N" ]; then
            order=$(python3 -c "import math; print(f'{math.log(float(\"$prev_L2\")/float(\"$L2\"))/math.log(float(\"$N\")/float(\"$prev_N\")):.2f}')" 2>/dev/null || echo "N/A")
        else
            order="--"
        fi

        # Print results
        printf "%-8s %-12s %-14s %-10s\n" "$N" "$dx" "$L2" "$order"

        # Save to master CSV
        echo "$test_name,$c0,$phi0,$L,$eps,$stretch,$N,$dx,$L2,$order" >> "$MASTER_CSV"

        prev_L2=$L2
        prev_N=$N
    done

    # Check final order
    if [ -n "$order" ] && [ "$order" != "--" ] && [ "$order" != "N/A" ]; then
        order_check=$(python3 -c "o=float('$order'); print('PASS' if 1.9 <= o <= 2.1 else 'WARN' if 1.5 <= o <= 2.5 else 'FAIL')")
        echo ""
        echo "Final convergence order: $order - $order_check (expected: 2.00)"
    fi
}

# Build solver first
echo "Building solver..."
make -j4
if [ $? -ne 0 ]; then
    echo "Build failed!"
    exit 1
fi
echo ""

# ========================================
# Test 1: Baseline (current default)
# ========================================
run_convergence_test "baseline" 0.001 50 100 12 1.0

# ========================================
# Test 2: Different concentrations
# ========================================
run_convergence_test "c0_0.0001M" 0.0001 50 500 12 1.0   # Very dilute, large Debye length
run_convergence_test "c0_0.01M" 0.01 50 100 12 1.0       # Moderate
run_convergence_test "c0_0.1M" 0.1 50 100 12 1.0         # Higher concentration

# ========================================
# Test 3: Different potentials
# ========================================
run_convergence_test "phi_25mV" 0.001 25 100 12 1.0      # Low potential (linear regime)
run_convergence_test "phi_100mV" 0.001 100 100 12 1.0    # Higher potential
run_convergence_test "phi_200mV" 0.001 200 100 12 1.0    # High potential (nonlinear regime)

# ========================================
# Test 4: Different domain lengths
# ========================================
run_convergence_test "L_50nm" 0.001 50 50 12 1.0         # Shorter domain
run_convergence_test "L_200nm" 0.001 50 200 12 1.0       # Longer domain

# ========================================
# Test 5: Non-uniform grid
# ========================================
run_convergence_test "stretch_2.0" 0.001 50 100 12 2.0   # Moderate stretching
run_convergence_test "stretch_3.0" 0.001 50 100 12 3.0   # Stronger stretching

# ========================================
# Test 6: Aqueous electrolyte (higher dielectric)
# Note: eps=80 gives λ_D ≈ 10 nm, so we need L ≥ 300 nm for L/λ_D ≥ 30
# ========================================
run_convergence_test "eps_80" 0.001 50 300 80 1.0        # Water-like, longer domain for larger λ_D

# ========================================
# Summary
# ========================================
echo ""
echo "========================================"
echo "  Summary"
echo "========================================"
echo ""
echo "All results saved to: $MASTER_CSV"
echo ""

# Count pass/warn/fail
echo "Analyzing convergence orders..."
python3 << 'PYEOF'
import csv
import sys

pass_count = 0
warn_count = 0
fail_count = 0
results = []

with open('results/convergence/all_convergence_data.csv', 'r') as f:
    reader = csv.DictReader(f)
    current_test = None
    last_order = None

    for row in reader:
        if row['test_name'] != current_test:
            if current_test is not None and last_order is not None:
                try:
                    o = float(last_order)
                    if 1.9 <= o <= 2.1:
                        status = "PASS"
                        pass_count += 1
                    elif 1.5 <= o <= 2.5:
                        status = "WARN"
                        warn_count += 1
                    else:
                        status = "FAIL"
                        fail_count += 1
                    results.append((current_test, o, status))
                except:
                    pass
            current_test = row['test_name']
        last_order = row['order']

    # Handle last test
    if current_test is not None and last_order is not None:
        try:
            o = float(last_order)
            if 1.9 <= o <= 2.1:
                status = "PASS"
                pass_count += 1
            elif 1.5 <= o <= 2.5:
                status = "WARN"
                warn_count += 1
            else:
                status = "FAIL"
                fail_count += 1
            results.append((current_test, o, status))
        except:
            pass

print(f"{'Test Name':<20} {'Order':>8} {'Status':>8}")
print("-" * 40)
for name, order, status in results:
    print(f"{name:<20} {order:>8.2f} {status:>8}")

print("-" * 40)
print(f"PASS: {pass_count}, WARN: {warn_count}, FAIL: {fail_count}")

if fail_count > 0:
    print("\n!!! SOME TESTS FAILED - INVESTIGATION NEEDED !!!")
    sys.exit(1)
elif warn_count > 0:
    print("\nSome tests showed warnings - may need attention.")
else:
    print("\nAll tests passed with 2nd-order convergence!")
PYEOF

echo ""
echo "Done."
