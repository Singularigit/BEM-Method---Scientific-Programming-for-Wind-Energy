"""Main execution script for wind turbine BEM model analysis."""
import os
import matplotlib.pyplot as plt


# Data loading
from piwe_bem_mouxtin_ae.data_loader import load_blade_geometry, load_airfoil_polars, load_operational_strategy

# BEM solver and performance curves
from piwe_bem_mouxtin_ae.bem_solver import solve_bem, compute_power_thrust_curves

# Plotting functions
from piwe_bem_mouxtin_ae.performance_curves import (
    plot_performance_curves,
    plot_spanwise_variables,
    plot_operational_strategy,
    plot_cl_cd_vs_alpha_all
)
from piwe_bem_mouxtin_ae.plots import plot_coords_data, plot_polar_data

def main():
    """
    Run full BEM analysis workflow:
      1. Plot raw coordinate & polar data
      2. Load geometry & operational data
      3. Generate induction & performance curves
      4. Run a single-case BEM solve and spanwise plot
      5. Save all figures to outputs/
    """
    # Paths and output directory setup
    base_path = os.path.join("inputs", "IEA-15-240-RWT")
    blade_file = os.path.join(base_path, "IEA-15-240-RWT_AeroDyn15_blade.dat")
    polar_folder = os.path.join(base_path, "Airfoils")
    operational_file = os.path.join(base_path, "IEA_15MW_RWT_Onshore.opt")
    output_dir = "outputs"
    os.makedirs(output_dir, exist_ok=True)


    # 1. Visualize coordinate and polar files
    # Plot all airfoil shapes
    plot_coords_data(base_path, output_dir, start=0, end=50)
    # Plot raw polar Cl, Cd, Cm for sample range
    plot_polar_data(base_path, output_dir, start=0, end=50)

    # 2. Load BEM inputs
    # Blade geometry: span (r), chord (c), twist (beta), airfoil IDs
    r, c, beta, af_id = load_blade_geometry(blade_file)
    # Airfoil polars dictionary: {index: (alpha, Cl, Cd)}
    polar_database = load_airfoil_polars(polar_folder)
    # Operational strategy: arrays of V0, pitch, RPM, plus reference power & thrust
    v0_array, pitch_array, rpm_array, power_ref, thrust_ref = load_operational_strategy(operational_file)

    # 3. Combined polar plot for all airfoils
    plot_cl_cd_vs_alpha_all(polar_database, output_dir)

    # 4. Plot optimal operational strategy
    fig_strategy = plot_operational_strategy(v0_array, pitch_array, rpm_array)
    fig_strategy.savefig(os.path.join(output_dir, "operational_strategy.jpg"), dpi=300, bbox_inches="tight")

    # 5. Compute full power & thrust curves
    v0_array, power_curve, thrust_curve, _ = compute_power_thrust_curves(
        (r, c, beta, af_id),
        (v0_array, pitch_array, rpm_array),
        polar_database
    )

    fig_perf = plot_performance_curves(v0_array, power_curve, thrust_curve, power_ref, thrust_ref)
    fig_perf.savefig(os.path.join(output_dir, "power_thrust_curves.jpg"), dpi=300, bbox_inches="tight")

    # 6. Single-case BEM solution (e.g., index 10)
    case_idx = 10
    T, M, P, a, a_prime = solve_bem(
        r, c, beta, af_id,
        v0_array[case_idx],
        pitch_array[case_idx],
        rpm_array[case_idx],
        polar_database
    )

    #Summary of results
    print("=== Turbine performance ===")
    print(f"Wind speed: {v0_array[case_idx]:.2f} m/s")
    print(f"Thrust: {T:.2f} N")
    print(f"Torque: {M:.2f} Nm")
    print(f"Power: {P / 1000:.2f} kW")

    # 7. Spanwise induction factors for the chosen case
    fig_span = plot_spanwise_variables(r, a, a_prime, v0_array[case_idx], rpm_array[case_idx])
    fig_span.savefig(os.path.join(output_dir, "spanwise_induction.jpg"), dpi=300, bbox_inches="tight")

    print("All figures saved in 'outputs/'")

if __name__ == "__main__":
    main()
