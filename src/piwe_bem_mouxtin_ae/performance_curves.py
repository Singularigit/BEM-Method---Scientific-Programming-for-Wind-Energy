import numpy as np
import matplotlib.pyplot as plt
import os

def plot_performance_curves(v0_array, power_curve, thrust_curve, power_ref=None, thrust_ref=None, figsize=(12, 8)):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    ax1.plot(v0_array, power_curve, 'b-', linewidth=2, label='BEM Prediction')
    if power_ref is not None:
        ax1.plot(v0_array, power_ref, 'r--', linewidth=2, label='Reference')
    ax1.set_xlabel('Wind Speed [m/s]')
    ax1.set_ylabel('Power [kW]')
    ax1.set_title('Power Curve')
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend()

    ax2.plot(v0_array, thrust_curve, 'g-', linewidth=2, label='BEM Prediction')
    if thrust_ref is not None:
        ax2.plot(v0_array, thrust_ref, 'r--', linewidth=2, label='Reference')
    ax2.set_xlabel('Wind Speed [m/s]')
    ax2.set_ylabel('Thrust [kN]')
    ax2.set_title('Thrust Curve')
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.legend()

    plt.tight_layout()
    return fig

def plot_cp_ct_surfaces(pitch_grid, tsr_grid, cp_surface, ct_surface, figsize=(14, 6)):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    cp_levels = np.linspace(0, np.max(cp_surface), 20)
    cp_contour = ax1.contourf(pitch_grid, tsr_grid, cp_surface, levels=cp_levels, cmap='viridis')
    ax1.set_xlabel('Pitch Angle [deg]')
    ax1.set_ylabel('Tip Speed Ratio [-]')
    ax1.set_title('Power Coefficient (Cp)')
    fig.colorbar(cp_contour, ax=ax1, label='Cp [-]')

    ct_levels = np.linspace(0, np.max(ct_surface), 20)
    ct_contour = ax2.contourf(pitch_grid, tsr_grid, ct_surface, levels=ct_levels, cmap='plasma')
    ax2.set_xlabel('Pitch Angle [deg]')
    ax2.set_ylabel('Tip Speed Ratio [-]')
    ax2.set_title('Thrust Coefficient (Ct)')
    fig.colorbar(ct_contour, ax=ax2, label='Ct [-]')

    plt.tight_layout()
    return fig

def plot_spanwise_variables(r, a, a_prime, v0, rpm, figsize=(10, 8)):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)

    ax1.plot(r, a, 'b-', linewidth=2)
    ax1.set_xlabel('Blade Span [m]')
    ax1.set_ylabel('Axial Induction Factor (a)')
    ax1.set_title(f'Spanwise Axial Induction Distribution (V0={v0} m/s, RPM={rpm})')
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.set_xlim([r.min(), r.max()])

    ax2.plot(r, a_prime, 'r-', linewidth=2)
    ax2.set_xlabel('Blade Span [m]')
    ax2.set_ylabel("Tangential Induction Factor (a')")
    ax2.set_title(f'Spanwise Tangential Induction Distribution (V0={v0} m/s, RPM={rpm})')
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.set_xlim([r.min(), r.max()])

    plt.tight_layout()
    return fig

def plot_spanwise_forces(r, dT, dM, v0, rpm, figsize=(10, 8)):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)

    ax1.plot(r, dT, 'b-', linewidth=2)
    ax1.set_xlabel('Blade Span [m]')
    ax1.set_ylabel('dT/dr [N/m]')
    ax1.set_title(f'Spanwise Thrust Distribution (V0={v0} m/s, RPM={rpm})')
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.set_xlim([r.min(), r.max()])

    ax2.plot(r, dM, 'r-', linewidth=2)
    ax2.set_xlabel('Blade Span [m]')
    ax2.set_ylabel('dM/dr [Nm/m]')
    ax2.set_title(f'Spanwise Torque Distribution (V0={v0} m/s, RPM={rpm})')
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.set_xlim([r.min(), r.max()])

    plt.tight_layout()
    return fig

def calculate_aep(v0_array, power_curve, wind_dist, rated_power=15000):
    wind_dist = wind_dist / np.sum(wind_dist)
    power_produced = np.array(power_curve) * wind_dist * 8760
    aep = np.sum(power_produced) / 1000  # GWh
    cf = aep / (rated_power * 8.760)
    return aep, cf

def plot_operational_strategy(v0_array, pitch_array, rpm_array):
    fig, ax1 = plt.subplots()
    ax1.plot(v0_array, pitch_array, 'b-', label="Pitch angle (θₚ)")
    ax1.set_xlabel("Wind speed (m/s)")
    ax1.set_ylabel("Pitch angle (deg)", color='b')

    ax2 = ax1.twinx()
    ax2.plot(v0_array, rpm_array, 'r--', label="RPM")
    ax2.set_ylabel("Rotational speed (RPM)", color='r')

    plt.title("Optimal Operational Strategy")
    fig.tight_layout()
    plt.grid(True)
    return fig

def plot_cl_cd_vs_alpha(polar_database, af_id_example):
    alpha, cl, cd = polar_database[af_id_example]
    fig, ax = plt.subplots()
    ax.plot(alpha, cl, label="Cl")
    ax.plot(alpha, cd, label="Cd")
    ax.set_xlabel("Angle of attack (°)")
    ax.set_ylabel("Coefficient")
    ax.set_title(f"Cl and Cd vs α (Airfoil #{af_id_example})")
    ax.legend()
    ax.grid(True)
    return fig

def plot_cl_cd_vs_alpha_all(polar_database, output_dir):
    import matplotlib.pyplot as plt
    os.makedirs(output_dir, exist_ok=True)

    plt.figure(figsize=(10, 6))
    for af_id, (alpha, cl, cd) in polar_database.items():
        if len(alpha) == len(cl):
            plt.plot(alpha, cl, label=f"Cl-{af_id:02d}", alpha=0.7)
        if len(alpha) == len(cd):
            plt.plot(alpha, cd, label=f"Cd-{af_id:02d}", linestyle='--', alpha=0.7)

    plt.xlabel("Angle of attack (°)")
    plt.ylabel("Coefficient")
    plt.title("Cl and Cd vs α for All Airfoils")
    plt.grid(True)
    plt.legend(fontsize='x-small', ncol=2)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "cl_cd_vs_alpha_all.jpg"), dpi=300)
    plt.close()
