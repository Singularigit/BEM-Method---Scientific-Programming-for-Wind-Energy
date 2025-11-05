"""
Blade Element Momentum (BEM) solver for wind turbine aerodynamic analysis.
Includes:
- solve_bem(): BEM calculation for a single operating point.
- compute_power_thrust_curves(): Generate power and thrust curves.
- compute_cp_ct_surfaces(): Compute Cp/Ct maps over pitch and TSR ranges.
"""

import numpy as np
from piwe_bem_mouxtin_ae.airfoil_tools import interpolate_airfoil_coefficients

def solve_bem(r, c, beta, af_id, v0, pitch, rpm, polar_database):
    """Solve BEM equations for thrust, torque, power, and induction factors."""
    rho = 1.225  # kg/m^3
    num_blades = 3
    rotor_radius = np.max(r)
    omega = rpm * 2 * np.pi / 60

    a = np.zeros_like(r)
    a_prime = np.zeros_like(r)
    dr = np.gradient(r)

    for i, radius in enumerate(r):
        if radius < 0.01 * rotor_radius:
            continue

        solidity = (num_blades * c[i]) / (2 * np.pi * radius)

        for _ in range(100):
            phi = np.arctan2((1 - a[i]) * v0, (1 + a_prime[i]) * omega * radius)

            sin_phi = np.sin(phi)
            if np.abs(sin_phi) < 1e-6:
                loss_factor = 0.99
            else:
                exponent = -((num_blades / 2) * (rotor_radius - radius) / (radius * sin_phi))
                tip_loss = (2 / np.pi) * np.arccos(np.exp(exponent))
                loss_factor = tip_loss

            alpha = np.rad2deg(phi) - (pitch + beta[i])
            if af_id[i] in polar_database:
                alpha_data, cl_data, cd_data = polar_database[af_id[i]]
                cl, cd = interpolate_airfoil_coefficients(alpha_data, cl_data, cd_data, alpha)
            else:
                cl, cd = 0.0001, 0.35

            cn = cl * np.cos(phi) + cd * np.sin(phi)
            ct = cl * np.sin(phi) - cd * np.cos(phi)

            if cn <= 0.01:
                new_a = 0
            else:
                term = (4 * loss_factor * sin_phi ** 2) / (solidity * cn)
                new_a = min(1 / (term + 1), 0.4)

            if ct <= 0.01 or np.abs(sin_phi) < 0.01 or np.abs(np.cos(phi)) < 0.01:
                new_a_prime = 0
            else:
                term = (4 * loss_factor * sin_phi * np.cos(phi)) / (solidity * ct)
                new_a_prime = 0 if term <= 1 else 1 / (term - 1)

            if np.abs(new_a - a[i]) < 1e-5 and np.abs(new_a_prime - a_prime[i]) < 1e-5:
                break

            relax = 0.5
            a[i] = relax * a[i] + (1 - relax) * new_a
            a_prime[i] = relax * a_prime[i] + (1 - relax) * new_a_prime

            a[i] = np.clip(a[i], 0, 0.5)
            a_prime[i] = np.clip(a_prime[i], -0.5, 0.5)

    d_thrust = np.zeros_like(r)
    d_torque = np.zeros_like(r)

    for i, radius in enumerate(r):
        if radius < 0.01 * rotor_radius:
            continue

        phi = np.arctan2((1 - a[i]) * v0, (1 + a_prime[i]) * omega * radius)
        sin_phi = np.sin(phi)
        if np.abs(sin_phi) < 1e-6:
            loss_factor = 0.99
        else:
            exponent = -((num_blades / 2) * (rotor_radius - radius) / (radius * sin_phi))
            loss_factor = (2 / np.pi) * np.arccos(np.exp(exponent))

        d_thrust[i] = (
            4 * np.pi * radius * rho * v0**2 * a[i] * (1 - a[i]) * loss_factor * dr[i]
        )
        d_torque[i] = (
            4 * np.pi * radius**3 * rho * v0 * omega * a_prime[i] * (1 - a[i])
            * loss_factor * dr[i]
        )

    thrust = np.sum(d_thrust)
    torque = np.sum(d_torque)
    power = torque * omega
    return thrust, torque, power, a, a_prime

def compute_power_thrust_curves(blade_data, operational_data, polar_database):
    """Compute power, thrust, and torque curves over wind speed range."""
    r, c, beta, af_id = blade_data
    v0_array, pitch_array, rpm_array = operational_data

    rated_power = 15000
    rated_wind_speed = 10.0

    for i in range(1, len(v0_array)):
        if pitch_array[i] > pitch_array[i - 1] + 1.0:
            rated_wind_speed = v0_array[i]
            break

    power_curve = []
    thrust_curve = []
    torque_curve = []

    for v0, pitch, rpm in zip(v0_array, pitch_array, rpm_array):
        thrust, torque, power, *_ = solve_bem(r, c, beta, af_id, v0, pitch, rpm, polar_database)

        if v0 > rated_wind_speed:
            power = min(power, rated_power * 1000)
            if v0 > rated_wind_speed + 2:
                thrust *= 0.85 * rated_wind_speed / v0

        power_curve.append(power / 1000)
        thrust_curve.append(thrust / 1000)
        torque_curve.append(torque)

    return v0_array, power_curve, thrust_curve, torque_curve

def compute_cp_ct_surfaces(blade_data, polar_database, pitch_range, tsr_range, v0=10.0):
    """Compute Cp and Ct surfaces as function of pitch and TSR."""
    r, c, beta, af_id = blade_data
    rotor_radius = np.max(r)
    rho = 1.225
    area = np.pi * rotor_radius**2

    pitch_grid, tsr_grid = np.meshgrid(pitch_range, tsr_range)
    cp_surface = np.zeros_like(pitch_grid)
    ct_surface = np.zeros_like(pitch_grid)

    for i, tsr_row in enumerate(tsr_grid):
        for j, tsr in enumerate(tsr_row):
            pitch = pitch_grid[i, j]
            rpm = (tsr * v0 * 60) / (2 * np.pi * rotor_radius)
            thrust, _, power, *_ = solve_bem(r, c, beta, af_id, v0, pitch, rpm, polar_database)

            cp = power / (0.5 * rho * area * v0**3)
            ct = thrust / (0.5 * rho * area * v0**2)
            cp_surface[i, j] = cp
            ct_surface[i, j] = ct

    return pitch_grid, tsr_grid, cp_surface, ct_surface
