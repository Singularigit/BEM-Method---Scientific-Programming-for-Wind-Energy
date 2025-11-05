# tests/test_bem.py

import os
import numpy as np
import pytest

# Use headless backend for matplotlib to avoid Tkinter errors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from piwe_bem_mouxtin_ae.data_loader import ( 
    load_blade_geometry,
    load_airfoil_polars,
    load_operational_strategy,
)
from piwe_bem_mouxtin_ae.bem_solver import solve_bem, compute_power_thrust_curves
from piwe_bem_mouxtin_ae.performance_curves import (
    plot_performance_curves,
    plot_operational_strategy,
    plot_cl_cd_vs_alpha,
    plot_spanwise_variables,
    calculate_aep,
)
from piwe_bem_mouxtin_ae.airfoil_tools import interpolate_airfoil_coefficients
from piwe_bem_mouxtin_ae.turbine_classes import GeneralWindTurbine, WindTurbine


@pytest.fixture(scope="session")
def paths_and_data():
    base = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    blade_file = os.path.join(base, "inputs", "IEA-15-240-RWT", "IEA-15-240-RWT_AeroDyn15_blade.dat")
    polar_folder = os.path.join(base, "inputs", "IEA-15-240-RWT", "Airfoils")
    opt_file = os.path.join(base, "inputs", "IEA-15-240-RWT", "IEA_15MW_RWT_Onshore.opt")

    r, c, beta, af_id = load_blade_geometry(blade_file)
    polar_db = load_airfoil_polars(polar_folder)

    return {
        "r": r,
        "c": c,
        "beta": beta,
        "af_id": af_id,
        "polar_db": polar_db,
        "opt_file": opt_file,
    }


def test_load_blade_and_airfoil(paths_and_data):
    r = paths_and_data["r"]
    c = paths_and_data["c"]
    beta = paths_and_data["beta"]
    af_id = paths_and_data["af_id"]
    polar_db = paths_and_data["polar_db"]

    assert isinstance(r, np.ndarray)
    assert isinstance(c, np.ndarray)
    assert isinstance(beta, np.ndarray)
    assert r.shape == c.shape == beta.shape
    assert isinstance(af_id, np.ndarray) and af_id.size > 0
    assert isinstance(polar_db, dict) and polar_db


def test_interpolate_airfoil_coefficients(paths_and_data):
    r = paths_and_data["r"]
    af_id = paths_and_data["af_id"]
    polar_db = paths_and_data["polar_db"]

    idx = len(r) // 2
    key = af_id[idx]
    alpha, cl_vals, cd_vals = polar_db[key]

    cl, cd = interpolate_airfoil_coefficients(alpha, cl_vals, cd_vals, 5.0)

    assert isinstance(cl, float)
    assert isinstance(cd, float)
    assert cl >= 0 and cd >= 0


def test_solve_bem_and_physical_results(paths_and_data):
    r = paths_and_data["r"]
    c = paths_and_data["c"]
    beta = paths_and_data["beta"]
    af_id = paths_and_data["af_id"]
    polar_db = paths_and_data["polar_db"]

    T, M, P, a, a_prime = solve_bem(r, c, beta, af_id, 10.0, 0.0, 7.5, polar_db)

    assert T > 0
    assert M > 0
    assert P > 0
    assert np.all(a >= 0) and np.all(a <= 1)
    assert np.all(a_prime >= 0) and np.all(a_prime <= 1)


def test_load_operational_strategy(paths_and_data):
    v0, pitch, rpm, power_ref, thrust_ref = load_operational_strategy(paths_and_data["opt_file"])
    assert len(v0) == len(pitch) == len(rpm)
    assert len(power_ref) == len(thrust_ref)


def test_compute_power_thrust_curves(paths_and_data):
    r = paths_and_data["r"]
    c = paths_and_data["c"]
    beta = paths_and_data["beta"]
    af_id = paths_and_data["af_id"]
    polar_db = paths_and_data["polar_db"]

    v0, pitch, rpm, power_ref, thrust_ref = load_operational_strategy(paths_and_data["opt_file"])
    v0_out, P_curve, T_curve, Q_curve = compute_power_thrust_curves(
        (r, c, beta, af_id), (v0, pitch, rpm), polar_db
    )

    assert isinstance(v0_out, (list, np.ndarray))
    L = len(v0_out)
    assert len(P_curve) == len(T_curve) == len(Q_curve) == L


def test_plotting_functions_do_not_raise(paths_and_data, tmp_path):
    polar_db = paths_and_data["polar_db"]
    af_id = paths_and_data["af_id"]

    fig1 = plot_cl_cd_vs_alpha(polar_db, af_id[0])
    fig1.savefig(tmp_path / "cl_cd.png")
    plt.close(fig1)

    v0, pitch, rpm, power_ref, thrust_ref = load_operational_strategy(paths_and_data["opt_file"])
    fig2 = plot_operational_strategy(v0, pitch, rpm)
    fig2.savefig(tmp_path / "ops.png")
    plt.close(fig2)

    v0_out, P_curve, T_curve, Q_curve = compute_power_thrust_curves(
        (paths_and_data["r"], paths_and_data["c"], paths_and_data["beta"], af_id),
        (v0, pitch, rpm), polar_db
    )
    fig3 = plot_performance_curves(v0_out, P_curve, T_curve, power_ref, thrust_ref)
    fig3.savefig(tmp_path / "perf.png")
    plt.close(fig3)

    assert (tmp_path / "cl_cd.png").exists()
    assert (tmp_path / "ops.png").exists()
    assert (tmp_path / "perf.png").exists()


def test_calculate_aep(paths_and_data):
    v0, pitch, rpm, power_ref, thrust_ref = load_operational_strategy(paths_and_data["opt_file"])
    wind_dist = np.ones_like(v0)
    aep, cf = calculate_aep(v0, power_ref, wind_dist)
    assert aep > 0
    assert 0 <= cf <= 1


def test_spanwise_plot(paths_and_data, tmp_path):
    r = paths_and_data["r"]
    a = np.linspace(0.1, 0.3, len(r))
    a_prime = np.linspace(0.01, 0.05, len(r))
    fig = plot_spanwise_variables(r, a, a_prime, v0=10.0, rpm=7.5)
    fig.savefig(tmp_path / "spanwise_plot.png")
    plt.close(fig)
    assert (tmp_path / "spanwise_plot.png").exists()


def test_general_windturbine_vectorized():
    turbine = GeneralWindTurbine(240, 150, 15000, 3, 10, 25)
    wind = np.array([2, 5, 10, 15, 30])
    power = turbine.get_power(wind)
    assert power.shape == wind.shape
    assert np.all(power >= 0)


def test_windturbine_interp():
    power_data = np.array([[3, 0], [10, 15000], [25, 0]])
    turbine = WindTurbine(240, 150, 15000, 3, 10, 25, power_data)
    assert turbine.get_power(5) > 0
    assert turbine.get_power(10) == 15000
    assert turbine.get_power(2) == 0

from piwe_bem_mouxtin_ae.performance_curves import plot_cl_cd_vs_alpha_all, plot_cp_ct_surfaces
from piwe_bem_mouxtin_ae.data_loader import load_airfoil_coordinates

def test_plot_cl_cd_vs_alpha_all(paths_and_data, tmp_path):
    polar_db = paths_and_data["polar_db"]
    output_dir = tmp_path / "polars"
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_cl_cd_vs_alpha_all(polar_db, str(output_dir))
    assert (output_dir / "cl_cd_vs_alpha_all.jpg").exists()


def test_plot_cp_ct_surfaces(paths_and_data):
    r = paths_and_data["r"]
    c = paths_and_data["c"]
    beta = paths_and_data["beta"]
    af_id = paths_and_data["af_id"]
    polar_db = paths_and_data["polar_db"]

    pitch_range = np.linspace(-2, 10, 3)
    tsr_range = np.linspace(5, 10, 3)

    pitch_grid, tsr_grid, cp_surface, ct_surface = plot_cp_ct_surfaces(
        (r, c, beta, af_id),
        polar_db,
        pitch_range,
        tsr_range,
        v0=10.0
    )

    assert cp_surface.shape == pitch_grid.shape == tsr_grid.shape
    assert ct_surface.shape == cp_surface.shape


def test_load_airfoil_coordinates(paths_and_data):
    base = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    airfoil_folder = os.path.join(base, "inputs", "IEA-15-240-RWT", "Airfoils")
    coords_db = load_airfoil_coordinates(airfoil_folder)
    assert isinstance(coords_db, dict)
    assert all(isinstance(coords_db[k], tuple) for k in coords_db)


def test_load_coords_from_plots_module(tmp_path):
    # minimal dummy test to improve coverage for plots.py load_coords
    from piwe_bem_mouxtin_ae.plots import load_coords
    dummy_file = tmp_path / "AF00_Coords.txt"
    dummy_file.write_text("!  x/c        y/c\n0.0 0.0\n1.0 0.0\n")
    coords = load_coords(str(dummy_file))
    assert isinstance(coords, tuple)
    assert len(coords[0]) == 2

def test_plot_cp_ct_surface_plot(tmp_path):
    from piwe_bem_mouxtin_ae.performance_curves import plot_cp_ct_surfaces

    # Create dummy pitch and TSR grids
    pitch_range = np.linspace(-2, 10, 3)
    tsr_range = np.linspace(5, 10, 3)
    pitch_grid, tsr_grid = np.meshgrid(pitch_range, tsr_range)

    # Create dummy Cp and Ct data
    cp_surface = np.random.uniform(0, 0.5, size=pitch_grid.shape)
    ct_surface = np.random.uniform(0.5, 1.5, size=pitch_grid.shape)

    # Call plotting function
    fig = plot_cp_ct_surfaces(pitch_grid, tsr_grid, cp_surface, ct_surface)

    # Save and verify
    output_file = tmp_path / "cp_ct_surface.png"
    fig.savefig(output_file)
    plt.close(fig)

    assert output_file.exists()

def test_plot_cp_ct_surfaces(tmp_path):
    from piwe_bem_mouxtin_ae.performance_curves import plot_cp_ct_surfaces

    # Create synthetic input grids and Cp/Ct surfaces
    pitch_range = np.linspace(-2, 10, 3)
    tsr_range = np.linspace(5, 10, 3)
    pitch_grid, tsr_grid = np.meshgrid(pitch_range, tsr_range)

    cp_surface = np.random.random(pitch_grid.shape)
    ct_surface = np.random.random(pitch_grid.shape)

    fig = plot_cp_ct_surfaces(pitch_grid, tsr_grid, cp_surface, ct_surface)

    # Save and assert image created
    output_path = tmp_path / "cp_ct_surface_plot.png"
    fig.savefig(output_path)
    plt.close(fig)

    assert output_path.exists()
    

from piwe_bem_mouxtin_ae.airfoil_tools import plot_airfoil_polars, plot_airfoil_shapes

def test_plot_airfoil_polars(tmp_path, paths_and_data):
    polar_db = paths_and_data["polar_db"]
    af_id = paths_and_data["af_id"]

    # FIX: wrap index in a list
    fig = plot_airfoil_polars(polar_db, [af_id[0]])
    path = tmp_path / "polar_plot.png"
    fig.savefig(path)
    plt.close(fig)

    assert path.exists()


def test_plot_airfoil_shapes(tmp_path):
    # Generate dummy airfoil coordinates
    dummy_coords = {
        1: (np.array([0.0, 0.5, 1.0]), np.array([0.0, 0.1, 0.0])),
        2: (np.array([0.0, 0.5, 1.0]), np.array([0.0, 0.2, 0.0])),
    }

    fig = plot_airfoil_shapes(dummy_coords)
    path = tmp_path / "airfoil_shapes.png"
    fig.savefig(path)
    plt.close(fig)

    assert path.exists()


def test_solve_bem_extreme_tip_loss(paths_and_data):
    from piwe_bem_mouxtin_ae.bem_solver import solve_bem

    r = paths_and_data["r"]
    c = paths_and_data["c"]
    beta = paths_and_data["beta"]
    af_id = paths_and_data["af_id"]
    polar_db = paths_and_data["polar_db"]

    # TSR ~ 0.1 to force extreme axial induction / tip loss behavior
    T, M, P, a, a_prime = solve_bem(r, c, beta, af_id, 10.0, 0.0, 0.1, polar_db)

    assert np.all(a >= 0) and np.all(a <= 1)
    assert np.all(a_prime >= 0)
    assert T >= 0 and M >= 0 and P >= 0


def test_compute_power_thrust_invalid_afid():
    from piwe_bem_mouxtin_ae.bem_solver import compute_power_thrust_curves

    # Airfoil ID missing from polar DB
    r = np.array([1.0, 2.0])
    c = np.array([0.5, 0.5])
    beta = np.array([0.0, 0.0])
    af_id = np.array([999, 999])  # invalid ID

    v0 = np.array([8.0])
    pitch = np.array([0.0])
    rpm = np.array([5.0])
    polar_db = {}  # empty polar DB

    v_out, P_curve, T_curve, Q_curve = compute_power_thrust_curves(
        (r, c, beta, af_id), (v0, pitch, rpm), polar_db
    )

    # Expect failure = outputs are all zero or nan
    assert np.allclose(P_curve, 0) or np.any(np.isnan(P_curve))


from piwe_bem_mouxtin_ae.plots import load_airfoil_data

def test_load_airfoil_data_sample():
    sample = "inputs/IEA-15-240-RWT/Airfoils/IEA-15-240-RWT_AeroDyn15_Polar_00.dat"
    aoa, cl, cd, cm = load_airfoil_data(sample)
    assert len(aoa) == len(cl) == len(cd)

from piwe_bem_mouxtin_ae.plots import load_coords

def test_load_coords_sample():
    sample = "inputs/IEA-15-240-RWT/Airfoils/IEA-15-240-RWT_AF00_Coords.txt"
    x, y = load_coords(sample)
    assert len(x) == len(y)

from piwe_bem_mouxtin_ae.plots import plot_polar_data

def test_plot_polar_data_runs(tmp_path):
    plot_polar_data(
        base_path="inputs/IEA-15-240-RWT",
        output_dir=tmp_path,
        start=0,
        end=1  # reduce time
    )
    assert (tmp_path / "cl_vs_aoa.jpg").exists()
    assert (tmp_path / "cd_vs_aoa.jpg").exists()
    # Cm plot may not exist if cm=None in data

from piwe_bem_mouxtin_ae.plots import plot_coords_data

def test_plot_coords_data_runs(tmp_path):
    plot_coords_data(
        base_path="inputs/IEA-15-240-RWT",
        output_dir=tmp_path,
        start=0,
        end=1
    )
    assert (tmp_path / "airfoil_shapes.jpg").exists()


import numpy as np
import pytest
from piwe_bem_mouxtin_ae.turbine_classes import GeneralWindTurbine, WindTurbine


def test_general_turbine_power_curve():
    turbine = GeneralWindTurbine(150, 100, 8000, 3, 12, 25, name="TestTurbine")
    wind_speeds = np.array([0, 2.5, 6, 12, 20, 30])
    power = turbine.get_power(wind_speeds)
    assert power[0] == 0  # Below cut-in
    assert power[-1] == 0  # Above cut-out
    assert power[2] > 0
    assert power[3] == 8000
    fig, ax = turbine.plot_power_curve()
    assert fig is not None and ax is not None


def test_wind_turbine_with_curve_data():
    power_curve = np.array([
        [0, 0],
        [3, 0],
        [6, 2000],
        [12, 8000],
        [25, 8000],
        [30, 0]
    ])
    turbine = WindTurbine(150, 100, 8000, 3, 12, 25, power_curve, name="DataTurbine")

    wind_speeds = np.array([0, 3, 6, 12, 25, 30])
    power = turbine.get_power(wind_speeds)
    assert np.allclose(power, [0, 0, 2000, 8000, 8000, 0])
    fig, ax = turbine.plot_power_curve()
    assert fig is not None and ax is not None


def test_get_power_scalar_behavior():
    power_curve = np.array([
        [0, 0],
        [6, 2000],
        [12, 8000],
        [25, 8000],
    ])
    turbine = WindTurbine(150, 100, 8000, 3, 12, 25, power_curve)

    assert turbine.get_power(0) == 0
    assert turbine.get_power(10) > 0
    assert turbine.get_power(30) == 0  # outside range

