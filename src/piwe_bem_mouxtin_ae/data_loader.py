"""
Data loading utilities for wind turbine BEM analysis.

Provides functions to load blade geometry, airfoil data, and operational strategy.
"""

import os
import glob
import numpy as np


def load_blade_geometry(filepath):
    """
    Load blade geometry from AeroDyn15 blade data file.

    Parameters:
    -----------
    filepath : str
        Path to the AeroDyn15 blade data file

    Returns:
    --------
    tuple
        (r, c, beta, af_id) arrays for:
        r : blade span positions [m]
        c : chord lengths [m]
        beta : twist angles [deg]
        af_id : airfoil indices
    """
    data = np.loadtxt(filepath, skiprows=8)
    r = data[:, 0]
    c = data[:, 5]
    beta = data[:, 4]
    af_id = data[:, 6].astype(int)
    return r, c, beta, af_id


def load_airfoil_polars(airfoil_folder):
    """
    Load all airfoil polar files into a dictionary.

    Parameters:
    -----------
    airfoil_folder : str
        Path to folder with airfoil polar files

    Returns:
    --------
    dict
        Keys are airfoil indices, values are (alpha, Cl, Cd) arrays
    """
    polar_database = {}
    polar_files = glob.glob(os.path.join(airfoil_folder, "*Polar*.dat"))

    for filepath in polar_files:
        filename = os.path.basename(filepath)
        try:
            af_idx_str = filename.split("_")[-1].split(".")[0]
            af_idx = int(af_idx_str) + 1
        except (IndexError, ValueError):
            continue

        with open(filepath, encoding="utf-8") as file:
            lines = file.readlines()

        data_lines = []
        in_table = False
        for line in lines:
            if "!    Alpha      Cl      Cd" in line:
                in_table = True
                continue
            if in_table and not line.strip().startswith("!"):
                data_lines.append(line)

        try:
            data = np.loadtxt(data_lines)
            alpha, cl, cd = data[:, 0], data[:, 1], data[:, 2]
            polar_database[af_idx] = (alpha, cl, cd)
        except (ValueError, IndexError):
            continue

    return polar_database


def load_airfoil_coordinates(airfoil_folder):
    """
    Load all airfoil coordinate files into a dictionary.

    Parameters:
    -----------
    airfoil_folder : str
        Path to folder with coordinate files

    Returns:
    --------
    dict
        Keys are airfoil indices, values are (x, y) coordinate arrays
    """
    coords_database = {}
    coord_files = glob.glob(os.path.join(airfoil_folder, "*Coords*.txt"))

    for filepath in coord_files:
        filename = os.path.basename(filepath)
        try:
            af_idx_str = filename.split("_AF")[1].split("_")[0]
            af_idx = int(af_idx_str) + 1
        except (IndexError, ValueError):
            continue

        with open(filepath, encoding="utf-8") as file:
            lines = file.readlines()

        coord_start = None
        for i, line in enumerate(lines):
            if "!  x/c        y/c" in line and i > 10:
                coord_start = i + 1
                break

        if coord_start is not None:
            try:
                data = np.loadtxt(lines[coord_start:])
                x, y = data[:, 0], data[:, 1]
                coords_database[af_idx] = (x, y)
            except (ValueError, IndexError):
                continue

    return coords_database


def load_operational_strategy(filepath):
    """
    Load wind speed, pitch, rpm, power and thrust data.

    Parameters:
    -----------
    filepath : str
        Path to the operational strategy file

    Returns:
    --------
    tuple
        v0 : array - wind speeds [m/s]
        pitch : array - pitch angles [deg]
        rpm : array - rotational speeds [rpm]
        power : array - power output [kW]
        thrust : array - thrust force [kN]
    """
    data = np.loadtxt(filepath, skiprows=1)
    v0 = data[:, 0]
    pitch = data[:, 1]
    rpm = data[:, 2]
    power = data[:, 3]
    thrust = data[:, 4]
    return v0, pitch, rpm, power, thrust
