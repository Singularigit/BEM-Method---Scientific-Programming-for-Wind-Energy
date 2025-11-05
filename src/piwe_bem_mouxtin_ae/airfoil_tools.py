"""
Tools for working with airfoil data in wind turbine BEM analysis.
Provides functions for interpolation and visualization of airfoil data.
"""

import numpy as np
import matplotlib.pyplot as plt


def interpolate_airfoil_coefficients(alpha_data, cl_data, cd_data, alpha_query):
    """
    Interpolates Cl and Cd for a given angle of attack alpha_query.

    Parameters:
    -----------
    alpha_data : array
        Array of angle of attack values [deg]
    cl_data : array
        Array of lift coefficient values
    cd_data : array
        Array of drag coefficient values
    alpha_query : float
        Angle of attack to interpolate for [deg]

    Returns:
    --------
    tuple
        (cl, cd) - Interpolated lift and drag coefficients
    """
    cl = np.interp(alpha_query, alpha_data, cl_data)
    cd = np.interp(alpha_query, alpha_data, cd_data)
    return cl, cd


def plot_airfoil_shapes(coords_database, num_airfoils=None, figsize=(12, 10)):
    """
    Plot airfoil shapes from the coordinate database.

    Parameters:
    -----------
    coords_database : dict
        Dictionary with airfoil indices as keys and (x, y) coordinate arrays as values
    num_airfoils : int, optional
        Number of airfoils to plot (None for all)
    figsize : tuple, optional
        Figure size

    Returns:
    --------
    fig : matplotlib.figure.Figure
        The figure with airfoil plots
    """
    fig, ax = plt.subplots(figsize=figsize)

    airfoil_indices = sorted(coords_database.keys())

    if num_airfoils is not None:
        airfoil_indices = airfoil_indices[:num_airfoils]

    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / len(airfoil_indices)) for i in range(len(airfoil_indices))]

    for i, af_idx in enumerate(airfoil_indices):
        x, y = coords_database[af_idx]
        ax.plot(x, y, color=colors[i], label=f"Airfoil {af_idx}")

    ax.set_xlabel('x/c [-]')
    ax.set_ylabel('y/c [-]')
    ax.set_title('Airfoil Shapes')
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.set_aspect('equal')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(min(airfoil_indices), max(airfoil_indices)))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Airfoil Index')

    plt.tight_layout()
    return fig


def plot_airfoil_polars(polar_database, airfoil_indices=None, figsize=(14, 8)):
    """
    Plot lift and drag coefficients for selected airfoils.

    Parameters:
    -----------
    polar_database : dict
        Dictionary with airfoil indices as keys and (alpha, cl, cd) arrays as values
    airfoil_indices : list, optional
        List of airfoil indices to plot (None for all)
    figsize : tuple, optional
        Figure size

    Returns:
    --------
    fig : matplotlib.figure.Figure
        The figure with polar plots
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    if airfoil_indices is None:
        airfoil_indices = sorted(polar_database.keys())

    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / len(airfoil_indices)) for i in range(len(airfoil_indices))]

    for i, af_idx in enumerate(airfoil_indices):
        alpha, cl, cd = polar_database[af_idx]
        ax1.plot(alpha, cl, color=colors[i], label=f"Airfoil {af_idx}")
        ax2.plot(alpha, cd, color=colors[i], label=f"Airfoil {af_idx}")

    ax1.set_xlabel('Angle of Attack [deg]')
    ax1.set_ylabel('Lift Coefficient (Cl) [-]')
    ax1.set_title('Lift Coefficient vs Angle of Attack')
    ax1.grid(True, linestyle='--', alpha=0.7)

    ax2.set_xlabel('Angle of Attack [deg]')
    ax2.set_ylabel('Drag Coefficient (Cd) [-]')
    ax2.set_title('Drag Coefficient vs Angle of Attack')
    ax2.grid(True, linestyle='--', alpha=0.7)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(min(airfoil_indices), max(airfoil_indices)))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax2)
    cbar.set_label('Airfoil Index')

    plt.tight_layout()
    return fig
