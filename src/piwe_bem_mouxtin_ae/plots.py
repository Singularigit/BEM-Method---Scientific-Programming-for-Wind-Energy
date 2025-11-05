"""
Visualization tools for wind turbine BEM analysis.
Functions for generating various plots from airfoil and polar data.
"""
import os
import matplotlib.pyplot as plt


def load_airfoil_data(filepath):
    """
    Load airfoil aerodynamic data from a polar file.
    
    Parameters:
    -----------
    filepath : str
        Path to the airfoil polar data file
    
    Returns:
    --------
    tuple
        (aoa, cl, cd, cm) - Arrays of angle of attack, lift, drag, and moment coefficients
    """
    aoa, cl, cd, cm = [], [], [], []
    with open(filepath, 'r', encoding='utf-8') as file:
        for line in file:
            line = line.strip()
            if line and not (line.startswith("#") or line.startswith("!")):
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        a = float(parts[0])
                        c = float(parts[1])
                        d = float(parts[2])
                        m = float(parts[3]) if len(parts) >= 4 else None
                        aoa.append(a)
                        cl.append(c)
                        cd.append(d)
                        cm.append(m)
                    except ValueError:
                        continue
    return aoa, cl, cd, cm


def load_coords(filepath):
    """
    Load airfoil coordinates from a coordinate file.
    
    Parameters:
    -----------
    filepath : str
        Path to the airfoil coordinate file
    
    Returns:
    --------
    tuple
        (x_coords, y_coords) - Arrays of x and y coordinates
    """
    x_coords, y_coords = [], []
    with open(filepath, 'r', encoding='utf-8') as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith("!"):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        x = float(parts[0])
                        y = float(parts[1])
                        x_coords.append(x)
                        y_coords.append(y)
                    except ValueError:
                        continue
    return x_coords, y_coords


def plot_polar_data(base_path, output_dir, polar_pattern="IEA-15-240-RWT_AeroDyn15_Polar_{:02d}.dat", start=0, end=50):
    """
    Plot polar data (lift, drag, and moment coefficients) for multiple airfoils.
    
    Parameters:
    -----------
    base_path : str
        Base path to the airfoil data directory
    output_dir : str
        Directory to save output plots
    polar_pattern : str, optional
        Format string pattern for polar filenames
    start : int, optional
        Starting airfoil index
    end : int, optional
        Ending airfoil index
    """
    os.makedirs(output_dir, exist_ok=True)

    # Lift coefficient plot
    plt.figure(figsize=(10, 6))
    for i in range(start, end):
        path = os.path.join(base_path, "Airfoils", polar_pattern.format(i))
        if os.path.exists(path):
            aoa, cl, _, _ = load_airfoil_data(path)
            if len(aoa) == len(cl):
                plt.plot(aoa, cl, label=f"{i:02d}")
    plt.xlabel("Angle of Attack (°)")
    plt.ylabel("Lift Coefficient (Cl)")
    plt.title("Cl vs AoA")
    plt.legend(fontsize='small', ncol=2)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "cl_vs_aoa.jpg"), dpi=300)
    plt.close()

    # Drag coefficient plot
    plt.figure(figsize=(10, 6))
    for i in range(start, end):
        path = os.path.join(base_path, "Airfoils", polar_pattern.format(i))
        if os.path.exists(path):
            aoa, _, cd, _ = load_airfoil_data(path)
            if len(aoa) == len(cd):
                plt.plot(aoa, cd, label=f"{i:02d}")
    plt.xlabel("Angle of Attack (°)")
    plt.ylabel("Drag Coefficient (Cd)")
    plt.title("Cd vs AoA")
    plt.legend(fontsize='small', ncol=2)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "cd_vs_aoa.jpg"), dpi=300)
    plt.close()

    # Moment coefficient plot
    plt.figure(figsize=(10, 6))
    for i in range(start, end):
        path = os.path.join(base_path, "Airfoils", polar_pattern.format(i))
        if os.path.exists(path):
            aoa, _, _, cm = load_airfoil_data(path)
            if len(aoa) == len(cm) and all(m is not None for m in cm):
                plt.plot(aoa, cm, label=f"{i:02d}")
    plt.xlabel("Angle of Attack (°)")
    plt.ylabel("Moment Coefficient (Cm)")
    plt.title("Cm vs AoA")
    plt.legend(fontsize='small', ncol=2)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "cm_vs_aoa.jpg"), dpi=300)
    plt.close()


def plot_coords_data(base_path, output_dir, coords_pattern="IEA-15-240-RWT_AF{:02d}_Coords.txt", start=0, end=50):
    """
    Plot airfoil coordinates for multiple airfoils.
    
    Parameters:
    -----------
    base_path : str
        Base path to the airfoil data directory
    output_dir : str
        Directory to save output plots
    coords_pattern : str, optional
        Format string pattern for coordinate filenames
    start : int, optional
        Starting airfoil index
    end : int, optional
        Ending airfoil index
    """
    os.makedirs(output_dir, exist_ok=True)
    plt.figure(figsize=(10, 6))
    for i in range(start, end):
        path = os.path.join(base_path, "Airfoils", coords_pattern.format(i))
        if os.path.exists(path):
            x, y = load_coords(path)
            if len(x) == len(y):
                plt.plot(x, y, label=f"{i:02d}")
    plt.xlabel("x/c")
    plt.ylabel("y/c")
    plt.title("Airfoil Shapes")
    plt.legend(fontsize='small', ncol=2)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "airfoil_shapes.jpg"), dpi=300)
    plt.close()