# Wind Turbine Modeling Using Blade‑Element‑Momentum (BEM) Theory

# **Team**: Mouxtin A.E

## **Overview**
This repository implements a steady‑state Blade‑Element‑Momentum (BEM) model to predict the aerodynamic performance of a reference IEA 15 MW offshore wind turbine. It computes power, thrust, and torque as functions of wind speed, rotor speed, and blade pitch. The model simulates he aerodynamic performance of a horizontal-axis wind turbine. The BEM formulation combines momentum theory and blade element theory to iteratively solve for the axial and tangential induction factors \( a \) and \( a' \), which modify the local inflow velocity experienced by the blade elements.

The aerodynamic forces are computed using lift and drag coefficients $c_l$ and $C_d$, interpolated from airfoil polars as functions of the angle of attack $\alpha$. These are used to calculate the normal and tangential force coefficients:

```math
C_n = C_l \cos\phi + C_d \sin\phi, \quad
C_t = C_l \sin\phi - C_d \cos\phi
```

With these, the induction factors are updated using:

```math
a = \frac{1}{\left( \frac{4 \sin^2 \phi}{\sigma(r) C_n} + 1 \right)}, \quad
a' = \frac{1}{\left( \frac{4 \sin \phi \cos \phi}{\sigma(r) C_t} - 1 \right)}
```

Once convergence is achieved, the distributed loads are integrated over the blade span to yield total thrust \( T \), torque \( M \), and power \( P \):

```math
P = M \omega, \quad
C_P = \frac{P}{\frac{1}{2} \rho A V_0^3}, \quad
C_T = \frac{T}{\frac{1}{2} \rho A V_0^2}
```

The goal is to provide a computational tool for evaluating performance metrics as functions of inflow wind speed $V_0$, rotor speed $\omega$, and blade pitch angle $\theta_p$, enabling aerodynamic analysis and operational optimization of large-scale wind turbines.


---

## **Quick Start Guide**
To set up and run the simulation, follow these steps:

### **1️. Clone the Repository somewhere**
```sh
git clone https://github.com/Singularigit/BEM-Method---Scientific-Programming-for-Wind-Energy.git
```
```sh
cd BEM-Method---Scientific-Programming-for-Wind-Energy
```
Ensure you have Python 3.11+:
```sh
python --version
```
<!-- ```sh
conda create -n windsim python=3.11 -y
conda activate windsim
pip install -r requirements.txt
``` -->

### **2. Set Up Virtual Environment (recommended)**
```sh
python -m venv .venv
```
Do not forget to switch to this beautiful new enviroment
```sh
source .venv/Scripts/activate # GitBash
```
```sh
source .venv/bin/activate   # Mac/Linux
```
```sh
.\.venv\Scripts\Activate.ps1 # Windows PowerShell
```

### **3. Install dependencies of the package**
```sh
pip install --upgrade pip
```
```sh
pip install -e .[dev]
```
Sit back, relax and wait for approximately 2 minutes, no need to worry, in my village they say ''The good thing, takes time''.
   

### **4. Run the package**
```sh
python examples/main.py
```
   This will load inputs/, run the BEM solver and save figures to outputs/.

---

## **Project Structure**
The package uses a simple, layered architecture: a **data-loading** layer (`data_loader.py`, `turbine_classes.py`), a **computation** layer (`bem_solver.py`, `performance_curves.py`) runs the BEM algorithm and assembles thrust, torque, and power curves.A **visualization** layer (`plots.py`, `airfoil_tools.py`) produces all of the figures. An `examples/main.py` script ties these layers together and a parallel `tests/` suite verifies the modules with > 80 % coverage.

```text
final-project-mouxtin-ae/
├── COLLABORATION.md             # team collaboration plan
├── LICENSE                      # project license
├── README.md                    # you are here
├── examples/                   
│   ├── LEANWIND_8MW_164_...
│   └── main.py                  # executable example
├── inputs/                      # provided turbine geometry, polars & strategy
│   └── IEA-15-240-RWT/
│       ├── Airfoils/            # airfoil coordinate & polar files
│       ├── IEA-15-240-...
│       ├── IEA_15MW_...
│       └── rotor_diagram.jpeg
├── outputs/                     # generated plots
│   ├── airfoil_shapes.jpg
│   ├── cd_vs_aoa.jpg
│   ├── cl_cd_vs_alpha_all.jpg
│   ├── cl_vs_aoa.jpg
│   ├── cm_vs_aoa.jpg
│   ├── operational_strategy.jpg
│   ├── power_thrust_curves.jpg
│   └── spanwise_induction.jpg
├── pyproject.toml               # project metadata & dependencies
├── src/
│   └── piwe_bem_mouxtin_ae/     # installable BEM package by our great team
│       ├── _init_.py
│       ├── airfoil_tools.py
│       ├── bem_solver.py
│       ├── data_loader.py
│       ├── performance_curves.py
│       ├── plots.py
│       └── turbine_classes.py
├── tests/                       # pytest scripts
│   └── test_bem.py
└── .gitignore                   # files and folders to ignore
```
---
## Classes

### `src/piwe_bem_mouxtin_ae/turbine_classes.py`

#### GeneralWindTurbine
Models a wind turbine using a theoretical “cubic” power curve.

**Constructor**  
`GeneralWindTurbine(rotor_diameter, hub_height, rated_power, v_in, v_rated, v_out, name=None)`

**Methods**  
- `get_power(wind_speed)` → returns power (kW) for scalar or array inputs  
- `plot_power_curve(v_range=None, label=None)` → returns a Matplotlib figure

---

#### WindTurbine
Subclass of `GeneralWindTurbine` that uses actual power-curve data via interpolation.

**Constructor**  
`WindTurbine(rotor_diameter, hub_height, rated_power, v_in, v_rated, v_out, power_curve_data, name=None)`

**Methods**  
- `get_power(wind_speed)` → interpolated power values

---
## Testing & Quality

- **Tests:** 26 tests executed, all passed. Code coverage is **89%**, exceeding the 80% requirement (`pytest --cov=src tests/`).
- **Linting:** `pylint src/` reports a score of **9.19/10**, well above the 8.0 threshold.

To verify on your own machine, and you installed in edditable mode with dev dependencies:
### run the full test to check coverage
```sh
pytest --cov=src tests/
```
### run linting to check code quality and style
```sh
pylint src/
```
---

## Generated Output (view in your local Markdown preview)

The package produces 8 figures demonstrating each core capability:

1. **Blade & Airfoil Data**  
   - **Data loading & parsing** (`src/piwe_bem_mouxtin_ae/data_loader.py`):  
     - `load_blade_geometry`  
     - `load_airfoil_polars`  
     - `load_airfoil_coordinates`  
     - `load_operational_strategy`

2. **Airfoil Shapes** (`src/piwe_bem_mouxtin_ae/plots.py`)  
   - `plot_coords_data(base_path, output_dir, …)`  
   <div align="center">
     <img src="outputs/airfoil_shapes.jpg" width="50%" />
   </div>

3. **Airfoil Polars** (`src/piwe_bem_mouxtin_ae/plots.py`)  
   - **`plot_polar_data(...)`** generates three plots:
     <table width="100%">
       <tr>
         <td align="center">
           <img src="outputs/cl_vs_aoa.jpg" width="100%" /><br>
           **Cl vs AoA**
         </td>
         <td align="center">
           <img src="outputs/cd_vs_aoa.jpg" width="100%" /><br>
           **Cd vs AoA**
         </td>
       </tr>
       <tr>
         <td colspan="2" align="center">
           <img src="outputs/cm_vs_aoa.jpg" width="50%" /><br>
           **Cm vs AoA**
         </td>
       </tr>
     </table>

   - **`plot_cl_cd_vs_alpha_all(...)`**  
     <div align="center">
       <img src="outputs/cl_cd_vs_alpha_all.jpg" width="60%" /><br>
       **Cl & Cd vs α for all airfoils**
     </div>

4. **Induction Factors** (`src/piwe_bem_mouxtin_ae/performance_curves.py`)  
   - `plot_spanwise_variables(r, a, a_prime, v0, rpm)`  
   <div align="center">
     <img src="outputs/spanwise_induction.jpg" width="50%" />
   </div>

5. **Operational Strategy** (`src/piwe_bem_mouxtin_ae/performance_curves.py`)  
   - `plot_operational_strategy(v0_array, pitch_array, rpm_array)`  
   <div align="center">
     <img src="outputs/operational_strategy.jpg" width="40%" />
   </div>

6. **Performance Curves** (`src/piwe_bem_mouxtin_ae/performance_curves.py`)  
   - `plot_performance_curves(v0_array, power_curve, thrust_curve, power_ref, thrust_ref)`  
   <div align="center">
     <img src="outputs/power_thrust_curves.jpg" width="60%" />
   </div>

---
## Team Collaboration

See [Collaboration.md](Collaboration.md) for detailed roles, workflow, and communication channels.

---

