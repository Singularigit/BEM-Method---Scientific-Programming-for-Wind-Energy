"""
Wind Turbine classes for modeling power output based on wind speed.
"""

import numpy as np
import matplotlib.pyplot as plt


class GeneralWindTurbine:
    """
    A general wind turbine class that calculates power output based on theoretical power curve.
    """

    def __init__(self, rotor_diameter, hub_height, rated_power, v_in, v_rated, v_out, name=None):
        """
        Initialize a general wind turbine with theoretical power curve.
        """
        self.rotor_diameter = rotor_diameter
        self.hub_height = hub_height
        self.rated_power = rated_power
        self.v_in = v_in
        self.v_rated = v_rated
        self.v_out = v_out
        self.name = name

    def get_power(self, wind_speed):
        """
        Calculate power output for a given wind speed.
        """
        if isinstance(wind_speed, (list, np.ndarray)):
            return np.array([self._calc_power(v) for v in wind_speed])
        return self._calc_power(wind_speed)

    def _calc_power(self, wind_speed):
        """
        Helper method to calculate power for a single wind speed value.
        """
        if wind_speed < self.v_in or wind_speed > self.v_out:
            return 0.0
        if wind_speed < self.v_rated:
            return self.rated_power * (wind_speed / self.v_rated) ** 3
        return self.rated_power

    def plot_power_curve(self, v_range=None, label=None):
        """
        Plot the power curve for this turbine.
        """
        if v_range is None:
            v_range = np.linspace(0, self.v_out + 2, 100)

        power = self.get_power(v_range)

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(v_range, power, label=label or self.name or "Theoretical")
        ax.set_xlabel("Wind Speed (m/s)")
        ax.set_ylabel("Power (kW)")
        ax.set_title("Wind Turbine Power Curve")
        ax.grid(True)

        if label or self.name:
            ax.legend()

        return fig, ax


class WindTurbine(GeneralWindTurbine):
    """
    A wind turbine class that uses actual power curve data for calculating power output.
    """

    def __init__(self, rotor_diameter, hub_height, rated_power, v_in, v_rated, v_out,
                 power_curve_data, name=None):
        """
        Initialize a wind turbine with actual power curve data.
        """
        super().__init__(rotor_diameter, hub_height, rated_power, v_in, v_rated, v_out, name)
        self.power_curve_data = np.array(power_curve_data)

    def get_power(self, wind_speed):
        """
        Calculate power output using interpolation on the power curve data.
        """
        wind_speeds = self.power_curve_data[:, 0]
        power_values = self.power_curve_data[:, 1]

        v_min = np.min(wind_speeds)
        v_max = np.max(wind_speeds)

        if isinstance(wind_speed, (list, np.ndarray)):
            wind_speed = np.array(wind_speed)
            result = np.zeros_like(wind_speed, dtype=float)
            mask = (wind_speed >= v_min) & (wind_speed <= v_max)
            result[mask] = np.interp(wind_speed[mask], wind_speeds, power_values)
            return result

        if wind_speed < v_min or wind_speed > v_max:
            return 0.0
        return np.interp(wind_speed, wind_speeds, power_values)
