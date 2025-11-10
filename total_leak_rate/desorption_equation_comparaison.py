import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, h


def theoretical_desorption(t, T, area, E_range):
    """Calculate desorption term with specified energy range"""
    E_min, E_max = E_range  # eV
    n_energies = 100
    E_values = np.linspace(E_min, E_max, n_energies) * 1.6e-19  # Convert to Joules

    # Pre-exponential factors
    A_factors = 1e13 * np.exp(-(E_values - E_values[0]) / (k * T))

    nu_0 = k * T / h
    desorption = np.zeros_like(t, dtype=np.float64)

    for A_i, E_i in zip(A_factors, E_values):
        desorption += (
            A_i * np.exp(-E_i / (k * T)) * np.exp(-nu_0 * t * np.exp(-E_i / (k * T)))
        )

    return desorption * area


def empirical_desorption(t, J_0, area):
    """Simple 1/t empirical model"""
    t_0 = 3600  # 1 hour in seconds
    return J_0 * area * (t_0 / t)


# Parameters
T = 293.15  # Temperature in Kelvin
area = 100  # Surface area
t = np.logspace(2, 8, 1000)  # Time points from 100s to 10^8s

# Different energy ranges to compare
energy_ranges = [
    (0.5, 0.7),  # Low energy range
    (0.8, 1.2),  # Medium energy range
    (1.3, 1.5),  # High energy range
]

plt.figure(figsize=(10, 6))

# Plot empirical model for reference
J_0 = 1e-6  # Reference flux
emp_des = empirical_desorption(t, J_0, area)
plt.loglog(t / 3600, emp_des, "k--", label="Empirical (1/t)", linewidth=2)

# Plot theoretical models for different energy ranges
colors = ["blue", "green", "red"]
for (E_min, E_max), color in zip(energy_ranges, colors):
    # Calculate theoretical desorption
    theo_des = theoretical_desorption(t, T, area, (E_min, E_max))
    label = f"Theoretical ({E_min}-{E_max} eV)"
    plt.loglog(t / 3600, theo_des, color=color, label=label, linewidth=2)

    # Plot individual energy contributions
    E_demo = np.array([E_min + 0.1, (E_min + E_max) / 2, E_max - 0.1]) * 1.6e-19
    A_demo = 1e13 * np.exp(-(E_demo - E_demo[0]) / (k * T))

    for E_i, A_i in zip(E_demo, A_demo):
        des_i = (
            A_i
            * np.exp(-E_i / (k * T))
            * np.exp(-(k * T / h) * t * np.exp(-E_i / (k * T)))
        )
        label = f"E = {E_i/1.6e-19:.2f} eV"
        plt.loglog(t / 3600, des_i * area, ":", color=color, alpha=0.3, linewidth=1.5)

plt.grid(True, which="both", ls="-", alpha=0.2)
plt.xlabel("Time (hours)")
plt.ylabel("Desorption Rate (mbar*L/s)")
plt.title("Desorption Rates for Different Activation Energy Ranges")
plt.ylim(1e-20, 1e10)
plt.legend()
plt.tight_layout()
plt.show()

# Calculate and print slopes at different time regions
print("\nSlope Analysis at Different Times:")
print("Energy Range | Early Slope | Late Slope")
print("-" * 45)

for E_min, E_max in energy_ranges:
    theo_des = theoretical_desorption(t, T, area, (E_min, E_max))

    # Calculate slopes in early and late time regions
    early_idx = slice(100, 200)
    late_idx = slice(-200, -100)

    early_slope = np.polyfit(np.log10(t[early_idx]), np.log10(theo_des[early_idx]), 1)[
        0
    ]
    late_slope = np.polyfit(np.log10(t[late_idx]), np.log10(theo_des[late_idx]), 1)[0]

    print(f"{E_min}-{E_max} eV  | {early_slope:.2f}      | {late_slope:.2f}")
