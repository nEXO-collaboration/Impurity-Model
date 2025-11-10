"""Comparison of theoretical and empirical outgassing / desorption models."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, eV

# Physical constants
R = 8.314462618 * 10  # Gas constant in mbarL/mol/K
NA = 6.02214076e23  # Avogadro's number
TAU_0 = 1e-13  # Vibrational frequency in s

# Default model parameters
DEFAULT_PARAMS = {
    # Diffusion parameters
    "c_0": 1e15,  # Initial concentration
    "D": 1e-9,  # Diffusion coefficient
    "d": 0.1,  # Material thickness
    "area": 100,  # Surface area
    # Desorption parameters
    "T": 293.15,  # Temperature in Kelvin
    "n_energies": 100,  # Number of energy levels
    "E_min": 0.8,  # Minimum energy in eV
    "E_max": 1.2,  # Maximum energy in eV
    # Empirical model parameters
    "alpha": 1.0,  # Power law decay parameter
    "t_0": 3600,  # Reference time (1 hour in seconds)
}

# Numerical controls
DIFFUSION_TERMS = 200  # Number of Fourier terms kept in the diffusion series


def build_desorption_terms(model_params, energies_ev=None):
    """Pre-compute prefactors and decay constants for the desorption model."""
    temperature = model_params["T"]
    if energies_ev is None:
        energies_ev = np.linspace(
            model_params["E_min"], model_params["E_max"], model_params["n_energies"]
        )
    energies_ev = np.asarray(energies_ev, dtype=float)
    energies_joules = energies_ev * eV

    per_level_concentration = model_params["c_0"] / model_params["n_energies"]
    boltzmann = np.exp(-energies_joules / (k * temperature))
    prefactor = per_level_concentration * R * temperature / (NA * TAU_0) * boltzmann
    lifetimes = TAU_0 * np.exp(energies_joules / (k * temperature))
    return energies_ev, prefactor, lifetimes


def _evaluate_desorption_components(time_values, prefactor, lifetimes):
    """Return the per-energy contributions (without the surface-area factor)."""
    times = np.asarray(time_values, dtype=float)
    scalar_input = times.ndim == 0
    if scalar_input:
        times = times.reshape(1)

    decay = np.exp(-np.outer(1.0 / lifetimes, times))
    return prefactor[:, None] * decay, scalar_input


def theoretical_diffusion(time_values, model_params, n_terms=DIFFUSION_TERMS):
    """Calculate the diffusion term using a configurable number of Fourier modes."""
    surface_concentration = model_params["c_0"]
    diffusion_coeff = model_params["D"]
    thickness = model_params["d"]

    times = np.asarray(time_values, dtype=float)
    scalar_input = times.ndim == 0
    if scalar_input:
        times = times.reshape(1)

    mode_indices = np.arange(n_terms)
    eigenvalues = ((np.pi * (2 * mode_indices + 1) / thickness) ** 2) * diffusion_coeff
    decay = np.exp(-np.outer(eigenvalues, times))
    diffusion = (4 * surface_concentration * diffusion_coeff / thickness) * decay.sum(
        axis=0
    )

    if scalar_input:
        return float(diffusion[0])
    return diffusion


def theoretical_desorption(
    time_values, model_params, prefactor=None, lifetimes=None
):
    """First-order desorption rate with Arrhenius kinetics.

    Rate constant: k_d = τ₀⁻¹ exp(-E/kT) ⇒ coverage ∝ exp(-t / (τ₀ exp(E/kT)))
    """
    if prefactor is None or lifetimes is None:
        _, prefactor, lifetimes = build_desorption_terms(model_params)

    components, scalar_input = _evaluate_desorption_components(
        time_values, prefactor, lifetimes
    )
    total = model_params["area"] * components.sum(axis=0)

    if scalar_input:
        return float(total[0])
    return total


def empirical_diffusion(time_values, model_params):
    """Calculate diffusion term of the empirical three-regime model."""
    surface_concentration = model_params["c_0"]
    diffusion_coeff = model_params["D"]
    thickness = model_params["d"]
    diffusion_time_constant = thickness**2 / (np.pi**2 * diffusion_coeff)

    diffusion = np.zeros_like(time_values)
    transition_short = 0.5 * diffusion_time_constant
    transition_long = 3.0 * diffusion_time_constant

    mask_short = time_values <= transition_short
    mask_mid = (time_values > transition_short) & (time_values <= transition_long)
    mask_long = time_values > transition_long

    diffusion[mask_short] = surface_concentration * np.sqrt(
        diffusion_coeff / (np.pi * time_values[mask_short])
    )
    diffusion[mask_mid] = (
        4 * surface_concentration * diffusion_coeff / thickness
        * np.exp(-time_values[mask_mid] / diffusion_time_constant)
    )

    value_at_transition = (
        4 * surface_concentration * diffusion_coeff / thickness
        * np.exp(-transition_long / diffusion_time_constant)
    )
    scale_long = value_at_transition * (transition_long**3)
    diffusion[mask_long] = scale_long / time_values[mask_long] ** 3

    return diffusion


def empirical_desorption(time_values, model_params, prefactor=None, lifetimes=None):
    """Empirical 1/t^alpha model matched to the theoretical rate at t=t_0."""
    reference_time = model_params["t_0"]
    alpha = model_params["alpha"]

    reference = theoretical_desorption(reference_time, model_params, prefactor, lifetimes)

    times = np.asarray(time_values, dtype=float)
    scalar_input = times.ndim == 0
    if scalar_input:
        times = times.reshape(1)

    empirical = reference * (reference_time / times) ** alpha
    if scalar_input:
        return float(empirical[0])
    return empirical


def diffusion_asymptotic(time_values, model_params, regime="short"):
    """Calculate asymptotic diffusion behavior for short/long times."""
    surface_concentration = model_params["c_0"]
    diffusion_coeff = model_params["D"]
    thickness = model_params["d"]
    times = np.asarray(time_values, dtype=float)

    if regime == "short":
        return surface_concentration * np.sqrt(diffusion_coeff / (np.pi * times))
    return (
        4 * surface_concentration * diffusion_coeff / thickness
        * np.exp(-((np.pi / thickness) ** 2) * diffusion_coeff * times)
    )


def compute_model_components(time_grid, model_params=None):
    """Build the diffusion/desorption components for plotting."""
    if model_params is None:
        model_params = DEFAULT_PARAMS

    _, prefactor, lifetimes = build_desorption_terms(model_params)
    theo_diff = theoretical_diffusion(
        time_grid, model_params, n_terms=DIFFUSION_TERMS
    )
    theo_des = theoretical_desorption(
        time_grid, model_params, prefactor, lifetimes
    )
    emp_diff = empirical_diffusion(time_grid, model_params)
    emp_des = empirical_desorption(time_grid, model_params, prefactor, lifetimes)

    critical_time = model_params["d"] ** 2 / (6 * model_params["D"])
    diffusion_tau = model_params["d"] ** 2 / (np.pi**2 * model_params["D"])
    short_mask = time_grid <= 0.1 * critical_time
    long_mask = time_grid >= 10 * critical_time

    diffusion_short = np.zeros_like(time_grid)
    diffusion_long = np.zeros_like(time_grid)
    diffusion_short[short_mask] = diffusion_asymptotic(
        time_grid[short_mask], model_params, "short"
    )
    diffusion_long[long_mask] = diffusion_asymptotic(
        time_grid[long_mask], model_params, "long"
    )

    return {
        "prefactor": prefactor,
        "lifetimes": lifetimes,
        "theo_diff": theo_diff,
        "theo_des": theo_des,
        "emp_diff": emp_diff,
        "emp_des": emp_des,
        "critical_time": critical_time,
        "diffusion_tau": diffusion_tau,
        "short_mask": short_mask,
        "long_mask": long_mask,
        "diffusion_short": diffusion_short,
        "diffusion_long": diffusion_long,
    }


def plot_diffusion_panel(axes_handle, time_grid, components):
    """Draw the diffusion comparison subplot."""
    theo_diff = components["theo_diff"]
    emp_diff = components["emp_diff"]
    diffusion_short = components["diffusion_short"]
    diffusion_long = components["diffusion_long"]
    short_mask = components["short_mask"]
    long_mask = components["long_mask"]
    diffusion_tau = components["diffusion_tau"]
    critical_time = components["critical_time"]

    axes_handle.loglog(
        time_grid / 3600,
        theo_diff,
        label="Theoretical",
        linestyle="-",
        linewidth=4,
        alpha=0.5,
    )
    axes_handle.loglog(
        time_grid / 3600,
        emp_diff,
        label="Empirical",
        linestyle="--",
        linewidth=4,
        alpha=0.5,
    )
    axes_handle.loglog(
        time_grid[short_mask] / 3600,
        diffusion_short[short_mask],
        label="Short-time approx.",
        linestyle=":",
        color="red",
        linewidth=4,
        alpha=0.6,
    )
    axes_handle.loglog(
        time_grid[long_mask] / 3600,
        diffusion_long[long_mask],
        label="Long-time approx.",
        linestyle=":",
        color="green",
        linewidth=4,
        alpha=0.6,
    )
    axes_handle.grid(True, which="both", ls="-", alpha=0.2)
    axes_handle.set_ylabel("Diffusion Rate (mbar*L/s)")
    axes_handle.set_title("Diffusion Components")
    axes_handle.axvline(
        x=0.5 * diffusion_tau / 3600,
        color="gray",
        linestyle=":",
        alpha=0.7,
        label="τ/2",
        linewidth=2,
    )
    axes_handle.axvline(
        x=3.0 * diffusion_tau / 3600,
        color="gray",
        linestyle=":",
        alpha=0.7,
        label="3τ",
        linewidth=2,
    )
    axes_handle.axvline(
        x=critical_time / 3600,
        color="gray",
        linestyle="-",
        alpha=0.7,
        label="t_c",
        linewidth=2,
    )
    axes_handle.legend()


def plot_desorption_panel(axes_handle, time_grid, components, model_params=None):
    """Draw the desorption subplot along with sample energy channels."""
    if model_params is None:
        model_params = DEFAULT_PARAMS

    theo_des = components["theo_des"]
    emp_des = components["emp_des"]

    axes_handle.loglog(
        time_grid / 3600,
        theo_des,
        label="Theoretical (sum)",
        linestyle="-",
        linewidth=4,
        alpha=0.5,
    )
    axes_handle.loglog(
        time_grid / 3600,
        emp_des,
        label="Empirical",
        linestyle="--",
        linewidth=4,
        alpha=0.5,
    )

    energies_demo_ev = np.array([0.85, 0.95, 1.05])
    _, demo_prefactor, demo_lifetimes = build_desorption_terms(
        model_params, energies_ev=energies_demo_ev
    )
    demo_components, _ = _evaluate_desorption_components(
        time_grid, demo_prefactor, demo_lifetimes
    )
    for energy, component in zip(energies_demo_ev, demo_components):
        axes_handle.loglog(
            time_grid / 3600,
            model_params["area"] * component,
            ":",
            color="0.4",
            linewidth=2,
            label=f"E = {energy:.2f} eV",
        )

    axes_handle.set_ylim(1e-20, 1e10)
    axes_handle.grid(True, which="both", ls="-", alpha=0.2)
    axes_handle.set_ylabel("Desorption Rate (mbar*L/s)")
    axes_handle.set_title("Desorption Components")
    axes_handle.legend()


def plot_total_panel(axes_handle, time_grid, components):
    """Draw the total outgassing subplot with component markers."""
    theo_diff = components["theo_diff"]
    theo_des = components["theo_des"]
    emp_diff = components["emp_diff"]
    emp_des = components["emp_des"]

    axes_handle.loglog(
        time_grid / 3600,
        theo_diff + theo_des,
        label="Theoretical Total",
        linestyle="-",
        linewidth=4,
        alpha=0.5,
    )
    axes_handle.loglog(
        time_grid / 3600,
        emp_diff + emp_des,
        label="Empirical Total",
        linestyle="--",
        linewidth=4,
        color="orange",
        alpha=0.5,
    )
    axes_handle.loglog(
        time_grid / 3600,
        theo_diff,
        label="Theoretical Diffusion",
        linestyle="None",
        marker="x",
        markevery=0.03,
        linewidth=2,
        color="blue",
        alpha=1,
    )
    axes_handle.loglog(
        time_grid / 3600,
        theo_des,
        label="Theoretical Desorption",
        marker="+",
        markevery=0.03,
        linewidth=2,
        color="blue",
        alpha=1,
        linestyle="None",
    )
    axes_handle.loglog(
        time_grid / 3600,
        emp_diff,
        label="Empirical Diffusion",
        marker="x",
        markevery=0.03,
        linewidth=2,
        color="orange",
        alpha=1,
        linestyle="None",
    )
    axes_handle.loglog(
        time_grid / 3600,
        emp_des,
        label="Empirical Desorption",
        marker="+",
        markevery=0.03,
        linewidth=2,
        color="orange",
        alpha=1,
        linestyle="None",
    )

    axes_handle.grid(True, which="both", ls="-", alpha=0.2)
    axes_handle.set_xlabel("Time (hours)")
    axes_handle.set_ylabel("Total Outgassing Rate (mbar*L/s)")
    axes_handle.set_title("Total Outgassing")
    axes_handle.legend()


def main():
    """Entry point for plotting the model comparison."""
    time_grid = np.logspace(2, 8, 1000)
    components = compute_model_components(time_grid, DEFAULT_PARAMS)

    fig = plt.figure(figsize=(9, 7))
    grid_spec = plt.GridSpec(2, 2, height_ratios=[1, 1.2])

    plot_diffusion_panel(fig.add_subplot(grid_spec[0, 0]), time_grid, components)
    plot_desorption_panel(
        fig.add_subplot(grid_spec[0, 1]), time_grid, components, DEFAULT_PARAMS
    )
    plot_total_panel(fig.add_subplot(grid_spec[1, :]), time_grid, components)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
