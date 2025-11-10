import numpy as np
from .constants import (
    BOLTZMANN_CONSTANT_EV,
    IDEAL_GAS_MOLAR_VOLUME,
    AVOGADRO_NUMBER,
    GXE_DENSITY,
)


def get_diff_temp(setup) -> None:
    """
    Calculate temperature-dependent diffusion constants for the setup.

    Args:
        setup: The outgassing setup object.

    Raises:
        ValueError: If diffusion or activation energy is not set.
    """
    if setup.diffusion is None or setup.activation_energy is None:
        raise ValueError("Diffusion or activation energy not set.")

    setup.diffusion_constants = [
        setup.diffusion
        * np.exp(
            setup.activation_energy
            / BOLTZMANN_CONSTANT_EV
            * ((1.0 / 293.15) - (1.0 / temp))
        )
        for temp in setup.temperatures
    ]


def get_initial_impurities(setup, units: str) -> None:
    """
    Calculate initial impurities for the setup based on the given units.

    Args:
        setup: The outgassing setup object.
        units: The units for impurity calculation ('ppm', 'ppb', 'ppt', or '#').

    Raises:
        ValueError: If required attributes are not set or if units are unsupported.
    """
    required_attributes = ["volume", "solubility", "abundance", "molar_mass", "xe_mass"]
    if any(getattr(setup, attr) is None for attr in required_attributes):
        raise ValueError("System attributes not fully set for impurity calculation.")

    impurity_volume = setup.volume * setup.solubility * setup.abundance
    impurity_mass = impurity_volume / IDEAL_GAS_MOLAR_VOLUME * setup.molar_mass

    if units.startswith("pp"):
        conversion_factors = {"ppm": 1e6, "ppb": 1e9, "ppt": 1e12}
        if units not in conversion_factors:
            raise ValueError("Unsupported unit. Use 'ppm', 'ppb', 'ppt', or '#'.")
        setup.initial_impurities = (
            impurity_mass / setup.xe_mass * conversion_factors[units]
        )
    elif units == "#":
        setup.initial_impurities = (impurity_mass / setup.molar_mass) * AVOGADRO_NUMBER
    else:
        raise ValueError("Unsupported unit. Use 'ppm', 'ppb', 'ppt', or '#'.")


def solve_diffusion_equation(
    time: float, diff: float, thickness: float, conc: float, n_terms: int = 1000
) -> float:
    """
    Solve the diffusion equation for a given time, diffusion constant, thickness, and concentration.

    Args:
        time: Time in seconds.
        diff: Diffusion constant.
        thickness: Material thickness.
        conc: Initial concentration.
        n_terms: Number of terms to use in the series expansion (default: 1000).

    Returns:
        float: The solution to the diffusion equation.
    """
    n_range = np.arange(n_terms)
    terms = (
        (1.0 / ((2.0 * n_range + 1.0) ** 2))
        * np.exp(-((np.pi * (2.0 * n_range + 1.0) / thickness) ** 2) * diff * time)
        * (conc * 8.0 * thickness)
        / (np.pi**2 * 2.0)
    )
    return np.sum(terms)


def solve_flow_rate(
    time: float,
    diff: float,
    thickness: float,
    conc: float,
    area: float,
    n_terms: int = 1000,
) -> float:
    """
    Calculate the flow rate for a given time, diffusion constant, thickness, concentration, and area.

    Args:
        time: Time in seconds.
        diff: Diffusion constant.
        thickness: Material thickness.
        conc: Initial concentration.
        area: Surface area.
        n_terms: Number of terms to use in the series expansion (default: 1000).

    Returns:
        float: The calculated flow rate.
    """
    n_range = np.arange(n_terms)
    terms = np.exp(-((np.pi * (2.0 * n_range + 1.0) / thickness) ** 2) * diff * time)
    return np.sum(terms) * (4.0 * conc * diff / thickness) * area


def solve_steel_flow_rate_vs_pumping_time(
    unbaked_flow_rate: float, area: float, initial_pumped_time: float, time: float
) -> float:
    """
    Calculate the steel flow rate versus pumping time.

    Args:
        unbaked_flow_rate: Unbaked flow rate.
        area: Surface area.
        initial_pumped_time: Initial pumped time in seconds.
        time: Current time in seconds.

    Returns:
        float: The calculated steel flow rate.
    """
    return unbaked_flow_rate * area * initial_pumped_time / time


def solve_electron_lifetime(
    initial_impurities: float,
    circulation_rate: float,
    purification_efficiency: float,
    out_diffusion: float,
    purifier_output: float,
    timestamp: float,
    xe_mass: float,
    field_factor: float,
) -> float:
    """
    Calculate the electron lifetime for given parameters.

    Args:
        initial_impurities: Initial impurity concentration.
        circulation_rate: Circulation rate.
        purification_efficiency: Purification efficiency.
        out_diffusion: Out-diffusion rate.
        purifier_output: Purifier output.
        timestamp: Current time.
        xe_mass: Xenon mass.
        field_factor: Field factor.

    Returns:
        float: The calculated electron lifetime.
    """
    factor_exp = np.exp(
        -GXE_DENSITY * purification_efficiency * circulation_rate * timestamp / xe_mass
    )
    denominator = initial_impurities * factor_exp + (
        (out_diffusion + purifier_output * circulation_rate)
        / (purification_efficiency * circulation_rate)
    ) * (1 - factor_exp)

    return float("inf") if denominator == 0 else field_factor / denominator
