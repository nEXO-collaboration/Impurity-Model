import numpy as np
from .constants import (
    BOLTZMANN_CONSTANT_EV,
    IDEAL_GAS_MOLAR_VOLUME,
    AVOGADRO_NUMBER,
    GXE_DENSITY,
)


def get_diff_temp(setup) -> None:
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
    if (
        setup.volume is None
        or setup.solubility is None
        or setup.abundance is None
        or setup.molar_mass is None
        or setup.xe_mass is None
    ):
        raise ValueError("System attributes not fully set for impurity calculation.")

    impurity_volume = setup.volume * setup.solubility * setup.abundance
    impurity_mass = impurity_volume / IDEAL_GAS_MOLAR_VOLUME * setup.molar_mass

    if "pp" in units:
        conversion_factor = (
            1e6
            if units == "ppm"
            else 1e9 if units == "ppb" else 1e12 if units == "ppt" else 1
        )
        setup.initial_impurities = impurity_mass / setup.xe_mass * conversion_factor
    elif units == "#":
        setup.initial_impurities = (impurity_mass / setup.molar_mass) * AVOGADRO_NUMBER
    else:
        raise ValueError("Unsupported unit. Use 'ppm', 'ppb', 'ppt', or '#'.")


def solve_diffusion_equation(
    time: float, diff: float, thickness: float, conc: float
) -> float:
    terms = [
        (1.0 / ((2.0 * n + 1.0) ** 2))
        * np.exp(-((np.pi * (2.0 * n + 1.0) / thickness) ** 2) * diff * time)
        * (conc * 8.0 * thickness)
        / (np.pi**2 * 2.0)
        for n in range(1000)
    ]
    return sum(terms)


def solve_flow_rate(
    time: float, diff: float, thickness: float, conc: float, area: float
) -> float:
    terms = [
        np.exp(-((np.pi * (2.0 * n + 1.0) / thickness) ** 2) * diff * time)
        * (4.0 * conc * diff)
        / thickness
        for n in range(1000)
    ]
    return sum(terms) * area


def solve_steel_flow_rate_vs_pumping_time(
    unbaked_flow_rate: float, area: float, initial_pumped_time: float, time: float
) -> float:
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
    factor_exp = np.exp(
        -GXE_DENSITY * purification_efficiency * circulation_rate * timestamp / xe_mass
    )
    denominator = initial_impurities * factor_exp + (
        (out_diffusion + purifier_output * circulation_rate)
        / (purification_efficiency * circulation_rate)
    ) * (1 - factor_exp)

    return float("inf") if denominator == 0 else field_factor / denominator
