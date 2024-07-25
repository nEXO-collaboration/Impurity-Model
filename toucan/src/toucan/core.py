import json
import numpy as np
from typing import List, Union, Optional, Dict, Tuple
from .utils import get_time_stamps, plot_data
from .constants import (
    IDEAL_GAS_MOLAR_VOLUME,
    AVOGADRO_NUMBER,
    BOLTZMANN_CONSTANT_EV,
    BOLTZMANN_CONSTANT_L,
    GXE_DENSITY,
)


class Outgassing_setup:
    """
    This class represents a system with various attributes and methods for outgassing analysis.

    Attributes:
    name (str): The name of the setup.
    material (str): The material of the system.
    solute (str): The solute in the system.
    version (str): The version of the system.
    temperatures (list): The temperatures of the system.
    time (list of list of floats): Timestamps of the system.
    diffusion (float): The diffusion constant of the solute in the material.
    solubility (float): The solubility of the solute in the material.
    activation_energy (float): The activation energy of the solute in the material.
    abundance (float): The abundance of the solute in the air.
    molar_mass (float): The molar mass of the solute.
    xe_mass (float): The mass of Xenon in the system.
    volume (float): The volume of the system.
    area (float): The area of the system.
    thickness (float): The thickness of the system.
    """

    def __init__(
        self,
        setup: str,
        material: Optional[str] = None,
        solute: Optional[str] = None,
        version: Optional[str] = None,
        data_file: str = "library.json",
    ):
        # Load data from the JSON file
        try:
            with open(data_file, "r") as file:
                data = json.load(file)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"The data file '{data_file}' was not found. Please ensure it's in the current working directory or provide the correct path."
            )

        self.name: str = setup
        self.material: Optional[str] = material
        self.solute: Optional[str] = solute
        self.version: Optional[str] = version

        self.temperatures: List[Union[int, float]] = []
        self.diffusion_constants: List[float] = []

        # Retrieve material properties safely from JSON data
        if material and solute:
            material_props = data["Material"].get(material, {}).get(solute, {})
            self.diffusion: Optional[float] = material_props.get("Diffusion Constant")
            self.solubility: Optional[float] = material_props.get("Solubility")
            self.activation_energy: Optional[float] = material_props.get(
                "Activation Energy"
            )

        # Retrieve system properties safely from JSON data
        if setup and material and version:
            system_props = (
                data["System"].get(setup, {}).get(material, {}).get(version, {})
            )
            self.volume: float = system_props.get("Volume")
            self.area: float = system_props.get("Area")
            self.thickness: float = system_props.get("Thickness")

        # Retrieve gas properties safely from JSON data
        if solute:
            gas_props = data["Gas"].get(solute, {})
            self.abundance: Optional[float] = gas_props.get("Abundance in Air")
            self.molar_mass: Optional[float] = gas_props.get("Molar Mass")

        # Retrieve Xenon Mass and Field Factor from JSON data
        self.xe_mass: Optional[float] = data["System"].get(setup, {}).get("Xenon Mass")
        self.field_factor: Optional[float] = (
            data["System"].get(setup, {}).get("Field Factor")
        )
        self.comment: Optional[str] = data["System"].get(setup, {}).get("Comment")

    def __str__(self) -> str:
        """
        Return a string representation of the system with attributes that are not None or empty.
        """
        attributes = vars(self)
        non_empty_attributes = {
            item: attributes[item]
            for item in attributes
            if attributes[item] not in [None, [], ""]
        }
        return "\n".join(
            f"{key}: {value}" for key, value in non_empty_attributes.items()
        )

    def get_diff_temp(self) -> None:
        """
        Calculate the diffusion constant at the instance's temperatures.
        """
        if self.diffusion is None or self.activation_energy is None:
            raise ValueError("Diffusion or activation energy not set.")

        self.diffusion_constants = [
            self.diffusion
            * np.exp(
                self.activation_energy
                / BOLTZMANN_CONSTANT_EV
                * ((1.0 / 293.15) - (1.0 / temp))
            )
            for temp in self.temperatures
        ]

    def get_initial_impurities(self, units: str) -> None:
        """
        Calculate the initial amount of impurities in various units.
        """
        if (
            self.volume is None
            or self.solubility is None
            or self.abundance is None
            or self.molar_mass is None
            or self.xe_mass is None
        ):
            raise ValueError(
                "System attributes not fully set for impurity calculation."
            )

        impurity_volume = self.volume * self.solubility * self.abundance
        impurity_mass = impurity_volume / IDEAL_GAS_MOLAR_VOLUME * self.molar_mass

        # Unit conversion and impurity_mass calculation
        if "pp" in units:
            conversion_factor = (
                1e6
                if units == "ppm"
                else 1e9 if units == "ppb" else 1e12 if units == "ppt" else 1
            )
            self.initial_impurities = (
                impurity_mass / self.xe_mass * conversion_factor
            )  # this is not correct
        elif units == "#":
            self.initial_impurities = (
                impurity_mass / self.molar_mass
            ) * AVOGADRO_NUMBER
        else:
            raise ValueError("Unsupported unit. Use 'ppm', 'ppb', 'ppt', or '#'.")

    def get_impurities_vs_time(self, time: List[List[float]]) -> List[List[float]]:
        """
        Calculate the impurities over time considering diffusion constants.
        """
        if not self.diffusion_constants:
            raise ValueError(
                "Diffusion constants are not initialized. Call get_diff_temp() first."
            )

        impurities = []
        # Use the same time segment for all temperatures if only one timestamp
        time_segments = (
            [time[0]] * len(self.diffusion_constants) if len(time) == 1 else time
        )

        for i, (diff_constant, time_segment) in enumerate(
            zip(self.diffusion_constants, time_segments)
        ):
            segment_impurities = []
            if (
                i == 0 or len(time) == 1
            ):  # For the first segment or if there's only one time segment
                current_impurity = self.initial_impurities
            else:
                # Start with the last impurity value of the previous segment
                current_impurity = impurities[i - 1][-1]

            for timestamp in time_segment:
                current_impurity = solve_diffusion_equation(
                    timestamp, diff_constant, self.thickness, current_impurity
                )
                segment_impurities.append(current_impurity)

            impurities.append(segment_impurities)

        return impurities

    def get_flow_rate_vs_time(
        self,
        time: List[List[float]],
        impurities: List[List[float]] = [],
        units: str = "#",
    ) -> List[List[float]]:
        """
        Calculate the flow rate over time considering diffusion constants and constraints.
        """
        if not impurities:
            raise ValueError("Impurities list is empty. Calculate impurities first.")

        initial_concentration = [x[0] / (self.volume * 1e3) for x in impurities]

        flow_rates = []
        time_segments = (
            [time[0]] * len(self.diffusion_constants) if len(time) == 1 else time
        )

        for diff_constant, temp, time_segment in zip(
            self.diffusion_constants, self.temperatures, time_segments
        ):
            segment_flow_rates = []
            for timestamp in time_segment:
                flow_rate = solve_flow_rate(
                    timestamp,
                    diff_constant,
                    self.thickness,
                    initial_concentration[0],
                    self.area,
                )

                if units == "mBar Liter":
                    flow_rate *= (BOLTZMANN_CONSTANT_L * temp) / AVOGADRO_NUMBER

                segment_flow_rates.append(flow_rate)

            flow_rates.append(segment_flow_rates)

        return flow_rates

    def get_steel_flow_rate_vs_pumping_time(
        self,
        time: List[List[float]],
        unbaked_flow_rate: float,
        initial_pumped_time: float,
    ) -> List[List[float]]:
        """
        Calculate the flow rate over time for steel considering pumping time.
        """
        if unbaked_flow_rate is None or initial_pumped_time is None:
            raise ValueError(
                "unbaked_flow_rate and initial_pumped_time must be provided."
            )
        flow_rates = []
        # Ensure time attribute is not empty or null
        if not time or time[0] is None:
            raise ValueError("Time attribute is not set for the system.")

        # Calculate flow rate for each time point
        flow_rate = [
            solve_steel_flow_rate_vs_pumping_time(
                unbaked_flow_rate, self.area, initial_pumped_time, t
            )
            for t in time[0]
        ]

        # Append calculated flow rate to the flow_rates attribute
        flow_rates.append(flow_rate)
        return flow_rates

    def get_electron_lifetime_vs_time(
        self,
        time: List[List[float]],
        initial_impurities: float,
        circulation_rate: float,
        purification_efficiency: float,
        out_diffusion: float,
        purifier_output: float = 0,
    ) -> Tuple[List[float], Dict[str, float]]:
        """
        Calculate the electron lifetime for each timestamp in the given list of timestamps.
        Returns a tuple with the list of electron lifetimes and a dictionary of parameters used.
        """
        if self.xe_mass is None or self.field_factor is None:
            raise ValueError("Xenon mass or Field Factor is not initialized.")

        electron_lifetimes = [
            solve_electron_lifetime(
                initial_impurities,
                circulation_rate,
                purification_efficiency,
                out_diffusion,
                purifier_output,
                timestamp,
                self.xe_mass / 1e3,
                self.field_factor,
            )
            for timestamp in time[0]
        ]

        # Dictionary of parameters used in the calculation
        params = {
            "initial_impurities": initial_impurities,
            "circulation_rate": circulation_rate,
            "purification_efficiency": purification_efficiency,
            "out_diffusion": out_diffusion,
            "purifier_output": purifier_output,
        }

        return electron_lifetimes, params


def solve_diffusion_equation(
    time: float, diff: float, thickness: float, conc: float
) -> float:
    """
    Solve the diffusion equation over a given time period.
    """
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
    """
    Calculate the flow rate of a substance through a medium.
    """
    terms = [
        np.exp(-((np.pi * (2.0 * n + 1.0) / thickness) ** 2) * diff * time)
        * (4.0 * conc * diff)
        / thickness
        for n in range(1000)
    ]
    return sum(terms) * area


def solve_steel_flow_rate_vs_pumping_time(
    unbaked_flow_rate: float,
    area: float,
    initial_pumped_time: float,
    time: float,
) -> float:
    """
    Calculate the steel flow rate for a given time.
    """
    flow_rate = unbaked_flow_rate * area * initial_pumped_time / time
    return flow_rate


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
    Calculate the electron lifetime for a given timestamp.
    """
    factor_exp = np.exp(
        -GXE_DENSITY * purification_efficiency * circulation_rate * timestamp / xe_mass
    )
    denominator = initial_impurities * factor_exp + (
        (out_diffusion + purifier_output * circulation_rate)
        / (purification_efficiency * circulation_rate)
    ) * (1 - factor_exp)

    if denominator == 0:
        return float("inf")
    else:
        return field_factor / denominator
