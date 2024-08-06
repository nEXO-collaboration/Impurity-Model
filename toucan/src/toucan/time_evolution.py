from typing import List, Dict, Tuple
from .calculations import (
    solve_diffusion_equation,
    solve_flow_rate,
    solve_steel_flow_rate_vs_pumping_time,
    solve_electron_lifetime,
)
from .constants import BOLTZMANN_CONSTANT_L, AVOGADRO_NUMBER


class TimeEvolution:
    @staticmethod
    def get_impurities_vs_time(setup, time: List[List[float]]) -> List[List[float]]:
        if not setup.diffusion_constants:
            raise ValueError(
                "Diffusion constants are not initialized. Call get_diff_temp() first."
            )

        impurities = []
        time_segments = (
            [time[0]] * len(setup.diffusion_constants) if len(time) == 1 else time
        )

        for i, (diff_constant, time_segment) in enumerate(
            zip(setup.diffusion_constants, time_segments)
        ):
            segment_impurities = []
            current_impurity = (
                setup.initial_impurities
                if i == 0 or len(time) == 1
                else impurities[i - 1][-1]
            )

            for timestamp in time_segment:
                current_impurity = solve_diffusion_equation(
                    timestamp, diff_constant, setup.thickness, current_impurity
                )
                segment_impurities.append(current_impurity)

            impurities.append(segment_impurities)

        return impurities

    @staticmethod
    def get_flow_rate_vs_time(
        setup,
        time: List[List[float]],
        impurities: List[List[float]] = [],
        units: str = "#",
    ) -> List[List[float]]:
        if not impurities:
            raise ValueError("Impurities list is empty. Calculate impurities first.")

        initial_concentration = [x[0] / (setup.volume * 1e3) for x in impurities]

        flow_rates = []
        time_segments = (
            [time[0]] * len(setup.diffusion_constants) if len(time) == 1 else time
        )

        for diff_constant, temp, time_segment in zip(
            setup.diffusion_constants, setup.temperatures, time_segments
        ):
            segment_flow_rates = []
            for timestamp in time_segment:
                flow_rate = solve_flow_rate(
                    timestamp,
                    diff_constant,
                    setup.thickness,
                    initial_concentration[0],
                    setup.area,
                )

                if units == "mBar Liter":
                    flow_rate *= (BOLTZMANN_CONSTANT_L * temp) / AVOGADRO_NUMBER

                segment_flow_rates.append(flow_rate)

            flow_rates.append(segment_flow_rates)

        return flow_rates

    @staticmethod
    def get_steel_flow_rate_vs_pumping_time(
        setup,
        time: List[List[float]],
        unbaked_flow_rate: float,
        initial_pumped_time: float,
    ) -> List[List[float]]:
        if unbaked_flow_rate is None or initial_pumped_time is None:
            raise ValueError(
                "unbaked_flow_rate and initial_pumped_time must be provided."
            )
        if not time or time[0] is None:
            raise ValueError("Time attribute is not set for the system.")

        flow_rate = [
            solve_steel_flow_rate_vs_pumping_time(
                unbaked_flow_rate, setup.area, initial_pumped_time, t
            )
            for t in time[0]
        ]

        return [flow_rate]

    @staticmethod
    def get_electron_lifetime_vs_time(
        setup,
        time: List[List[float]],
        initial_impurities: float,
        circulation_rate: float,
        purification_efficiency: float,
        out_diffusion: float,
        purifier_output: float = 0,
    ) -> Tuple[List[float], Dict[str, float]]:
        if setup.xe_mass is None or setup.field_factor is None:
            raise ValueError("Xenon mass or Field Factor is not initialized.")

        electron_lifetimes = [
            solve_electron_lifetime(
                initial_impurities,
                circulation_rate,
                purification_efficiency,
                out_diffusion,
                purifier_output,
                timestamp,
                setup.xe_mass / 1e3,
                setup.field_factor,
            )
            for timestamp in time[0]
        ]

        params = {
            "initial_impurities": initial_impurities,
            "circulation_rate": circulation_rate,
            "purification_efficiency": purification_efficiency,
            "out_diffusion": out_diffusion,
            "purifier_output": purifier_output,
        }

        return electron_lifetimes, params
