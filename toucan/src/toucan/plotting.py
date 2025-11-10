import matplotlib.pyplot as plt
from typing import List, Tuple, Optional, Dict, Union
from .utils import get_time_stamps
from .calculations import get_diff_temp, get_initial_impurities
from .time_evolution import TimeEvolution
from .setup import OutgassingSetup
import numpy as np


def plot_metals(
    setups: Union[OutgassingSetup, List[OutgassingSetup]],
    time_range: Tuple[float, float] = (1e-1, 1e4),
    time_scale: str = "Hours",
    rate_unit: str = "mbar·L/s·cm²",
    unbaked_flow_rates: Union[float, List[float]] = 1e-8,
    initial_pumped_times: Union[float, List[float]] = 2600,
    fig_size: Tuple[int, int] = (10, 6),
    y_scale: str = "log",
    x_scale: str = "log",
) -> plt.Figure:
    """
    Plot the outgassing rate for one or multiple metal setups. Works best for the first 10h of pumping.

    Args:
        setups (Union[OutgassingSetup, List[OutgassingSetup]]): The setup(s) to plot.
        time_range (Tuple[float, float]): The time range to plot (start, end).
        time_scale (str): The time scale to use ("Hours", "Days", or "Seconds").
        rate_unit (str): The unit for the outgassing rate.
        unbaked_flow_rates (Union[float, List[float]]): The unbaked flow rate(s).
        initial_pumped_times (Union[float, List[float]]): The initial pumped time(s) in seconds.
        fig_size (Tuple[int, int]): The figure size (width, height).
        y_scale (str): Scale for y-axis ("log" or "linear").
        x_scale (str): Scale for x-axis ("log" or "linear").

    Returns:
        plt.Figure: The generated matplotlib figure.
    """
    # Convert single setup to list for consistent processing
    if not isinstance(setups, list):
        setups = [setups]

    # Ensure unbaked_flow_rates and initial_pumped_times are lists
    if not isinstance(unbaked_flow_rates, list):
        unbaked_flow_rates = [unbaked_flow_rates] * len(setups)
    if not isinstance(initial_pumped_times, list):
        initial_pumped_times = [initial_pumped_times] * len(setups)

    # Create figure and axis
    fig, ax = plt.subplots(figsize=fig_size)

    # Convert time range to seconds
    conversion_factor = {"Seconds": 1, "Hours": 3600, "Days": 86400}[time_scale]
    time_range_seconds = (
        time_range[0] * conversion_factor,
        time_range[1] * conversion_factor,
    )

    # Generate time array
    time_array = np.logspace(
        np.log10(time_range_seconds[0]), np.log10(time_range_seconds[1]), num=1000
    )

    # Plot for each setup
    for setup, unbaked_flow_rate, initial_pumped_time in zip(
        setups, unbaked_flow_rates, initial_pumped_times
    ):
        # Calculate outgassing rate
        outgassing_rate = (
            unbaked_flow_rate * setup.area * initial_pumped_time
        ) / time_array

        # Plot
        ax.plot(
            time_array / conversion_factor,
            outgassing_rate / setup.area,
            label=f"{setup.version}, Area={setup.area:.2f} cm²",
        )

    # Set labels and scales
    ax.set_xlabel(f"Time ({time_scale})")
    ax.set_ylabel(f"Outgassing Rate ({rate_unit})")
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)

    # Set grid and legend
    ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.legend()

    # Set title
    plt.title("Metal Outgassing Rate vs Time")

    # Adjust layout
    fig.tight_layout()

    return fig


def plot_polymers(
    setups: Union[List[OutgassingSetup], OutgassingSetup],
    time_range: Union[Tuple[float, float], Dict[float, Tuple[float, float]]],
    time_scale: str = "Seconds",
    plot_type: str = "impurities",
    fig_size: Tuple[int, int] = (10, 6),
    temperatures: Optional[List[float]] = None,
    activation_energy: Optional[float] = None,
    x_scale: str = "linear",
    spacing: float = 1,
) -> plt.Figure:
    if not isinstance(setups, list):
        setups = [setups]

    y_label = (
        "Total Number of Impurities"
        if plot_type == "impurities"
        else r"Outgassing Rate [mBar$\,\cdot\,$Liter/s]"
    )

    fig, ax = plt.subplots(figsize=fig_size)
    ax.set_xlabel(f"Time [{time_scale}]")
    ax.set_ylabel(y_label)
    ax.set_xscale(x_scale)
    ax.set_yscale("log")
    ax.grid(True, which="both", ls="--", alpha=0.5)

    conversion_factor = {"Seconds": 1, "Hours": 3600, "Days": 86400}[time_scale]

    temperature_changes_plotted = False

    for setup in setups:
        if activation_energy is not None:
            setup.activation_energy = activation_energy
        if temperatures is not None:
            setup.temperatures = temperatures

        original_temperatures = setup.temperatures.copy()
        setup_time_ranges = (
            time_range
            if isinstance(time_range, dict)
            else {temp: time_range for temp in setup.temperatures}
        )

        all_timestamps = []
        all_data = []
        last_impurity = None

        for i, temp in enumerate(original_temperatures):
            setup.temperatures = [temp]
            get_diff_temp(setup)

            if i == 0 or not isinstance(time_range, dict):
                get_initial_impurities(setup, units="#")
            elif last_impurity is not None:
                setup.initial_impurities = last_impurity

            start, end = setup_time_ranges[temp]
            time_range_seconds = (start * conversion_factor, end * conversion_factor)
            timestamps = get_time_stamps(
                points=time_range_seconds, spacing=spacing, time_scale="Seconds"
            )

            if plot_type == "impurities":
                data = TimeEvolution.get_impurities_vs_time(setup, time=timestamps)
            else:
                impurities = TimeEvolution.get_impurities_vs_time(
                    setup, time=timestamps
                )
                data = TimeEvolution.get_flow_rate_vs_time(
                    setup, time=timestamps, impurities=impurities, units="mBar Liter"
                )

            plot_timestamps = [t / conversion_factor for t in timestamps[0]]

            if not isinstance(time_range, dict):
                label = f"{setup.version}, {'Thickness' if plot_type == 'impurities' else 'Area'} (cm{'²' if plot_type == 'outgassing' else ''})={setup.thickness if plot_type == 'impurities' else setup.area}, Temperature (K)={temp}"
                ax.plot(plot_timestamps, data[0], label=label)
            else:
                all_timestamps.extend(plot_timestamps)
                all_data.extend(data[0])
                last_impurity = data[0][-1]

        if isinstance(time_range, dict):
            label = f"{setup.version}, {'Thickness' if plot_type == 'impurities' else 'Area'} (cm{'²' if plot_type == 'outgassing' else ''})={setup.thickness if plot_type == 'impurities' else setup.area}"
            ax.plot(all_timestamps, all_data, label=label)

            if not temperature_changes_plotted:
                for i, temp in enumerate(original_temperatures[1:], 1):
                    change_time = setup_time_ranges[temp][0]
                    ax.axvline(x=change_time, color="r", linestyle="--", alpha=0.5)
                    ax.text(
                        change_time,
                        ax.get_ylim()[1],
                        f"{original_temperatures[i-1]}K → {temp}K",
                        rotation=90,
                        va="top",
                        ha="right",
                    )
                temperature_changes_plotted = True

        setup.temperatures = original_temperatures

    ax.legend()
    fig.tight_layout()
    return fig


def plot_electron_lifetime(
    setup: OutgassingSetup,
    time_range: Tuple[float, float],
    out_diffusion_values: List[float],
    initial_impurities: float = 1,
    circulation_rate: float = 200 / 60,
    purification_efficiency: float = 1,
    time_scale: str = "Days",
    fig_size: Tuple[int, int] = (10, 6),
    y_scale: str = "log",
    spacing: float = 1 / 24,
) -> plt.Figure:
    """
    Calculate and plot electron lifetime results for different out-diffusion rates.

    Args:
        setup (OutgassingSetup): The setup to plot electron lifetime for.
        time_range (Tuple[float, float]): The time range to plot (start, end).
        out_diffusion_values (List[float]): List of out-diffusion rates to calculate for.
        initial_impurities (float): Initial impurities concentration.
        circulation_rate (float): Circulation rate.
        purification_efficiency (float): Purification efficiency.
        time_scale (str): The time scale to use ("Days", "Hours", or "Seconds").
        fig_size (Tuple[int, int]): The figure size (width, height).
        y_scale (str): The y-axis scale ("linear" or "log").
        spacing (float): The spacing between data points.

    Returns:
        plt.Figure: The generated matplotlib figure.
    """
    lifetime_timestamps = get_time_stamps(
        points=list(time_range), spacing=spacing, time_scale=time_scale
    )

    fig, ax = plt.subplots(figsize=fig_size)
    ax.set_xlabel(f"Time ({time_scale})")
    ax.set_ylabel("Electron Lifetime (ms)")
    ax.set_yscale(y_scale)
    ax.grid(True, which="both", ls="--", alpha=0.5)

    for out_diffusion_value in out_diffusion_values:
        electron_lifetimes, params = TimeEvolution.get_electron_lifetime_vs_time(
            setup,
            time=lifetime_timestamps,
            initial_impurities=initial_impurities,
            circulation_rate=circulation_rate,
            purification_efficiency=purification_efficiency,
            out_diffusion=out_diffusion_value,
        )

        label = (
            f"Out Diffusion rate (ppb.l/s)={params['out_diffusion']}, "
            f"Flow rate (SLPM)={params['circulation_rate']*60}, "
            f"Efficiency={params['purification_efficiency']}"
        )
        ax.plot(lifetime_timestamps[0], electron_lifetimes, label=label)

    ax.legend()
    fig.tight_layout()
    return fig
