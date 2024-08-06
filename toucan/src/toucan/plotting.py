import matplotlib.pyplot as plt
from typing import List, Tuple, Optional, Dict, Union
from .utils import get_time_stamps
from .calculations import get_diff_temp, get_initial_impurities
from .time_evolution import TimeEvolution
from .setup import OutgassingSetup


def _internal_plot(
    data_sets: List[Tuple[List[float], List[float], str]],
    x_label: str,
    y_label: str,
    fig_size: Tuple[int, int] = (10, 6),
    x_scale: Optional[str] = None,
    y_scale: Optional[str] = None,
) -> None:
    plt.figure(figsize=fig_size)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # Set the x and y scales if specified
    if x_scale:
        plt.xscale(x_scale)
    if y_scale:
        plt.yscale(y_scale)

    plt.grid(which="both", linestyle="--", linewidth=0.5, color="gray")
    plt.minorticks_on()

    for x_data, y_data, label in data_sets:
        plt.plot(x_data, y_data, label=label)

    plt.legend()
    plt.show()


def plot_outgassing_rate(
    setup,
    time_range: Tuple[float, float] = (1e-1, 1e4),
    time_scale: str = "Hours",
    rate_unit: str = "mbar.l/s",
    unbaked_flow_rate: float = 1e-8,
    initial_pumped_time: float = 2600,
    fig_size: Tuple[int, int] = (10, 6),
) -> None:
    """
    Plot the outgassing rate for a given setup.
    """
    time_stamps = get_time_stamps(
        points=list(time_range), spacing=1, time_scale=time_scale
    )

    flow_rates = TimeEvolution.get_steel_flow_rate_vs_pumping_time(
        setup,
        time=time_stamps,
        unbaked_flow_rate=unbaked_flow_rate,
        initial_pumped_time=initial_pumped_time,
    )

    data_sets = [
        (time_stamps[0], flow_rates[0], f"{setup.version}, Area (cm²)={setup.area}")
    ]

    _internal_plot(
        data_sets=data_sets,
        x_label=f"Time ({time_scale})",
        y_label=f"Outgassing Rate ({rate_unit})",
        fig_size=fig_size,
        x_scale="log",
        y_scale="log",
    )


def plot_polymers(
    setups: Union[List[OutgassingSetup], OutgassingSetup],
    time_range: Union[Tuple[float, float], Dict[float, Tuple[float, float]]],
    time_scale: str = "Seconds",
    plot_type: str = "impurities",
    fig_size: Tuple[int, int] = (10, 6),
    temperatures: Optional[List[float]] = None,
    activation_energy: Optional[float] = None,
    x_scale: str = "linear",
    spacing: float = 1,  # Default spacing of 1 second
) -> plt.Figure:
    """
    Plot impurities or outgassing rates for one or multiple polymer setups.
    """
    if not isinstance(setups, list):
        setups = [setups]

    # Define y_label
    if plot_type == "impurities":
        y_label = "Total Number of Impurities"
    elif plot_type == "outgassing":
        y_label = r"Outgassing Rate [mBar$\,\cdot\,$Liter/s]"
    else:
        raise ValueError("Invalid plot_type. Use 'impurities' or 'outgassing'.")

    fig, ax = plt.subplots(figsize=fig_size)
    ax.set_xlabel(f"Time [{time_scale}]")
    ax.set_ylabel(y_label)
    ax.set_xscale(x_scale)
    ax.set_yscale("log")
    ax.grid(which="both", linestyle="--", linewidth=0.5, color="gray")

    # Convert time to seconds if necessary
    time_conversion = {"Seconds": 1, "Hours": 3600, "Days": 86400}
    conversion_factor = time_conversion.get(time_scale, 1)

    for setup in setups:
        if activation_energy is not None:
            setup.activation_energy = activation_energy

        if temperatures is not None:
            setup.temperatures = temperatures

        original_temperatures = setup.temperatures.copy()

        # Handle time ranges
        if isinstance(time_range, tuple):
            # Single time range for all temperatures
            setup_time_ranges = {temp: time_range for temp in setup.temperatures}
            single_range = True
        else:
            # Multiple time ranges for different temperatures
            setup_time_ranges = time_range
            single_range = False

        all_timestamps = []
        all_data = []
        last_impurity = None

        for i, temp in enumerate(original_temperatures):
            setup.temperatures = [temp]  # Set single temperature for calculation
            get_diff_temp(setup)  # Recalculate diffusion constants for this temperature

            if i == 0 or single_range:
                get_initial_impurities(setup, units="#")
            elif last_impurity is not None:
                setup.initial_impurities = last_impurity

            start, end = setup_time_ranges[temp]
            time_range_seconds = (start * conversion_factor, end * conversion_factor)

            # Generate timestamps
            timestamps = get_time_stamps(
                points=time_range_seconds,
                spacing=spacing,
                time_scale="Seconds",  # Always use seconds internally
            )

            if plot_type == "impurities":
                data = TimeEvolution.get_impurities_vs_time(setup, time=timestamps)
            else:  # plot_type == "outgassing"
                impurities = TimeEvolution.get_impurities_vs_time(
                    setup, time=timestamps
                )
                data = TimeEvolution.get_flow_rate_vs_time(
                    setup, time=timestamps, impurities=impurities, units="mBar Liter"
                )

            plot_timestamps = [t / conversion_factor for t in timestamps[0]]

            if single_range:
                label = (
                    f"{setup.version}, "
                    f"{'Thickness' if plot_type == 'impurities' else 'Area'} "
                    f"(cm{'²' if plot_type == 'outgassing' else ''})="
                    f"{setup.thickness if plot_type == 'impurities' else setup.area}, "
                    f"Temperature (K)={temp}"
                )
                ax.plot(plot_timestamps, data[0], label=label)
            else:
                all_timestamps.extend(plot_timestamps)
                all_data.extend(data[0])
                last_impurity = data[0][-1]

        if not single_range:
            label = (
                f"{setup.version}, "
                f"{'Thickness' if plot_type == 'impurities' else 'Area'} "
                f"(cm{'²' if plot_type == 'outgassing' else ''})="
                f"{setup.thickness if plot_type == 'impurities' else setup.area}"
            )
            ax.plot(all_timestamps, all_data, label=label)

            # Add vertical lines and text to indicate temperature changes
            for i, temp in enumerate(original_temperatures[1:], 1):
                change_time = setup_time_ranges[temp][0]
                ax.axvline(x=change_time, color="r", linestyle="--", alpha=0.5)
                ax.text(
                    change_time,
                    ax.get_ylim()[1],
                    f"{temp}K",
                    rotation=90,
                    va="top",
                    ha="right",
                )

        # Restore original temperatures
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
    x_label: Optional[str] = None,
    y_label: str = "Electron Lifetime (ms)",
    y_scale: str = "log",
    spacing: float = 1 / 24,  # Default spacing of 1 hour
) -> None:
    """
    Calculate and plot electron lifetime results for different out-diffusion rates.

    :param setup: OutgassingSetup instance
    :param time_range: Tuple of (start_time, end_time)
    :param out_diffusion_values: List of out-diffusion rates to calculate for
    :param initial_impurities: Initial impurities concentration
    :param circulation_rate: Circulation rate
    :param purification_efficiency: Purification efficiency
    :param time_scale: Unit for time axis (default: "Days")
    :param fig_size: Figure size as (width, height)
    :param x_label: Label for x-axis (default: "Time ({time_scale})")
    :param y_label: Label for y-axis
    :param y_scale: Scale for y-axis ('linear' or 'log')
    :param spacing: Time spacing for data points
    """
    # Generate timestamps
    lifetime_timestamps = get_time_stamps(
        points=list(time_range), spacing=spacing, time_scale=time_scale
    )

    # Calculate electron lifetimes for each out-diffusion value
    electron_lifetime_results = []
    params_results = []
    for out_diffusion_value in out_diffusion_values:
        electron_lifetimes, params = TimeEvolution.get_electron_lifetime_vs_time(
            setup,
            time=lifetime_timestamps,
            initial_impurities=initial_impurities,
            circulation_rate=circulation_rate,
            purification_efficiency=purification_efficiency,
            out_diffusion=out_diffusion_value,
        )
        electron_lifetime_results.append(electron_lifetimes)
        params_results.append(params)

    # Plotting
    plt.figure(figsize=fig_size)
    if x_label is None:
        x_label = f"Time ({time_scale})"
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.yscale(y_scale)
    plt.grid(which="both", linestyle="--", linewidth=0.5, color="gray")
    plt.minorticks_on()

    for lifetimes, params in zip(electron_lifetime_results, params_results):
        label = (
            f"Out Diffusion rate (ppb.l/s)={params['out_diffusion']}, "
            f"Flow rate (SLPM)={params['circulation_rate']*60}, "
            f"Efficiency={params['purification_efficiency']}"
        )
        plt.plot(lifetime_timestamps[0], lifetimes, label=label)

    plt.legend()
    plt.tight_layout()
    plt.show()
