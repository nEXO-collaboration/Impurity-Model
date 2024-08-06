import matplotlib.pyplot as plt
from typing import List, Tuple, Optional, Union


def get_time_stamps(
    points: List[Union[int, float]],
    spacing: Union[int, float],
    time_scale: str = "Seconds",
) -> List[List[float]]:
    """
    Generate a list of lists containing timestamps in seconds.
    Each sublist represents the range from one point to the next,
    incremented by the given spacing value.
    """
    scale_factors = {"Days": 86400, "Hours": 3600, "Seconds": 1}
    if time_scale not in scale_factors:
        raise ValueError(
            "Unsupported time scale. Choose from 'Days', 'Hours', or 'Seconds'."
        )
    converted_points = [point * scale_factors[time_scale] for point in points]
    timestamps = []
    for start_point, end_point in zip(converted_points, converted_points[1:]):
        current_time = start_point
        time_segment = []
        while current_time < end_point:
            time_segment.append(current_time)
            current_time += spacing * scale_factors[time_scale]
        timestamps.append(time_segment)
    return timestamps


def plot_data(
    fig_size: Tuple[int, int],
    x_label: str,
    y_label: str,
    data_sets: List[Tuple[List[float], List[float], str]],
    x_scale: Optional[str] = None,
    y_scale: Optional[str] = None,
    time_unit: str = "Seconds",
) -> None:
    """
    Plot data sets with customizable labels and scales.
    """
    plt.figure(figsize=fig_size)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if x_scale:
        plt.xscale(x_scale)
    if y_scale:
        plt.yscale(y_scale)
    plt.grid(which="both", linestyle="--", linewidth=0.5, color="gray")
    plt.minorticks_on()

    # Convert time to the specified unit
    time_conversion_factor = 1
    if time_unit == "Days":
        time_conversion_factor = 86400  # Seconds in a day
    elif time_unit == "Hours":
        time_conversion_factor = 3600  # Seconds in an hour

    for time_data, data, label in data_sets:
        adjusted_time_data = [t / time_conversion_factor for t in time_data]
        plt.plot(adjusted_time_data, data, label=label)
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

    :param setup: Outgassing_setup instance
    :param time_range: Tuple of (min_time, max_time)
    :param time_scale: Unit for time axis (default: "Hours")
    :param rate_unit: Unit for rate axis (default: "mbar.l/s")
    :param unbaked_flow_rate: Unbaked flow rate (default: 1e-8)
    :param initial_pumped_time: Initial pumped time in seconds (default: 2600)
    :param fig_size: Figure size as (width, height) (default: (10, 6))
    """
    # Generate time stamps
    time_stamps = get_time_stamps(
        points=list(time_range), spacing=1, time_scale=time_scale
    )

    # Calculate flow rates
    flow_rates = setup.get_steel_flow_rate_vs_pumping_time(
        time=time_stamps,
        unbaked_flow_rate=unbaked_flow_rate,
        initial_pumped_time=initial_pumped_time,
    )

    # Prepare data for plotting
    data_sets = [
        (time_stamps[0], flow_rates[0], f"{setup.version}, Area (cm²)={setup.area}")
    ]

    # Plot the data
    plot_data(
        fig_size=fig_size,
        x_label=f"Time ({time_scale})",
        y_label=f"Outgassing Rate ({rate_unit})",
        data_sets=data_sets,
        x_scale="log",
        y_scale="log",
        time_unit=time_scale,
    )


def convert_time(time: List[float], from_unit: str, to_unit: str) -> List[float]:
    """
    Convert time from one unit to another.

    :param time: List of time points
    :param from_unit: Original time unit
    :param to_unit: Desired time unit
    :return: List of converted time points
    """
    conversion_factors = {"Seconds": 1, "Hours": 3600, "Days": 86400}
    factor = conversion_factors[from_unit] / conversion_factors[to_unit]
    return [t * factor for t in time]
