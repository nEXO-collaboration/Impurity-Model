import matplotlib.pyplot as plt
import numpy as np
from typing import List, Tuple, Optional, Union

def get_time_stamps(
    points: List[Union[int, float]],
    spacing: Union[int, float],
    time_scale: str = "Seconds"
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
    time_unit: str = "Seconds"
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

    for time_data, data, label in data_sets:
        adjusted_time_data = [t / time_conversion_factor for t in time_data]
        plt.plot(adjusted_time_data, data, label=label)

    plt.legend()
    plt.show()