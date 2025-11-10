from typing import List, Union


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
