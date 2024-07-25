# TOUCAN: Toolkit for Outgassing Contamination Analysis

TOUCAN is a Python package designed for analyzing outgassing and impurities in gas systems, with a focus on xenon gas experiments. It provides tools for calculating diffusion rates, impurity levels, and electron lifetimes in various experimental setups.

## Installation

You can install TOUCAN using pip:

```
pip install toucan-outgassing
```

## Usage

To use TOUCAN, you need to provide a JSON file (`library.json`) containing the necessary data for your analysis. This file should be in your working directory or its path should be specified when creating an `Outgassing_setup` instance.

Here's a quick example of how to use TOUCAN in a Jupyter notebook:

```python
import toucan as tc

# Create an Outgassing_setup instance
# Assuming 'library.json' is in the current working directory
setup = tc.Outgassing_setup("EXO-200", "Teflon", "Oxygen")

# If 'library.json' is in a different location, specify the path:
# setup = tc.Outgassing_setup("EXO-200", "Teflon", "Oxygen", data_file="/path/to/library.json")

# Generate timestamps
timestamps = tc.get_time_stamps([0, 1, 2], 0.1, "Days")

# Calculate impurities over time
impurities = setup.get_impurities_vs_time(timestamps)

# Plot the results
tc.plot_data((10, 6), "Time (Days)", "Impurities", 
             [(timestamps[0], impurities[0], "Teflon Oxygen Impurities")])

# Use constants
print(f"Avogadro's number: {tc.AVOGADRO_NUMBER}")
print(f"Ideal gas molar volume: {tc.IDEAL_GAS_MOLAR_VOLUME} liter/mol")
```

For more detailed usage instructions, please refer to the documentation.

## Data File

The `library.json` file should contain the necessary data for your analysis, including material properties, system configurations, and gas properties. Ensure this file is in your working directory or provide its path when initializing the `Outgassing_setup` class.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.