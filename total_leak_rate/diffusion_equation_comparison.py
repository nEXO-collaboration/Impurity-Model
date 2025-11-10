import numpy as np
import matplotlib.pyplot as plt


def J1(t, c0, D, d):
    """First equation"""
    n = np.arange(0, 10000)  # Sum to 10000 terms for approximation
    sum_terms = np.sum(np.exp(-((np.pi * (2 * n + 1) / d) ** 2 * D * t)))
    return (4 * c0 * D / d) * sum_terms


def J2(t, c0, D, d):
    """Second equation"""
    n = np.arange(1, 10000)  # Sum to 10000 terms for approximation
    sum_terms = np.sum((-1) ** n * np.exp(-(n**2 * d**2) / (4 * D * t)))
    return c0 * np.sqrt(D / (np.pi * t)) * (1 + 2 * sum_terms)


# Parameters
c0 = 1  # Initial concentration
D = 1e-5  # Diffusion coefficient (cm^2/s)
d = 0.1  # Thickness (cm)

# Time array (log scale), now stopping at 10^4
t = np.logspace(-2, 4, 10000)

# Calculate J1 and J2
J1_values = [J1(ti, c0, D, d) for ti in t]
J2_values = [J2(ti, c0, D, d) for ti in t]

# Plotting
plt.figure(figsize=(10, 6))
plt.loglog(t, J1_values, label="J1")
plt.loglog(t, J2_values, label="J2")
plt.xlabel("Time (s)")
plt.ylabel("J (cm^-2 s^-1)")
plt.title("Comparison of Diffusion Equations")
plt.legend()
plt.grid(True)
plt.show()

# Calculate and print the maximum relative difference
rel_diff = np.abs((np.array(J1_values) - np.array(J2_values)) / np.array(J1_values))
max_rel_diff = np.max(rel_diff)
print(f"Maximum relative difference: {max_rel_diff:.2%}")

# Print specific values for analysis
print("\nSpecific values for analysis:")
print("Time (s) \t J1 (cm^-2 s^-1) \t J2 (cm^-2 s^-1) \t Relative Difference")
for ti in [0.01, 0.1, 1, 10, 100, 1000, 10000]:
    j1 = J1(ti, c0, D, d)
    j2 = J2(ti, c0, D, d)
    diff = abs(j1 - j2) / j1
    print(f"{ti:.2e} \t {j1:.4e} \t {j2:.4e} \t {diff:.2%}")

# Calculate characteristic time
t_char = d**2 / D
print(f"\nCharacteristic diffusion time (d^2/D): {t_char:.2e} s")
