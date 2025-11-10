import numpy as np
import matplotlib.pyplot as plt

# Parameters
D = 1e-5  # Diffusion coefficient [cm²/s]
d = 0.5  # Thickness [cm]
c0 = 1  # Initial concentration

# Time scales
tc_1 = d**2 / (6 * D)  # First critical time
tc_2 = 0.5 * d**2 / (np.pi**2 * D)  # Second critical time (from three-step model)

# Time array (log scale)
t = np.logspace(-2, 5, 1000)  # From 0.01s to 10⁵s


def full_solution(t):
    """Full series solution of diffusion equation"""
    n_terms = 100  # Number of terms in the sum
    J = np.zeros_like(t)

    for n in range(n_terms):
        # Using the equation from the chapter
        J += np.exp(-((np.pi * (2 * n + 1) / d) ** 2) * D * t)

    return 4 * c0 * D / d * J


def short_time_approx(t):
    """Short time approximation"""
    return c0 * np.sqrt(D / (np.pi * t))


def long_time_approx(t):
    """Long time approximation (first term of series)"""
    return 4 * c0 * D / d * np.exp(-np.pi**2 * D * t / d**2)


# Calculate solutions
J_full = full_solution(t)
J_short = short_time_approx(t)
J_long = long_time_approx(t)

# Plotting
plt.figure(figsize=(12, 8))
plt.loglog(t, J_full, "k-", label="Full solution")
plt.loglog(t, J_short, "--", label="Short-time approximation")
plt.loglog(t, J_long, "--", label="Long-time approximation")

# Add vertical lines for critical times
plt.axvline(tc_1, color="r", linestyle=":", label=f"tc_1 = d²/6D = {tc_1:.1f}s")
plt.axvline(tc_2, color="g", linestyle=":", label=f"tc_2 = 0.5d²/π²D = {tc_2:.1f}s")

plt.grid(True, which="both", ls="-", alpha=0.2)
plt.xlabel("Time [s]")
plt.ylabel("Outgassing rate [a.u.]")
plt.title("Comparison of Diffusion Solutions and Critical Times")
plt.legend()


# Calculate relative errors at both critical times
def relative_error(t_crit):
    """Calculate relative errors between approximations and full solution at t_crit"""
    idx = np.abs(t - t_crit).argmin()
    err_short = abs(J_short[idx] - J_full[idx]) / J_full[idx]
    err_long = abs(J_long[idx] - J_full[idx]) / J_full[idx]
    return err_short, err_long


err_short_1, err_long_1 = relative_error(tc_1)
err_short_2, err_long_2 = relative_error(tc_2)

print(f"\nRelative errors at tc_1 = {tc_1:.1f}s:")
print(f"Short-time approximation: {err_short_1:.1%}")
print(f"Long-time approximation: {err_long_1:.1%}")

print(f"\nRelative errors at tc_2 = {tc_2:.1f}s:")
print(f"Short-time approximation: {err_short_2:.1%}")
print(f"Long-time approximation: {err_long_2:.1%}")

plt.show()
