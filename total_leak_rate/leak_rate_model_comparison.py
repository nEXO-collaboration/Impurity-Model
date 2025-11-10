import numpy as np
import matplotlib.pyplot as plt

# Constants
R = 8.314  # J/(mol·K)
T = 165  # K
V = 22.4e-3  # m³/mol
k_S = 1e11  # L/(mol·s)

# Variables
eta = 1  # purifier efficiency
F = np.linspace(50, 650, 100)  # SLPM
F_m3s = F * 1.66667e-5  # convert SLPM to m³/s
tau = np.logspace(-3, -2, 100)  # s, range from 1 ms to 10 ms
n_p = 0  # 1e-12  # mol/L


def leak_rate_limit_simple(r, tau):
    """Calculate leak rate limit using the simple model from Section 1.1"""
    return 2.72e-10 * (r / tau)


def leak_rate_limit_dynamic(eta, F, k_S, tau, n_p):
    """Calculate leak rate limit using the dynamic model from Section 1.2"""
    return (R * T / V) * ((eta * F) / (k_S * tau) - n_p * F)


# 1. Comparison window (Comparing models)
plt.figure("Model Comparisons", figsize=(18, 6))

# 1.1 Compare models across flow rates
plt.subplot(1, 3, 1)  # Create subplot 1 (1 row, 3 columns, 1st plot)
J_simple_flow = leak_rate_limit_simple(F, 5e-3)
J_dynamic_flow = leak_rate_limit_dynamic(eta, F_m3s, k_S, 5e-3, n_p)
plt.plot(F, J_simple_flow, label="Simple Model")
plt.plot(F, J_dynamic_flow, label="Dynamic Model")
plt.xlabel("Flow Rate (SLPM)")
plt.ylabel("Leak Rate Limit (mbar·L/s)")
plt.title("Leak Rate Limit vs Flow Rate")
plt.legend()
plt.grid(True)
plt.yscale("log")

# 1.2 Compare models across electron lifetimes
plt.subplot(1, 3, 2)  # Create subplot 2 (1 row, 3 columns, 2nd plot)
J_simple_tau = leak_rate_limit_simple(350, tau)
J_dynamic_tau = leak_rate_limit_dynamic(eta, 350 * 1.66667e-5, k_S, tau, n_p)
plt.plot(tau * 1000, J_simple_tau, label="Simple Model")
plt.plot(tau * 1000, J_dynamic_tau, label="Dynamic Model")
plt.xlabel("Electron Lifetime (ms)")
plt.ylabel("Leak Rate Limit (mbar·L/s)")
plt.title("Leak Rate Limit vs Electron Lifetime")
plt.legend()
plt.grid(True)
plt.xscale("log")
plt.yscale("log")

# 1.3 Model Ratio Analysis (Simple/Dynamic)
plt.subplot(1, 3, 3)  # Create subplot 3 (1 row, 3 columns, 3rd plot)
ratio_flow = J_simple_flow / J_dynamic_flow  # Inverted ratio
ratio_tau = J_simple_tau / J_dynamic_tau  # Inverted ratio
plt.plot(F, ratio_flow, label="Ratio (Simple/Dynamic) vs Flow Rate")
plt.plot(tau * 1000, ratio_tau, label="Ratio (Simple/Dynamic) vs Electron Lifetime")
plt.xlabel("Flow Rate (SLPM) / Electron Lifetime (ms)")
plt.ylabel("Ratio (Simple/Dynamic)")
plt.title("Model Ratio Analysis (Inverted)")
plt.legend()
plt.grid(True)

plt.tight_layout()

# 2. Dynamic behavior window (Efficiency and Concentration sensitivity)
plt.figure("Dynamic Behavior", figsize=(12, 6))

# 2.1 Sensitivity to purifier efficiency
plt.subplot(1, 2, 1)  # Create subplot 1 (1 row, 2 columns, 1st plot)
eta_range = np.linspace(0.9, 0.9999, 100)
J_dynamic_eta = leak_rate_limit_dynamic(eta_range, 350 * 1.66667e-5, k_S, 5e-3, n_p)
plt.plot(eta_range, J_dynamic_eta)
plt.xlabel("Purifier Efficiency")
plt.ylabel("Leak Rate Limit (mbar·L/s)")
plt.title("Sensitivity to Purifier Efficiency")
plt.grid(True)

# 2.2 Sensitivity to purifier output concentration
plt.subplot(1, 2, 2)  # Create subplot 2 (1 row, 2 columns, 2nd plot)
n_p_range = np.logspace(-14, -10, 100)
J_dynamic_n_p = leak_rate_limit_dynamic(eta, 350 * 1.66667e-5, k_S, 5e-3, n_p_range)
plt.semilogx(n_p_range, J_dynamic_n_p)
plt.xlabel("Purifier Output Concentration (mol/L)")
plt.ylabel("Leak Rate Limit (mbar·L/s)")
plt.title("Sensitivity to Purifier Output Concentration")
plt.grid(True)

plt.tight_layout()

# Show both windows
plt.show()

# Print key results
print(
    f"Average ratio (Simple/Dynamic) vs Flow Rate: {np.mean(ratio_flow):.4f} ± {np.std(ratio_flow):.4f}"
)
print(
    f"Average ratio (Simple/Dynamic) vs Electron Lifetime: {np.mean(ratio_tau):.4f} ± {np.std(ratio_tau):.4f}"
)
