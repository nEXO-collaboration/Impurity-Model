from .setup import OutgassingSetup
from .calculations import get_diff_temp, get_initial_impurities
from .time_evolution import TimeEvolution
from .utils import get_time_stamps, convert_time
from .plotting import (
    plot_outgassing_rate,
    plot_polymers,
    plot_electron_lifetime,
)
from .constants import (
    IDEAL_GAS_MOLAR_VOLUME,
    AVOGADRO_NUMBER,
    BOLTZMANN_CONSTANT_EV,
    BOLTZMANN_CONSTANT_L,
    GXE_DENSITY,
)

__all__ = [
    "OutgassingSetup",
    "get_diff_temp",
    "get_initial_impurities",
    "TimeEvolution",
    "get_time_stamps",
    "convert_time",
    "plot_outgassing_rate",
    "plot_polymers",
    "plot_electron_lifetime",
    "IDEAL_GAS_MOLAR_VOLUME",
    "AVOGADRO_NUMBER",
    "BOLTZMANN_CONSTANT_EV",
    "BOLTZMANN_CONSTANT_L",
    "GXE_DENSITY",
]

__version__ = "0.2.0"
