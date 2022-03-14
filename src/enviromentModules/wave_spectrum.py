"""
License
This file is part of UiS-Aqua.

UiS-Aqua is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

UiS-Aqua is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with UiS-Aqua. If not, see <https://www.gnu.org/licenses/>.

-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Any questions about this code,
please email: hui.cheng@uis.no \n
"""

import numpy as np
from numpy import pi


def pierson_moskowitz_spectra(omega, hs, tp):
    """
    Generate Pierson-Moskowitz wave spectrum \n
    ref: DNVGL-RP-C205 ver.2018,pp64
    :param omega: numpy.ndarray | Array of frequencies
    :param hs: float  |  significant wave height [m]
    :param tp: float  |  peak wave period [s]
    :return: numpy.ndarray  |     Array of shape omega with wave energy densities
    """
    omega_p = float(2 * pi / tp)
    spectra = 5.0 / 16.0 * pow(hs, 2) * pow(omega_p, 4) * pow(omega, -5) * np.exp(-1.25 * pow(omega / omega_p, -4))
    return spectra


def jonswap_spectra(omega, hs, tp, gamma=3.3, gamma_auto=False):
    """
    Generate JONSWAP spectrum \n
    The Jonswap wave spectrum is expected to be a reasonable model for:
    3.6 < Tp/sqrt(hs) < 5 \n
    ref: DNVGL-RP-C205 ver.2018,pp64
    :param omega: numpy.ndarray | Array of frequencies
    :param hs: float  |  significant wave height [m]
    :param tp: float  |  peak wave period [s]
    :param gamma: float | peak shape parameter (default: 3.3)
    :param gamma_auto: Boolean (True or False)  |  The value of gamma will be calculated automatically
    :return: numpy.ndarray  |     Array of shape omega with wave energy densities
    """
    # default values
    sigma_low = 0.07
    sigma_high = 0.09
    if gamma_auto:
        if tp / np.sqrt(hs) <= 3.6:
            gamma = 5.0
        elif tp / np.sqrt(hs) >= 5:
            gamma = 1.0
        else:
            gamma = np.exp(5.75 - 1.15 * tp / np.sqrt(hs))

    # Pierson-Moskowitz
    omega_p = float(2 * pi / tp)
    # pm_spectra = 5.0 / 16.0 * pow(hs, 2) * pow(omega_p, 4) * pow(omega, -5) * np.exp(-1.25 * pow(omega / omega_p, -4))
    pm_spectra = pierson_moskowitz_spectra(omega, hs, tp)
    # JONSWAP
    a_gamma = 1 - 0.287 * np.log(gamma)
    sigma = np.ones(omega.shape) * sigma_low
    sigma[omega > omega_p] = sigma_high
    spectra = a_gamma * pm_spectra * pow(gamma, np.exp(-0.5 * pow((omega - omega_p) / (sigma * omega_p), 2)))
    return spectra


if __name__ == "__main__":
    sp1 = jonswap_spectra(np.linspace(0.01, 3, 1000), 5, 8, 3.3)
    print(sp1)
