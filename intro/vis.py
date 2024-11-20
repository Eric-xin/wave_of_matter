import numpy as np
import matplotlib.pyplot as plt

# Constants
h = 6.62607015e-34  # Planck's constant (J·s)
c = 3.0e8           # Speed of light (m/s)
k_B = 1.380649e-23  # Boltzmann's constant (J/K)

# De Broglie wavelength
def de_broglie_wavelength(momentum):
    return h / momentum

# Planck's black body radiation
def planck_law(wavelength, temperature):
    return (2 * h * c**2) / (wavelength**5 * (np.exp((h * c) / (wavelength * k_B * temperature)) - 1))

# Plotting de Broglie wavelength
momentum = np.linspace(1e-24, 1e-20, 500)  # Momentum range (kg·m/s)
wavelength = de_broglie_wavelength(momentum)

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(momentum, wavelength)
plt.yscale("log")
plt.title("De Broglie Wavelength")
plt.xlabel("Momentum (kg·m/s)")
plt.ylabel("Wavelength (m), log scale")
plt.grid(True)

# Plotting Planck's black body radiation
wavelengths = np.linspace(1e-9, 3e-6, 500)  # Wavelength range (m)
temperature = 5800  # Temperature (K)
intensity = planck_law(wavelengths, temperature)

plt.subplot(1, 2, 2)
plt.plot(wavelengths * 1e9, intensity)  # Convert wavelength to nm for better visualization
plt.title("Planck's Black Body Radiation")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Intensity")
plt.grid(True)

plt.tight_layout()
plt.show()