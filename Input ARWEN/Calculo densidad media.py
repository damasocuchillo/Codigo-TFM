import numpy as np

NA = 6.022e23

# Long[m]
R = 600e-6 
H= 400e-6
e= 8e-6
L= 35e-6

# Densidad [kg/m^3]
rho_CH = 0.9
rho_Cu_LD = 9
rho_Cu_HD = 8960


S_total = R*H
S_rod = e*L

S_Cu_20 = 20 * S_rod
S_Cu_50 = 50 * S_rod

S_CH_20 = S_total - S_Cu_20
S_CH_50 = S_total - S_Cu_50

rho_A = (S_CH_20 * rho_CH + S_Cu_20 * rho_Cu_LD) / S_total
rho_B = (S_CH_50 * rho_CH + S_Cu_50 * rho_Cu_LD) / S_total

rho_C = (S_CH_20 * rho_CH + S_Cu_20 * rho_Cu_HD) / S_total
rho_D = (S_CH_50 * rho_CH + S_Cu_50 * rho_Cu_HD) / S_total


print(f'Caso A rho = {rho_A} kg/m^3')
print(f'Caso B rho = {rho_B} kg/m^3')
print(f'Caso C rho = {rho_C} kg/m^3')
print(f'Caso D rho = {rho_D} kg/m^3')

import numpy as np
import matplotlib.pyplot as plt


# -------------------------
# Datos
# -------------------------
rho = np.array([rho_A, rho_B, rho_C, rho_D])

sigma_t = np.array([
    2.8e-6,
    4e-6,
    9.0e-6,
    1.755e-5
])

labels = ["Caso A", "Caso B", "Caso C", "Caso D"]

# -------------------------
# Plot
# -------------------------
plt.figure(figsize=(6,5))

plt.scatter(rho, sigma_t, s=90)

for i, txt in enumerate(labels):
    plt.annotate(
        txt,
        (rho[i], sigma_t[i]),
        textcoords="offset points",
        xytext=(6,6)
    )

plt.xlabel(r"Densidad $\overline{\rho}$ [kg/m$^3$]")
plt.ylabel(r"$\sigma_{RT} \cdot \tau_{\mathrm{choque}}$")

plt.title("Crecimiento efectivo de RTI vs densidad")

plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.show()