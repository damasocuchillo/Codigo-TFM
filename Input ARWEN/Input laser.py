import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import brentq

N = 1000  # Numero de puntos
R = 300e-6  # Radio de la lente [m]
d = 200e-6  # Distancia de enfoque [m]
sigma = 0.5  # Ancho pulso [m]
m = 4  # Exponente Gaussiano

# HAZ LASER con direccion definida
x_cell = np.zeros(N + 1)  # Bordes de celda
x = np.zeros(N)  # Puntos x iniciales del haz
y = np.full(N, -d)  # Puntos y iniciales del haz

# Añadimos desviación angular aleatoria de hasta ±5 grados
alpha= np.pi / 2 # DIrección del haz respeto al eje x
desviacion_grados = 5
desviaciones = np.random.uniform(-desviacion_grados, desviacion_grados, N)
desviaciones_rad = np.radians(desviaciones)

kx = np.cos(alpha - desviaciones_rad)  # componente x
ky = np.sin(alpha - desviaciones_rad)  # componente y

def Integral(r):
    return np.exp((-r ** m)/ (sigma ** m)) * r

def g(x, target):
    val, _ = quad(Integral, 0, x)
    return val - target

# Integrar de 0 a R
A, _ = quad(Integral, 0, R)

# Asignar primer borde como 0
x_cell[0] = 0.0

# Calcular los bordes de celda, empezando desde i=1
for i in range(1, N + 1):
    target = A * i / N
    a = x_cell[i - 1] + 1e-15
    b = R
    root = brentq(g, a, b, args=(target,))
    x_cell[i] = root

print('x    y   kx  ky')
for i in range(len(x)):
    x[i] = (x_cell[i] + x_cell[i + 1]) / 2
    print(f'{i+1} {x[i]:.6e} {y[i]:.6e} {kx[i]:.6f} {ky[i]:.6f}')


ruta= r'C:\Users\34620\Documents\TFM\Input ARWEN\Python\input_laser.txt'

# Guardar en archivo
with open(ruta, "w") as f:
    f.write("# Input Laser ARWEN \n")
    f.write(f"# Puntos = {N}\n")
    f.write(f"# Radio de la lente R = {R} m\n")
    f.write(f"# Distancia de enfoque d= {d} m\n")
    f.write(f"# Ancho de pulso sigma = {sigma}\n")
    f.write(f"# Exponente Gaussiano m = {m}\n")
    f.write(f"# Desviación angular aleatoria de hasta ±{desviacion_grados} grados\n")
    f.write("# x    y   kx  ky\n")
    for i in range(len(x)):
        x[i] = (x_cell[i] + x_cell[i + 1]) / 2
        f.write(f"{x[i]:.6e} {y[i]:.6e} {kx[i]:.6f} {ky[i]:.6f}\n")


# r values
r = np.linspace(0, R, 1000)
I_vals = np.exp(-(np.abs(r) / sigma)**m)

r= r * 1e6
plt.figure()
plt.plot(r, I_vals)
plt.ylim(0, 1)

plt.xlabel("r [μm]")
plt.ylabel(r'$\frac{I(r)}{I_0}$')
plt.title(r'Distribución de intensidad del láser $\frac{I(r)}{I_0}$')
plt.show()

print(I_vals[0])




