import numpy as np
import matplotlib.pyplot as plt

# Definir las e
sigma = 1.0
delta = 1.0
limite = delta + 2**(1/6)*sigma

# Definir el rango de la distancia
r = np.linspace(0, limite, 1000)

# Calcular el potencial
U_rc = np.zeros_like(r)
U_rc[r < limite] = 4*((sigma/(r[r < limite] - delta))**12 - (sigma/(r[r < limite] - delta))**6 + 1/4)

d = np.linspace(0.0, 3, 100)

# Graficar el potencial
plt.plot(r, U_rc)

# Agregar línea vertical y texto en el eje x
plt.axvline(limite, color='black', linestyle='dotted')
plt.axvline(delta, color='black', linestyle='dotted')

plt.text(2.3, -2, r'$\Delta_{ij} + 2^{1/6}\sigma$', ha='center', va='center', fontsize=12)
plt.text(1.05, -2, r'$\Delta_{ij}$', ha='center', va='center', fontsize=12)
plt.text(1.77, 10, r'$U^{rc}(r_{ij})$', ha='center', va='center', fontsize=12)


# Configurar el tamaño de la figura
fig = plt.gcf()
fig.set_size_inches(8, 6)

# Configurar los límites de los ejes
plt.xlim([0.8, 3])
plt.ylim([-5, 20])

y0 = 0*d
plt.plot(d, y0, color='black')

# Configurar los nombres de los ejes
plt.xlabel("$r_{ij}$", fontsize=15)
plt.ylabel("$U^{rc}(r_{ij})$", fontsize=15)


plt.xticks([])
plt.yticks([])

# Agregar la marca del 0 al eje y
yticks = plt.gca().get_yticks().tolist()
if 0 not in yticks:
    yticks.append(0)
plt.gca().set_yticks(yticks)


# Mostrar la gráfica
plt.show()
