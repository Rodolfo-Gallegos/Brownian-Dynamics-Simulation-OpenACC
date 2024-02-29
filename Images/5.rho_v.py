import matplotlib.pyplot as plt
import numpy as np

# Definir el radio y el incremento de radio
r = 5
dr = 1

# Definir el centro de los círculos
center = (0, 0)

# Crear una figura
fig, ax = plt.subplots()


# Graficar el círculo interior
circle1 = plt.Circle(center, r-4*dr, color='black', fill=False, linewidth=2)
ax.add_artist(circle1)

# Graficar el círculo interior
circle1 = plt.Circle(center, r-3*dr, color='black', fill=False, linewidth=2)
ax.add_artist(circle1)

# Graficar el círculo interior
circle1 = plt.Circle(center, r-2*dr, color='black', fill=False, linewidth=2)
ax.add_artist(circle1)

# Graficar el círculo interior
circle1 = plt.Circle(center, r-dr, color='black', fill=False, linewidth=2)
ax.add_artist(circle1)

# Graficar el círculo medio
circle2 = plt.Circle(center, r, color='black', fill=False, linewidth=2)
ax.add_artist(circle2)


# Calcular el punto medio entre los dos círculos
midpoint = (center[0], center[1]+(r+dr)/2)

# Marcar el punto medio con una etiqueta "V"
# ax.annotate('V', midpoint, fontsize=20, color='green', ha='center', va='center')

# Configurar las etiquetas de los ejes
ax.set_xlabel('')
ax.set_ylabel('')

plt.plot([0, 5], [0, 0], color='black', linewidth=1)
plt.plot([0, .84], [0, .48], color='black', linewidth=1)
plt.plot([0, 3.25], [0, -2.27], color='black', linewidth=1)

plt.xlim([-7, 7])
plt.ylim([-7, 7])

plt.scatter(-1.12, -1, color='r')
plt.scatter(-1.12, 1, color='r')
plt.scatter(1.12, 1, color='r')
plt.scatter(1.12, -1, color='r')

plt.scatter(-4.5, 0, color='b')
plt.scatter(0, -3.5, color='b')
plt.scatter(0, 3.5, color='b')
plt.scatter(-3.5, 0, color='b')
plt.scatter(-3.28, -3.07, color='b')
plt.scatter(3.28, -3.07, color='b')
plt.scatter(3.28, 3.07, color='b')
plt.scatter(-2.97, 3.07, color='b')
plt.scatter(1.8, 1.8, color='b')
plt.scatter(-1.8, 1.8, color='b')
plt.scatter(-1.8, -1.8, color='b')
plt.scatter(1.8, -1.8, color='b')

plt.xticks([])
plt.yticks([])

# Mostrar los valores de r y dr en la gráfica
plt.text(.10, 0.25, '$dr$', fontsize=12, rotation=35)
plt.text(3.15, .2, '$5dr$', fontsize=12)
plt.text(1.8, -1.35, '$4dr$', fontsize=12, rotation=-32)

# Mostrar la gráfica
plt.show()
