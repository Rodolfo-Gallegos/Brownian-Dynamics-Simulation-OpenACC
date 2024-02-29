import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Definir la función de costo
def cost_function(x, y):
    return x**2 + y**2

# Definir la derivada parcial de la función de costo con respecto a x
def partial_x(x, y):
    return 2*x

# Definir la derivada parcial de la función de costo con respecto a y
def partial_y(x, y):
    return 2*y

# Definir el tamaño del paso del descenso de gradiente
alpha = 0.1

# Definir el punto inicial
x0, y0 = 5, 5

# Crear una figura 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Crear una malla de puntos en los ejes x e y
x_range = np.linspace(-5, 5, 50)
y_range = np.linspace(-5, 5, 50)
X, Y = np.meshgrid(x_range, y_range)

# Calcular el valor de la función de costo en cada punto de la malla
Z = cost_function(X, Y)

# Graficar la función de costo en 3D
ax.plot_surface(X, Y, Z, cmap='coolwarm')

# Definir la trayectoria del descenso de gradiente
x = [x0]
y = [y0]
for i in range(100):
    # Calcular la dirección del descenso de gradiente
    dx = -alpha * partial_x(x[-1], y[-1])
    dy = -alpha * partial_y(x[-1], y[-1])
    
    # Calcular el nuevo punto en la trayectoria
    x_new = x[-1] + dx
    y_new = y[-1] + dy
    
    # Agregar el nuevo punto a la trayectoria
    x.append(x_new)
    y.append(y_new)


# Configurar los límites de los ejes
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])
ax.set_zlim([0, 50])

# Configurar los nombres de los ejes
ax.w_xaxis.set_ticklabels([])
ax.w_yaxis.set_ticklabels([])
ax.w_zaxis.set_ticklabels([])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_zlabel('')

# Mostrar la gráfica
plt.show()
