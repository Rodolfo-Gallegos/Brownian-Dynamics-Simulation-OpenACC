import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Definimos los parámetros
num_points = 100   # Número de puntos a generar
cube_size = 19.0  # Tamaño del cubo
arrow_scale = 1 # Escala de las flechas
point_size = 50   # Tamaño de los puntos

# Generamos una matriz de posiciones aleatorias dentro del cubo
positions = np.random.uniform(-cube_size/2, cube_size/2, size=(num_points, 3))

# Generamos una matriz de direcciones aleatorias
directions = np.random.randn(num_points, 3)

# Normalizamos las direcciones
directions /= np.linalg.norm(directions, axis=1)[:, np.newaxis]

# Graficamos los puntos y las flechas
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-cube_size/2, cube_size/2)
ax.set_ylim(-cube_size/2, cube_size/2)
ax.set_zlim(-cube_size/2, cube_size/2)
ax.set_axis_off()  # Oculta los ejes de coordenadas
ax.set_facecolor('white')  # Cambia el color de fondo a azul
ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], s=point_size, c='black')  # Cambia el color de los puntos a negro
for i in range(num_points):
    ax.quiver(positions[i, 0], positions[i, 1], positions[i, 2],
              directions[i, 0], directions[i, 1], directions[i, 2],
              length=arrow_scale, normalize=True, color='gray', arrow_length_ratio=0.5)  # Cambia el color de las flechas a gris
# Agregamos las aristas del cubo
vertices = np.array([[-10, -10, -10],
                     [-10, -10,  10],
                     [-10,  10,  10],
                     [-10,  10, -10],
                     [ 10, -10, -10],
                     [ 10, -10,  10],
                     [ 10,  10,  10],
                     [ 10,  10, -10]])
for i in range(4):
    ax.plot(*zip(vertices[i], vertices[(i+1)%4]), color='black')
    ax.plot(*zip(vertices[i+4], vertices[((i+1)%4)+4]), color='black')
    ax.plot(*zip(vertices[i], vertices[i+4]), color='black')
plt.show()
