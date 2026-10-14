# Grafica varias salidas de un codigo que resuelve una sola ecuación de
# conservación, como la ec. de advección o la de Burgers, de una por una,
# usando el modo "interactivo" de matplotlib.

import sys
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------------------------

# Especificar los números mínimo y máximo de salida a graficar
# Supone que los archivos de datos se llaman "Lax_xx.txt", donde
# xx es el número de salida (con ceros de relleno a la izquierda)
min_out = 0
max_out = 100

# Plantilla para los nombrs de archivo
fname_base = fname = "Lax_{:02d}.txt"

# ------------------------------------------------------------------------------

# Crea una figura nueva -- esto lo hacemos una sola vez
plt.figure(figsize=(8,6))

# Activar modo "interactivo" (que nos permite alterar una figura)
plt.ion()

for nout in range(min_out, max_out+1):

  # Generar nombre de archivo de datos
  fname = fname_base.format(nout)

  # Cargamos con numpy los datos, separamos las columnas
  A = np.loadtxt(fname)
  x = A[:, 0]
  var = A[:, 1]

  # Comentar esta instrucción se quieren graficar todas las salidas
  # encimadas en el mismo plot
  plt.clf()

  # Grafica var vs. x usando el línea continua y puntos pequeños
  # en los valores, en color azul
  plt.plot(x, var, ls="-", marker=".", color="blue")

  # Título y etiquetas de los ejes
  plt.title(fname)
  plt.xlabel("Posición x")
  plt.ylabel("Variable dependiente u")

  # Límites horizontales de la gráfica
  plt.xlim(0, 1)

  # Fijamos los límites verticales a los de la condición inicial
  if nout == min_out:
    ylims = plt.ylim()
  else:
    plt.ylim(ylims)
  # plt.ylim(-0.05, 2.05)

  # Ajusta los márgenes
  plt.tight_layout()

  # Opcional: guardar la figura a disco:
  #plt.savefig(fname.replace(".txt", ".png"))

  # Agregar una pausa
  plt.pause(0.1)
