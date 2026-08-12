# Grafica la salida de un código que resuelve una sola ecuación de
# conservación, como la ec. de advección o la de Burgers

import sys
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------------------------

# Nombre del archivo de datos a cargar, obtenido de la línea de comandos
# También se puede especificar manualmente aquí
fname = sys.argv[1]
# fname = "prob1_20.txt"

# ------------------------------------------------------------------------------

# Crea una figura nueva
plt.figure(figsize=(8,6))

# Cargamos con numpy los datos, separamos las columnas
A = np.loadtxt(fname)
x = A[:, 0]
var = A[:, 1]

# Comentar esta instrucción se quieren graficar todas las salidas
# encimadas en el mismo plot
plt.clf()

# Grafica var vs. x usando el estilo "." (puntos pequeños) y en color azul
plt.plot(x, var, "-o", color="blue")

# Título y etiquetas de los ejes
plt.title(fname)
plt.xlabel("Posición x")
plt.ylabel("Variable dependiente U")

# Ajusta los márgenes
plt.tight_layout()

# Mostrar la figura
plt.show()
