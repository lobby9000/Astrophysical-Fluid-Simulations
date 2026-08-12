# ==============================================================================
# Lax_advection.py -- VERSION COMPLETA
#
# por J.C. Toledo-Roy
# 23 / mar / 2021
#
# Este programa resuelve la ecuacion lineal de advección empleando el metodo de
# diferencias finitas de Lax.
#
# El problema a resolver consiste en la condicion inicial:
#  u = 1
# Las condiciones de frontera oscilatorias:
#   izquierda: 1+sin(2pit)/0.5
#   derecha: salida libre
# Se utiliza 1 celda fantasma en cada extremo del arreglo de datos.
# ==============================================================================
import time
import numpy as np

# ==============================================================================

# CONSTANTES Y VARIABLES GLOBALES

# Parámetros constantes de la simulación
NX = 500         # Tamaño de la malla
XL = 0.0         # Coordenada física del extremo izquierdo
XR = 1.0         # Coordenada física del extremo derecho
TFIN = 1.0       # Tiempo final de integración
CFL = 0.9        # Parametro de Courant
DTOUT = 0.05     # Intervalo de tiempo para escribir a disco

# Para tubo de choque: coord de la separación entre estados iniciales
X0 = 0.2

# Para ecuación de advección: velocidad de ondas ('a' en las notas)
A = 1.0

# Constantes derivadas de las anteriores
DX = (XR-XL)/NX      # Espaciamiento de la malla

# Variables globales (sus valors serán asignados en initflow)
U  = [0]*(NX+2)    # Variables conservadas actuales
# Alternativa: U = np.zeros(NX+2)
# La diferencia: La primera es una lista, la segunda es un arrgelo de numpy
# que soportan operaciones matemáticas.
UP = [0]*(NX+2)    # Variables conservadas "avanzadas"
F  = [0]*(NX+2)    # Flujos físicos
dt = 0             # Paso de tiempo
t = 0              # Tiempo actual
it = 0             # Iteración actual
start_time = 0     # Tiempo de inicio
next_tout = 0       # Tiempo para el siguiente output
nout = 0           # Número de salida

# ==============================================================================

# Impone las condiciones iniciales
def initflow(U):

  global t, it, next_tout, nout

  # Inicializar los valores de U en todo el dominio
  # Nótese que también llenamos las celdas fantasma
  
  for i in range(0, NX+2):

   U[i]= 1

  # Inicializar otras variables
  t = 0
  it = 0
  nout = 0
  next_tout = 0

# ==============================================================================

# Escribe a disco el estado de la simulación
def output(U):

  global nout, next_tout
  
  # Generar el nombre del archivo de salida
  fname = "temp_4/Lax_{:02d}.txt".format(nout)

  # Abrir el archivo para escritura
  fout = open(fname, "w")

  # Escribir los valores de U al archivo (sólo celdas físicas)
  for i in range(1, NX+1):
    x = XL + i*DX
    fout.write("{} {}\n".format(x, U[i]))

  # Cerrar archivo
  fout.close()

  print("Se escribió salida", nout)

  # Avanzar variables de output
  nout = nout + 1
  next_tout = DTOUT * nout

# ==============================================================================

# Aplicar condiciones de frontera a celdas fantasma
# El arreglo pasado es al que aplicaremos las BCs
def boundary(U):

  # Asignamos los valores de U[0] y U[NX+1], según las BCs deseadas
  
  # En este caso, al querérlas periódicas la izquierda debe ser la de U[nx] y
  # la derecha, la de U[1]:

  U[0] = 1+np.sin((2*np.pi*t)/0.5)
  U[NX+1] = U[NX]
  
  pass

# ==============================================================================

# Calcular los flujos físicos F
def fluxes(U, F):
  for i in range(0, NX+2):
     F[i]=U[i]
        
  pass

# ==============================================================================

# Calcula el paso de tiempo resultante de la condición CFL
def timestep(U):

  # Para la ec. de advección, max_u es simplemente A
  max_speed = abs(A)

  # Para otros casos, debemos calcular el valor máximo de |velocidad|
  # max_speed = 0
  # for i in range(1, NX+1):
  #   if abs(U[i]) > max_speed:
  #     max_speed = abs(U[i])
  
  dt = CFL * DX / max_speed

  return dt

# ==============================================================================

# Aplica el método de Lax para obtener las UP a partir de las U
# Supone que los flujos F ya fueron actualizados
def Lax(U, F, UP):

  for i in range(1, NX+1):
      
      UP[i] = (U[i+1]+U[i-1])/2-(dt/(2*DX))*(F[i+1]-F[i-1])
  
  pass

# ==============================================================================

# Hace un paso de tiempo, volcando las UPs sobre las Us y avanzando variables
def step(U, UP):

  global t, it

  for i in range(0, NX+2):
    U[i] = UP[i]

  t = t + dt
  it = it + 1

# ==============================================================================

# PROGRAMA PRINCIPAL

# Condición inicial e inicializaciones
initflow(U)

# Escribir condición inicial a disco
output(U)

# Tiempo de inicio de la simulación
start_time = time.time()

# Bucle temporal
while (t <= TFIN):

  # Actualizar el paso de tiempo
  dt = timestep(U)

  # Actualizar flujos físicos
  fluxes(U, F)

  # Aplicar método de Lax para actualizar las UP
  Lax(U, F, UP)

  # Aplicar condiciones de frontera a las UP
  boundary(UP)

  # Avanzar en el tiempo
  step(U, UP)

  # Escribir a disco
  if (t >= next_tout):
    output(U)

# Terminar
elapsed = time.time() - start_time
print("\nSe calcularon {} iteraciones {:.1f} s\n".format(it, elapsed))

# ==============================================================================

# Grafica varias salidas de un codigo que resuelve una sola ecuación de
# conservación, como la ec. de advección o la de Burgers, de una por una,
# usando el modo "interactivo" de matplotlib.

import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------

# Especificar los números mínimo y máximo de salida a graficar
# Supone que los archivos de datos se llaman "Lax_xx.txt", donde
# xx es el número de salida (con ceros de relleno a la izquierda)
min_out = 0
max_out = 100

# Plantilla para los nombrs de archivo
fname_base = fname = "temp_4/Lax_{:02d}.txt"

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
  #if nout == min_out:
    #ylims = plt.ylim()
  #else:
    #plt.ylim(ylims)
  plt.ylim(-0.05, 2.05)

  # Ajusta los márgenes
  plt.tight_layout()

  # Opcional: guardar la figura a disco:
  plt.savefig(fname.replace(".txt", ".png"))

  # Agregar una pausa
  plt.pause(0.1)

