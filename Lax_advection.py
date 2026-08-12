# ==============================================================================
# Lax_advection.py -- VERSION INCOMPLETA
#
# por J.C. Toledo-Roy
# 23 / mar / 2021
#
# Este programa resuelve la ecuacion lineal de advección empleando el metodo de
# diferencias finitas de Lax.
#
# El problema a resolver consiste en la condicion inicial:
#  u = 2 para x <= 0.2
#   u = 1 para x > 0.2
# Las condiciones de frontera son:
#   izquierda: inflow con valor u = 2,
#   derecha: salida libre
# Se utiliza 1 celda fantasma en cada extremo del arreglo de datos.
# ==============================================================================

import time

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

  # ************** #
  # TU CÓDIGO AQUI #
  # ************** #

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
  fname = "temp/Lax_{:02d}.txt".format(nout)

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

  # ************** #
  # TU CÓDIGO AQUI #
  # ************** #
  pass

# ==============================================================================

# Calcular los flujos físicos F
def fluxes(U, F):

  # ************** #
  # TU CÓDIGO AQUI #
  # ************** #
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

  # ************** #
  # TU CÓDIGO AQUI #
  # ************** #
  pass

# ==============================================================================

# Hace un paso de tiempo, volcando las UPs sobre las Us y avanzando variables
def step(U, UP):

  global t, it

  for i in range(1, NX+1):
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