#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 12:29:40 2023

@author: looby9000
"""


# ==============================================================================
# Ecuaciones de Euler 
# Este programa resuelve la ecuacion lineal de advección empleando el metodo de
# Macormack o Lax para condiciones iniciales y de frontera dadas.
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
TFIN = 0.2      # Tiempo final de integración
CFL = 0.9        # Parametro de Courant
DTOUT = TFIN/20     # Intervalo de tiempo para escribir a disco
GAMMA = 1.4      # Parámetro del gas.


ETA = 0.37


# Para tubo de choque: coord de la separación entre estados iniciales
X0 = 0.5

# Para ecuación de advección: velocidad de ondas ('a' en las notas)
A = 1.0

# Constantes derivadas de las anteriores
DX = (XR-XL)/NX      # Espaciamiento de la malla



# Variables globales (sus valores serán asignados en initflow)

NEQ = 3                  # Dimensión extra, número de ecuaciones. 
PRIM = np.zeros((NX+2,NEQ))  # Arreglo de primitivas.
U  = np.zeros((NX+2,NEQ))    # Variables conservadas actuales
UP = np.zeros((NX+2,NEQ))    # Variables conservadas "avanzadas"
UT = np.zeros((NX+2,NEQ))    # Variable para Macormack.
F  = np.zeros((NX+2,NEQ))    # Flujos físicos


dt = 0             # Paso de tiempo
t = 0              # Tiempo actual
it = 0             # Iteración actual
start_time = 0     # Tiempo de inicio
next_tout = 0       # Tiempo para el siguiente output
nout = 0           # Número de salida

# ==============================================================================

# Impone las condiciones iniciales
def initflow(U, rhol, ul, pl, rhor, ur, pr):

  global t, it, next_tout, nout

  u1l = rhol
  u2l = rhol*ul
  u3l = 0.5*rhol*(ul**2)+pl/(GAMMA-1)
  u1r = rhor
  u2r = rhor*ur
  u3r = 0.5*rhor*(ur**2)+pr/(GAMMA-1)
  
  # Inicializar los valores de U en todo el dominio
  # Nótese que también llenamos las celdas fantasma

  for i in range(0,NX+2):
          
          if i*DX <= X0:
              U[i][0] = u1l
              U[i][1] = u2l
              U[i][2] = u3l
             
              
          if i*DX > X0:
              U[i][0] = u1r
              U[i][1] = u2r
              U[i][2] = u3r
              
          
  # Inicializar otras variables
  t = 0
  it = 0
  nout = 0
  next_tout = 0

# ==============================================================================


# ==============================================================================
def update_prim(U, PRIM):
    
    for i in range(0, NX+2):
        PRIM[i,0] = U[i,0]
        PRIM[i,1] = U[i,1]/U[i,0]
        PRIM[i,2] = (GAMMA-1)*(U[i,2]-0.5*(U[i,1]*U[i,1])/(U[i,0]))


# ==============================================================================
# Escribe a disco el estado de la simulación
def output(PRIM):

  global nout, next_tout
  
  # Generar el nombre del archivo de salida
  fname = "eulerMac_1/Mac_{:02d}.txt".format(nout)

  # Abrir el archivo para escritura
  fout = open(fname, "w")

  # Escribir los valores de U al archivo (sólo celdas físicas)
  for i in range(1, NX+1):
   x = XL + i*DX
   fout.write("{:.5f} {:.5f} {:.5f} {:.5f}\n".format(x, PRIM[i,0], PRIM[i,1], PRIM[i,2]))

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

  for e in range(0,NEQ):
        U[0][e] = U[1][e]
        U[NX+1][e] = U[NX][e]
        
  pass



# ==============================================================================

# Calcular los flujos físicos de Euler 
def fluxes_Euler(PRIM, F):
  for i in range(0, NX+2):
     F[i,0] = PRIM[i,0]*PRIM[i,1]
     F[i,1] = PRIM[i,0]*PRIM[i,1]**2 + PRIM[i,2]
     F[i,2] = PRIM[i,1]*(0.5*PRIM[i,0]*PRIM[i,1]**2+(GAMMA/(GAMMA-1))*PRIM[i,2])
        
  pass

# ==============================================================================

# Calcula el paso de tiempo resultante de la condición CFL
def timestep(PRIM):
    
  max_speed = 0
  for i in range(1, NX+1): 
      cs = np.sqrt(GAMMA*PRIM[i,2]/PRIM[i,0])
      if abs(PRIM[i,1])+cs> max_speed:
       max_speed = abs(PRIM[i,1])+cs
  
  dt = CFL * DX / max_speed

  return dt

# ==============================================================================

def Macormack(UT, U, F, UP):
    
    for i in range(1,NX+1):
        
        for e in range(0,NEQ):
            
         UT[i,e] = U[i,e]- (dt/DX)*(F[i+1,e]-F[i,e])
            
    boundary(UT)
    update_prim(UT, PRIM)
    fluxes_Euler(PRIM, F)
        
    for i in range(1, NX+1):
        
        for e in range(0,NEQ):
            UP[i,e] = (U[i,e]+UT[i,e])/2-(dt/(2*DX))*(F[i,e]-F[i-1,e])
    
    pass

# ==============================================================================

# Hace un paso de tiempo, volcando las UPs sobre las Us y avanzando variables
def step_viscoso(U, UP):

  global t, it
  
  for e in range(0, NEQ):
      U[0,e] = UP[0,e]
      U[NX+1,e] = UP[NX+1,e]
  
  for i in range(0, NX+1):
      for e in range(0, NEQ):
          if (UP[i+1,e]-UP[i,e])*(UP[i,e]-UP[i-1,e])<0:
              U[i,e] = UP[i,e] + ETA*(UP[i+1,e]+UP[i-1,e]-2*UP[i,e])
          else:
              U[i,e]=UP[i,e]

  t = t + dt
  it = it + 1
  
# ==============================================================================

# PROGRAMA PRINCIPAL

# Condición inicial e inicializaciones
initflow(U, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1)
# Actualizar Primitivas
  
update_prim(U, PRIM)

# Escribir condición inicial a disco
output(PRIM)

# Tiempo de inicio de la simulación
start_time = time.time()

# Bucle temporal
while (t <= TFIN):
    

  # Actualizar el paso de tiempo
  dt = timestep(PRIM)

  # Actualizar flujos físicos
  fluxes_Euler(PRIM, F)

  # Aplicar método de Lax para actualizar las UP
  Macormack(UT, U, F, UP)

  # Aplicar condiciones de frontera a las UP
  boundary(UP)

  # Avanzar en el tiempo
  step_viscoso(U, UP)
  
  # Actualizar Primitivas
    
  update_prim(U, PRIM)

  # Escribir a disco
  if (t >= next_tout):
    output(PRIM)


# Terminar
elapsed = time.time() - start_time
print("\nSe calcularon {} iteraciones {:.1f} s\n".format(it, elapsed))

# ==============================================================================

# Este programa grafica la salida del código que resuelve las ecuaciones
# de Euler en 1D, en formato ASCII.

# En el archivo de datos, en formato ASCII, cada fila debe tener 4 columnas
# separadas por espacios:
# coord x
# densidad
# velocidad
# presión

# ==============================================================================

import matplotlib.pyplot as plt
import numpy as np

# ==============================================================================

# Nombre del archivo de datos a leer (en formato ASCII)
fname = "eulerMac_1/Mac_20.txt"


# Leemos los datos ASCII del archivo que se cargan como una matriz que
# tiene 4 columnas (posición x, densidad, velocidad, presión) y tantas
# filas como celdas en la simulación.
datos = np.loadtxt(fname)

# Para más comodido le damos definimos variables individuales para las
# columnas
x = datos[:, 0]
dens = datos[:, 1]
vel = datos[:, 2]
pres = datos[:, 3]

# También calculamos la temperatura; nótese que pres y dens son vectores
# (arreglos de numpy), así que esta operación se realiza sobre todos los
# elementos en todo el vector
temp = pres/dens

# Crear nueva figura de tamaño 10 x 8 "pulgadas"
plt.figure(figsize=(10,8))

# Con subplot(filas, columnas, posición) especificamos que queremos un arreglo
# de sub-gráficas con ese número de filas y columnas, y que actualmente queremos
# graficar en esa posición (la numeración empieza en 1 y corre de izquierda a
# derecha, y de arriba hacia abajo)
plt.subplot(2,2,1)
# El comando plot() recibe un vector con los valores de x y otro con los valores
# de y como los primeros dos argumentos. Los demás argumentos son opciones,
# como el estilo de los marcadores y el color.
plt.plot(x, dens, marker=".", color="blue")
# Título de esta (sub-)gráfica
plt.title("Densidad")
# Rango del eje x que queremos mostrar
plt.xlim(0, 1)
# Rango del eje y que queremos mostrar
plt.ylim(0, 1.05)
# Cuadrícula
plt.grid(ls=":")

# Los siguientes 3 grupos grafican diferentes variables en subplots distintos

plt.subplot(2,2,2)
plt.plot(x, vel, marker=".", color="green")
plt.title("Velocidad")
plt.xlim(0, 1)
plt.ylim(-0.05, 1.2)
plt.grid(ls=":")

plt.subplot(2,2,3)
plt.plot(x, pres, marker=".", color="darkorange")
plt.title("Presión")
plt.xlim(0, 1)
plt.ylim(0, 1.2)
plt.grid(ls=":")

plt.subplot(2,2,4)
plt.plot(x, temp, marker=".", color="red")
plt.title("Temperatura")
plt.xlim(0, 1)
plt.ylim(0.7, 1.2)
plt.grid(ls=":")

# Ajusta los márgenes entre subplots para que queden más estrechos pero sin
# sobreponer subplots
plt.tight_layout()

# Descomentar para guardar la gráfica
plt.savefig("Shock1M_20.png")

# Muestra la figura en pantalla
plt.show()