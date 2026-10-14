# ==============================================================================

# Parámetros, condición inicial y condición de frontera para simular un
# jet astrofísico

# ==============================================================================

import numpy as np

# Constantes astrofísicas y fundamentales
AMU  = 1.660538782e-24     # Unidad de masa atómica (g)
KB   = 1.380650400e-16     # Constante de Boltzmann (erg/K)
PC   = 3.085677588e+18     # Parsec (cm)
AU   = 1.495978707e+13     # Unidad astronómica (cm)
YR   = 3.155673600e+7      # Año sideral terrestre (s)
KYR  = 3.155673600e+10     # Mil años en segundos (s)
MSUN = 1.988920000e+33     # Masa solar (g)
KPS  = 1.0e5               # Un km/s en cm/s

# ==============================================================================
# PARÁMETROS DE LA SIMULACIÓN
# ==============================================================================
# Todas las cantidades físicas en unidades cgs

# Parámetros de la simulación
NEQ = 4              # Número de ecuaciones
NX = 400             # Número de celdas en X
NY = 150             # Número de celdas en Y

X1 = 0               # Coordenada física del extremo izquierdo
X2 = 4000*AU         # Coordenada física del extremo derecho
Y1 = -750*AU               # Coordenada física del extremo inferior
Y2 = 750*AU         # Coordenada física del extremo superior

TFIN = 100*YR        # Tiempo final de integración
DTOUT = 1*YR         # Intervalo para escribir a disco

CFL = 0.1            # Parametro de Courant
ETA = 5e-3           # Viscosidad artficial (para Macormack)

GAMMA = 5.0/3.0      # Razón de capacidades caloríficas
mu0 = 1.3            # Masa promedio por partícula del gas neutro
mui = 0.61           # Masa promedio por partícula del gas ionizado
ion_thresh = 1e4;    # Temperatura de ionización del gas

DX = (X2-X1)/(NX+1)      # Espaciamiento en X de la malla
DY = (Y2-Y1)/(NY+1)      # Espaciamiento en Y de la malla

# Solucionador numérico
# SOLVER = "SOLVER_LAX"
SOLVER = "SOLVER_MACORMACK"

# Enfriamiento radiativo
USE_COOLING = 1

# Condiciones de frontera en cada lado de la caja computacional
BC_LEFT   = "BC_FREEFLOW"
BC_RIGHT  = "BC_FREEFLOW"
BC_BOTTOM = "BC_FREEFLOW"
BC_TOP    = "BC_FREEFLOW"

# Directorio donde escribir las salidas (usar "./" para dir actual)
OUT_DIR = "Proyecto_2/"

# La "plantilla" para los nombres de archivos de salida, incluyendo
# el formato para el número de salida
# Por ejemplo, "output_%02i.txt" usará dos cifras (con un 0 si necesario)
OUT_FNAME = "output_%04i.txt"

# Formato de salidas
OUTPUT_FORMAT = "FORMAT_BINARY"

# Otras constantes
AXISYMMETRIC = False

# Arreglos de coordenadas de las celdas
coords_X = X1 + (np.arange(NX+2)+0.5)*DX
coords_Y = Y1 + (np.arange(NY+2)+0.5)*DY
XX, YY = np.meshgrid(coords_X, coords_Y, indexing='ij')

# ==============================================================================
# PARÁMETROS DEL JET
# ==============================================================================
# Todas las cantidades físicas en unidades cgs

# Parámetros del medio interestelar
n_ism = 1.0e3       # Densidad numérica
T_ism = 1000.0      # Temperatura
u_ism = v_ism = 0.0         # Componentes de la velocidad
mu_ism = 1.3        # Masa promedio por partícula [amu]

# Parámetros del jet
n_jet = 1.0e6       # Densidad numérica
T_jet = 1000        # Temperatura
v0_jet = 150*KPS    # Velocidad promedio
tau_jet = 40*YR     # Periodo de variabilidad
y_jet = 0           # Posición vertical de la región de inyección
r_jet = 25*AU       # Ancho de la región de inyección
l_jet = 30*AU       # Largo de la región de inyección
dv_jet = 0.15       # Porcentaje de variación de la velocidad

# Parámetros calculados del jet
omega = 2*np.pi/tau_jet
rho_jet = n_jet * mu_ism * AMU
P_jet = n_jet * KB * T_jet
y1_jet = y_jet - r_jet
y2_jet = y_jet + r_jet

# ==============================================================================
# CONDICIÓN INICIAL
# ==============================================================================

# La condición inicial simplemente impone el medio interestelar
# El jet es impuesto por la condición de frontera
def IC_custom(UU):

  rho_ism = n_ism * mu_ism * AMU
  P_ism = n_ism * KB * T_ism
  E_ism = 0.5*rho_ism*(u_ism**2 + v_ism**2) + P_ism/(GAMMA-1)

  # Llenamos las Us con medio interestelar en todo el dominio
  UU[0] = rho_ism
  UU[1] = rho_ism * u_ism
  UU[2] = rho_ism * v_ism
  UU[3] = E_ism

# ==============================================================================
# CONDICIÓN DE FRONTERA CUSTOM
# ==============================================================================

# Impone el jet dentro del área de inyeción en cada paso de tiempo
def BC_custom(UU, t):

  vj = v0_jet * (1 + dv_jet*np.cos(omega*t))

  inside = (XX <= l_jet) & (YY >= y1_jet) & (YY <= y2_jet)

  UU[0][inside] = rho_jet
  UU[1][inside] = rho_jet * vj
  UU[2][inside] = 0
  UU[3][inside] = 0.5*rho_jet*vj**2 + P_jet/(GAMMA-1)

# ------------------------------------------------------------------------------