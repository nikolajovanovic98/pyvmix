import numpy as np
import matplotlib.pyplot as plt
import pyicon as pyic
from netCDF4 import Dataset
import sys
import xarray as xr

# Initialize default parameters
# -----------------------------
exec(open('./pyvmix_defaults.py').read())

# Modify parameters
# -----------------
run = __file__.split('/')[-1][:-3]
path_data = f'/Users/nbruegge/work/pyvmix/{run}/'

savefig = False
path_fig = '../pics_nac/'
nnf=0

nz = 500
dz = 0.5*np.ones(nz)

fcor = 2*np.pi/86400 * np.sin(np.pi/180. * 45)

c_k = 0.1

kvAv = 'tke'

deltaT = 1800.
nt = 5*86400//deltaT
lsave = 3*3600//deltaT

#d_ek = np.sqrt(2.*Av0/fcor)
#T_ek = np.sqrt( (taux0/fcor)**2 + (tauy0/fcor)**2 )
#tke_min = 1e-33

# Initialize grid and variables
# -----------------------------
exec(open('./pyvmix_initialization.py').read())

# Modify initial conditions
# -------------------------
#b = 5e-3*np.exp(zt/200.)
b = 0.*zt
#b = 1e-4*zt/100.

# Modify forcing
# --------------
# wind forcing
#taux0 = 2e-5
taux0 = 1.5e-4
Qsurf0 = 500.

b0 = 1.*b
lam_b = 0.*1./(1.*86400.)

#ylim_hov = [-25, 0]

def wind_forcing(time):
  taux = taux0
  tauy = 0.
  return taux, tauy

def buoyancy_forcing(time):
  #Qsurf = Qsurf0
  Qsurf = Qsurf0*np.sin(2*np.pi*time/86400.)
  if Qsurf<-250.:
    Qsurf = -250.
  return Qsurf

# Run the model
# -------------
exec(open('./pyvmix_main.py').read())

# Visualize results
# -----------------
plt.close('all')

exec(open('./plot_defaults.py').read())

plt.show()
