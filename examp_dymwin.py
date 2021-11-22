import numpy as np
import matplotlib.pyplot as plt
import pyicon as pyic
from netCDF4 import Dataset
import sys
import xarray as xr

import pyvmix

# Initialize settings
# -------------------
S = pyvmix.Model()

# Modify default settings
# -----------------------
S.run = __file__.split('/')[-1][:-3]
S.path_data = f'/Users/nbruegge/work/pyvmix/{S.run}/'

S.savefig = False
S.path_fig = '../pics_nac/'
S.nnf=0

S.nz = 500
S.dz = 0.5*np.ones(S.nz)
#S.nz = 125
#S.dz = 2.*np.ones(S.nz)

S.fcor = 2*np.pi/86400 * np.sin(np.pi/180. * 45)

S.cu = 0.1
S.tke_min = 1e-6

S.kvAv = 'tke'

S.deltaT = 1800.
S.nt = 5*86400//S.deltaT
S.lsave = 1*3600//S.deltaT

# Initialize grid and variables
# -----------------------------
S.initialize()

# Modify initial conditions
# -------------------------
#S.b = 5e-3*np.exp(S.zt/200.)
S.b = 0.*S.zt
#S.b = 1e-4*S.zt/100.

# Modify forcing
# --------------
# wind forcing
#S.taux0 = 1.5e-4
S.Qsurf0 = 800.
S.uair = 4.
S.vair = 0.

S.b0 = 1.*S.b
S.lam_b = 0.*1./(1.*86400.)

S.ylim_hov = [-25, 0]

def wind_forcing(S, time):
  #taux = S.taux0
  #tauy = 0.
  usurf = S.uvel[0]
  vsurf = S.vvel[0]
  taux = (
    S.rho_air * S.cdrag / S.rho0
    * np.sqrt((S.uair-usurf)**2+(S.vair-vsurf)**2)
    * (S.uair-usurf)
    )
  tauy = (
    S.rho_air * S.cdrag / S.rho0
    * np.sqrt((S.uair-usurf)**2+(S.vair-vsurf)**2)
    * (S.vair-vsurf)
    )
  return taux, tauy
S.wind_forcing = wind_forcing

def buoyancy_forcing(S, time):
  #Qsurf = S.Qsurf0
  Qsurf = S.Qsurf0*np.sin(2*np.pi*time/86400.)
  Qcut = -100.
  if Qsurf<Qcut:
    Qsurf = Qcut
  return Qsurf
S.buoyancy_forcing = buoyancy_forcing

# Run the model
# -------------
S.run_model()

# Visualize results
# -----------------
plt.close('all')

exec(open('./plot_defaults.py').read())

plt.show()
