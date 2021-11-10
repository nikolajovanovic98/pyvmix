# Vertical grid
# -------------
zu = np.concatenate(([0],-dz.cumsum()))
zt = 0.5*(zu[1:]+zu[:-1])
dzt = np.concatenate(([-zt[0]],zt[:-1]-zt[1:],[-zu[-1]+zt[-1]]))

# Allocate initial variables
# --------------------------
uvel = np.zeros((nz))
vvel = np.zeros((nz))
b = np.zeros((nz))
tke = np.zeros((nz+1))
kv = np.zeros((nz+1))
Av = np.zeros((nz+1))

# --- mixing profile (for kvAv = 'profile')
kv_prof = np.zeros((nz+1))
Av_prof = np.zeros((nz+1))

Tu_wnd = np.zeros((nz))
Tv_wnd = np.zeros((nz))
Tu_bot = np.zeros((nz))
Tv_bot = np.zeros((nz))
Tb_sfl = np.zeros((nz))

# Allocate forcing variables
# --------------------------
wvel = np.zeros((nz+1))
dpdx = np.zeros((nz))
dpdy = np.zeros((nz))

def wind_forcing(time):
  taux = 0.
  tauy = 0.
  return taux, tauy

def buoyancy_forcing(time):
  Qsurf = 0.
  return Qsurf

# Allocate output variables
# -------------------------
nsave = int(nt/lsave)

# --- time series of profiles
uvel_s = np.zeros((nsave,nz))
vvel_s = np.zeros((nsave,nz))
b_s = np.zeros((nsave,nz))
tke_s = np.zeros((nsave,nz+1))
kv_s = np.zeros((nsave,nz+1))
Av_s = np.zeros((nsave,nz+1))
Lmix_s = np.zeros((nsave,nz+1))
N2_s = np.zeros((nsave,nz+1))

# ---- time series vint MKE equation
Tke_cor = np.zeros((nsave)) 
Tke_hpr = np.zeros((nsave)) 
Tke_wnd = np.zeros((nsave)) 
Tke_bot = np.zeros((nsave)) 
Tke_vdf = np.zeros((nsave)) 
Tke_vds = np.zeros((nsave)) 
Tke_vfl = np.zeros((nsave)) 
Tke_tot = np.zeros((nsave)) 

# --- time series vint TKE equation
Ttke_tot = np.zeros((nsave)) 
Ttke_bpr = np.zeros((nsave)) 
Ttke_spr = np.zeros((nsave)) 
Ttke_dis = np.zeros((nsave)) 
Ttke_vdf = np.zeros((nsave)) 
Ttke_bck = np.zeros((nsave)) 

# --- time series of the forcing
Qsurf_ts = np.zeros((nsave))
taux_ts = np.zeros((nsave))
tauy_ts = np.zeros((nsave))
