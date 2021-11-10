
# Settings 
# --------
# mixing type: choose between 'tke' and 'profile'
kvAv = 'tke'

# --- convective adjustment
conv_adj = False
kv_conv = 1e-1

# --- background mixing
kv_back = 0.
Av_back = 0.

# Physical constants
# ------------------
grav = 9.81
tAlpha = 2e-4
rho0 = 1024.
cp = 4.18e3
#fcor = 1e-4
fcor = 2.*np.pi/86400.

# Time stepping
# --------------
nt = 24*2 * 20
deltaT = 1800.
lsave = 2*3

# TKE parameter
# -------------
cu = 0.1
cd = 3.75
alpha = 30.
ceps = 0.7
tke_min = 1e-6
Lmix_min = 1e-8
# bottom drag 0.001-0.003
bottomDragQuadratic = 0.*0.002

# Default plotting
# ----------------
ylim_hov = [None, None]

