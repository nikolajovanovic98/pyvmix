import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import xarray as xr
import sys

## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver(ac, bc, cc, dc):
  '''
  TDMA solver, a b c d can be NumPy array type or Python list type.
  refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
  '''
  #nf = len(a)     # number of equations
  nf = ac.shape[0]
  #ac, bc, cc, dc = map(np.array, (a, b, c, d))     # copy the array
  for it in range(1, nf):
    mc = ac[it]/bc[it-1]
    bc[it] = bc[it] - mc*cc[it-1] 
    dc[it] = dc[it] - mc*dc[it-1]

  xc = ac
  xc[-1] = dc[-1]/bc[-1]

  for il in range(nf-2, -1, -1):
    xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]
  #del bc, cc, dc  # delete variables from memory
  return xc

nt = int(nt)
lsave = int(lsave)

ls = 0
for l in range(nt):
  time = l*deltaT
  #if (l%10)==0:
  #  print('l = %d/%d' % (l, nt))
  #taux = taux0#*(1.-np.heaviside(time-10.*86400., 0.) )
  #tauy = tauy0#*(1.-np.heaviside(time-10.*86400., 0.) )
  ##Qsurf = Qsurf0#*(1.-np.heaviside(time-10.*86400., 0.) )
  #Qsurf = Qsurf0*np.sin(2*np.pi*time/86400.)
  taux, tauy = wind_forcing(time)
  Qsurf = buoyancy_forcing(time)

  # save variables of previous time step
  uvel_old = 1.*uvel
  vvel_old = 1.*vvel
  b_old    = 1.*b
  tke_old  = 1.*tke

  # vert derivatives
  N2 = np.zeros((nz+1))
  uz = np.zeros((nz+1))
  vz = np.zeros((nz+1))
  N2[1:-1] = (b[:-1]-b[1:])/dzt[1:-1]
  uz[1:-1] = (uvel[:-1]-uvel[1:])/dzt[1:-1]
  vz[1:-1] = (vvel[:-1]-vvel[1:])/dzt[1:-1]

  # paramter
  Ri = N2/(uz+1e-33)
  Pr = np.maximum(1, np.minimum(10., 6.6*Ri))
  cb = cu/Pr

  # mixing length
  Lmix = np.sqrt(2*tke/np.abs(N2+1e-33))
  #Lmix[1:] = np.minimum(Lmix[1:], Lmix[:-1]+dzt[:-1])
  Lmix[0] = 0.
  Lmix[-1] = 0.
  for k in range(1,nz):
    Lmix[k] = np.minimum(Lmix[k], Lmix[k-1]+dzt[k-1])
  Lmix[nz] = np.minimum(Lmix[nz], Lmix_min+dzt[nz])
  for k in range(nz-1,1,-1):
    Lmix[k] = np.minimum(Lmix[k], Lmix[k+1]+dzt[k])
  Lmix = np.maximum(Lmix, Lmix_min)

  # diffusivities
  if kvAv=='tke':
    kv = cu * tke**0.5*Lmix + kv_back
    Av = cb * tke**0.5*Lmix + Av_back
  elif kvAv=='profile':
    kv = 1.*kv_prof + kv_back
    Av = 1.*Av_prof + Av_back
  if conv_adj:
    kv[N2<0] = kv_conv
  ktke = alpha*0.5*(Av[1:]+Av[:-1])

  # ==========
  # tke
  # ==========
  Tt_bpr = -kv*N2
  Tt_spr = Av*(uz**2+vz**2)

  # tke diffusion
  """
  dz K dz p = (K(k-1)*(p(k)-p(k-1))/dzw(k-1) - K(k)*(p(k+1)-p(k))/dzw(k))/dzt(k)
    with a(k) = K(k-1)/(dzw(k-1)*dzt(k))
    and  c(k) = K(k)/(dzw(k)*dzt(k))
    and  b(k) = a(k)+c(k)
  dz K dz p = -a(k) p(k-1) + b(k) p(k) - c(k) p(k+1)
  upper BC: a(ks) = 0; b(ks) = c(ks)
  lower BC: c(ks) = 0; b(ks) = a(ks)
  """
  am = np.zeros((nz+1))
  cm = np.zeros((nz+1))
  bm = np.zeros((nz+1))
  am[1:]  = ktke/(dz*dzt[1:])
  cm[:-1] = ktke/(dz*dzt[:-1])
  bm[1:-1] = am[1:-1]+cm[1:-1]
  # ----- upper and lower boundary conditions
  bm[0]    = cm[0]
  bm[-1]   = am[-1]

  implicite_tke = True
  tke_bim = 1.*tke
  if implicite_tke:
    # ----- solve tri-diagonal diffusion matrix
    atr = -deltaT*am
    btr = 1+deltaT*bm
    ctr = -deltaT*cm
    dtr = 1.*tke
    dtr += deltaT*(Tt_bpr+Tt_spr)
    # --- adding dissipation
    btr[1:-1] += deltaT*ceps*tke[1:-1]**0.5/(Lmix[1:-1]+1e-33)
    tke = TDMAsolver(1.*atr, 1.*btr, 1.*ctr, 1.*dtr)
  #else:
  #  flux = np.zeros((nz+1))
  #  flux = ktke*(tke[:-1]-tke[1:])/dz
  #  Tt_vdf = np.zeros((nz+1))
  #  Tt_vdf[1:-1] = (flux[:-1]-flux[1:])/dzt[1:-1]
  #  Tt_vdf[0]  = (0.       - flux[0])/dzt[0]
  #  Tt_vdf[-1] = (flux[-1] - 0.     )/dzt[-1]

  #dz K dz p = -a(k) p(k-1) + b(k) p(k) - c(k) p(k+1)
  Tt_vdf = np.zeros((nz+1))
  Tt_vdf[1:-1] = am[1:-1]*tke[:-2] - bm[1:-1]*tke[1:-1] + cm[1:-1]*tke[2:]
  Tt_vdf[0]  = -bm[0] *tke[0]  + cm[0] *tke[1]
  Tt_vdf[-1] =  am[-1]*tke[-2] - bm[-1]*tke[-1]
  Tt_dis = np.zeros((nz+1))
  Tt_dis[1:-1] = - ceps*tke_bim[1:-1]**0.5*tke[1:-1]/(Lmix[1:-1]+1e-33)

  Tt     = Tt_bpr + Tt_spr + Tt_dis + Tt_vdf

  # ==========
  # hor. vel.
  # ==========
  # velocity tendencies
  impliciteCoriolis = True
  if not impliciteCoriolis:
    Tu_cor = +fcor*vvel
    Tv_cor = -fcor*uvel
  else:
    """
    Implicite Coriolis term: 
    u(l+1)-u(l) =  f dt v(l+1)
    v(l+1)-v(l) = -f dt u(l+1)
    
    u(l+1) - f dt v(l+1) = u(l)
    v(l+1) + f dt u(l+1) = v(l)
    
    |1     - f dt | (u(l+1))   (u(l))
    |f dt    1    | (v(l+1)) = (v(l))

    |u(l+1)|                   |1      f dt| |u(l)|
    |      | = 1/(f^2 dt^2+1)  |           | |    |
    |v(l+1)|                   |-f dt  1   | |v(l)|
    """
    tmp = 1./(fcor**2*deltaT**2+1.)
    ulp1 = (uvel + fcor*deltaT*vvel)/(fcor**2*deltaT**2+1.)
    vlp1 = (vvel - fcor*deltaT*uvel)/(fcor**2*deltaT**2+1.)
    Tu_cor = (ulp1-uvel)/deltaT
    Tv_cor = (vlp1-vvel)/deltaT

  Tu_hpr = - dpdx
  Tv_hpr = - dpdy
  Tu_wnd[0] = taux/dz[0]
  Tv_wnd[0] = tauy/dz[0]
  Tu_bot[-1] = -bottomDragQuadratic/dz[-1]*np.sqrt(uvel[-1]**2+vvel[-1]**2)*uvel[-1]
  Tv_bot[-1] = -bottomDragQuadratic/dz[-1]*np.sqrt(uvel[-1]**2+vvel[-1]**2)*vvel[-1]
  Tu = Tu_cor + Tu_hpr + Tu_wnd + Tu_bot #+ Tu_vdf
  Tv = Tv_cor + Tv_hpr + Tv_wnd + Tv_bot #+ Tv_vdf

  # diffusion
  """
  T_diff = d/dz k d/dz p
         = 1/dz[k] {    A[k]  ( p[k-1]-p[k]   ) / dzt[k] 
                     -  A[k+1]( p[k]  -p[k+1] ) / dzt[k+1] }
         =    A[k]/(dz[k]*dzt[k]) * p[k-1] 
           - {A[k]/(dz[k]*dzt[k])+A[k+1]/(dz[k]*dzt[k+1])} * p[k]
           +  A[k+1]/(dz[k]*dzt[k+1]) * p[k+1]
         = a[k] * p[k-1] - (a+c)[k] * p[k] + c[k] * p[k+1]
  Neuman:
  k=1:  T_diff = 1/dz[k] { -  A[k+1]( p[k]  -p[k+1] ) / dzt[k+1] } 
  k=nz: T_diff = 1/dz[k] {    A[k]  ( p[k-1]-p[k]   ) / dzt[k]   }
  Dirichlet:
  k=1:  T_diff = -(2*a+b)[k]p[k] + c[k]*p[k+1]
  k=nz: T_diff = a[k]*p[k-1] - (a+2*c)[k]p[k]
  """
  am = Av[:-1]/(dz*dzt[:-1])
  cm = Av[1:] /(dz*dzt[1:])
  bm = am+cm
  bm[0] = cm[0]
  #bm[0] = 2.*am[0]+cm[0] # no-slip
  freeSlipBottom = False
  if freeSlipBottom:   # free-slip
    bm[-1] = am[-1]
  else:                # no-slip
    bm[-1] = am[-1]+2.*cm[-1]

  implicite = True
  if implicite:
    # ----- solve tri-diagonal diffusion matrix
    atr = -deltaT*am
    btr = 1+deltaT*bm
    ctr = -deltaT*cm
    # ----- for u
    dtru = 1.*uvel
    dtru += deltaT*(Tu)
    uvel = TDMAsolver(1.*atr, 1.*btr, 1.*ctr, 1.*dtru)
    # ----- for v
    dtrv = 1.*vvel
    dtrv += deltaT*(Tv)
    vvel = TDMAsolver(1.*atr, 1.*btr, 1.*ctr, 1.*dtrv)

  Tu_vdf = np.zeros((nz))
  Tu_vdf[1:-1] = am[1:-1]*uvel[:-2] - bm[1:-1]*uvel[1:-1] + cm[1:-1]*uvel[2:]
  Tu_vdf[0]  = -bm[0] *uvel[0]  + cm[0] *uvel[1]
  #Tu_vdf[-1] =  am[-1]*uvel[-2] - bm[-1]*uvel[-1]
  # always add no-slip contribution to bottom tendency
  Tu_vdf[-1] =  am[-1]*uvel[-2] - am[-1]*uvel[-1]
  if not freeSlipBottom:
    Tu_bot[-1] +=  -2.*cm[-1]*uvel[-1]

  Tv_vdf = np.zeros((nz))
  Tv_vdf[1:-1] = am[1:-1]*vvel[:-2] - bm[1:-1]*vvel[1:-1] + cm[1:-1]*vvel[2:]
  Tv_vdf[0]  = -bm[0] *vvel[0]  + cm[0] *vvel[1]
  #Tv_vdf[-1] =  am[-1]*vvel[-2] - bm[-1]*vvel[-1]
  # always add no-slip contribution to bottom tendency
  Tv_vdf[-1] =  am[-1]*vvel[-2] - am[-1]*vvel[-1] 
  if not freeSlipBottom:
    Tv_bot[-1] +=  -2.*cm[-1]*vvel[-1]

  Tu += Tu_vdf
  Tv += Tv_vdf

  # explicite diffusion:  d/dz(k d/dz phi)
  #flux = Av*uz
  #Tu_vdf = (flux[:-1]-flux[1:])/dz
  #flux = Av*vz
  #Tv_vdf = (flux[:-1]-flux[1:])/dz
  #flux = Av*vz
  #Tv_vdf = (flux[:-1]-flux[1:])/dz

  # ==========
  # buoycancy
  # ==========
  # buoyancy tendencies
  bi = np.zeros((nz+1))
  bi[1:-1] = 0.5*(b[:-1]+b[1:]) # FIXME: not true for unequal grid spacing
  flux = wvel*bi
  Tb_vad = -(flux[:-1]-flux[1:])/dz
  Tb_res = lam_b*(b0-b)
  #flux = kv*N2
  #Tb_vdf = (flux[:-1]-flux[1:])/dz
  Tb_sfl[0] = grav*tAlpha/(rho0*cp)/dz[0]*Qsurf  # for [Qsurf]=W/m^2
  Tb = Tb_vad + Tb_sfl + Tb_res #+ Tb_vdf

  # diffusion
  am = kv[:-1]/(dz*dzt[:-1])
  cm = kv[1:] /(dz*dzt[1:])
  bm = am+cm
  bm[0] = cm[0]
  bm[-1] = am[-1]

  if implicite:
    # ----- solve tri-diagonal diffusion matrix
    atr = -deltaT*am
    btr = 1+deltaT*bm
    ctr = -deltaT*cm
    # ----- for b
    dtrb = 1.*b
    dtrb += deltaT*(Tb)
    b = TDMAsolver(1.*atr, 1.*btr, 1.*ctr, 1.*dtrb)

  Tb_vdf = np.zeros((nz))
  Tb_vdf[1:-1] = am[1:-1]*b[:-2] - bm[1:-1]*b[1:-1] + cm[1:-1]*b[2:]
  Tb_vdf[0]  = -bm[0] *b[0]  + cm[0] *b[1]
  Tb_vdf[-1] =  am[-1]*b[-2] - bm[-1]*b[-1]

  Tb += Tb_vdf

  # ==========
  # Adams-Bashforth timestepping
  # ==========
  if l==0:
    Tu_old = 1.*Tu
    Tv_old = 1.*Tv
    Tb_old = 1.*Tb
    Tt_old = 1.*Tt

  epsab = 0.01
  if not implicite:
    uvel = uvel + deltaT*((1.5+epsab)*Tu - (0.5+epsab)*Tu_old)
    vvel = vvel + deltaT*((1.5+epsab)*Tv - (0.5+epsab)*Tv_old)
    b    = b    + deltaT*((1.5+epsab)*Tb - (0.5+epsab)*Tb_old)
  if not implicite_tke:
    tke  = tke  + deltaT*((1.5+epsab)*Tt - (0.5+epsab)*Tt_old)

  tke_before_bck = 1.*tke
  tke[tke<tke_min] = tke_min
  Tt_bck = (tke-tke_before_bck)/deltaT

  Tu_old = 1.*Tu
  Tv_old = 1.*Tv
  Tb_old = 1.*Tb
  Tt_old = 1.*Tt

  Tu_tot = (uvel-uvel_old)/deltaT
  Tv_tot = (vvel-vvel_old)/deltaT
  Tb_tot = (b-b_old)/deltaT
  Tt_tot = (tke-tke_old)/deltaT

  # saving variables
  if (l%lsave)==0:
    print('saving at l = %d / %d: %.1f'%(l, nt, 100.*l/nt))
    b_s[ls,:] = b
    uvel_s[ls,:] = uvel
    vvel_s[ls,:] = vvel
    tke_s[ls,:] = tke
    kv_s[ls,:] = kv
    Av_s[ls,:] = Av
    Lmix_s[ls,:] = Lmix
    N2_s[ls,:] = N2

    Tke_cor[ls] = ((ulp1*Tu_cor + vlp1*Tv_cor)*dz).sum()
    Tke_hpr[ls] = ((uvel*Tu_hpr + vvel*Tv_hpr)*dz).sum()
    Tke_wnd[ls] = ((uvel*Tu_wnd + vvel*Tv_wnd)*dz).sum()
    Tke_bot[ls] = ((uvel*Tu_bot + vvel*Tv_bot)*dz).sum()
    Tke_vdf[ls] = ((uvel*Tu_vdf + vvel*Tv_vdf)*dz).sum()
    Tke_tot[ls] = ((uvel*Tu_tot + vvel*Tv_tot)*dz).sum()

    Ttke_tot[ls] = (Tt_tot*dzt).sum()
    Ttke_bpr[ls] = (Tt_bpr*dzt).sum()
    Ttke_spr[ls] = (Tt_spr*dzt).sum()
    Ttke_dis[ls] = (Tt_dis*dzt).sum()
    Ttke_vdf[ls] = (Tt_vdf*dzt).sum()
    Ttke_bck[ls] = (Tt_bck*dzt).sum()

    Qsurf_ts[ls] = Qsurf
    taux_ts[ls] = taux
    tauy_ts[ls] = tauy

    """
    u d/dz Av d/dz u = d/dz u Av d/dz u - av (d/dz u)^2

    """
    tmp = -Av*(uz)**2 - Av*(vz)**2
    Tke_vds[ls] = ((0.5*(tmp[1:]+tmp[:-1]))*dz).sum()
    Tke_vfl[ls] = ((   uvel*Tu_vdf + vvel*Tv_vdf 
                     - 0.5*(tmp[1:]+tmp[:-1])   )*dz).sum()


    ls += 1

Tt_err = Tt_tot - (Tt_bpr+Tt_spr+Tt_dis+Tt_vdf+Tt_bck)
Tb_err = Tb_tot - (Tb_vad+Tb_res+Tb_vdf+Tb_sfl)
Tu_err = Tu_tot - (Tu_cor+Tu_hpr+Tu_wnd+Tu_bot+Tu_vdf)
Tv_err = Tv_tot - (Tv_cor+Tv_hpr+Tv_wnd+Tv_bot+Tv_vdf)

mke_ts = (0.5*(uvel_s**2+vvel_s**2)*dz[np.newaxis,:]).sum(axis=1)
tke_ts = (tke_s*dzt[np.newaxis,:]).sum(axis=1)
uvel_ts = (uvel_s*dz[np.newaxis,:]).sum(axis=1)
vvel_ts = (vvel_s*dz[np.newaxis,:]).sum(axis=1)
#ape = (0.5*b_s**2/(N2i+1e-33))*dz[np.newaxis,:]).sum(axis=1)
b2 = (b_s**2*dz[np.newaxis,:]).sum(axis=1)
time = deltaT*lsave * np.arange(nsave)
