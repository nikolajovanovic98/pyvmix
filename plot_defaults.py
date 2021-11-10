# ---
hca, hcb = pyic.arrange_axes(6,2, plot_cb=False, asp=2., fig_size_fac=2.,
                            sharex=False, sharey=True, xlabel="", ylabel="depth [m]")
ii=-1

nclev=nsave; cols = plt.cm.get_cmap('jet')(np.arange(nclev)/(float(nclev-1)))

ii+=1; ax=hca[ii]; cax=hcb[ii]
for l in range(nsave):
  ax.plot(b_s[l,:], zt, color=cols[l,:])
ax.set_title('b')

ii+=1; ax=hca[ii]; cax=hcb[ii]
for l in range(nsave):
  ax.plot(uvel_s[l,:]/1e-2, zt, color=cols[l,:])
ax.set_title('uvel [cm/s]')

ii+=1; ax=hca[ii]; cax=hcb[ii]
for l in range(nsave):
  ax.plot(vvel_s[l,:]/1e-2, zt, color=cols[l,:])
ax.set_title('vvel [cm/s]')

ii+=1; ax=hca[ii]; cax=hcb[ii]
for l in range(nsave):
  ax.plot(tke_s[l,:]/1e-3, zu, color=cols[l,:])
ax.set_title('tke/1e-3 [m$^2$/s$^2$]')

ii+=1; ax=hca[ii]; cax=hcb[ii]
for l in range(nsave):
  ax.plot(kv_s[l,:], zu, color=cols[l,:])
ax.set_title('kv')

ii+=1; ax=hca[ii]; cax=hcb[ii]
for l in range(nsave):
  ax.plot(Av_s[l,:], zu, color=cols[l,:], marker='.')
ax.set_title('Av')

ii+=1; ax=hca[ii]; cax=hcb[ii]
for l in range(nsave):
  ax.plot(Lmix_s[l,:], zu, color=cols[l,:])
ax.set_title('Lmix')

ii+=1; ax=hca[ii]; cax=hcb[ii]
for l in range(nsave):
  ax.plot(N2_s[l,:], zu, color=cols[l,:])
ax.set_title('N2')

nnf+=1
if savefig:
  print('save figure: %s_%02d.pdf' % (__file__.split('/')[-1][:-3], nnf))
  plt.savefig("%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], nnf))

# --- tendencies
hca, hcb = pyic.arrange_axes(2,2, plot_cb=False, asp=1.5, fig_size_fac=2.,
                            sharex=False, sharey=False, xlabel="", ylabel="")
ii=-1

marker='.'

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(Tt_bpr, zu, label='bpr', marker=marker)
ax.plot(Tt_spr, zu, label='spr', marker=marker)
ax.plot(Tt_dis, zu, label='dis', marker=marker)
ax.plot(Tt_vdf, zu, label='vdf', marker=marker)
ax.plot(Tt_bck, zu, label='bck', marker=marker)
ax.plot(Tt_tot, zu, label='tot', marker=marker)
ax.plot(Tt_err, zu, label='err', marker=marker)
ax.legend(loc='best')
ax.set_title('tke tendencies')

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(Tb_vad, zt, label='vad', marker=marker)
ax.plot(Tb_res, zt, label='res', marker=marker)
ax.plot(Tb_vdf, zt, label='vdf', marker=marker)
ax.plot(Tb_sfl, zt, label='sfl', marker=marker)
ax.plot(Tb_tot, zt, label='tot', marker=marker)
ax.plot(Tb_err, zt, label='err', marker=marker)
ax.legend(loc='best')
ax.set_title('buoyancy tendencies')

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(Tu_cor, zt, label='cor', marker=marker)
ax.plot(Tu_hpr, zt, label='hpr', marker=marker)
ax.plot(Tu_wnd, zt, label='wnd', marker=marker)
ax.plot(Tu_bot, zt, label='bot', marker=marker)
ax.plot(Tu_vdf, zt, label='vdf', marker=marker)
ax.plot(Tu_tot, zt, label='tot', marker=marker)
ax.plot(Tu_err, zt, label='err', marker=marker)
ax.legend(loc='best')
ax.set_title('uvel tendencies')

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(Tv_cor, zt, label='cor', marker=marker)
ax.plot(Tv_hpr, zt, label='hpr', marker=marker)
ax.plot(Tv_wnd, zt, label='wnd', marker=marker)
ax.plot(Tv_bot, zt, label='bot', marker=marker)
ax.plot(Tv_vdf, zt, label='vdf', marker=marker)
ax.plot(Tv_tot, zt, label='tot', marker=marker)
ax.plot(Tv_err, zt, label='err', marker=marker)
ax.legend(loc='best')
ax.set_title('vvel tendencies')

nnf+=1
if savefig:
  print('save figure: %s_%02d.pdf' % (__file__.split('/')[-1][:-3], nnf))
  plt.savefig("%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], nnf))

# --- time series
hca, hcb = pyic.arrange_axes(4,2, plot_cb=False, asp=1., fig_size_fac=2.,
                            sharex=False, sharey=False, xlabel="", ylabel="")
ii=-1

tpl = time/86400
tstr = 'time [days]'

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(tpl, mke_ts)
ax.set_title('mke_ts')

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(tpl, tke_ts)
ax.set_title('tke_ts')

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(tpl, uvel_ts)
ax.plot(tpl, vvel_ts)
ax.set_title('uvel_ts and vvel_ts')

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(tpl, Tke_cor, label='cor')
ax.plot(tpl, Tke_hpr, label='hpr')
ax.plot(tpl, Tke_wnd, label='wnd')
ax.plot(tpl, Tke_bot, label='bot')
ax.plot(tpl, Tke_vdf, label='vdf')
ax.plot(tpl, Tke_vds, label='vds')
ax.plot(tpl, Tke_vfl, label='vfl')
ax.plot(tpl, Tke_tot, label='tot')
ax.legend(loc='best')
ax.set_title('ke budget')

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(tpl, Ttke_bpr, label='bpr')
ax.plot(tpl, Ttke_spr, label='spr')
ax.plot(tpl, Ttke_dis, label='dis')
ax.plot(tpl, Ttke_vdf, label='vdf')
ax.plot(tpl, Ttke_bck, label='bck')
ax.plot(tpl, Ttke_tot, label='tot')
ax.legend(loc='best')
ax.set_title('tke budget')

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(tpl, Qsurf_ts, label='Qsurf')
ax.set_title('Qsurf')

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(tpl, taux_ts, label='taux')
ax.plot(tpl, tauy_ts, label='tauy')
ax.set_title('wind stress')
ax.legend()

for ax in hca:
  ax.set_xlabel(tstr)

nnf+=1
if savefig:
  print('save figure: %s_%02d.pdf' % (__file__.split('/')[-1][:-3], nnf))
  plt.savefig("%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], nnf))

# --- time series
hca, hcb = pyic.arrange_axes(4,2, plot_cb=True, asp=0.5, fig_size_fac=1.5,
                            sharex=False, sharey=False, 
                            xlabel=tstr, ylabel="depth [m]")
ii=-1

ii+=1; ax=hca[ii]; cax=hcb[ii]
pyic.shade(tpl, zt, b_s.transpose(), ax=ax, cax=cax, clim='auto') 
ax.set_title('buoyancy')

ii+=1; ax=hca[ii]; cax=hcb[ii]
pyic.shade(tpl, zu, N2_s.transpose(), ax=ax, cax=cax, clim='auto') 
ax.set_title('stratifiaction N2')

ii+=1; ax=hca[ii]; cax=hcb[ii]
pyic.shade(tpl, zu, tke_s.transpose(), ax=ax, cax=cax, clim=[-6,0], logplot=True) 
ax.set_title('TKE')

ii+=1; ax=hca[ii]; cax=hcb[ii]
pyic.shade(tpl, zt, uvel_s.transpose(), ax=ax, cax=cax, clim='sym') 
ax.set_title('zon. velocity')

ii+=1; ax=hca[ii]; cax=hcb[ii]
pyic.shade(tpl, zt, vvel_s.transpose(), ax=ax, cax=cax, clim='sym') 
ax.set_title('mer. velocity')

ii+=1; ax=hca[ii]; cax=hcb[ii]
pyic.shade(tpl, zu, np.sqrt(uvel_s**2+vvel_s**2).transpose(), ax=ax, cax=cax, clim=0.2) 
ax.set_title('\sqrt(u^2+v^2)')

for ax in hca:
  ax.set_ylim(ylim_hov)

nnf+=1
if savefig:
  print('save figure: %s_%02d.pdf' % (__file__.split('/')[-1][:-3], nnf))
  plt.savefig("%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], nnf))
