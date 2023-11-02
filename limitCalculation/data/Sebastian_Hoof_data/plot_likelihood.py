import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import griddata

npoints = 300
lgm, lgg, ll = np.genfromtxt("cast2017_loglike.dat", unpack=True)
lgm = np.log10(lgm)
lgg = np.log10(lgg)
ll -= np.max(ll)

x0 = np.linspace(lgm.min(), lgm.max(), npoints)
y0 = np.linspace(lgg.min(), lgg.max(), npoints)
z0 = griddata((lgm, lgg), ll, (x0[None,:], y0[:,None]), method='linear')

fig, ax = plt.subplots()

plt.contourf(x0, y0, np.exp(z0), cmap='YlOrRd', levels=25, zorder=-9)

cs = plt.colorbar(ticks=[0.2*i for i in range(6)], pad=0.025)
cs.ax.set_ylabel('Profile likelihood', labelpad=12, rotation=270)
cs.ax.tick_params(labelsize=6, length=2, pad=2)
cs.ax.set_ylim(0,1)

ax.contour(x0, y0, -2.0*z0, colors='k', linestyles='-', linewidths=1, levels=[1, 4, 9])
ax.set_rasterization_zorder(-1)

ax.set_xlabel(r"ALP mass $\log_{10}(m_a/\mathrm{eV})$")
ax.set_ylabel(r"ALP-photon coupling $\log_{10}(g_{a\gamma}/\mathrm{GeV}^{-1})$")

ax.set_xlim([lgm.min(), lgm.max()])
ax.set_ylim([lgg.min(), lgg.max()])

fig.tight_layout(pad=0.3)
plt.savefig("cast2017_loglike.pdf")#, backend='pgf')
plt.show()