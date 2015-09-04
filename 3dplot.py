import os
os.chdir(os.path.expanduser("~/Desktop/reflectivity/dynXRD/"))
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import pylab as pl
import numpy as np



Q = pl.loadtxt('SrTiO3_Q.dat', unpack = True)
Energy= pl.loadtxt('SrTiO3_Energy.dat', unpack = True)
refl_matrix_fit = pl.loadtxt('SrTiO3_fit_matrix.dat', unpack = True)
refl_matrix_data = pl.loadtxt('SrTiO3_data_matrix.dat', unpack = True)


fig=pl.figure()
ax=fig.gca(projection='3d')
Energy1, Q1=np.meshgrid(Energy, Q)




ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_zscale('log')
ax.autoscale(enable=True, axis='z', tight=True)


#surf_fit=ax.plot_surface(Q1, Energy1, refl_matrix_fit, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#surf_data=ax.plot_surface(Q1, Energy1, refl_matrix_data, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# fig.colorbar(surf_fit, shrink=0.5, aspect=5)
# fig.colorbar(surf_data, shrink=0.5, aspect=5)

# ax.plot_wireframe(Q1, Energy1, refl_matrix_data, rstride=10, cstride=10, color='blue')
# ax.plot_wireframe(Q1, Energy1, refl_matrix_fit, rstride=10, cstride=10, color='red')

ax.plot_surface(Q1, Energy1, refl_matrix_data, rstride=9, cstride=9, alpha=0.3, color='blue')
ax.plot_surface(Q1, Energy1, refl_matrix_fit, rstride=9, cstride=9, alpha=0.3, color='red')

pl.show()