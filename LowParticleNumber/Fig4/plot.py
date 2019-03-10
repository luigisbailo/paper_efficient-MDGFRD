import matplotlib
import numpy as np
import matplotlib.pyplot as pl


pl.style.use(['style1'])

pl.figure(figsize=(3.37, 2.5))

alpha=8.4

dt=0.1

data = np.loadtxt ('out.txt')

pl.plot ( np.arange(0,10,0.001), np.sqrt(np.arange(0,10,0.001)/1000)*np.sqrt(dt)*alpha,linewidth=0.5,linestyle='--')
pl.plot ( data[:,0]*1000, data [:,1], 'o')

pl.xlabel (r'$D\thinspace [\frac{\mu m^2}{s}]$', fontsize=10) 
pl.ylabel (r'$\rho\thinspace [nm$]',fontsize=10)
pl.xlim (0,10.5)	
#pl.yticks([0.1,0.2,0.3])
pl.tick_params('both',labelsize=10)

pl.tight_layout()

# pl.show()
pl.savefig('GF-BF.eps')



# f, (ax1,ax2) = pl.subplots (2, sharex='col') 

# dt=0.1

# data = np.loadtxt ('dt' + str(dt) +'.txt')

# ax1.plot ( np.arange(0,10,0.001), np.sqrt(np.arange(0,10,0.001)/1000)*np.sqrt(dt)*alpha)
# ax1.plot ( data[:,0]*1000, data [:,1], 'o')

# ax1.set_ylabel (r'$\rho\thinspace [nm$]')#, fontsize=20 )
# ax1.set_xlim (0,10)	
# ax1.set_yticks([0.1,0.2,0.3])


# dt = 0.01

# data = np.loadtxt ('dt' + str(dt) +'.txt')

# ax2.plot ( np.arange(0,10,0.001), np.sqrt(np.arange(0,10,0.001)/1000)*np.sqrt(dt)*alpha)
# ax2.plot ( data[:,0]*1000, data [:,1], 'o')

# ax2.set_xlabel (r'$D\thinspace [\frac{\mu m}{s}]$')#, fontsize=20 )
# ax2.set_ylabel (r'$\rho\thinspace [nm]$')#, fontsize=20 )

# ax2.set_xlim (0,10)
# ax2.set_yticks([0.025,0.05,0.075,0.1])

# pl.tight_layout()

# #pl.show()
# pl.savefig('GF-BF.pdf')
