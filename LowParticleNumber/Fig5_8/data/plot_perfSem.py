import matplotlib
import numpy as np
import matplotlib.pyplot as pl
import os

pl.style.use(['style1'])

f, ((ax1)) = pl.subplots(1,1,sharey=True)

f.set_figheight(3.5)
f.set_figwidth(3.7)
# f.set_figwidth(3.37)

filename = "1.perf.txt"
data = np.loadtxt(filename)
filename = "1.data.txt"
param = np.loadtxt (filename)
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.errorbar(100000000/data[:,0]/data[:,0]/data[:,0]/6.022,data[:,5],yerr=data[:,11],color='#B86637',marker='s',markersize=4,linestyle='-',linewidth=1, label='MD-GFRD')
ax1.errorbar(100000000/data[:,0]/data[:,0]/data[:,0]/6.022,data[:,2],yerr=data[:,8],color='#377EB8',marker='H',markersize=4,linestyle='-',linewidth=1, label='New Scheme')
ax1.plot(100000000/data[:,0]/data[:,0]/data[:,0]/6.022,data[:,6],color='#E41A1C',linestyle='--',markersize=4, label='BD')

text = r"$D = $" + str(int(1000*param[0])) + r"$\thinspace\frac{\mu m^2}{s}$"
box=dict(facecolor='white', edgecolor='grey', boxstyle='round')

ax1.text(0.168,0.81,text,fontsize=8,color='grey',bbox=box,horizontalalignment='center',verticalalignment='center',transform = ax1.transAxes)

text = r"$D_1 = $" + str(int(1000*param[0])) + r"$\thinspace\frac{\mu m^2}{s}$"+"\n" +r"$D_2 = $" + str(int(1000*param[1])) + r"$\thinspace\frac{\mu m^2}{s}$"

ax1.legend(loc=4,fontsize=8)

ax1.set_xlim(0.0005,1500)
ax1.set_ylim (0.0001,150)

ax1.set_yticks([0.001,0.01,0.1,1,10,100])
ax1.set_xticks([0.001,0.01,0.1,1,10,100,1000])
ax1.tick_params('both',labelsize=10)

ax1.set_ylabel(r'$CPU\ time\ [s]$',fontsize=10)
# ax2.set_ylabel(r'$CPU\ time\ [s]$',fontsize=18)
ax1.set_xlabel(r'$Molar\ concentration\ [\mu M]$',fontsize=10)


pl.tight_layout()

f.subplots_adjust(wspace=0.05)
f.subplots_adjust(hspace=0.05)

#pl.show()
pl.savefig('perfSem.eps', dpi=600)
