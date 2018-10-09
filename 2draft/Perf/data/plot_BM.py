import matplotlib
import numpy as np
import matplotlib.pyplot as pl
import os

pl.style.use(['style1'])

f, ((ax1),(ax2)) = pl.subplots(1,2,sharey=True)

f.set_figheight(3.5)
f.set_figwidth(6.69)

filename = "1.BM.txt"
data = np.loadtxt(filename)
filename = "1.data.txt"
param = np.loadtxt (filename)
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.errorbar(100000000/data[:,0]/data[:,0]/data[:,0]/6.022,data[:,5],yerr=data[:,10],color='#9637B8',marker='^',markersize=4,linestyle='-',linewidth=1, label='MD-GFRD\ 1')
ax1.errorbar(100000000/data[:,0]/data[:,0]/data[:,0]/6.022,data[:,4],yerr=data[:,9],color='#B86637',marker='v',markersize=4,linestyle='-',linewidth=1, label='MD-GFRD\ 2')	
ax1.errorbar(100000000/data[:,0]/data[:,0]/data[:,0]/6.022,data[:,3],yerr=data[:,8],color='#B8375B',marker='s',markersize=4,linestyle='-.',linewidth=1,label='Hybrid scheme')
ax1.errorbar(100000000/data[:,0]/data[:,0]/data[:,0]/6.022,data[:,2],yerr=data[:,7],color='#B8A737',marker='H',markersize=4,linestyle='-',linewidth=1, label='New Scheme\ 1')	
ax1.errorbar(100000000/data[:,0]/data[:,0]/data[:,0]/6.022,data[:,1],yerr=data[:,6],color='#377EB8',marker='o',markersize=4,linestyle='-',linewidth=1,label='New scheme\ 2')

text = r"$D = $" + str(int(1000*param[0])) + r"$\thinspace\frac{\mu m^2}{s}$"
box=dict(facecolor='white', edgecolor='grey', boxstyle='round')
ax1.text(0.168,0.8,text,fontsize=8,color='grey',bbox=box,horizontalalignment='center',verticalalignment='center',transform = ax1.transAxes)

ax1.text(0.05,0.91,'(a)',fontsize=10,transform = ax1.transAxes)
ax2.text(0.05,0.91,'(b)',fontsize=10,transform = ax2.transAxes)

# ax1.legend(loc=4,fontsize=8)


filename = "2.BM.txt"
data = np.loadtxt(filename)
filename = "2.data.txt"
param = np.loadtxt (filename)
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.errorbar(100000000/data[:,0]/data[:,0]/data[:,0]/6.022,data[:,5],yerr=data[:,10],color='#9637B8',marker='^',markersize=4,linestyle='-',linewidth=1, label='MD-GFRD\ 1')
ax2.errorbar(100000000/data[:,0]/data[:,0]/data[:,0]/6.022,data[:,4],yerr=data[:,9],color='#B86637',marker='v',markersize=4,linestyle='-',linewidth=1, label='MD-GFRD\ 2')	
ax2.errorbar(100000000/data[:,0]/data[:,0]/data[:,0]/6.022,data[:,3],yerr=data[:,8],color='#B8375B',marker='s',markersize=4,linestyle='-.',linewidth=1,label='Hybrid scheme')
ax2.errorbar(100000000/data[:,0]/data[:,0]/data[:,0]/6.022,data[:,2],yerr=data[:,7],color='#B8A737',marker='H',markersize=4,linestyle='-',linewidth=1, label='New Scheme\ 1')	
ax2.errorbar(100000000/data[:,0]/data[:,0]/data[:,0]/6.022,data[:,1],yerr=data[:,6],color='#377EB8',marker='o',markersize=4,linestyle='-',linewidth=1,label='New scheme\ 2')

#text = r"$\thinspace[\mu m^2/s]$" + "D_2 = " + str(1000*param[1]) + r"$\thinspace[\mu m^2/s]$" + r"$t = $" + str(param[2]) + r"$\thinspace[ns]$"
#	ax.text (0.004,70,r"$D_1 = $" + str(1000*param[0]) + r"$\thinspace[\frac{\mu m^2}{s}]$",fontsize=7)
#	ax.text (2,70, r"$D_2 = $" + str(1000*param[1]) + r"$\thinspace[\frac{\mu m^2}{s}]$",fontsize=7)

text = r"$D_1 = $" + str(int(1000*param[0])) + r"$\thinspace\frac{\mu m^2}{s}$"+"\n" +r"$D_2 = $" + str(int(1000*param[1])) + r"$\thinspace\frac{\mu m^2}{s}$"
box=dict(facecolor='white', edgecolor='grey', boxstyle='round')
ax2.text(0.175,0.75,text,fontsize=8,color='grey',bbox=box,horizontalalignment='center',verticalalignment='center',transform = ax2.transAxes)

ax1.legend(loc=4,fontsize=8)

ax1.set_xlim(0.0005,1500)
ax2.set_xlim(0.0005,1500)

ax1.set_ylim(2,200000000)
ax2.set_ylim(2,200000000)

# ax1.set_yticks([0.1,1,10,100,1000])
ax1.set_yticks([10,100,1000,10000,100000,1000000,10000000,100000000])
ax2.set_yticks([10,100,1000,10000,100000,1000000,10000000,100000000])
# ax2.set_xticks([0.001,0.01,0.1,1,10,100])
ax1.tick_params('y',labelsize=10)
ax2.tick_params('both',labelsize=10)

ax1.set_ylabel(r'$Total\ number\ of\ brute-force \ steps$',fontsize=10)
ax1.set_xlabel(r'$Molar\ concentration\ [\mu M]$',fontsize=10)
ax2.set_xlabel(r'$Molar\ concentration\ [\mu M]$',fontsize=10)

pl.tight_layout()

f.subplots_adjust(wspace=0.05)
f.subplots_adjust(hspace=0.05)

#pl.show()
pl.savefig('BM.eps',dpi=600)
