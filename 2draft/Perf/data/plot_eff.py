import matplotlib
import numpy as np
import matplotlib.pyplot as pl

pl.style.use(['style1'])

# pl.figure(figsize=(5, 3))
f, ((ax1),(ax2)) = pl.subplots(1,2,sharey=True)

f.set_figheight(3.5)
f.set_figwidth(6.69)


data_burst = np.loadtxt("1.burst.txt")
data_GF = np.loadtxt("1.GF.txt")
filename = "1.data.txt"
param = np.loadtxt (filename)
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.plot(100000000/data_burst[:,0]/data_burst[:,0]/data_burst[:,0]/6.022,data_burst[:,5]/data_GF[:,5],color='#9637B8',marker='^',markersize=4,linestyle='-',linewidth=1, label='MD-GFRD\ 1')
ax1.plot(100000000/data_burst[:,0]/data_burst[:,0]/data_burst[:,0]/6.022,data_burst[:,4]/data_GF[:,4],color='#B86637',marker='v',markersize=4,linestyle='-',linewidth=1, label='MD-GFRD\ 2')	
ax1.plot(100000000/data_burst[:,0]/data_burst[:,0]/data_burst[:,0]/6.022,data_burst[:,3]/data_GF[:,3],color='#B8375B',marker='s',markersize=4,linestyle='-.',linewidth=1,label='Hybrid scheme')
ax1.plot(100000000/data_burst[:,0]/data_burst[:,0]/data_burst[:,0]/6.022,data_burst[:,2]/data_GF[:,2],color='#B8A737',marker='H',markersize=4,linestyle='-',linewidth=1, label='New Scheme\ 1')	
ax1.plot(100000000/data_burst[:,0]/data_burst[:,0]/data_burst[:,0]/6.022,data_burst[:,1]/data_GF[:,1],color='#377EB8',marker='o',markersize=4,linestyle='-',linewidth=1,label='New scheme\ 2')

text = r"$D = $" + str(int(1000*param[0])) + r"$\thinspace\frac{\mu m^2}{s}$"
box=dict(facecolor='white', edgecolor='grey', boxstyle='round')
ax1.text(0.168,0.8,text,fontsize=8,color='grey',bbox=box,horizontalalignment='center',verticalalignment='center',transform = ax1.transAxes)

ax1.text(0.05,0.91,'(a)',fontsize=10,transform = ax1.transAxes)
ax2.text(0.05,0.91,'(b)',fontsize=10,transform = ax2.transAxes)

# ax1.legend(loc=4,fontsize=10)


data_burst = np.loadtxt("2.burst.txt")
data_GF = np.loadtxt("2.GF.txt")
filename = "2.data.txt"
param = np.loadtxt (filename)
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.plot(100000000/data_burst[:,0]/data_burst[:,0]/data_burst[:,0]/6.022,data_burst[:,5]/data_GF[:,5],color='#9637B8',marker='^',markersize=4,linestyle='-',linewidth=1, label='MD-GFRD\ 1')
ax2.plot(100000000/data_burst[:,0]/data_burst[:,0]/data_burst[:,0]/6.022,data_burst[:,4]/data_GF[:,4],color='#B86637',marker='v',markersize=4,linestyle='-',linewidth=1, label='MD-GFRD\ 2')	
ax2.plot(100000000/data_burst[:,0]/data_burst[:,0]/data_burst[:,0]/6.022,data_burst[:,3]/data_GF[:,3],color='#B8375B',marker='s',markersize=4,linestyle='-.',linewidth=1,label='Hybrid scheme')
ax2.plot(100000000/data_burst[:,0]/data_burst[:,0]/data_burst[:,0]/6.022,data_burst[:,2]/data_GF[:,2],color='#B8A737',marker='H',markersize=4,linestyle='-',linewidth=1, label='New Scheme\ 1')	
ax2.plot(100000000/data_burst[:,0]/data_burst[:,0]/data_burst[:,0]/6.022,data_burst[:,1]/data_GF[:,1],color='#377EB8',marker='o',markersize=4,linestyle='-',linewidth=1,label='New scheme\ 2')

text = r"$D_1 = $" + str(int(1000*param[0])) + r"$\thinspace\frac{\mu m^2}{s}$"+"\n" +r"$D_2 = $" + str(int(1000*param[1])) + r"$\thinspace\frac{\mu m^2}{s}$"
box=dict(facecolor='white', edgecolor='grey', boxstyle='round')
ax2.text(0.175,0.75,text,fontsize=8,color='grey',bbox=box,horizontalalignment='center',verticalalignment='center',transform = ax2.transAxes)

ax1.legend(loc=4,fontsize=8)

ax1.set_xlim(0.0005,150)
ax2.set_xlim(0.0005,150)


ax1.set_ylim(0.002,0.5)
ax2.set_ylim(0.002,0.5)

ax1.set_yticks([0.01,0.1])
ax2.set_xticks([0.001,0.01,0.1,1,10,100])
ax1.tick_params('y',labelsize=10)
ax2.tick_params('both',labelsize=10)

ax1.set_ylabel(r'$Bursting\ probability$', fontsize=10)
ax1.set_xlabel(r'$Molar\ concentration\ [\mu M]$',fontsize=10)
ax2.set_xlabel(r'$Molar\ concentration\ [\mu M]$',fontsize=10)

pl.tight_layout()

f.subplots_adjust(wspace=0.05)
f.subplots_adjust(hspace=0.05)

#pl.show()
pl.savefig('eff.eps',dpi=600)
