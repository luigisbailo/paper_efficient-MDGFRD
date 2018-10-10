import matplotlib
import numpy as np
import matplotlib.pyplot as pl
import math

pl.style.use(['style1'])

f, ((ax1),(ax2)) = pl.subplots(2,1,sharex=True)

f.set_figheight(5)
f.set_figwidth(3.37)

data = np.loadtxt('diff1.out')


D=0.01

ax1.plot ( np.arange(0,600000,1000)/1000, np.arange(0,600000,1000)*6*D/pow(10,6),linewidth=0.5,linestyle='--',color='#E41A1C')
ax1.errorbar(data[:,0]/1000,data[:,4]/pow(10,6),yerr=data[:,8]/pow(10,6),marker='x',markersize=8,color='#B8A737',linestyle='None',label='BD')
# ax1.errorbar(data[:,0]/1000,data[:,3]/pow(10,6),yerr=data[:,7]/pow(10,6),marker='^',markersize=4,color='#9337B8',linestyle='None',label="MD-GFRD\ 1")
ax1.errorbar(data[:,0]/1000,data[:,2]/pow(10,6),yerr=data[:,6]/pow(10,6),marker='s',markersize=4,color='#B86637',linestyle='None',label="MD-GFRD")
ax1.errorbar(data[:,0]/1000,data[:,1]/pow(10,6),yerr=data[:,5]/pow(10,6),marker='o',markersize=3,color='#377EB8',linestyle='None',label='New scheme')

ax1.set_ylabel (r'$MSD\thinspace [\mu m^2]$', fontsize=10)

text = r"$D = 10$"  + r"$\thinspace\frac{\mu m^2}{s}$"
box=dict(facecolor='white', edgecolor='grey', boxstyle='round')
ax1.text(0.18,0.8,text,fontsize=8,color='grey',bbox=box,horizontalalignment='center',verticalalignment='center',transform = ax1.transAxes)

ax1.set_yticks([0,0.01,0.02,0.03])
ax2.set_yticks([0,0.005,0.01,0.015])


ax1.text(0.05,0.91,'(a)',fontsize=10,transform = ax1.transAxes)
ax2.text(0.05,0.91,'(b)',fontsize=10,transform = ax2.transAxes)



data = np.loadtxt('diff2.out')

D1=0.01
D2=0.001

ax2.plot ( np.arange(0,600000,1000)/1000, np.arange(0,600000,1000)*6*(D1+D2)/2/pow(10,6),linewidth=0.5,linestyle='--',color='#E41A1C')
ax2.errorbar(data[:,0]/1000,data[:,4]/pow(10,6),yerr=data[:,8]/pow(10,6),marker='x',markersize=8,color='#B8A737',linestyle='None',label='BD')
# ax2.errorbar(data[:,0]/1000,data[:,3]/pow(10,6),yerr=data[:,7]/pow(10,6),marker='^',markersize=4,color='#9337B8',linestyle='None',label="MD-GFRD\ 1")
ax2.errorbar(data[:,0]/1000,data[:,2]/pow(10,6),yerr=data[:,6]/pow(10,6),marker='s',markersize=4,color='#B86637',linestyle='None',label="MD-GFRD")
ax2.errorbar(data[:,0]/1000,data[:,1]/pow(10,6),yerr=data[:,5]/pow(10,6),marker='o',markersize=3,color='#377EB8',linestyle='None',label='New scheme')


ax2.set_xlabel (r'$t\thinspace [\mu s]$', fontsize=10 )
ax2.set_ylabel (r'$MSD\thinspace [\mu m^2]$', fontsize=10)


pl.tick_params('both',labelsize=10)

text = r"$D_1 = 10$" +  r"$\thinspace\frac{\mu m^2}{s}$"+"\n" +r"$D_2 = 1$" + r"$\thinspace\frac{\mu m^2}{s}$"
box=dict(facecolor='white', edgecolor='grey', boxstyle='round')
ax2.text(0.187,0.75,text,fontsize=8,color='grey',bbox=box,horizontalalignment='center',verticalalignment='center',transform = ax2.transAxes)

# ax2.set_xticks([0,10,20,30,40,50,60])

# ax1.set_ylim(0,0.0035)

pl.legend(loc=4,fontsize=8)

ax1.tick_params('y',labelsize=10)
ax2.tick_params('both',labelsize=10)


pl.tight_layout()

f.subplots_adjust(wspace=0.05)
f.subplots_adjust(hspace=0.05)


pl.savefig('diff.eps')
