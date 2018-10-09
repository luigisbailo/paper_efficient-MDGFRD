import matplotlib
import numpy as np
import matplotlib.pyplot as pl

pl.style.use(['style1'])
 

f  = pl.figure()

f.set_figheight(3.5)
f.set_figwidth(6.69)

ax2 = f.add_subplot(3,2,1)
ax2.tick_params(labelbottom='off')
ax3 = f.add_subplot(3,2,2,sharey = ax2)
ax3.tick_params(labelleft='off')
ax3.tick_params(labelbottom='off')
ax4 = f.add_subplot(3,2,3,sharex = ax2)
ax4.tick_params(labelbottom='off')
ax5 = f.add_subplot(3,2,4,sharey = ax4)
ax5.tick_params(labelleft='off')
ax5.tick_params(labelbottom='off')
ax5.set_yticks([1,1.04,1.08])
ax6 = f.add_subplot(3,2,5,sharex = ax4)
ax7 = f.add_subplot(3,2,6,sharey = ax6)
ax7.tick_params(labelleft='off')
ax3.set_yticks([1,1.04,1.08])

ax2.text(0.08,0.82,'(a)',fontsize=10,transform = ax2.transAxes)
ax3.text(0.08,0.82,'(b)',fontsize=10,transform = ax3.transAxes)
ax4.text(0.08,0.82,'(c)',fontsize=10,transform = ax4.transAxes)
ax5.text(0.08,0.82,'(d)',fontsize=10,transform = ax5.transAxes)
ax6.text(0.08,0.82,'(e)',fontsize=10,transform = ax6.transAxes)
ax7.text(0.08,0.82,'(f)',fontsize=10,transform = ax7.transAxes)


ax6.set_xlabel(r"$\alpha$",fontsize=10)
ax7.set_xlabel(r"$\alpha$",fontsize=10)
ax4.set_ylabel(r"$Relative\ CPU\ time$",fontsize=10)

ax6.set_xticks([6,8,10,12,14,16])


xlimit = [5,17]
ax2.set_xlim(xlimit)
ax3.set_xlim(xlimit)
ax4.set_xlim(xlimit)
ax5.set_xlim(xlimit)
ax6.set_xlim(xlimit)
ax7.set_xlim(xlimit)

ax2.set_ylim(0.995,1.04)
ax2.set_yticks ([1,1.01,1.02,1.03])

ax4.set_ylim(0.99,1.08)
ax4.set_yticks ([1,1.02,1.04,1.06])

ax6.set_ylim(0.97,1.2)
ax6.set_yticks ([1,1.05,1.1,1.15])

# ax6.set_ylim(0.98,1.14)
# ax6.set_yticks ([1,1.04,1.08,1.12])

f.subplots_adjust(wspace=0.05)
f.subplots_adjust(hspace=0.05)

text = 'text'

for i in range (2,8):

	if (i==2):
		ax = ax2
		text = r'$c = 10^{-2} \mu M$'
	if (i==3):
		ax = ax3
		text = r'$c = 10^{-1} \mu M$'
	if (i==4):
		ax = ax4
		text = r'$c = 10^{0} \mu M$'
	if (i==5):
		ax = ax5
		text = r'$c = 10^{1} \mu M$'		
	if (i==6):
		ax = ax6
		text = r'$c = 10^{2} \mu M$'
	if (i==7):
		ax = ax7
		text = r'$c = 10^{3} \mu M$'

	if (i==8):
		ax = ax8

	ax2.tick_params('y',labelsize=10)
	ax4.tick_params('y',labelsize=10)
	ax6.tick_params('both',labelsize=10)
	ax7.tick_params('x',labelsize=10)

	filename = str(i) + ".out"
	data=np.loadtxt(filename)	
	param = np.loadtxt (str(i) + ".data")
	ax.plot (data[:,0],data[:,1]/np.min(data[:,1]),color='#377EB8',marker='.',linestyle='None')
	box=dict(facecolor='white', edgecolor='grey', boxstyle='round')
	ax.text(0.365,0.82,text,fontsize=8,color='grey',bbox=box,horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
	


pl.savefig ('alphaPerf.eps')