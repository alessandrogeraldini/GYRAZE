import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as scintpol
import scipy.integrate as scint
import matplotlib
import os
#plt.rcParams.update({'font.size': 10}) 
#plt.rc('xtick', labelsize=18)  
#plt.rc('ytick', labelsize=18) 


if (os.path.isfile("inputfile.txt")):
	alphadeg, gammae, num_spec, niovne, TiovTe, miovme, set_current, jorphi = np.loadtxt('inputfile.txt', unpack = True)
elif (os.path.isfile("input_physparams.txt")):
	file = open("input_physparams.txt")
	gammaflag = 0
	counter = 0
	listline = []
	for line in file:
		listlineold = listline
		listline = line.split()
		if (len(listlineold) > 1) and ( (listlineold[1] == "gamma_ref") or (counter == 2) ):
			gamma = float(listline[0])
		if (len(listlineold) > 1) and ( (listlineold[1] == "alphadeg") or ( (counter == 0) and (str.isnumeric(line[0]) == True) ) ):
			alphadeg = float(listline[0])
		if (len(listlineold) > 1) and ( (listlineold[1] == "TioverTe") or (counter == 5) ):
			TiovTe = float(listline[0])
		if (len(listlineold) > 0) and (str.isnumeric(line[0]) == True):
			counter += 1

		##print(line)
else:
	prinf("ERROR: no input file found. Exit code now")
	exit(1)

xx, phi, ni, ne = np.loadtxt('phi_n_MP.txt', unpack = True)

gg = [ np.sqrt(x) for x in xx] 
chargedensitymp = [ nnii - nnee for (nnee, nnii) in zip(ne, ni) ]

phipp = [0.0]
phippp = []
phip = []

for i in range(len(phi)-2):
	phip.append((phi[i+1] - phi[i])/(xx[i+1]-xx[i]))
	phipp.append(2.0*( (phi[i+2] - phi[i+1])/(xx[i+2] - xx[i+1]) - (phi[i+1] - phi[i])/(xx[i+1] - xx[i]) ) / (xx[i+2] - xx[i]))

for i in range(len(phi)-2):
	phippp.append((phipp[i+1] - phipp[i])/(xx[i+1]-xx[i]))


phip.append(0.0)
phip.append(0.0)
phippp.append(0.0)
phippp.append(0.0)
phipp.append(phipp[len(phipp)-1])
phipp[0] = phipp[1]

phipar = [- 0.25*0.5*(xval - 3.5)**2 - 0.4 for xval in xx]
#phipar = [- 0.5*(xval - 2.0)**2 - 0.8 for xval in xx]
#plt.plot(xx, phippp)
plt.plot(xx, phipp)
plt.plot(xx, phip)
plt.plot(xx, phi)
plt.plot(xx, phipar)
plt.xlim(0.0, 12.0)
plt.ylim(-3.0, 3.0)
#plt.show()
plt.clf()

#flag = 0
#xx, phi = np.loadtxt('phispline.txt', unpack = True)
#x_orbit = []
#vx_orbit = []
#vxsq_orbit = [ 0.0 - 0.5*(xval - 2.5)**2 - phival for (xval, phival) in zip(xx, phi) ]
#for i in range(len(xx)):
#	xval = xx[i]
#	vxsqval = vxsq_orbit[i]
#	if (vxsqval > 0.0):
#		vx_orbit.append(np.sqrt(vxsqval))
#		x_orbit.append(xval)
#		if (flag == 0):
#			xbottom = xval - (xval-xx[i-1])*vxsqval/(vxsqval - vxsq_orbit[i-1])
#			x_orbit.insert(0, xbottom)
#			vx_orbit.insert(0, 0.0)
#			flag = 1
#	else:
#		if (flag == 1):
#			xtop = xval - (xval-xx[i-1])*vxsqval/(vxsqval - vxsq_orbit[i-1])
#			x_orbit.append(xtop)
#			vx_orbit.append(0.0)
#		flag = 0
#
#minusvx_orbit = [-vxval for vxval in vx_orbit]
#mu = (1.0/np.pi)*np.trapz(vx_orbit, x_orbit)
#print(mu)
#plt.plot(x_orbit, vx_orbit)
#plt.plot(x_orbit, minusvx_orbit)
#plt.show()
#plt.clf()
#exit()

#xx, phi = np.loadtxt('phi_MP_4deg.txt', unpack = True)
#xxs, phis = np.loadtxt('phi_DS_4deg.txt', unpack = True)
fMP = scintpol.interp1d(xx, phi);
niMP = scintpol.interp1d(xx, ni);
neMP = scintpol.interp1d(xx, ne);

resol = 10.0
flagMP = 0
plt.plot(xx, phi)
plt.plot(xx, phi, '.')
plt.xlabel(r'$x/\rho_B$', fontsize = 20)
plt.ylabel(r'$e\phi_{mp}/T_e$', fontsize = 20)
plt.xlim(0.0, 10.0)
plt.tight_layout()
plt.savefig('phi_MP.pdf', dpi = 200)
plt.clf()
plt.plot(gg, phi)
plt.xlabel(r'$\sqrt{x/\rho_B}$', fontsize = 20)
plt.ylabel(r'$e\phi_{mp}/T_e$', fontsize = 20)
plt.tight_layout()
plt.savefig('phi_sqrtx_MP.pdf', dpi = 200)
plt.clf()
plt.plot(gg, ni, '-')
plt.plot(gg, ne, '--')
plt.plot(gg, chargedensitymp, '-.')
plt.xlabel(r'$x/\rho_B$', fontsize = 20)
plt.ylabel(r'$n_{mp}/n_{\infty}$', fontsize = 20)
plt.ylim(-0.05, 1.0)
plt.tight_layout()
plt.savefig('ng_MP.pdf', dpi = 200)
plt.clf()
plt.plot(xx, ni, '-')
plt.plot(xx, ni, '.')
plt.plot(xx, ne, '--')
plt.plot(xx, chargedensitymp, '-.')
plt.xlabel(r'$x/\rho_B$', fontsize = 20)
plt.ylabel(r'$n_{mp}/n_{\infty}$', fontsize = 20)
plt.ylim(-0.05, 1.0)
plt.tight_layout()
plt.savefig('nx_MP.pdf', dpi = 200)
plt.clf()
plt.plot(phi, ni)
plt.plot(phi, ne, '--')
plt.xlabel(r'$n_{mp}/n_{\infty}$', fontsize = 20)
plt.ylabel(r'$e\phi_{mp}/T_e$', fontsize = 20)
plt.tight_layout()
plt.savefig('nphi_MP.pdf', dpi = 200)
plt.clf()

fig, ax = plt.subplots(2)
ax[0].plot(xx, phi)
ax[0].set_xlabel(r'$x/\rho_B$')
ax[0].set_ylabel(r'$e\phi_{mp}/T_e$')
#ax[0].set_xlim(0.0, 9.0)

ax[1].plot(xx, ni, label = r'$n_{i,mp}$')
ax[1].plot(xx, ne, '--', label = r'$n_{e,mp}$')
ax[1].plot(xx, chargedensitymp, '-.', label = r'$n_{i}-n_{e}$')
ax[1].set_xlabel(r'$x/\rho_B$')
ax[1].set_ylabel('$n_{mp}/n_{\infty}$')
#ax[1].set_xlim(0.0, 9.0)
ax[1].set_ylim(-0.05, 1.05)
ax[1].legend(loc = 'best', fontsize = 8.5)
plt.tight_layout()
plt.savefig('phi_n_mp_2.pdf', dpi = 200)
plt.clf()


Te = 1.0
Ts = 1.0 + Te
plt.plot(xx, phi)

def xifp(x, phi):
	y = 1.0/np.sqrt( -3.0 -2.0*phi + 4.0*np.exp(phi) - np.exp(2.0*phi) )
	return y	

numb = 4.0
#phi0 = np.log(1.12*alpha)
phi0 = -2.15
xarr = np.arange(0.0, 20.0, 0.01)
phiarr = np.arange(phi0, -0.000000000001, 0.0001)
phiarrf = [p + (1.0/numb)*np.exp(numb*phi0)*(np.exp(-numb*p)-1.0) for p in phiarr]
#phif = scint.odeint(phip, np.log(alpha), xarr)
#xf = scint.odeint(xip, -0.5/np.sqrt(2.0*np.log(alpha)), phiarr)
#phif = scint.odeint(phip, np.log(alpha), xarr)
xif = scint.odeint(xifp, 0.0, phiarrf)
xif  = np.resize(xif, len(xif))
f1 = scintpol.interp1d(xif, phiarr)
phi1 = f1(xarr)
#xarr /= np.sqrt(-2.0*np.log(alpha))
plt.plot(xarr, phi1, ls = '--', lw = 1.0)
plt.xlim(0.0, 20.0)
#plt.ylim(-1.0, -0.0)
plt.savefig('fig-phiChodura.pdf', dpi = 200)
plt.clf()

