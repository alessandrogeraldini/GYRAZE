import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as scint
import matplotlib
import os

#plt.rcParams.update({'font.size': 10}) 
#plt.rc('xtick', labelsize=18)  
#plt.rc('ytick', labelsize=18) 

if (os.path.isfile("inputfile.txt")):
	alphadeg, gamma, num_spec, niovne, TiovTe, miovme, set_current, jorphi = np.loadtxt('inputfile.txt', unpack = True)
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

lambdaDoverrhoi= 0.1
#gg, phi = np.loadtxt('phidata_compare.txt', unpack = True)
xx, phi, ni, ne = np.loadtxt('phi_n_MP.txt', unpack = True)
xxs, phis, nis, nes = np.loadtxt('phi_n_DS.txt', unpack = True)

phispp = [ -(1/gamma**2)*(phis[i+2] - 2*phis[i+1] + phis[i])/(xxs[1]**2) for i in range(len(phis)-2) ]
phispp.insert(0, phispp[0])
phispp.append(0.0)

chargedensity = [ nnii - nnee for (nnee, nnii) in zip(nes, nis) ]
chargedensitymp = [ nnii - nnee for (nnee, nnii) in zip(ne, ni) ]

#xx, phi = np.loadtxt('phi_MP_4deg.txt', unpack = True)
#xxs, phis = np.loadtxt('phi_DS_4deg.txt', unpack = True)
xxs_D = [gamma*xval for xval in xxs]
phis_scaled = [phi[0] - p for p in phis]
xxs_scaled = [lambdaDoverrhoi*x for x in xxs]
nis_scaled = [n for n in nis]
nes_scaled = [n for n in nes]
xxs_scaled.append(100.0)
nis_scaled.append(1.0)
nes_scaled.append(1.0)
phisapp = np.append(phis, 0.0)
EW = (phis[1] - phis[0])/(xxs_D[1] - xxs_D[0])
phismodel = [ phis[0]*np.exp(EW*x/phis[0]) for x in xxs_D ]

fMP = scint.interp1d(xx, phi);
fDS = scint.interp1d(xxs_scaled, phisapp);
niMP = scint.interp1d(xx, ni);
niDS = scint.interp1d(xxs_scaled, nis_scaled);
neMP = scint.interp1d(xx, ne);
neDS = scint.interp1d(xxs_scaled, nes_scaled);

xx_match = []
phi_match = []
ni_match =[]
ne_match =[]
resol = 10.0
flagDS = 0
flagMP = 0
niDSvalold = 0.0
neDSvalold = 0.0
for i in range(int(resol*10.0/lambdaDoverrhoi)):
	x = i*lambdaDoverrhoi/resol
	neDSval = neDS(x)
	niMPval = niMP(x)
	if (i != 0):
		if ( (abs(neDSval-neDSvalold) > 0.1) or (flagDS==1)):
			 neDSval = 1.0
			 flagDS = 1
		if ( (abs(niMPval-niMPvalold) > 0.1) or (flagMP == 1)):
			 niMPval = 1.0
			 flagMP = 1
	xx_match.append(x)
	phi_match.append(fMP(x) + fDS(x))
	ni_match.append(niDS(x)*niMPval)
	ne_match.append(neDSval*neMP(x))
	#print(ne_match[i], neMP(x))
	#print(ni_match[i])
	niMPvalold = niMPval
	neDSvalold = neDSval
	

chargedensity_match = [ nnii - nnee for (nnee, nnii) in zip(ne_match, ni_match) ]
chargedensity_match = [ 0.0 for (nnee, nnii) in zip(ne_match, ni_match) ]

plt.plot(xx_match, ni_match)
plt.plot(xx_match, ne_match)
#plt.xlim(0.0, 5.0)
plt.xlabel(r'$x/\rho_i$', fontsize = 20)
plt.ylabel(r'$n/n_{\infty}$', fontsize = 20)
plt.xlim(0.0, 9.0)
plt.tight_layout()
plt.savefig('n_MS.pdf', dpi = 200)
plt.clf()

plt.plot(xx_match, phi_match)
#plt.xlim(0.0, 5.0)
plt.xlabel(r'$x/\rho_i$', fontsize = 20)
plt.ylabel(r'$e(\phi_{ds} + \phi_{mp})/T_e$', fontsize = 20)
plt.tight_layout()
plt.savefig('phi_MS.pdf', dpi = 200)
plt.clf()

gg = [x**0.5 for x in xx]
#print(xx, phi)
phipinfMP = (0.025138  - 0.018802)/(9.225348- 8.775355)
#phipinfMP = phi
phi0 = 0.022144*9.0**4
#print(phi0)

#xan = np.arange(0.0, 15.0)
#phian = [-phi0/(x)**4 for x in xan]

plt.plot(xxs_D, phis)
plt.plot(xxs_D, phismodel, '--')
plt.xlabel(r'$x/\lambda_D$', fontsize = 20)
plt.ylabel(r'$e\phi_{ds}/T_e$', fontsize = 20)
plt.tight_layout()
plt.savefig('phi_DS.pdf', dpi = 200)
plt.clf()
plt.plot(xxs_D, phis)
plt.xlabel(r'$x/\lambda_D$', fontsize = 20)
plt.ylabel(r'$e\phi_{ds}/T_e$', fontsize = 20)
plt.ylim(-0.1, 0.0)
plt.tight_layout()
plt.savefig('phi_DSzoom.pdf', dpi = 200)
plt.clf()
plt.plot(xx, phi)
plt.xlabel(r'$x/\rho_B$', fontsize = 20)
plt.ylabel(r'$e\phi_{mp}/T_e$', fontsize = 20)
plt.tight_layout()
plt.savefig('phi_MP.pdf', dpi = 200)
plt.clf()
plt.plot(gg, phi, 'o')
plt.plot(gg, phi)
plt.xlim(0.0, 0.5)
plt.ylim(-2.1, -1.6)
plt.xlabel(r'$\sqrt{x/\rho_B}$', fontsize = 20)
plt.ylabel(r'$e\phi_{mp}/T_e$', fontsize = 20)
plt.tight_layout()
plt.savefig('phi_sqrtx_MP.pdf', dpi = 200)
plt.clf()
plt.plot(xxs_D, nis, label = r'$n_{i}$')
plt.plot(xxs_D, nes, '--', label = r'$n_{e}$')
plt.plot(xxs_D, chargedensity, '-.', label = r'$n_{i}-n_{e}$')
plt.plot(xxs_D, phispp, ':', label = r'$-\partial_{xx} \phi$')
plt.ylabel(r'$n_{ds}/n_{ds\infty}$', fontsize = 20)
plt.xlabel(r'$x/\lambda_D$', fontsize = 20)
plt.ylim(-0.05, 1.0)
plt.legend(loc = 'best', fontsize = 14)
plt.tight_layout()
plt.savefig('nx_DS.pdf', dpi = 200)
plt.clf()
plt.plot(gg, ni, label = r'$n_{i}$')
plt.plot(gg, ni, 'o')
plt.plot(gg, ne, '--', label = r'$n_e$')
plt.plot(gg, chargedensitymp, '-.', label = r'$n_i - n_e$')
plt.xlabel(r'$x/\rho_B$', fontsize = 20)
plt.ylabel(r'$n_{mp}/n_{\infty}$', fontsize = 20)
plt.ylim(-0.05, 1.0)
plt.tight_layout()
plt.savefig('ng_MP.pdf', dpi = 200)
plt.clf()
plt.plot(xx, ni, label = r'$n_{i}$')
plt.plot(xx, ne, '--', label = r'$n_e$')
plt.plot(xx, chargedensitymp, '-.', label = r'$n_i - n_e$')
plt.xlabel(r'$x/\rho_B$', fontsize = 20)
plt.ylabel(r'$n_{mp}/n_{\infty}$', fontsize = 20)
plt.ylim(-0.05, 1.0)
plt.tight_layout()
plt.savefig('nx_MP.pdf', dpi = 200)
plt.clf()
plt.plot(phis, nis)
plt.plot(phis, nes, '--')
plt.xlabel(r'$n_{ds}/n_{\infty}$', fontsize = 20)
plt.ylabel(r'$e\phi_{ds}/T_e$', fontsize = 20)
plt.tight_layout()
plt.savefig('nphi_DS.pdf', dpi = 200)
plt.clf()
plt.plot(phi, ni)
plt.plot(phi, ne, '--')
plt.xlabel(r'$n_{mp}/n_{\infty}$', fontsize = 20)
plt.ylabel(r'$e\phi_{mp}/T_e$', fontsize = 20)
plt.tight_layout()
plt.savefig('nphi_MP.pdf', dpi = 200)
plt.clf()

fig, ax = plt.subplots(2, 2)
ax[0][0].plot(xxs_D, phis)
ax[0][0].set_xlabel(r'$x/\lambda_D$')
ax[0][0].set_ylabel(r'$e\phi_{ds}/T_e$')
ax[0][0].set_xlim(0.0, 9.0)
ax[0][0].set_ylim(phi_match[0], 0.0)
ax[0][1].plot(xx, phi)
ax[0][1].set_xlabel(r'$x/\rho_B$')
ax[0][1].set_ylabel(r'$e\phi_{mp}/T_e$')
ax[0][1].set_xlim(0.0, 9.0)
ax[0][1].set_ylim(phi_match[0], 0.0)

#ax[2][0].plot([x**0.5 for x in xx_match], phi_match)
#ax[0][2].plot(xx_match, phi_match)
#ax[0][2].set_xlabel(r'$x/\rho_B$')
#ax[0][2].set_ylabel(r'$e\phi_{ms}/T_e$')
#ax[0][2].set_xlim(0.0, 9.0)
#ax[0][2].set_ylim(phi_match[0], 0.0)
##ax[2][0].set_xlim(0.0, 3.0)
#ax[1][2].plot(xx_match, ni_match, label = r'$n_i$')
#ax[1][2].plot(xx_match, ne_match, '--', label = r'$n_e$')
#ax[1][2].plot(xx_match, chargedensity_match, '-.', label = r'$n_{i}-n_{e}$')
#ax[1][2].set_xlabel(r'$x/\rho_B$')
#ax[1][2].set_ylabel(r'$n_{ms}/n_{\infty}$')
#ax[1][2].set_xlim(0.0, 9.0)
ax[1][0].plot(xxs_D, nis, label = r'$n_{\rm i,ds}/n_{\rm i,ds}(\infty)$')
ax[1][0].plot(xxs_D, nes, '--', label = r'$n_{\rm e, ds}/n_{\rm i,ds}(\infty)$')
ax[1][0].plot(xxs_D, chargedensity, '-.', label = r'$(n_{\rm i,ds}-n_{\rm e,ds})/n_{\rm i,ds}(\infty)$')
ax[1][0].plot(xxs_D, phispp, ':', label = r'$-\epsilon_0 \partial_{xx} \phi_{\rm ds} / (en_{\rm i,ds}(\infty))$')
ax[1][0].set_ylabel('term in Poisson\'s equation')
ax[1][0].set_xlabel(r'$x/\lambda_D$')
ax[1][0].set_xlim(0.0, 9.0)
ax[1][1].set_ylim(-0.05, 1.05)
ax[1][1].plot(xx, ni, label = r'$n_{\rm i,mp}/n_{\rm i,mp}(\infty)$')
ax[1][1].plot(xx, ne, '--', label = r'$n_{\rm e,mp}/n_{\rm i,mp}(\infty)$')
ax[1][1].plot(xx, chargedensitymp, '-.', label = r'$(n_{\rm i,mp}-n_{\rm e,mp})/n_{\rm i,mp}(\infty)$')
ax[1][1].set_xlabel(r'$x/\rho_{\rm B}$')
ax[1][1].set_ylabel('term in quasineutrality')
ax[1][1].set_xlim(0.0, 9.0)
ax[1][1].set_ylim(-0.05, 1.05)
ax[1][0].set_ylim(-0.05, 1.05)
ax[1][0].legend(loc = 'best', fontsize = 8.5)
ax[1][1].legend(loc = (0.22, 0.07), fontsize = 8.5)
#ax[1][1].legend(loc = 'best', fontsize = 8.5)
#ax[1][2].legend(loc = 'best', fontsize = 8.5)
plt.tight_layout()
plt.savefig('phi_n_mpds_2x2.pdf', dpi = 200)
plt.clf()

matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 
matplotlib.rc('axes', labelsize=12) 
fig, ax = plt.subplots(2)
ax[0].plot(xx_match, phi_match, lw = 2.0)
ax[0].set_xlabel(r'$x/\rho_B$')
ax[0].set_ylabel(r'$e\phi/T_e$')
ax[0].set_xlim(0.0, 3.0)
ax[0].set_ylim(phi_match[0], 0.0)
ax[1].plot(xx_match, ni_match, label = r'$n_i$', lw = 2.0)
ax[1].plot(xx_match, ne_match, '--', label = r'$n_e$', lw = 2.0)
ax[1].plot(xx_match, chargedensity_match, '-.', lw = 2.0, label = r'$n_{i}-n_{e}$')
ax[1].set_xlabel(r'$x/\rho_B$')
ax[1].set_ylabel(r'$n/n_{\infty}$')
ax[1].set_xlim(0.0, 3.0)
ax[1].legend(loc = 'best', fontsize = 12)
plt.tight_layout()
plt.savefig('phi_n_ms_2x1.pdf', dpi = 200)
plt.clf()
#plt.plot(xx, phi)
#phipar = [-27.0/(x)**4 for x in xx]
#plt.plot(xx, phipar)
#plt.show()
