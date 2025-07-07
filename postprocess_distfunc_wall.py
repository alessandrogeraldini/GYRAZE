### Author: Alessandro Geraldini
### LARGE GYRO-ORBIT MODEL ION VELOCITY DISTRIBUTION IN PLASMA
### AT A SOLID TARGET IN A SHALLOW-ANGLE MAGNETIC FIELD

import numpy as np
import scipy.integrate as scinteg
import scipy.interpolate as scintpol
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.special as scsp
import h5py

plotvx=1
plotspeedangle= 1
plotEthetaphi = 1
wall = 1

alphadeg, gammae, num_spec, niovne, TiovTe, miovme, set_current, jorphi = np.loadtxt('inputfile.txt', unpack = True)
JW, phiW, Phie, Phii, qe, qi = np.loadtxt('misc_output.txt', unpack = True)
alpha = alphadeg*np.pi/180
xMP, phiMP, niMP, neMP = np.loadtxt("phi_n_MP.txt", unpack = True)
vyDSE, muDSE, chiM, DeltaM = np.loadtxt('Fi_W.txt', unpack = True)
F = np.loadtxt("Fi_mpe.txt")
print(F[0][1], F[1][0], "if second value is zero, F is function of (mu, Uminmu), as expected")
Umufile = open("Fi_mpe_args.txt", 'r')
Umu = []
for line in Umufile:
	val = line.split()
	Umu.append([float(v) for v in val])	

Uminmu = Umu[1]
mu = Umu[0]

vx0DSE = [ np.sqrt(2.0*(chiMnum - 0.5*vy*vy - phiMP[0]/TiovTe)) for (chiMnum, vy) in zip(chiM, vyDSE) ]
print(vx0DSE)
print(DeltaM)

###########################################################
########## DEFINE INCOMING DISTRIBUTION FUNCTIONS #########
###########################################################

distfunc = scintpol.RegularGridInterpolator((mu, Uminmu), F, bounds_error=False, fill_value = 0.0)
mufunc = scintpol.RegularGridInterpolator((vyDSE,), muDSE, bounds_error=False, fill_value = 0.0)
invchiMfunc = scintpol.RegularGridInterpolator((chiM,), vyDSE, bounds_error=False, fill_value=0.0)
chiMfunc = scintpol.RegularGridInterpolator((vyDSE,), chiM, bounds_error=False, fill_value=0.0)
vx0DSEfunc = scintpol.RegularGridInterpolator((vyDSE,), vx0DSE, bounds_error=False, fill_value=0.0)
DeltaMfunc = scintpol.RegularGridInterpolator((vyDSE,), DeltaM, bounds_error=False, fill_value = 0.0)

print(distfunc([0.0, 1.0]), distfunc([1.0, 0.0]), "same as before: if the second value is non-zero, distfunc is a function of (mu, Uminmu)")

def tophat(x, h1, h2):
	if (h1<x and x<h2):
		yt = 1.0
	else: 
		yt = 0.0
	return yt

def fDSE(vx, vy, vz):
	TOP = tophat(vx, vx0DSEfunc([vy]), np.sqrt(vx0DSEfunc([vy])**2.0 + 2.0*alpha*vz*DeltaMfunc([vy]))) 
	arg1 = mufunc([vy])[0]
	arg2 = chiMfunc([vy])[0] + 0.5*vz**2 - mufunc([vy])[0]
	y = distfunc([[arg1, arg2]])
	yval = y[0]*4.0*TOP
	if (vy < vyDSE[0]): # otherwise mufunc(vy<vyDSE) would be assigned a value of zero, and the distribution function F(mu, Uminmu) ends up being larger non-zero
		yval = 0.0  # anyway vy > vyDSE is enforced in loops later on
	return yval

phiDS0 = -phiW - phiMP[0]
#phiDS0 = -2.7 - phiMP[0] ## for special case test?

def fWall(vx, vy, vz, phiDSval):
	TOP = tophat(vx, np.sqrt(vx0DSEfunc([vy])[0]**2 - 2.0*phiDSval), np.sqrt(vx0DSEfunc([vy])[0]**2.0 - 2.0*phiDSval + 2.0*alpha*vz*DeltaMfunc([vy])[0])) 
	arg1 = mufunc([vy])[0]
	arg2 = chiMfunc([vy])[0] + 0.5*vz**2 - mufunc([vy])[0]
	#arg2 = arg2[0]
	y = distfunc([[arg1, arg2]])
	yval = y[0]*4.0*TOP
	if (vy < vyDSE[0]): # otherwise mufunc(vy<vyDSE) would be assigned a value of zero, and the distribution function F(mu, Uminmu) ends up being larger non-zero
		yval = 0.0  # anyway vy > vyDSE is enforced in loops later on
	return yval


cbcolor = 'inferno'


#plt.rc('text', usetex=True)
#plt.rc('font',**{'family':'serif','serif':['Palatino']})


########################################################
#################    PLOT    ###########################
########################################################

if plotvx==1:
	vxarr =   [ 0.1*i for i in range(60) ]

	fxarr = []
	fvxarr = []
	fxbarvxarr = []
	for vx in vxarr:
		print("vx = ", vx)
		fyarr = []
		fxbararr = []
		fxbarxbararr = []
		vyarr = [ vyDSE[0] + 0.2*i for i in range(25) ]
		for vy in vyarr:
			fzarr = []
			vzarr =   [ 0.2*i for i in range(25) ]
			for vz in vzarr:
				if wall == 1:
					fzarr.append(fWall(vx, vy, vz, phiDS0/TiovTe))
				else:	
					fzarr.append(fDSE(vx, vy, vz))
				
			fyarr.append(fzarr)
			fxbararr.append(scinteg.trapezoid(fzarr, vzarr))

		fxarr.append(fyarr)
		fvxarr.append(scinteg.trapezoid(fxbararr, vyarr))

	fvxarr = np.array(fvxarr)
	vxarr = np.array(vxarr)
	density = np.trapezoid(fvxarr, vxarr)
	flux = np.trapezoid(fvxarr*vxarr, vxarr)
	flow = flux/density
	Bohm = np.trapezoid(fvxarr/(vxarr**2), vxarr)/density
	threevxmin4moment = (1.0/(2*TiovTe)**2)*np.trapezoid(3*fvxarr/(vxarr**4), vxarr)/density

	np.savetxt('data_f0x_simplemodel_alpha'+str(alphadeg)+'TiovTe'+str(TiovTe), (vxarr, fvxarr))
	plt.xlim(0.0, 4.0)
	fvxarr /= max(fvxarr)
	plt.plot(vxarr, fvxarr, c= 'g')
	plt.axvline(vx0DSE[0], lw = 0.5, ls = '--')
	plt.text(-0.04, 0.34, r'$\sim \hat{v}_{x,E}$', fontsize = 14)
	plt.text(0.27, 0.48, r'$\sim \sqrt{\hat{v}_{x,E}^2 + D \alpha }$', fontsize = 14)
	plt.annotate("",
		    xy=(0.0, 0.4),
		    xytext=(0.45, 0.4),
		    arrowprops=dict(arrowstyle="<->",
				    connectionstyle="arc3", color='k', lw=1),
		    )
	plt.annotate("",
		    xy=(0.0, 0.43),
		    xytext=(1.57, 0.43),
		    arrowprops=dict(arrowstyle="<->",
				    connectionstyle="arc3", color='k', lw=1),
		    )
	plt.title(r'$(b)$', fontsize = 14)
	plt.xlabel(r'$v_x/v_{\rm t,i}$', fontsize = 14)
	plt.ylabel(r'$f_x(v_x)/f_{x,\rm max}$', fontsize = 14)
	plt.savefig('fig-f0x_alpha'+str(alphadeg)+'TiovTe'+str(TiovTe)+'.pdf', dpi = 200)
	with h5py.File('fxyz_Wall.h5', 'w') as f:
		f.create_dataset('distribution function', data=fxarr)
		f.create_dataset('vx', data=vxarr)
		f.create_dataset('vy', data=vyarr)
		f.create_dataset('vz', data=vzarr)


if plotspeedangle==1:
	#f_energyangle = open('f-energyangle-alpha'+str(alpha)+'-TiovTe'+str(TiovTe)+'.txt', 'w') 
	#Earr = np.arange(0.0, 9.0*(1.0+1.0/TiovTe), 0.1)
	Emax = - phiDS0/TiovTe - phiMP[0]/TiovTe + 10.0
	Emin = - phiDS0/TiovTe - phiMP[0]/TiovTe 
	dtheta = np.pi/180.0
	deltaE = 0.4
	Earr = np.arange(Emin, Emax, deltaE)
	fW = []
	fK = []
	phiangarr = []
	thetaarr = []
	Earrmesh = []
	fWfull = []
	i = 0
	for E in Earr:
		fWfull.append([])
		phiangarr.append([])
		print("energy = ", E)
		#vyarr = np.arange(vyDSE[0] , invchiMfunc([E + phiDS0/TiovTe + phiMP[0]/TiovTe])[0], 0.05)
		j=0
		
		thetaarr.append(np.linspace(np.arcsin(vyDSE[0]/np.sqrt(2*E)), np.pi/2, 40))
		Earrmesh.append([])
		for theta in thetaarr[i]:
			fWfull[i].append([])
			Earrmesh[i].append(Earr[i])
			phiangarr[i].append([])
			vx = np.sqrt(2.0*E)*np.cos(theta)
			#print("energy = ", E, "angle = ", theta)
			#if (np.sqrt(2*E)*np.sin(theta) > vyDSE[0]):
			#vyarr = np.arange(vyDSE[0], np.sqrt(2*E)*np.sin(theta), 0.05)
			phiangarr[i][j] = np.linspace(np.arctan(vyDSE[0]/np.sqrt(0.0000001+ 2.0*E - vyDSE[0]**2 - vx**2)), np.pi/2, 40)
			intgrd = []
			#for vy in vyarr:
			for phiang in phiangarr[i][j]:
				#vz = np.sqrt(2.0*E - 2.0*chiMfunc([vy])[0] + 2.0*phiMP[0]/TiovTe)
				#vz = np.sqrt(2.0*E - vy**2 - vx**2)
				vz = np.sqrt(2.0*E)*np.sin(theta)*np.cos(phiang)
				vy = np.sqrt(2.0*E)*np.sin(theta)*np.sin(phiang)
				#phiang = np.arctan(vy/vz)
				#phiangarr[i][j].append(phiang)
				J = np.sqrt(2.0*E)*np.sin(theta)
				#J = np.sqrt(2.0*E)*np.sin(theta)/vz
				if (wall == 1):
					fval = fWall(vx, vy, vz, phiDS0)*J
					fWfull[i][j].append(fval)
					intgrd.append(fval)
				else:
					fval = fDSE(vx, vy, vz)*J
					fWfull[i][j].append(fval)
					intgrd.append(fval)

			#print(intgrd)
			#intgrd = np.array(intgrd)
			#intgrd = [intgr for intgr in intgrd]
			#fW1 = scinteg.trapezoid(intgrd, vyarr)
			fW1 = scinteg.trapezoid(intgrd, phiangarr[i][j])
			#else:
				#fW1 = 0.0

			fW.append(fW1)
			#f_energyangle.write(str(fW1)+' ')
			j+=1
				

		#f_energyangle.write('\n\n')
		i+=1

	#f_energyangle.close()
	Earr = np.array(Earr)
	Earrmesh = np.array(Earrmesh)
	thetaarr = np.array(thetaarr)
	fW = np.array(fW)
	fWfull = np.array(fWfull)
	fWfull /= np.max(fWfull)
	print(fW.shape)
	#fW = [np.log(f) for f in fW]
	fW = np.reshape(fW, (len(Earr), len(thetaarr[0])))
	#print(fW)
	#cp = plt.contourf(thetaarr, Earr, fW)
	fig1, ax2 = plt.subplots()
	#CS = ax2.contourf(X, Y, Z, 10, cmap=plt.cm.bone, origin=origin)
	maxfW = np.max(fW)
	cmap1 = plt.get_cmap(cbcolor)
	norm1 = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
	print(thetaarr.shape, Earrmesh.shape, fW.shape)
	thetaarrdeg = [t*180/np.pi for t in thetaarr]
	cp = ax2.contourf(thetaarrdeg, Earrmesh, fW/maxfW, 10, cmap = cmap1, norm = norm1)
	ax2.patch.set_facecolor(cmap1(1.0/20))
	#cp.set_clim(vmin = 0.0, vmax = 1.0)
	plt.colorbar(cp)
	#plt.clim(0.0, 1.0)
	plt.title('energy('+r'$E$'+')-angle('+r'$\theta$'+') distribution', fontsize = 20)
	plt.ylabel(r'$E / T_{\rm i}$', fontsize = 20)
	plt.xlabel(r'$\theta (^{\circ})$', fontsize = 20)
	plt.xticks(fontsize = 16)
	plt.yticks(fontsize = 16)
	#plt.xlim(0.0, 90.0)
	plt.ylim(0.0, Emax)
	plt.tight_layout()
	plt.savefig("fig-energyangle_alpha="+str(alphadeg)+"_TiovTe="+str(TiovTe)+".pdf", dpi = 100)
	np.savetxt('data_fspeedangle_'+str(alphadeg)+'_TiovTe'+str(TiovTe), fW)
	phiangarr = np.array(phiangarr)
	print(phiangarr.shape)
	with h5py.File('fEthetaphiv2.h5', 'w') as f:
		f.create_dataset('distribution function', data=fWfull)
		f.create_dataset('E', data=Earr)
		f.create_dataset('theta', data=thetaarr)
		f.create_dataset('phi', data=phiangarr)
	

if plotEthetaphi==1:
	Emax = - phiDS0/TiovTe - phiMP[0]/TiovTe + 10.0
	Emin = - phiDS0/TiovTe - phiMP[0]/TiovTe 
	dtheta = np.pi/90.0
	thetaarr = np.arange(0.0, np.pi/2.0, dtheta) 
	phiarr = np.arange(0.0, np.pi/2.0, dtheta) 
	deltaE = 0.4
	Earr = np.arange(Emin, Emax, deltaE)
	fE = []
	for E in Earr:
		print("energy = ", E)
		ftheta = []
		for theta in thetaarr:
			fphi = []
			for phiang in phiarr:
				vz = np.sqrt(2.0*E)*np.sin(theta)*np.cos(phiang)
				vy = np.sqrt(2.0*E)*np.sin(theta)*np.sin(phiang)
				vx = np.sqrt(2.0*E)*np.cos(theta)
				J  = np.sqrt(2.0*E)*np.sin(theta)
				if (wall == 1):
					fphi.append(fWall(vx, vy, vz, phiDS0)*J)
				else:
					fphi.append(fDSE(vx, vy, vz)*J)

			ftheta.append(fphi)

		fE.append(ftheta)

	with h5py.File('fEthetaphi.h5', 'w') as f:
		f.create_dataset('distribution function', data=fE)
		f.create_dataset('E', data=Earr)
		f.create_dataset('theta', data=thetaarr)
		f.create_dataset('phi', data=phiarr)


