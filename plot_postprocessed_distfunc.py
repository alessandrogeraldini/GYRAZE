import numpy as np
import h5py
import scipy.interpolate as scintpol
import matplotlib.pyplot as plt
import matplotlib as mpl
cbcolor = 'inferno'

version_angle = 1
energy_input = float(input("input the desired energy: "))

f = h5py.File('fxyz_Wall.h5', 'r')
print(f)
distfunc = f['distribution function']
vx = f['vx']
vy = f['vy']
vz = f['vz']
print(vx[30], vy[3], vz[15])
Eval = 0.5*(vx[30]**2 + vy[3]**2 + vz[15]**2) 
thetaval = np.atan(np.sqrt(vy[3]**2 + vz[15]**2)/vx[30])
phival = np.atan(vy[3]/vz[15])
print(Eval, thetaval, phival)
print(distfunc[30][3][15])

fangle = h5py.File('fEthetaphi.h5', 'r')
distfuncangle = fangle['distribution function']
E = fangle['E']
theta = fangle['theta']
phi = fangle['phi']
k=0
while (E[k] < energy_input):
	k+=1
print(E[k])
distfuncsliced = np.array(distfuncangle[k][:][:])
maxf = np.max(distfuncsliced)
distfuncsliced /= maxf
fig, ax = plt.subplots()
cmap1 = plt.get_cmap(cbcolor)
norm1 = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
phideg = phi[:]*180/np.pi
thetadeg = phi[:]*180/np.pi
cp = ax.contourf(phideg, thetadeg, distfuncsliced, 10, cmap = cmap1, norm = norm1)
ax.patch.set_facecolor(cmap1(1.0/20))
#cp.set_clim(vmin = 0.0, vmax = 1.0)
plt.colorbar(cp)
plt.title('angle('+r'$\theta$-$\phi$'+') distribution, energy ='+r'$E = '+str(E[k]), fontsize = 20)
plt.ylabel(r'$theta (^{\circ})$', fontsize = 20)
plt.xlabel(r'$\phi (^{\circ})$', fontsize = 20)
plt.xlim(0.0, 90.0)
plt.ylim(0.0, 90.0)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.tight_layout()
plt.savefig("fig-angles_v1.pdf", dpi = 100)

i=0
while (phi[i] < phival):
	i+=1
j=0
while (theta[j] < thetaval):
	j+=1
i-=1
j-=1
k=0
while (E[k] < Eval):
	k+=1
print(E[k])
print(theta[j])
print(phi[i])
print(distfuncangle[k][j][i]/(np.sqrt(2*E[k])*np.sin(theta[j])), distfunc[30][3][15], "are they the same?")

Fenangfull = []
for i in range(len(E)):
	Fenang = []
	for j in range(len(theta)):
		Fone = distfuncangle[i][j][:]
		Fenang.append(np.trapezoid(Fone, phi))
	Fenangfull.append(Fenang)
		
fig, ax = plt.subplots()
cmap1 = plt.get_cmap(cbcolor)
norm1 = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
maxf = np.max(Fenangfull)
Fenangfull /= maxf
thetadeg = [t*180/np.pi for t in theta]
cp = ax.contourf(thetadeg, E, Fenangfull, 10, cmap = cmap1, norm = norm1)
ax.patch.set_facecolor(cmap1(1.0/20))
plt.colorbar(cp)
plt.title('energy-angle('+r'$E$-$\theta$'+') distribution', fontsize = 20)
plt.ylabel(r'$E $', fontsize = 20)
plt.xlabel(r'$\theta (^{\circ})$', fontsize = 20)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.ylim(0.0, E[len(E)-1])
plt.xlim(0.0, 90)
plt.tight_layout()
plt.savefig("fig-enangle_v1.pdf")

fanglev2 = h5py.File('fEthetaphiv2.h5', 'r')
distfuncanglev2 = fanglev2['distribution function']
Ev2 = fanglev2['E']
thetav2 = fanglev2['theta']
phiv2 = fanglev2['phi']
k=0
while (Ev2[k] < energy_input):
	k+=1
print(Ev2[k])
distfuncv2sliced = np.array(distfuncanglev2[k][:][:])
phiv2sliced = phiv2[k][:][:]*180/np.pi
thetav2forphi = thetav2[k][:]*180/np.pi
thetav2forphi = [thetav2forphi for t in thetav2forphi]
maxfv2 = np.max(distfuncv2sliced)
distfuncv2sliced /= maxfv2
fig, ax = plt.subplots()
cmap1 = plt.get_cmap(cbcolor)
norm1 = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
print(thetav2forphi, phiv2sliced)
phiv2sliced = np.transpose(phiv2sliced)
distfuncv2sliced = np.transpose(distfuncv2sliced)
cp = ax.contourf(phiv2sliced, thetav2forphi, distfuncv2sliced, 10, cmap = cmap1, norm = norm1)
ax.patch.set_facecolor(cmap1(1.0/20))
#cp.set_clim(vmin = 0.0, vmax = 1.0)
plt.colorbar(cp)
plt.title('angle('+r'$\theta$-$\phi$'+') distribution, energy ='+r'$E = '+str(Ev2[k]), fontsize = 20)
plt.ylabel(r'$theta (^{\circ})$', fontsize = 20)
plt.xlabel(r'$\phi (^{\circ})$', fontsize = 20)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.xlim(0.0, 90.0)
plt.ylim(0.0, 90.0)
plt.tight_layout()
plt.savefig("fig-angles_v2.pdf", dpi = 100)

Fenangfull = []
for i in range(len(Ev2)):
	Fenang = []
	for j in range(len(thetav2[i])):
		Fone = distfuncanglev2[i][j][:]
		phione = phiv2[i][j][:]
		Fenang.append(np.trapezoid(Fone, phione))
	Fenangfull.append(Fenang)
		
fig, ax = plt.subplots()
cmap1 = plt.get_cmap(cbcolor)
norm1 = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
maxf = np.max(Fenangfull)
Fenangfull /= maxf
thetav2deg = [t*180/np.pi for t in thetav2]
Efortheta = [Ev2 for t in thetav2deg[0]]
np.array(Efortheta)
Efortheta = np.transpose(Efortheta)
print(np.shape(Efortheta))
print(np.shape(thetav2deg))
cp = ax.contourf(thetav2deg, Efortheta, Fenangfull, 10, cmap = cmap1, norm = norm1)
ax.patch.set_facecolor(cmap1(1.0/20))
#cp.set_clim(vmin = 0.0, vmax = 1.0)
plt.colorbar(cp)
#plt.clim(0.0, 1.0)
plt.title('energy-angle('+r'$E$-$\theta$'+') distribution', fontsize = 20)
plt.ylabel(r'$E $', fontsize = 20)
plt.xlabel(r'$\theta (^{\circ})$', fontsize = 20)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.ylim(0.0, E[len(E)-1])
plt.xlim(0.0, 90)
plt.tight_layout()
plt.savefig("fig-enangle_v2.pdf")
