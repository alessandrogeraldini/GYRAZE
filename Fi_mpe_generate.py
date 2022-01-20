import numpy as np
import scipy.special as sp

alpha, rhoeoverlambdaD, n_species, tau, mioverme, fix_current, current = np.loadtxt('inputfile.txt')
Te = 1.0/(tau)
print(alpha, Te)
Ts = 1.0 + Te

maxE = 20.0 + 1.0*Te
maxvpar = 5.0
maxmu = 20.0
maxvperp = 5.0
U1 = 1.0
dvpar = 0.05
dvperp = 0.05
dmu = 0.05
dE = 0.1
#sizemu = int(maxmu/dmu)
sizemu = int(maxvperp/dvperp)

sizeUminmu= int(maxvpar/dvpar)
#Te = 0.0000001
if Te > 0.999: 
	#sizeUminmu= int(maxE/dE)
	coldelectrons = 0
elif (Te > 0.00001) and (Te < 0.999):
	#dEnarrow = dE*Te
	#sizeUminmu= int(maxE/dE)  + int(U1/dEnarrow) - int(U1/dE)
	coldelectrons = 0
else:
	coldelectrons = 1
	sizeUminmu= int(maxE/dE)


u = 0.0
condition = 2.0/Te 
flow = "flow not calculated yet"
if coldelectrons == 0:
	if Te>0.999:
		chodura = 9999999.0
		u = 0.01
		while (chodura > condition) or (chodura < condition - 0.05):
			if (chodura > condition):		
				u += 0.0001
			elif (chodura < condition):
				u -= 0.0001
			normalization =  (1 + sp.erf(u))*(1.0+2.0*u**2) + (2.0/np.sqrt(np.pi))*u*np.exp(-u**2) 
			chodura = 2.0*(1+sp.erf(u))/(normalization)
			flow = ( 0.5*(1+u**2)*np.exp(-u**2) + (1.5 + u**2)*(np.sqrt(np.pi)*u/2.0)*(1 + sp.erf(u)))*4/np.sqrt(np.pi)/normalization 
			print(u, flow, chodura, condition)
	else:
		chodura = 0.0
		while (chodura > condition) or (chodura < condition - 0.1):
			if (chodura > condition):		
				u -= 0.01
			elif (chodura < condition):
				u += 0.01
			normalization =  (np.sqrt(np.pi*u) - np.pi*np.exp(1/u)*(1-sp.erf(1/np.sqrt(u))))/(2.0*u**(1.5))
			print(u)
			chodura = np.pi*np.exp(1/u)*(1-sp.erf(1/np.sqrt(u)))/(2.0*np.sqrt(u)*normalization)
			print(u, flow, chodura, condition)

fp = open('Fi_mpe.txt', 'w')

#u=0.0
#normalization =  (1 + sp.erf(u))*(1.0+2.0*u**2) + (2.0/np.sqrt(np.pi))*u*np.exp(-u**2) 
if coldelectrons == 0:
	if Te>0.999:
		for i in range(sizemu):
			#mu = i*dmu
			mu = i*dvperp*i*dvperp
			for j in range(sizeUminmu):
				Uminmu = j*dvpar*j*dvpar
				#Uminmu = j*dE
				FF = (1.0/np.pi**(1.5))*(1.0/normalization)*(Uminmu)*np.exp(- Uminmu - mu + 2.0*u*np.sqrt(Uminmu) - u**2)

				if j==sizeUminmu-1:
					fp.write(str(FF)+'\n')
				else:
					fp.write(str(FF)+' ')

		fp.close()

		fp2 = open('Fi_mpe_args.txt', 'w')
		for j in range(sizemu):
			#mu = j*dmu
			mu = j*dvperp*j*dvperp
			if j == sizemu-1:
				fp2.write(str(mu)+'\n')
				
			else:
				fp2.write(str(mu)+' ')

		for i in range(sizeUminmu):
			Uminmu = i*dvpar*i*dvpar
			if i == sizeUminmu-1:
				fp2.write(str(Uminmu)+'\n')
				
			else:
				fp2.write(str(Uminmu)+' ')
	else:
		for i in range(sizemu):
			mu = i*dvperp*i*dvperp
			for j in range(sizeUminmu):
				#Uminmu = j*dEnarrow
				Uminmu = j*dvpar*j*dvpar
			#	else: 
			#		Uminmu = (j - int(U1/dEnarrow) + int(U1/dE))*dE

				FF = (1.0/(4.0*np.pi))*(1.0/normalization)*((Uminmu)/(1+u*(Uminmu)))*np.exp(- Uminmu - mu )

				if j==sizeUminmu-1:
					fp.write(str(FF)+'\n')
				else:
					fp.write(str(FF)+' ')

		fp.close()

		fp2 = open('Umufile.txt', 'w')
		for j in range(sizemu):
			mu = j*dvperp*j*dvperp
			if j == sizemu-1:
				fp2.write(str(mu)+'\n')
				
			else:
				fp2.write(str(mu)+' ')

		for i in range(sizeUminmu):
			Uminmu = i*dvpar*i*dvpar
			#if (i < int(U1/dEnarrow)):
			#	Uminmu = i*dEnarrow
			#else: 
			#	Uminmu = (i - int(U1/dEnarrow) + int(U1/dE))*dE

			if i == sizeUminmu-1:
				fp2.write(str(Uminmu)+'\n')
				
			else:
				fp2.write(str(Uminmu)+' ')

if coldelectrons == 1:
	for i in range(sizemu):
		mu = i*dE
		for j in range(sizeUminmu):
			Uminmu = j*dE
			FF = (1.0/(2.0*(np.pi)**(1.5)))*np.exp(- Uminmu - mu )
			if j==sizeUminmu-1:
				fp.write(str(FF)+'\n')
			else:
				fp.write(str(FF)+' ')
	fp.close()

	fp2 = open('Umufile.txt', 'w')
	for j in range(sizemu):
		mu = j*dE
		if j == sizemu-1:
			fp2.write(str(mu)+'\n')
			
		else:
			fp2.write(str(mu)+' ')

	for i in range(sizeUminmu):
		Uminmu = i*dE
		if i == sizeUminmu-1:
			fp2.write(str(Uminmu)+'\n')
			
		else:
			fp2.write(str(Uminmu)+' ')


fp2.close()

