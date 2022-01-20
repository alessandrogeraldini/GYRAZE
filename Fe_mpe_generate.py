import numpy as np
import scipy.special as sp

v_max = 7.0
#sizevpar = 35
#dmu = 0.2
dvpar = 0.2
dvperp = 0.2
sizemu = int(v_max/dvperp)
sizevpar = int(v_max/dvpar)
#sizemu = int(18.0/dmu)
print(sizemu)

fp = open('Fe_mpe.txt', 'w')
for i in range(sizemu):
	mu = 0.5*i*dvperp*i*dvperp
	#mu = i*dmu
	for j in range(sizevpar):
		vpar = j*dvpar
		Uminmu = 0.5*vpar**2
		FF = (1.0/(2.0*(np.pi)**(1.5)))*np.exp(- Uminmu - mu )
		if j==sizevpar-1:
			fp.write(str(FF)+'\n')
		else:
			fp.write(str(FF)+' ')
fp.close()

fp2 = open('Fe_mpe_args.txt', 'w')
for j in range(sizemu):
	mu = 0.5*j*dvperp*j*dvperp
	#mu = j*dmu
	if j == sizemu-1:
		fp2.write(str(mu)+'\n')
		
	else:
		fp2.write(str(mu)+' ')

for i in range(sizevpar):
	vpar = i*dvpar
	if i == sizevpar-1:
		fp2.write(str(vpar)+'\n')
		
	else:
		fp2.write(str(vpar)+' ')


fp2.close()

