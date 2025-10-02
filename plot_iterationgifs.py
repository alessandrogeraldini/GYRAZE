import numpy as np
import matplotlib.pyplot as plt
import imageio
import os.path

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

phixMP_pnglist=[]
nxMP_pnglist=[]

N=0
phixDS_pnglist=[]
ExDS_pnglist=[]
phippxMP_pnglist =[]
nxDS_pnglist=[]
errorDS_pnglist=[]
errorMP_pnglist=[]

while os.path.isdir('iteration'+str(N)):
	#MP
	xx, phi, ni, ne = np.loadtxt('iteration'+str(N)+'/phi_n_MP.txt', unpack = True)
	rho = [nin-nen for (nen,nin) in zip(ne, ni)]
	normrho = [1.0-nen/nin for (nen,nin) in zip(ne, ni)]
	i=len(ni)-1
	while (ni[i] == 0.0):
		i-=1
	sizenMP = i+1
	phipp = [ 2.0*(-(phi[i+2] - phi[i+1])/(xx[i+2] - xx[i+1]) + (phi[i+1] - phi[i])/(xx[i+1]-xx[i]))/(xx[i+2]-xx[i]) for i in range(len(phi)-2) ]
	plt.plot(xx, phi, lw = 2.0)
	plt.ylim(-3.0, 0.0)
	plt.xlim(0.0, 10.0)
	plt.title('potential at iteration # = '+str(N))
	plt.savefig('iteration'+str(N)+'/phix_MP.png')
	plt.close()
	phixMP_pnglist.append('iteration'+str(N)+'/phix_MP.png')
	xpp = [(xx[i]+xx[i+2])/4.0 + xx[i+1]/2.0 for i in range(len(xx)-2)]
	plt.plot(xpp, phipp, lw = 2.0)
	plt.ylim(0.0, 10.0)
	plt.xlim(0.0, 10.0)
	plt.title('electric field gradient at iteration # = '+str(N))
	plt.savefig('iteration'+str(N)+'/phippx_MP.png')
	phippxMP_pnglist.append('iteration'+str(N)+'/phippx_MP.png')
	plt.close()
	plt.plot(xx[0:sizenMP], ni[0:sizenMP], '-', lw = 2.0)
	plt.plot(xx[0:sizenMP], ne[0:sizenMP], '-', lw = 2.0)
	plt.ylim(0.0, 1.0)
	plt.xlim(0.0, 10.0)
	plt.title('density at iteration # = '+str(N))
	plt.savefig('iteration'+str(N)+'/nx_MP.png')
	plt.close()
	nxMP_pnglist.append('iteration'+str(N)+'/nx_MP.png')
	plt.plot(xx[0:sizenMP], normrho[0:sizenMP], '-', lw = 2.0)
	plt.xlim(0.0, 10.0)
	plt.title('error at iteration # = '+str(N))
	plt.savefig('iteration'+str(N)+'/error_MP.png')
	plt.close()
	errorMP_pnglist.append('iteration'+str(N)+'/error_MP.png')
	#DS
	if os.path.isfile('iteration'+str(N)+'/phi_n_DS.txt') == True:
		xx, phi, ni, ne = np.loadtxt('iteration'+str(N)+'/phi_n_DS.txt', unpack = True)
		xxs_D = [gamma*xval for xval in xx]
		i=len(ne)-1
		while (ne[i] == 0.0):
			i-=1
		sizenDS = i+1
		rho = [nin-nen for (nen,nin) in zip(ne, ni)]
		phipp = [ -(1/gamma**2)*(phi[i+2] - 2*phi[i+1] + phi[i])/(xx[1]**2) for i in range(len(phi)-2) ]
		#phipp.insert(0, phipp[0])
		phipp.insert(0, phipp[0] - (phipp[1] - phipp[0]))
		phipp.append(0.0)
		plt.plot(xxs_D, phi, lw = 2.0)
		plt.ylim(-0.7, 0.0)
		plt.xlim(0.0, 10.0)
		plt.title('potential at iteration # = '+str(N))
		plt.savefig('iteration'+str(N)+'/phix_DS.png')
		plt.close()
		phixDS_pnglist.append('iteration'+str(N)+'/phix_DS.png')
		phip = [(phi[i+1]-phi[i])/(xx[i+1]- xx[i]) for i in range(len(phi)-1)]
		#for p in phip:
		#	print(p)
		#	if (p< 0.0):
		#		print("Ex is negative")
		#		
		phip.append(phip[len(phip)-1])
		plt.plot(xx, phip, lw = 2.0)
		plt.title('E-field at iteration # = '+str(N))
		plt.savefig('iteration'+str(N)+'/Ex_DS.png')
		plt.close()
		ExDS_pnglist.append('iteration'+str(N)+'/Ex_DS.png')
		plt.plot(xxs_D[1:sizenDS], ni[1:sizenDS], '-', lw = 2.0)
		plt.plot(xxs_D[1:sizenDS], ne[1:sizenDS], '-', lw = 2.0)
		plt.plot(xxs_D[1:sizenDS], rho[1:sizenDS], '-', lw = 2.0)
		plt.plot(xxs_D[1:sizenDS], phipp[1:sizenDS], '-', lw = 2.0)
		plt.ylim(-0.2, 1.0)
		plt.xlim(0.0, 10.0)
		plt.title('density at iteration # = '+str(N))
		plt.savefig('iteration'+str(N)+'/nx_DS.png')
		plt.close()
		nxDS_pnglist.append('iteration'+str(N)+'/nx_DS.png')
		normerror = [1.0 - ne[i]/ni[i] - phipp[i]/ni[i] for i in range(len(phipp))]
		plt.plot(xx[1:sizenDS], normerror[1:sizenDS], '-', lw = 2.0)
		plt.xlim(0.0, 10.0)
		plt.title('error at iteration # = '+str(N))
		plt.savefig('iteration'+str(N)+'/error_DS.png')
		plt.close()
		errorDS_pnglist.append('iteration'+str(N)+'/error_DS.png')
	N+=1

with imageio.get_writer('phix_MP_iterations.gif', mode='I') as writer:
    for filename in phixMP_pnglist:
        image = imageio.imread(filename)
        writer.append_data(image)

with imageio.get_writer('phippx_MP_iterations.gif', mode='I') as writer:
    for filename in phippxMP_pnglist:
        image = imageio.imread(filename)
        writer.append_data(image)

with imageio.get_writer('nx_MP_iterations.gif', mode='I') as writer:
    for filename in nxMP_pnglist:
        image = imageio.imread(filename)
        writer.append_data(image)

with imageio.get_writer('phix_DS_iterations.gif', mode='I') as writer:
    for filename in phixDS_pnglist:
        image = imageio.imread(filename)
        writer.append_data(image)

with imageio.get_writer('Ex_DS_iterations.gif', mode='I') as writer:
    for filename in ExDS_pnglist:
        image = imageio.imread(filename)
        writer.append_data(image)

with imageio.get_writer('nx_DS_iterations.gif', mode='I') as writer:
    for filename in nxDS_pnglist:
        image = imageio.imread(filename)
        writer.append_data(image)

with imageio.get_writer('error_MP_iterations.gif', mode='I') as writer:
    for filename in errorMP_pnglist:
        image = imageio.imread(filename)
        writer.append_data(image)

with imageio.get_writer('error_DS_iterations.gif', mode='I') as writer:
    for filename in errorDS_pnglist:
        image = imageio.imread(filename)
        writer.append_data(image)

