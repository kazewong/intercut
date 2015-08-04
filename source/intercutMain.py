import numpy as np
import lineClass as lc
from astropy.io import fits
from scipy import interpolate
from scipy import optimize
import math
import LuminosityFunction
from CuttingFunction import *

fileio = np.genfromtxt('fileiov6.dat',dtype = "S200",delimiter = ';')
noiseA = []
lineA = []
field = []
parameter = [0,0,0,0]
zd = np.array([])
flux = np.array([])
noise = np.array([])
photocut = ''
flags = [1,1]# the first flag is radial noise,second is random fluctuation
#***********************************IO SECTION********************************
#this part read the variables from the program parameters file into the program and convert those variable into objects and local variable

for linenumber in range(fileio.size/fileio.ndim):
	if fileio[linenumber][0] == 'catalog':
		cata = fits.open(fileio[linenumber][1])
		table = cata[1].data
#***********************predefined vairables******************************************
#field correspoinding to the field in the fits file , please follow the convention below if addition variable wanted to be predefined
		locallist = [[3,4,5,6,7,8,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61],[0,0,0,0,0,0,0,0,26.5,26.1,25.9,25.1,0,0,0,0,0,0,0,0,0,0,0]]
		for i in range(len(locallist[0])):
			field.append(lc.photocut(cata[1].header[locallist[0][i]],locallist[1][i]))
			field[i].mag = np.array(table.field(locallist[0][i]))
			field[i].origin = np.array([-99.9 for loop in range(cata[1].header['naxis2'])])
			field[i].origin[np.where(field[i].mag > 0)] = (-0.4*((field[i].mag[np.where(field[i].mag >0)])+48.6))
			field[i].origin[np.where(field[i].mag >0)]= np.power(10,field[i].origin[np.where(field[i].mag >0)])
			field[i].flux = field[i].origin.copy()
		field = np.array(field)
		z = field[0].mag.copy()
		ebv = field[1].mag.copy()		
		dmz = field[2].mag.copy()
		mnyv = field[3].mag.copy()	
		Mstar = field[4].mag.copy()
		radius = field[5].mag.copy()
		mb = field[6].mag
		mv = field[7].mag
		mg = field[8].mag		
		mr = field[9].mag
		mi = field[10].mag
		mz = field[11].mag
		pi = 3.1415926
	if fileio[linenumber][0] == 'noise':
		noiseA.append(np.genfromtxt(fileio[linenumber][1]).transpose())# remember to change the noise data format
	if fileio[linenumber][0] == 'lines':
		datfile = np.genfromtxt(fileio[linenumber][1],dtype = "S20,float",delimiter = ';')
		lineA = ([0 for i in range(datfile.size/datfile.ndim)])
		lineDict = {}
		for i in range(datfile.size/datfile.ndim):
			lineA[i] = lc.line(datfile[i][0],datfile[i][1])
			lineA[i].flux = np.array(table.field('Flux_'+lineA[i].name))
			lineA[i].Lambda = np.array(table.field('Lambda_'+lineA[i].name))
			if flags[0] == 1:
				row = []
				column = []
				for linenumber2 in range(fileio.size/fileio.ndim):	
					if fileio[linenumber2][0] == 'noisepara':
						localpara = eval(fileio[linenumber2][1])
						row = np.ogrid[localpara[0]:localpara[1]:localpara[2]]
						column = np.ogrid[localpara[3]:localpara[4]:localpara[5]]	
				localinterp = interpolate.interp2d(row,column,noiseA[0],fill_value = 1.)
				localx = (((lineA[i].wavelength)*(1+z)))
				localy = np.array(radius*0.03)	
				lineA[i].noise = np.array([0. for j in range(lineA[i].flux.size)])
				for j in range(lineA[i].noise.size):
					lineA[i].noise[j] = localinterp(localx[j],localy[j])[0]
			else:
					lineA[i].noise = np.interp(lineA[i].wavelength*(1+z),noiseA[0][0],noiseA[0][1])#the firsxt column is the observed wavelength and the second column is the noise
	
#				if i == 0:lineA[i].noise = np.interp(z,noiseA[1][0],noiseA[1][1])
#				else :lineA[i].noise = np.interp(lineA[i].wavelength*(1+z),(noiseA[0][0]+1)*lineA[0].wavelength,noiseA[0][1])
			for linenumber2 in range(fileio.size/fileio.ndim):	
				if fileio[linenumber2][0] == 'condition':
					condition = np.genfromtxt(fileio[linenumber2][1],dtype = "S20,S500",delimiter = ';')
					for i2 in range(condition.size/condition.ndim):
						if condition[i2][0] == lineA[i].name:
							lineA[i].condition = np.append(lineA[i].condition,condition[i2][1])
					lineA[i].condition = np.array(lineA[i].condition)
			lineDict[lineA[i].name]=lineA[i]
	if fileio[linenumber][0] == 'loopNumber':
		parameter[0] = eval(fileio[linenumber][1])
	if fileio[linenumber][0] == 'binNumber':
		parameter[1] = eval(fileio[linenumber][1])
	if fileio[linenumber][0] == 'lowerLim':
		parameter[2] = eval(fileio[linenumber][1])
	if fileio[linenumber][0] == 'upperLim':
		parameter[3] = eval(fileio[linenumber][1])
	if fileio[linenumber][0] == 'output':
		outputname = fileio[linenumber][1]
	if fileio[linenumber][0] == 'photocut':
		photocut = fileio[linenumber][1]
	if fileio[linenumber][0] == 'flags':
		flags = eval(fileio[linenumber][1])

#*****************************END OF IO SECTION**************************************
print 'the output file is going to be '+str(outputname)
print 'the simulation number is '+str(parameter[0])
print 'the lines input are ' + str(len(lineA))
print 'the entries number is '+str(z.size)

#****************************************CALIBRATING SECTION********************************
#this part calibrate the catalogue 'CMC081211',only Ha,O3a and O3b need calibration 

# name = [0,0,0]
# for namematch in range(len(lineA)):
# 	if lineA[namematch].name == 'Ha':name[0] = namematch
# 	if lineA[namematch].name == 'OIIIa':name[1] = namematch
# 	if lineA[namematch].name == 'OIIIb':name[2] = namematch
# 
# 
# it = np.array(np.where(lineA[name[0]].flux > 0.0))
# dl = [0 for x in range((it[0].size/it[0].ndim))]
# lumHa =[0 for x in range((it[0].size/it[0].ndim))]
# lumO3a =[0 for x in range((it[0].size/it[0].ndim))]
# lumO3b =[0 for x in range((it[0].size/it[0].ndim))]
# 
# 
# 
# for i in range((it[0].size/it[0].ndim)):
# 	dl[i] = (10.**((0.2*dmz[it[0][i]])+1))*(3.09e18)
# 	lumHa[i] = (lineA[name[0]].flux[it[0][i]])*4*pi*dl[i]**2
# 	lumO3a[i] = (lineA[name[1]].flux[it[0][i]])*4*pi*dl[i]**2
# 	lumO3b[i] = (lineA[name[2]].flux[it[0][i]])*4*pi*dl[i]**2
# 
# lumHa2 = [0. for x in range((it[0].size/it[0].ndim))]
# lumHa2 = LuminosityFunction.CalibrationHa(z,lumHa,it)
# lumHa2 = np.power(10,lumHa2)#*1.4
# lineA[name[0]].flux = [-9999 for x in range((z.size/z.ndim))]
# 
# for i in range((it[0].size/it[0].ndim)):
#         lineA[name[0]].flux[it[0][i]] = (float(lumHa2[i])/(4*pi*(dl[i]**2)))
# 
# lumO3b2 = [0 for x in range((it[0].size/it[0].ndim))]
# lumO3b2 = LuminosityFunction.CalibrationO3b(z,lumO3b,it)
# lumO3b2 = np.power(10,lumO3b2)
# lineA[name[2]].flux = [-9999 for x in range((z.size/z.ndim))]
# for i in range((it[0].size/it[0].ndim-1)):
#         lineA[name[2]].flux[it[0][i]] = lumO3b2[i]/(4*pi*(dl[i]**2))
# 
# kO3b=4.05+2.659*(-2.156+1.509/.5007-0.198/(.5007**2)+0.011/(.5007**3))
# kO3a=4.05+2.659*(-2.156+1.509/.4959-0.198/(.4959**2)+0.011/(.4959**3))
# 
# for i in range (len(lineA[2].flux)):
#         lineA[name[1]].flux[i]=((lineA[name[2]].flux[i]*(10.**(0.4*ebv[i]*kO3b))/3.)*(10.**(-0.4*ebv[i]*kO3a)))
# 

#****************************END OF CALIBRATION******************************

#****************************RANDOM FLUCTUATION SECTION**********************
#this section perturb the flux of the spectral lines and photometric cut
for i in range(len(lineA)):
	lineA[i].flux = np.array(lineA[i].flux)
	lineA[i].origin = lineA[i].flux.copy()
if flags[1] == 1:
	ferr = np.random.normal(size = [parameter[0],len(lineA),lineA[0].flux.size])
	photoerr = np.random.normal(size = [len(field),z.size])
	for i in range(6,len(field)):
		fcut = np.power(10,-0.4*(56+1.+48.6))
		for j in range(field[i].mag.size):
			if (field[i].origin[j] > fcut):
				field[i].flux[j] = field[i].origin[j]+field[i].err*photoerr[i][j]
				if (field[i].flux[j] < fcut):field[i].flux[j] = fcut
			else:
				field[i].flux[j] = fcut
			field[i].mag[j] = -48.6-2.5*np.log10(field[i].flux[j])



for i in range(parameter[0]):
	print 'the loop number is '+str(i)
	for lineObject in lineDict:
		lineDict[lineObject].flux = lineDict[lineObject].origin.copy()
		lineDict[lineObject].intf = np.arange(lineDict[lineObject].flux.size)
	if flags[1]==1:	
		lineFluxPertubation(lineDict,1e-5)
#********************************END OF RANDOM FLUCTUATION SECTION*******************



#*******************************SECOND LINE IDENTIFICATION SECTION*******************
#this part evaluates the second line cutting condition and return the result as histogram

	secondLineCutting(lineDict,field)
	getFraction(lineDict,z,(0.8,2.4))

for lineObject in lineDict:
	lineDict[lineObject].fraction = np.array(lineDict[lineObject].fraction)
	lineDict[lineObject].fraction.reshape(lineDict[lineObject].fraction.size/parameter[1],parameter[1])
	lineDict[lineObject].getStat()

#*************************************END OF SECOND LINE IDENTIFICATION SECTION******************************

#***********************************************PLOTTING SECTION*********************************************
output = open(outputname,'w')
for i in range(len(lineA)):
	for j in range(parameter[1]):
		output.write(str(lineA[i].mean[j])+'\t')
	output.write('\n')
	
	for j in range(parameter[1]):
		output.write(str(lineA[i].sd[j])+'\t')
	output.write('\n')
	

output.close()


