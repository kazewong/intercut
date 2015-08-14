import numpy as np
import lineClass as lc
from astropy.io import fits
from scipy import interpolate
from scipy import optimize
import math
from LuminosityFunction import *
from CuttingFunction import *

fileio = np.genfromtxt('fileiov6.dat',dtype = "S200",delimiter = ';')
field = []
lineDict={}
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
		noise = np.genfromtxt(fileio[linenumber][1]).transpose()# remember to change the noise data format
	if fileio[linenumber][0] == 'lines':
		lineDict[eval(fileio[linenumber][1])[0]] = lc.line(eval(fileio[linenumber][1])[0],eval(fileio[linenumber][1])[1])
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
		
for line in lineDict:
	lineDict[line].flux = np.array(table.field('Flux_'+lineDict[line].name))
	lineDict[line].Lambda = np.array(table.field('Lambda_'+lineDict[line].name))
	for io in fileio:
		if io[0] == 'noisepara' and flags[0] == 1:
			localpara = eval(io[1])
			dataPoint = np.array([(lineDict[line].wavelength*(1+z)),np.array(radius*0.03)])
			grid = np.meshgrid(np.ogrid[localpara[0]:localpara[1]:localpara[2]],np.ogrid[localpara[3]:localpara[4]:localpara[5]])[0]
			lineDict[line].noise = interpolate.griddata(grid,noise,datapoint,fill_value = 1.)
		if io[0] == 'noisepara' and flags[0] == 0:
			lineDict[line].noise = np.interp(lineDict[line].wavelength*(1+z),noise[0],noise[1])
		if io[0] == 'condition':
			conditions = np.genfromtxt(io[1],dtype = 'S20,S500',delimiter = ';')
			for condition in conditions:
				if condition[0]==lineDict[line].name:
					lineDict[line].condition = np.append(lineDict[line].condition,condition[1])
			lineDict[line].condition = np.array(lineDict[line].condition)
			
O2noise = np.genfromtxt('noiseO2a2.dat').T
lineDict['OII'].noise = np.interp(lineDict['OII'].wavelength*(1+z),O2noise[0],O2noise[1])
	
#*****************************END OF IO SECTION**************************************
print 'the output file is going to be '+str(outputname)
print 'the simulation number is '+str(parameter[0])
print 'the lines input are ' + str(len(lineDict))
print 'the entries number is '+str(z.size)
#****************************************CALIBRATING SECTION********************************
#this part calibrate the catalogue 'CMC081211',only Ha,O3a and O3b need calibration 

kO3b=4.05+2.659*(-2.156+1.509/.5007-0.198/(.5007**2)+0.011/(.5007**3))
kO3a=4.05+2.659*(-2.156+1.509/.4959-0.198/(.4959**2)+0.011/(.4959**3))
it = np.array(np.where(lineDict['Ha'].flux>0)[0])
dl =np.array(10.**(0.2*dmz[it]+1)*3.09e18,dtype = 'float64')
lumHa = lineDict['Ha'].flux[it]*4*pi*dl**2
lumO3a = lineDict['OIIIa'].flux[it]*4*pi*dl**2
lumO3b = lineDict['OIIIb'].flux[it]*4*pi*dl**2
lumHa2 = np.array(CalibrationHa((1+z)*lineDict['OII']/lineDict['Ha']-1,lumHa,it))
lumO3b2 = np.array(CalibrationO3b((1+z)*lineDict['OII']/lineDict['OIIIb']-1,lumO3b,it))
lineDict['Ha'].flux =  np.array([-9999. for x in range(z.size)])
lineDict['OIIIb'].flux =  np.array([-9999. for x in range(z.size)])
lineDict['Ha'].flux[it] = lumHa2.astype('float64')/(4.*pi*dl**2.) 
lineDict['OIIIb'].flux[it] = lumO3b2.astype('float')/(4*pi*dl**2) 
lineDict['OIIIa'].flux[it] = (lineDict['OIIIb'].flux*(10.**(0.4*ebv*kO3b))/3.)*(10.**(-0.4*ebv*kO3a))

#****************************END OF CALIBRATION******************************

#****************************RANDOM FLUCTUATION SECTION**********************
#this section perturb the flux of the spectral lines and photometric cut
for i in lineDict:
	lineDict[i].flux = np.array(lineDict[i].flux)
	lineDict[i].origin = lineDict[i].flux.copy()
if flags[1] == 1:
	ferr = np.random.normal(size = [parameter[0],len(lineDict),lineDict[0].flux.size])
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
		if photocut != '':lineDict[lineObject].selection(eval(photocut))
	if flags[1]==1:	
		lineFluxPertubation(lineDict,1e-5)
#********************************END OF RANDOM FLUCTUATION SECTION*******************



#*******************************SECOND LINE IDENTIFICATION SECTION*******************
#this part evaluates the second line cutting condition and return the result as histogram

	secondLineCutting(lineDict,field)
	getFraction(lineDict,'OII',z,(1,2))

for lineObject in lineDict:
	lineDict[lineObject].fraction = np.array(lineDict[lineObject].fraction)
	lineDict[lineObject].fraction = lineDict[lineObject].fraction.reshape(lineDict[lineObject].fraction.size/parameter[1],parameter[1])
	lineDict[lineObject].getStat()

#*************************************END OF SECOND LINE IDENTIFICATION SECTION******************************

#***********************************************PLOTTING SECTION*********************************************
output = open(outputname,'w')
for i in lineDict:
	for j in range(parameter[1]):
		output.write(str(lineDict[i].mean[j])+'\t')
	output.write('\n')
	
	for j in range(parameter[1]):
		output.write(str(lineDict[i].sd[j])+'\t')
	output.write('\n')
	

output.close()


