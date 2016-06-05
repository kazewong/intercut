import numpy as np
import lineClass as lc
from astropy.io import fits
from scipy import interpolate
from scipy import optimize
import math
import LuminosityFunction as LF
import mpmath
from CuttingFunction import *

print 'Program start, reading input files..........'
fileio = np.genfromtxt('fileioWFIRST.dat',dtype = "S200",delimiter = ';')
field = []
lineDict={}
parameter = [0,0,0,0]
zd = np.array([])
flux = np.array([])
noise = np.array([])
photocut = ''
mainline = ''
cuttype= ''

flags = [1,1]# the first flag is radial noise,second is random fluctuation
#***********************************IO SECTION********************************
#this part read the variables from the program parameters file into the program and convert those variable into objects and local variable

for linenumber in range(fileio.size/fileio.ndim):
	if fileio[linenumber][0] == 'catalog':
		cata = fits.open(fileio[linenumber][1])
		table = cata[1].data
#***********************predefined vairables******************************************
#field correspoinding to the field in the fits file , please follow the convention below if addition variable wanted to be predefined
		locallist = [[3,4,5,6,7,8,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,108],[0,0,0,0,0,0,0,0,27.5,27.7,27.0,26.2,26.9,26.8,0,0,0,0,0,0,0,0,0,26.8],[0,0,0,0,0,0,0,0,0.39,0.39,0.39,0.39,0.12,0.14,0,0,0,0,0,0,0,0,0,0.12]]
		for i in range(len(locallist[0])):
			field.append(lc.photocut(cata[1].header[locallist[0][i]],locallist[1][i],locallist[2][i]))
			field[i].mag = np.array(table.field(locallist[0][i]))
			field[i].origin = np.array([-99.9 for loop in range(cata[1].header['naxis2'])])
			field[i].origin[np.where(field[i].mag > 0)] = np.power(10,(-0.4*((field[i].mag[np.where(field[i].mag >0)])+48.6)))
			field[i].flux = field[i].origin.copy()
		field = np.array(field)
		z = field[0].mag.copy()
		ebv = field[1].mag.copy()		
		dmz = field[2].mag.copy()
		mnyv = field[3].mag.copy()	
		Mstar = field[4].mag.copy()
		radius = field[5].mag.copy()*0.03
		mb = field[6].mag
		mv = field[7].mag
		mg = field[8].mag		
		mr = field[9].mag
		mi = field[10].mag
		mz = field[11].mag
		mj = field[12].mag
		mh = field[13].mag
		my = field[23].mag
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
	if fileio[linenumber][0] == 'mainline':
		mainline = lineDict[fileio[linenumber][1]].name 
	if fileio[linenumber][0] == 'cuttype':
		cuttype = fileio[linenumber][1]

if flags[0] ==1:
	for object in field:
		if object.dep != 0:
			object.dep = np.array(object.dep-1.25*np.log10(1+np.power((radius/object.r0),2)))
		else:
			object.dep = np.array([object.dep for i in range(object.mag.size)])
		object.err = np.power(10,(-0.4*(object.dep-23.9)))*1.e-29*0.2	
else:
	for object in field:
		object.dep = np.array([object.dep for i in range(object.mag.size)])
		object.err = np.power(10,(-0.4*(object.dep-23.9)))*1.e-29*0.2	
for line in lineDict:
	lineDict[line].flux = np.array(table.field('Flux_'+lineDict[line].name))
	lineDict[line].Lambda = np.array(table.field('Lambda_'+lineDict[line].name))
	for io in fileio:
		if io[0] == 'noisepara' and flags[0] == 1:
			localpara = eval(io[1])
			dataPoint = np.array([np.array(radius),(lineDict[line].wavelength*(1+z))]).T
			grid = [np.ogrid[localpara[3]:localpara[4]:localpara[5]],np.ogrid[localpara[0]:localpara[1]:localpara[2]]]
			lineDict[line].noise = interpolate.interpn(grid,noise,dataPoint,bounds_error=False,fill_value = 1.)	
		if io[0] == 'noisepara' and flags[0] == 0:
			lineDict[line].noise = np.interp(lineDict[line].wavelength*(1+z),noise[0],noise[1])
		if io[0] == 'condition':
			conditions = np.genfromtxt(io[1],dtype = 'S20,S5000',delimiter = ';')
			for condition in conditions:
				if condition[0]==lineDict[line].name:
					lineDict[line].condition = np.append(lineDict[line].condition,condition[1])
			lineDict[line].condition = np.array(lineDict[line].condition)
			
#O2noise = np.genfromtxt('noiseO2a2.dat').T
#lineDict['OII'].noise = np.interp(lineDict['OII'].wavelength*(1+z),O2noise[0],O2noise[1])
	
#*****************************END OF IO SECTION**************************************
print 'The output file is going to be '+str(outputname)
print 'The simulation number is '+str(parameter[0])
print str(len(lineDict)) + ' interloping emission line are considered'
print str(z.size) + ' galaxies read'
#****************************************CALIBRATING SECTION********************************
#this part calibrate the catalogue 'CMC081211',only Ha,O3a and O3b need calibration 

#lineDict["Ha"].flux = np.genfromtxt('CalibratedHa.LF')
#lineDict["OIIIb"].flux = np.genfromtxt('CalibratedO3b.LF')
kO3b=4.05+2.659*(-2.156+1.509/.5007-0.198/(.5007**2)+0.011/(.5007**3))
kO3a=4.05+2.659*(-2.156+1.509/.4959-0.198/(.4959**2)+0.011/(.4959**3))

lineDict["OIIIa"].flux=((lineDict["OIIIb"].flux*(10.**(0.4*ebv*kO3b))/3.)*(10.**(-0.4*ebv*kO3a)))

#****************************END OF CALIBRATION******************************

#****************************RANDOM FLUCTUATION SECTION**********************
#this section perturb the flux of the spectral lines and photometric cut
for i in lineDict:
	lineDict[i].flux = np.array(lineDict[i].flux)
	lineDict[i].origin = lineDict[i].flux.copy()
	
fcut = np.power(10.,-0.4*(56.+1.+48.6))

print 'Simulation Start..........'

for i in range(parameter[0]):
	for lineObject in lineDict:
		lineDict[lineObject].flux = lineDict[lineObject].origin.copy()
		lineDict[lineObject].intf = np.arange(lineDict[lineObject].flux.size)
		if photocut != '':lineDict[lineObject].selection(eval(photocut))
	if flags[1]==1:	
		ferr = np.random.normal(size = [len(lineDict),lineDict[mainline].flux.size])
		photoerr = np.random.normal(size = [len(field),z.size])
		for i in range(6,len(field)):
			field[i].flux[field[i].origin > fcut] = field[i].origin[field[i].origin > fcut]+field[i].err[field[i].origin > fcut]*photoerr[i][field[i].origin > fcut]
			field[i].flux[field[i].flux < fcut] = fcut
			field[i].mag = -48.6-2.5*np.log10(field[i].flux)
		lineFluxPertubation(lineDict,1e-5)
#********************************END OF RANDOM FLUCTUATION SECTION*******************



#*******************************SECOND LINE IDENTIFICATION SECTION*******************
#this part evaluates the second line cutting condition and return the result as histogram

	secondLineCutting(lineDict,field)
	getFraction(lineDict,mainline,z,(parameter[2],parameter[3]),bin=parameter[1])

print 'Writing result to output file..........'
for lineObject in lineDict:
	lineDict[lineObject].fraction = np.array(lineDict[lineObject].fraction)
	lineDict[lineObject].fraction = lineDict[lineObject].fraction.reshape(lineDict[lineObject].fraction.size/parameter[1],parameter[1]).T
	
if cuttype == 'total':getTotalIntf(lineDict)
else:getSelfIntf(lineDict,mainline)
#*************************************END OF SECOND LINE IDENTIFICATION SECTION******************************

#***********************************************PLOTTING SECTION*********************************************
output = open(outputname,'w')
output.write('The redshift of bins are '+str(np.linspace(parameter[2],parameter[3],parameter[1])+'\n'))
for i in lineDict:
	output.write(lineDict[i].name+'\n')
	for j in range(parameter[1]):
		output.write(str(lineDict[i].mean[j])+'\t')
	output.write('\n')
	
	for j in range(parameter[1]):
		output.write(str(lineDict[i].sd[j])+'\t')
	output.write('\n')
	

output.close()


