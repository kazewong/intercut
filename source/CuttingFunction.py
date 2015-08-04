import numpy as np

def secondLineCutting(lines,field):
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
	flux = np.array([0.0 for i in range(z.size)])
	noise = np.array([0.0 for i in range(z.size)])
	zd = np.array([0.0 for i in range(z.size)])
	wavelength = float
	for lineObject in lines:
		print 'cutting '+ lines[lineObject].name+' line, '+str(lines[lineObject].condition.size)+' condition found.'
		for lineObjectCondition in lines[lineObject].condition:
			if eval(lineObjectCondition)[0] == 'normal':
				flux = lines[eval(lineObjectCondition)[1]].flux
				noise = lines[eval(lineObjectCondition)[1]].noise
				wavelength = lines[eval(lineObjectCondition)[1]].wavelength
				lines[lineObject].selection(eval(lineObjectCondition)[2])
			if eval(lineObjectCondition)[0] == 'conditional':
				zd = (lines[lineObject].wavelength/lines[eval(lineObjectCondition)[1]].wavelength)*(1+z)-1
				lines[lineObject].conditionalSelection(eval(lineObjectCondition)[2],eval(lineObjectCondition)[3])

def lineFluxPertubation(lines,threshold):
	for lineObject in lines:
		ferr = np.random.normal(size = lines[lineObject][lines[lineObject].noise<threshold].size)
		lines[lineObject][lines[lineObject].noise<threshold].flux = lines[lineObject][lines[lineObject].noise<threshold].origin + lines[lineObject][lines[lineObject].noise<threshold].noise*ferr
		
def getFraction(lines,z,ranges,bin = 10):
	mainLine = lines[lines.items()[0][0]]
	hmain = np.array(np.histogram(z[mainLine.intf],bins = bin,range = ranges)[0],dtype = float)
	for lineObject in lines:
		if lines[lineObject] !=[]:
			hlocal = np.array(np.histogram((lines[lineObject].wavelength/lines[lines.items()[0][0]].wavelength)*(1+(z[lines[lineObject].intf]))-1)[0])
		else:
			hlocal = np.array(np.histogram(0,bins =bin,range = ranges)[0])
		lines[lineObject].fraction = np.append(lines[lineObject].fraction,np.divide(hlocal,hmain))
		
