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
	mj = field[12].mag
	mh = field[13].mag
	my = field[23].mag
	flux = np.array([0.0 for i in range(z.size)])
	noise = np.array([0.0 for i in range(z.size)])
	zd = np.array([0.0 for i in range(z.size)])
	wavelength = float
	for lineObject in lines:
		print 'cutting '+ lines[lineObject].name+' line, '+str(lines[lineObject].condition.size)+' condition found.'
		for lineObjectCondition in lines[lineObject].condition:
			flux = lines[eval(lineObjectCondition)[1]].flux
			noise = lines[eval(lineObjectCondition)[1]].noise
			wavelength = lines[eval(lineObjectCondition)[1]].wavelength
			zd = (lines[lineObject].wavelength/lines[eval(lineObjectCondition)[1]].wavelength)*(1+z)-1
			if eval(lineObjectCondition)[0] == 'normal':
				lines[lineObject].selection(eval(lineObjectCondition)[2])
			if eval(lineObjectCondition)[0] == 'conditional':
				print zd[0:10]
				lines[lineObject].conditionalSelection(eval(lineObjectCondition)[2],eval(lineObjectCondition)[3])

def lineFluxPertubation(lines,threshold):
	for lineObject in lines:
		ferr = np.random.normal(size = lines[lineObject].noise[lines[lineObject].noise<threshold].size)
		lines[lineObject].flux[lines[lineObject].noise<threshold] = lines[lineObject].origin[lines[lineObject].noise<threshold] + lines[lineObject].noise[lines[lineObject].noise<threshold]*ferr
		
def getFraction(lines,mainLineName,z,ranges,bin = 10):
	mainLine = lines[mainLineName]
	print mainLine.name, mainLine.intf
	hmain = np.array(np.histogram(z[mainLine.intf],bins = bin,range = ranges)[0],dtype = float)
	print mainLine.name,hmain
	for lineObject in lines:
		if lines[lineObject] !=[]:
			hlocal = np.array(np.histogram((lines[lineObject].wavelength/mainLine.wavelength)*(1+(z[lines[lineObject].intf]))-1,bins = bin,range = ranges)[0])
			print lines[lineObject].name,hlocal
		else:
			hlocal = np.array(np.histogram(0,bins =bin,range = ranges)[0])
		lines[lineObject].fraction = np.append(lines[lineObject].fraction,np.divide(hlocal,hmain))
		
def getTotalIntf(lines):
	fractionAll = np.array([[0. for i in range(lines.items()[0][1].fraction.shape[1])] for i in range(lines.items()[0][1].fraction.shape[0])])
	for line in lines:
		fractionAll = fractionAll + lines[line].fraction
		lines[line].mean = np.array([0. for i in range(lines.items()[0][1].fraction.shape[0])])
		lines[line].sd = np.array([0. for i in range(lines.items()[0][1].fraction.shape[0])])
	for line in lines:
		local = np.divide(lines[line].fraction,fractionAll)
		for bin in range(lines[line].fraction.shape[0]):
			lines[line].mean[bin] = np.mean(local[bin])
			lines[line].sd[bin] = np.std(local[bin])
	
def getSelfIntf(lines,mainline):
	for line in lines:
		lines[line].mean = np.array([0. for i in range(lines.items()[0][1].fraction.shape[0])])
		lines[line].sd = np.array([0. for i in range(lines.items()[0][1].fraction.shape[0])])
		local = np.divide(lines[line].fraction,lines[line].fraction+lines[mainline].fraction)
		for bin in range(lines[line].fraction.shape[0]):
			lines[line].mean[bin] = np.mean(local[bin])
			lines[line].sd[bin] = np.std(local[bin])