import numpy as np

def secondLineCutting(lines,z):
	for lineObject in lines:
		print 'cutting '+ lines[lineObject].name+' line, '+str(lines[lineObject].condition.size)+' condition found.'
		for lineObjectCondition in lines[lineObject].condition:
			if eval(lineObjectCondition[0]) == 'normal':
				flux = lines[eval(lineObjectCondition)[1]].flux
				noise = lines[eval(lineObjectCondition)[1]].noise
				wavelength = lines[eval(lineObjectCondition)[1]].wavelength
				lines[lineObject].selection(eval(lineObjectCondition)[2])
			if eval(lines.items()[lineObjectNumber].condition[lineObjectCondition])[0] == 'conditional':
				zd = (lines[lineObject].wavelength/lines[eval(lineObjectCondition)[1]].wavelength)*(1+z)-1
				lines[lineObject].conditionalSelection(eval(lineObjectCondition)[2],eval(lineObjectCondition)[3])w

def lineFluxPertubation(lines,threshold):
	for lineObject in lines:
		ferr = np.random.normal(size = lines[lineObject][lines[lineObject].noise<threshold].size)
		lines[lineObject][lines[lineObject].noise<threshold].flux = lines[lineObject][lines[lineObject].noise<threshold].origin + lines[lineObject][lines[lineObject].noise<threshold].noise*ferr
		
def getFraction(z,lines,ranges,bin = 10):
	mainLine = lines.items()[0]
	hmain = np.array(np.histogram(z[mainLine.intf],bins = bin,range = ranges)[0],dtype = float)
	for lineObject in lines:
		if lines[lineObject] !=[]:
			hlocal = np.array(np.histogram((lines[lineObject].wavelength/lines.item()[0].wavelength)*(1+(z[lines[lineObject].intf]))-1)[0])
		else:
			hlocal = np.array(np.histogram(0,bins =bin,range = ranges)[0])		
		lines[lineObject].fraction.append(np.divide(hlocal,hmain))

		