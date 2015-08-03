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
		
def  