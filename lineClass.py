import numpy as np	


class line:
	
	def __init__(self,name,wavelength):
		self.name = name
		self.wavelength = wavelength
	
	flux = np.array([],dtype=float)
	Lambda = np.array([],dtype=float)
	noise = np.array([],dtype=float)
	intf = np.array([],dtype=float)
	fraction = np.array([[]],dtype=float)
	origin = np.array([],dtype=float)
	mean = np.array([],dtype=float)
	sd = np.array([],dtype=float)
	condition = np.array([],dtype=float)
	photomCutForInt = np.array([],dtype=str)	
	def selection(self,condition):
		local = np.where(condition)[0]
		self.intf = np.intersect1d(local,self.intf)

	def conditionalSelection(self,condition,selection):
		print self.intf.size
		leftover = np.intersect1d(self.intf,np.where(condition*selection)[0])
		print leftover.size
		self.intf = np.intersect1d(self.intf,np.where(~(condition))[0])
		print self.intf.size
		self.intf = np.append(self.intf,leftover)
		self.intf.sort() 
		print self.intf.size
	
	
class photocut:
	def __init__(self,name,dep):
		self.name = name
		self.dep = dep
		self.err = np.power(10,(-0.4*(self.dep-23.9)))*1.e-29*0.2*2
	flux = np.array([],dtype=float)
	origin = np.array([],dtype=float)
	mag = np.array([],dtype=float)
