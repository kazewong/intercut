import numpy as np	


class line:
	
	def __init__(self,name,wavelength):
		self.name = name
		self.wavelength = wavelength
	
	flux = np.array([],dtype=float)
	Lambda = np.array([],dtype=float)
	noise = np.array([],dtype=float)
	intf = np.array([],dtype=float)
	fraction = []
	origin = np.array([],dtype=float)
	mean = np.array([],dtype=float)
	sd = np.array([],dtype=float)
	condition = np.array([],dtype=float)
	
	def selection(self,condition):
		local = np.where(condition)[0]
		self.intf = np.intersect1d(local,self.intf)

	def conditionalSelection(self,condition,selection):
		leftover = np.intersect1d(self.intf,np.where(condition*selection)[0])
		self.intf = np.intersect1d(self.intf,np.where(~(condition))[0])
		self.intf = np.append(self.intf,leftover)
		self.intf.sort() 
		
	def getStat(self):
		self.fraction = self.fraction.T
		self.mean = np.array([0. for i in range(self.fraction.shape[0])])
		self.sd = np.array([0. for i in range(self.fraction.shape[0])])
		self.fraction = np.divide(self.fraction,(self.fraction+1.))
		for i in range(self.fraction.shape[0]):
			self.mean[i] = np.mean(self.fraction[i])
			self.sd[i] = np.std(self.fraction[i])
	
	
class photocut:
	def __init__(self,name,dep):
		self.name = name
		self.dep = dep
		self.err = np.power(10,(-0.4*(self.dep-23.9)))*1.e-29*0.2*2
	flux = np.array([],dtype=float)
	origin = np.array([],dtype=float)
	mag = np.array([],dtype=float)
