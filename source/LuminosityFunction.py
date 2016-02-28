import numpy as np
import mpmath
from scipy import interpolate
from scipy import special
from astropy.io import fits
from scipy import optimize

#Here is the content of this script:
#luminosityParamter take the redshift of the object and interpolate those luminosity paramters according the mock catalogue.

#LuminosityFunction takes those luminosity parameters to construct correspoind luminosity function

#rootfinding part gives the new luminosity,i.e. inversing the function

#CalibrationHa and CalibrationO3b are the desired luminosity functions.

def luminosityParameterHa(z):
	bin1 = np.array([[0.6],[-2.51,41.72,-1.27]])
	bin2 = np.array([[1.2],[-2.7,42.18,-1.43]])
	output = np.array([[0. for i in range(z.size)] for i in range(3)])
	binsize = bin2[0][0] - bin1[0][0]
	for i in range(3):
		slope = (bin2[1][i]-bin1[1][i])/binsize
		output[i][np.where(z<bin1[0][0]-binsize/2)[0]] = bin1[1][i] - (slope *binsize/2)
		output[i][np.where((z>=bin1[0][0]-binsize/2)*(z<= bin2[0][0]+binsize/2))[0]] = bin1[1][i]+slope*(z-(bin1[0][0]-binsize/2))
		output[i][np.where(z>bin2[0][0]+binsize/2)[0]] = bin2[1][i]+slope*binsize/2
	output[0] = np.power(10,output[0])
	output[1] = np.power(10,output[1])
	return output
	
def luminosityParameterO3b(z):
	bin1 = np.array([[0.9],[-3.28,42.39,-1.5]])
	bin2 = np.array([[1.8],[-3.60,42.83,-1.5]])
	output = np.array([[0. for i in range(z.size)] for i in range(3)])
	binsize = bin2[0][0] - bin1[0][0]
	for i in range(3):
		slope = (bin2[1][i]-bin1[1][i])/binsize
		output[i][np.where(z<bin1[0][0]-binsize/2)[0]] = bin1[1][i] - (slope *binsize/2)
		output[i][np.where((z>=bin1[0][0]-binsize/2)*(z<= bin2[0][0]+binsize/2))[0]] = bin1[1][i]+slope*(z-(bin1[0][0]-binsize/2))
		output[i][np.where(z>bin2[0][0]+binsize/2)[0]] = bin2[1][i]+slope*binsize/2
	output[0] = np.power(10,output[0])
	output[1] = np.power(10,output[1])
	return output	

def LuminosityFunction(luminosityParameter,luminosity):
	return np.log10(luminosityParameter[0])+float(mpmath.log10(mpmath.gammainc(luminosityParameter[2]+1,luminosity/luminosityParameter[1])))

def rootfinding(luminosity,luminosityParameter,nold):
	return np.log10(luminosityParameter[0][0])+float(mpmath.log10(mpmath.gammainc(luminosityParameter[2][0]+1,luminosity/luminosityParameter[1][0])))-float(nold)

def CalibrationHa(z,oldLum,it):
	newLum = [0. for x in range((it.shape[0]))]
	oldpara = np.array([[0. for i in range(3)] for i in range(it.shape[0])])
	oldLF = np.array([0. for i in range(it.shape[0])])
	interpolatingSize =100
	for i in range(oldpara.shape[0]):
	        oldpara[i][0] = 1.37e-3
	        oldpara[i][2] = -1.35
	        if z[it[i]] <1.3 :oldpara[i][1] = 5.1e41*np.power((1+z[it[i]]),3.1)
	        if z[it[i]] >=1.3:oldpara[i][1] =6.8e42
	        oldLF[i] = LuminosityFunction(oldpara[i],oldLum[i])
	rangea = oldLF.copy()
	rangea.sort()
	rangea = rangea[rangea>0]
	grid = np.array([[0. for i in range(interpolatingSize)] for i in range(2)])
	mpgrid = mpmath.arange(mpmath.log10(rangea.min()),mpmath.log10(rangea.max()),(np.log10(rangea.max())-np.log10(rangea.min()))/interpolatingSize)
	grid[0] = np.ogrid[np.log10(rangea.min()):np.log10(oldLF.max()):interpolatingSize*1j]
	grid[1] = np.ogrid[0:6:interpolatingSize*1j]
	newLF = np.array([[0.for i in range(interpolatingSize)] for i in range(interpolatingSize)])
	for i in range(interpolatingSize):
		for j in range(interpolatingSize):
			newLF[i][j]= optimize.brentq(rootfinding,1e30,1e50,args=(luminosityParameterHa(grid[1][i]),mpgrid[j]))
	newLum = interpolate.interpn([grid[0],grid[1]],newLF,np.array([np.log10(oldLF),z[it]]).T,bounds_error=False,fill_value=0.0)
	newLum = np.power(10,newLum)
	return newLum

def CalibrationO3b(z,oldLum,it):
	newLum = [0. for x in range((it.shape[0]))]
	oldpara = np.array([[0. for i in range(3)] for i in range(it.size/it.ndim)])
	oldLF = np.array([0. for i in range(it.size/it.ndim)])
	interpolatingSize = 100
	for i in range(oldpara.shape[0]):
	        oldpara[i][0] = 1.37e-3
	        oldpara[i][2] = -1.35
	        if z[it[i]] <1.3 :oldpara[i][1] = 5.1e41*np.power((1+z[it[i]]),3.1)
	        if z[it[i]] >=1.3:oldpara[i][1] =6.8e42
	        oldLF[i] = LuminosityFunction(oldpara[i],oldLum[i])
	rangea = oldLF.copy()
	rangea.sort()
	rangea = rangea[rangea>1e-11]
	grid = np.array([[0. for i in range(interpolatingSize)] for i in range(2)])
	grid[0] = np.ogrid[np.log10(rangea.min()):np.log10(oldLF.max()):interpolatingSize*1j]
	grid[1] = np.ogrid[0:6:interpolatingSize*1j]
	newLF = np.array([[0.for i in range(interpolatingSize)] for i in range(interpolatingSize)])
	for i in range(interpolatingSize):
		for j in range(interpolatingSize):
	      	  newLF[i][j]= optimize.brentq(rootfinding,30,50,args=(luminosityParameterO3b(grid[1][i]),np.power(10,grid[0][j])))
	localinterp = interpolate.interp2d(grid[0],grid[1],newLF,fill_value=0.0)
	for i in range(oldLF.size):
		newLum[i]=float(localinterp(np.log10(oldLF[i]),z[it[i]]))
	newLum= np.power(10,newLum)
	return newLum
