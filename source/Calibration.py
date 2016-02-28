import numpy as np
import LuminosityFunction as LF
from astropy.io import fits
from scipy import optimize
from scipy import interpolate
import mpmath

cata = fits.open('CMC081211_all.fits')
table = cata[1].data
z=table.field(3)
dmz=table.field(5)
ebv=table.field(4)
Ha = table.field(85)
it = np.array(np.where(Ha>0)[0])
dl=np.array(10.**(0.2*dmz[it]+1)*3.09e18,dtype='float64')
lumHa = Ha[it]*4*np.pi*dl**2

logl1 = 35 +np.arange(134.)*0.1
y = np.genfromtxt('shiftl_ha_colbert.dat')
lumHa2 = [0 for x in range(it.shape[0])]
for i in range((it.shape[0])		):	
        iz = int(np.around(10*z[it[i]]))

        if iz <= 60:
                        lumHa2[i] = np.interp(np.log10(lumHa[i]),logl1,y[iz])     
        else:
                        lumHa2[i] = np.log10(lumHa[i])


for i in range (0,len(lumHa2)):
        lumHa2[i] = 10.**lumHa2[i]/(4*np.pi*(dl[i]**2))



newLum = [0. for x in range((it.shape[0]))]
oldpara = np.array([[0. for i in range(3)] for i in range(it.shape[0])])
oldLF = np.array([0. for i in range(it.shape[0])])
interpolatingSize =100
for i in range(oldpara.shape[0]):
        oldpara[i][0] = 1.37e-3
        oldpara[i][2] = -1.35
        if z[it[i]] <1.3 :oldpara[i][1] = 5.1e41*np.power((1+z[it[i]]),3.1)
        if z[it[i]] >=1.3:oldpara[i][1] =6.8e42
        oldLF[i] = LF.LuminosityFunction(oldpara[i],lumHa[i])
grid = np.array([[0. for i in range(interpolatingSize)] for i in range(2)])
mpgrid = mpmath.arange(oldLF[oldLF<-5].max(),oldLF.max(),(oldLF.max()-oldLF[oldLF<-5].max())/interpolatingSize)
grid[0] = np.ogrid[oldLF.min():oldLF.max():interpolatingSize*1j]
grid[1] = np.ogrid[0:6:interpolatingSize*1j]
newLF = np.array([[0.for i in range(interpolatingSize)] for i in range(interpolatingSize)])
for i in range(interpolatingSize):
	for j in range(interpolatingSize):
		newLF[i][j]= optimize.brentq(LF.rootfinding,1e30,1e50,args=(LF.luminosityParameterHa(grid[1][i]),mpgrid[j]))
newLum = interpolate.interpn([grid[1],mpgrid],newLF,np.array([z[it],oldLF]).T,bounds_error=False,fill_value=0.0)
for i in np.where(oldLF<-5)[0]:
	newLum[i] = optimize.brentq(LF.rootfinding,1e30,1e50,args=(LF.luminosityParameterHa(z[it][i]),oldLF[i]))
