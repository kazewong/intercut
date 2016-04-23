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
OIIIa=table.field(77)
OIIIb=table.field(81)
it = np.array(np.where(Ha>0)[0])
dl=np.array(10.**(0.2*dmz[it]+1)*3.09e18,dtype='float64')
lumHa = Ha[it]*4*np.pi*dl**2
lumO3a = OIIIa[it]*4*np.pi*dl**2
lumO3b = OIIIb[it]*4*np.pi*dl**2


# logl1 = 35 +np.arange(134.)*0.1
# y = np.genfromtxt('shiftl_ha_colbert.dat')
# lumHa2 = [0 for x in range(it.shape[0])]
# for i in range((it.shape[0])		):	
#         iz = int(np.around(10*z[it[i]]))
# 
#         if iz <= 60:
#                         lumHa2[i] = np.interp(np.log10(lumHa[i]),logl1,y[iz])     
#         else:
#                         lumHa2[i] = np.log10(lumHa[i])
# 
# 
# for i in range (0,len(lumHa2)):
#         lumHa2[i] = 10.**lumHa2[i]/(4*np.pi*(dl[i]**2))
# 
# 
# 
# newLum = [0. for x in range((it.shape[0]))]
# oldpara = np.array([[0. for i in range(3)] for i in range(it.shape[0])])
# oldLF = np.array([0. for i in range(it.shape[0])])
# interpolatingSize =100
# for i in range(oldpara.shape[0]):
#         oldpara[i][0] = 1.37e-3
#         oldpara[i][2] = -1.35
#         if z[it[i]] <1.3 :oldpara[i][1] = 5.1e41*np.power((1+z[it[i]]),3.1)
#         if z[it[i]] >=1.3:oldpara[i][1] =6.8e42
#         oldLF[i] = LF.LuminosityFunction(oldpara[i],lumHa[i])
# grid = np.array([[0. for i in range(interpolatingSize)] for i in range(2)])
# mpgrid = mpmath.linspace(oldLF[oldLF<-5].max(),oldLF.max(),interpolatingSize)
# grid[0] = np.ogrid[oldLF.min():oldLF.max():interpolatingSize*1j]
# grid[1] = np.ogrid[0:6:interpolatingSize*1j]
# newLF = np.array([[0.for i in range(interpolatingSize)] for i in range(interpolatingSize)])
# for i in range(interpolatingSize):
# 	for j in range(interpolatingSize):
# 		newLF[i][j]= optimize.brentq(LF.rootfinding,1e30,1e50,args=(LF.luminosityParameterHa(grid[1][i]),mpgrid[j]))
# newLum = interpolate.interpn([grid[1],mpgrid],newLF,np.array([z[it],oldLF]).T,bounds_error=False,fill_value=0.0)
# for i in np.where(oldLF<-5)[0]:
# 	newLum[i] = optimize.brentq(LF.rootfinding,1e30,1e50,args=(LF.luminosityParameterHa(z[it][i]),oldLF[i]))

logl1 = 34 +np.arange(141.)*0.1
y = np.genfromtxt('shiftl_o3_colbert_ha_geach.dat')
lumO3b2 = [0. for x in range((it.shape[0]))]
lumO3a2 = [0. for x in range((it.shape[0]))]

for i in range(it.shape[0]):	
        iz = int(np.around(10*z[it[i]]))
        if iz <= 60:
                        lumO3b2[i] = np.interp(np.log10(lumHa[i]),logl1,y[iz])
        else:
                        lumO3b2[i] = np.log10(lumO3b[i])

for i in range(0,len(lumO3b2)):
        lumO3b2[i] = 10.**lumO3b2[i]/(4*np.pi*(dl[i]**2))

kO3b=4.05+2.659*(-2.156+1.509/.5007-0.198/(.5007**2)+0.011/(.5007**3))
kO3a=4.05+2.659*(-2.156+1.509/.4959-0.198/(.4959**2)+0.011/(.4959**3))


lumO3a2 = np.array(lumO3a2)
lumO3b2 = np.array(lumO3b2)
for i in range (len(lumO3b2)):
        lumO3a2[i]= (lumO3b2[i]*(10.**(0.4*ebv[i]*kO3b))/3.)*(10.**(-0.4*ebv[i]*kO3a))
        
newLumO3b = [0. for x in range((it.shape[0]))]
oldpara = np.array([[0. for i in range(3)] for i in range(it.shape[0])])
oldLF = np.array([0. for i in range(it.shape[0])])
interpolatingSize =100
for i in range(oldpara.shape[0]):
        oldpara[i][0] = 1.37e-3
        oldpara[i][2] = -1.35
        if z[it[i]] <1.3 :oldpara[i][1] = 5.1e41*np.power((1+z[it[i]]),3.1)
        if z[it[i]] >=1.3:oldpara[i][1] =6.8e42
        oldLF[i] = LF.LuminosityFunction(oldpara[i],lumO3b[i])
grid = np.array([[0. for i in range(interpolatingSize)] for i in range(2)])
mpgrid = mpmath.linspace(oldLF[oldLF<-5].max(),oldLF.max(),interpolatingSize)
grid[0] = np.ogrid[oldLF.min():oldLF.max():interpolatingSize*1j]
grid[1] = np.ogrid[0:6:interpolatingSize*1j]
newLF = np.array([[0.for i in range(interpolatingSize)] for i in range(interpolatingSize)])
for i in range(interpolatingSize):
	for j in range(interpolatingSize):
		newLF[i][j]= optimize.brentq(LF.rootfinding,1e30,1e50,args=(LF.luminosityParameterO3b(grid[1][i]),mpgrid[j]))
newLumO3b = interpolate.interpn([grid[1],mpgrid],newLF,np.array([z[it],oldLF]).T,bounds_error=False,fill_value=0.0)
for i in np.where(oldLF<-5)[0]:
	newLumO3b[i] = optimize.brentq(LF.rootfinding,1e30,1e50,args=(LF.luminosityParameterO3b((z[it][i]+1)*(6563./5007.)-1),oldLF[i]))