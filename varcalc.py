import h5py as hdf
import numpy as np
from utils.iocode import getvar


def getdz(fil):
    height = getvar(fil, 'z_coords')
    dz = np.zeros_like(height)
#    for i in range(1,len(dz)-1):
#        dz[i]=.5*(height[i+1]-height[i-1])
    
    dz[:-1]=np.diff(height)
    dz[-1]=dz[-2]
    return dz

def getrho(fil):
    press = getvar(fil, 'press')
    tempk = getvar(fil, 'tempk')
    rho = (press*100) / (287*tempk)
    return rho

def meanprof(fil, varname, info = False):
    if varname == 'rho':
        var = getrho(fil)
    else:
        var = getvar(fil,varname)
    if info == True:
        print(varname,np.max(var),np.min(var))
    varout = np.mean(np.mean(var,1),1)
    return varout
