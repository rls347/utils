import numpy as np
from utils.iocode import getvar

def cloudtop(fil):
    cond = getvar(fil,'total_cond')
    z = getvar(fil,'z_coords')
    cond[cond<0.01]=0.0
    cond[cond>0.0]=1.0
    cond = cond*z[:,None,None]
    cloudtop = np.max(cond,0)
    return cloudtop


def findequal(x,y,z):
    zclose = np.argmin(np.abs(x-y))
    if np.max(x-y) <= 0 or np.min(x-y) >= 0:
        zm = z[zclose]
    else:
        if y[zclose] < x[zclose]:
            if y[zclose+1] < x[zclose+1]:
                z2 = zclose
                z1 = zclose-1
            else:
                z1 = zclose
                z2 = zclose+1
        else:
            if y[zclose+1] > x[zclose+1]:
                z2 = zclose
                z1 = zclose-1
            else:
                z1 = zclose
                z2 = zclose+1

        dydz = (y[z2]-y[z1])/(z[z2]-z[z1])
        dxdz = (x[z2]-x[z1])/(z[z2]-z[z1])

        zm = ( y[z1]-x[z1] + z[z1]*(dxdz-dydz) ) / (dxdz-dydz)

    return zm

