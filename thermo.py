import h5py as hdf
import numpy as np
from utils.iocode import getvar
from utils.varcalc import getdz
from utils.basiccalc import findequal

G = 9.8
CP = 1004.
RD = 287.
EPS = .622
P0 = 1000.
LV = 2.5e6
CV = CP-RD
GAMMAD = G/CP
KAPPA = RD/CP

def satmixratio(tempk,press):
    '''expects temperature in kelvin, pressure in millibars'''
    tc = tempk-273.15
    es = 6.112 * np.exp( (17.67*tc)/(tc+243.5) )  #in mb
    rs = (EPS *es) / (press-es)
    return rs

def find_lcl_sfc(height,tempk,q,press):  
    '''expects height, temperature in kelvin, q in g/kg, pressure in mb'''
    tprev = tempk[1]
    qsfc = q[1]/1000.
#    for z in range(2,len(q)):
#        p=press[z]
#        t = tprev - ((height[z]-height[z-1]) * GAMMAD)
#        rs = satmixratio(t,p)
#        tprev = t
#        if rs <= qsfc:
#            lcl = height[z]
#            tlcl = t
#            break
#    indlcl = z
    
    qline = np.ones(len(height))*qsfc
    qsline = np.zeros(len(height))
    tprev = tempk[1]
    for z in range(2,len(height)):
        p=press[z]
        t = tprev - ((height[z]-height[z-1]) * GAMMAD)
        qsline[z] = satmixratio(t,p)
        tprev = t

    testz = findequal(qline,qsline,height)
    #print(testz,lcl,height[indlcl-1],height[indlcl+1])

    return testz  # indlcl, lcl, tlcl   ###TODO decide whether to still return Temp LCL

def find_lcl_ml100(height,tempk,q,press):
    '''expects height, temperature in kelvin, q in g/kg, pressure in mb'''
    psfc = press[1]
    lev100 = np.argmin(np.abs(press-(psfc-100)))
    #lev100 = np.max(np.where(height<200.))

    rho = (press*100)/(RD*tempk)

    meanz = np.mean(height[1:lev100+1]*rho[1:lev100+1])/np.mean(rho[1:lev100+1])
    meant = np.mean(tempk[1:lev100+1]*rho[1:lev100+1])/np.mean(rho[1:lev100+1])
    meanq = np.mean(q[1:lev100+1]*rho[1:lev100+1])/np.mean(rho[1:lev100+1])

    zstart = np.argmin(np.abs(height-meanz))
    print(zstart)
    tprev = meant
    qsfc = meanq/1000.
#    for z in range(zstart,len(q)):
#        p=press[z]
#        t = tprev - ((height[z]-height[z-1]) * GAMMAD)
#        rs = satmixratio(t,p)
#        tprev = t
#        if rs <= qsfc:
#            lcl = height[z]
#            tlcl = t
#            break
#    indlcl = z

    qline = np.ones(len(height))*qsfc
    qsline = np.zeros(len(height))
    tprev = meant
    for z in range(2,len(height)):
        p = press[z]
        if z<zstart:
            t = tprev
        else:
            t = tprev - ((height[z]-height[z-1]) * GAMMAD)
        qsline[z] = satmixratio(t,p)
        tprev = t

    testz = findequal(qline,qsline,height)
 #   print(testz,lcl,height[indlcl-1],height[indlcl+1])



    return testz #indlcl, lcl, tlcl

def find_lcl_mostunstable(height,tempk,q,press):
    '''expects height, temperature in kelvin, q in g/kg, pressure in mb'''
    rho = (press*100)/(RD*tempk)
    theta = tempk*(P0/press)**(KAPPA)
    thetae = theta * np.exp((LV*(q/1000.))/(CP*tempk))
    thetae = thetae[height<5000]

    zstart = np.argmax(thetae)
    if zstart ==0:
        zstart =1
    tprev = tempk[zstart] 
    qsfc = q[zstart]/1000.
#    for z in range(zstart,len(q)):
#        p=press[z]
#        t = tprev - ((height[z]-height[z-1]) * GAMMAD)
#        rs = satmixratio(t,p)
#        tprev = t
#        if rs <= qsfc:
#            lcl = height[z]
#            tlcl = t
#            break
#    indlcl = z

    qline = np.ones(len(height))*qsfc
    qsline = np.zeros(len(height))
    tprev = tempk[zstart]
    for z in range(2,len(height)):
        p = press[z]
        if z<zstart:
            t = tprev
        else:
            t = tprev - ((height[z]-height[z-1]) * GAMMAD)
        qsline[z] = satmixratio(t,p)
        tprev = t

    testz = findequal(qline,qsline,height)
#    print(testz,lcl,height[indlcl-1],height[indlcl+1])

    return testz #indlcl, lcl, tlcl #, zstart

def moistadiabat(height,press,zlcl,tlcl,tempk,q):
    parcel = np.zeros_like(press)
    parcel[zlcl]=tlcl#tempk[zlcl]
    for z in range(zlcl+1,len(press)):
        t=parcel[z-1]
        #qv=q[z-1]/1000.
        qv = satmixratio(t,press[z-1])
        lqrt = (LV*qv)/(RD*t)
        lqrt2 = ((LV**2)*qv*(EPS+qv)) / (RD*(t**2))
        gammaPA = G * ( ((1+qv)*(1+lqrt)) / (CP + (qv*CV) + lqrt2) )
        parcel[z]= parcel[z-1] - ((height[z]-height[z-1]) * gammaPA)
    for z in range(zlcl-1,-1,-1):
        first = parcel[z+1] + ( (height[z+1]-height[z]) * GAMMAD)
        parcel[z]=first
    return parcel

###TODO fix this for the more specific LCLs
def capecalc(height,tempk,q,press,flag):
    if flag == 'sfc':
        zlcl, lcl, tlcl = find_lcl_sfc(height,tempk,q,press)
        zstart = 1
    elif flag == 'mu':
        zlcl, lcl, tlcl, zstart = find_lcl_mostunstable(height,tempk,q,press)
    else:
        zlcl, lcl, tlcl = find_lcl_ml100(height,tempk,q,press)
        zstart = np.max(np.where(height<200.))

    parcel = moistadiabat(height,press,zlcl,tlcl,tempk,q)
    parceltv = parcel*(1+(.61*satmixratio(parcel,press)))
    parceltv[:zlcl+1] = parcel[:zlcl+1]*(1+(.61*(q[0]/1000.)))
    Tv = tempk*(1+(.61*(q/1000.)))
    parceltv[:zstart+1] = Tv[:zstart+1]


    posarea = np.where(parceltv>Tv)
    intvar = (parceltv-Tv)*(G/Tv)

    if len(posarea[0])>1:
        lfc = np.min(np.asarray(posarea))
        if lfc<zlcl:
            lfc=zlcl
        el = np.max(np.asarray(posarea))
        intvarpos=intvar[lfc:el+1]
        zpos = height[lfc:el+1]
        intvarneg=intvar[0:lfc+1]
        zneg = height[0:lfc+1]

#    cape1 = np.sum(intvar[posarea])
        cape = np.trapz(intvarpos, dx=np.diff(zpos))
        cin = np.trapz(intvarneg, dx=np.diff(zneg))
        if el>10:
            xxparcel = np.linspace(parceltv[el],parceltv[el+1],height[el+1]-height[el])
            xxenv = np.linspace(Tv[el],Tv[el+1],height[el+1]-height[el])
            xxheight = np.linspace(height[el],height[el+1],height[el+1]-height[el])
            diffmin = np.argmin(np.abs(xxenv-xxparcel))
            lnb = xxheight[diffmin]
        else:
            lnb=0.
    else:
        cape=0
        cin= -5000
    return cape,cin

def getcapeparcel(filename,flag=None):
    try:
        fil = hdf.File(filename, 'r')
        height = getvar(fil, 'z_coords')
        if height.size<3:
            height = getvar('/nobackup/rstorer/convperts/revu/aug11-control/aug11-control-revu-001.h5','z_coords')
        dz = np.zeros_like(height)
        dz[0:]=np.diff(height)
        dz[-1]=dz[-2]
        tempk = getvar(fil, 'tempk')[:,10,10]
        q = getvar(fil, 'vapor')[:,10,10]
        press = getvar(fil, 'press')[:,10,10]
        fil.close()
    except:
        print('error reading variables')

    if flag == 'sfc':
        zlcl, lcl, tlcl = find_lcl_sfc(height,tempk,q,press)
    elif flag == 'mu':
        zlcl, lcl, tlcl, zstart  = find_lcl_mostunstable(height,tempk,q,press)
    else:
        zlcl, lcl, tlcl = find_lcl_ml100(height,tempk,q,press)


    parcel = moistadiabat(height,press,zlcl,tlcl,tempk,q)
    parceltv = parcel*(1+(.61*satmixratio(parcel,press)))
    parceltv[:zlcl+1] = parcel[:zlcl+1]*(1+(.61*q[0]))
    Tv = tempk*(1+(.61*(q/1000.)))
    posarea = np.where(parceltv>Tv)
    intvar = (parceltv-Tv)*(G/Tv)
    lfc = np.min(np.asarray(posarea))
    el = np.max(np.asarray(posarea))

    if lfc<zlcl:
        lfc=zlcl

    intvarpos=intvar[lfc:el+1]
    zpos = height[lfc:el+1]

    cape = np.trapz(intvarpos, dx=np.diff(zpos))
    return cape, parceltv-Tv
 



def get_cape(filename,flag=None,coords = None):
    fil = hdf.File(filename, 'r')
    height = getvar(fil, 'z_coords')
    if height.size <3:
        height = getvar('/nobackup/rstorer/convperts/revu/aug11-control/aug11-control-revu-001.h5','z_coords')
        dz = getdz('/nobackup/rstorer/convperts/revu/aug11-control/aug11-control-revu-001.h5')
    else:
        dz = getdz(fil)
    tempk = getvar(fil, 'tempk')
    q = getvar(fil, 'vapor')
    press = getvar(fil, 'press')
    fil.close()
#    except:
#        print 'error reading variables'

    if type(coords) is tuple:
        c0,c1 = coords
        capevar,cinvar = capecalc(height,tempk[:,c0,c1],q[:,c0,c1],press[:,c0,c1],flag)
    else:
        capevar = np.zeros((tempk.shape[1],tempk.shape[2]))
        cinvar = np.zeros_like(capevar)
        for c0 in range(capevar.shape[0]):
            for c1 in range(capevar.shape[1]):
                capevar[c0,c1],cinvar[c0,c1] = capecalc(height,tempk[:,c0,c1],q[:,c0,c1],press[:,c0,c1],flag)

    return capevar, cinvar
