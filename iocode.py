import numpy as np
import h5py as hdf
import hdf5plugin
from netCDF4 import Dataset
import xarray as xr
from scipy.interpolate import interp1d

def getvar(fil, varname):
    try:
        var = np.squeeze(fil[varname][:])
    except TypeError:
        if fil[-3:] == '.h5':
            filey = hdf.File(fil, 'r')
            var = np.squeeze(filey[varname])
            filey.close()
        if fil[-3:] == '.nc':
            filey = Dataset(fil, mode = "r")
            var = np.squeeze(np.asarray(filey.variables[varname]))
    return var
    
def regrid(varin, z, newz):
    try:
        func = interp1d(z,varin,axis=1)
        varout = func(newz)
    except IndexError:
        func = interp1d(z,varin)
        varout = func(newz)

    return varout
    
def writehdf(filename, vars):
    with hdf.File(filename,'w') as hf:
        for varname in vars.keys():
            var = vars[varname]
            a = hf.create_dataset(varname,data=var)
    return filename

def writenc(filename, vars):
    vardic = {}
    print('this code is incomplete')
    ds = xr.Dataset(vardic)
    ds.to_netcdf(filename)
    return filename

def ncinfo(f1):
    f=Dataset(f1,'r')
    for a in f.variables.keys():
        print(a,getvar(f1,a).shape)
    f.close()
    return

