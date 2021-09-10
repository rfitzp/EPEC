# -*-Python-*-
# Created by fitzpatrickr on 7 Sep 2021

# Script consolidates NETCDF file output from PHASE

import netCDF4 as nc
import numpy as np
import glob
import os

ncfiles = sorted(glob.glob('Outputs/ncFiles/*.nc'))

filenames = {}
keyno = -1
iresold = 0

for f in ncfiles:
    ds = nc.Dataset(f)
    mpol = np.asarray(ds['m_pol'])
    ires = mpol.size

    if iresold != ires:
        keyno = keyno + 1
        filenames[keyno] = [f]
    else:
        filenames[keyno].append(f)

    iresold = ires    
    
for k in filenames.keys():

    lock = 1
    for f in filenames[k]:
        
        ds = nc.Dataset(f)
        mpol = np.asarray(ds['m_pol'])
        psin = np.asarray(ds['PsiN_res'])

        if lock:
            time = np.asarray(ds['time'])
            irmp = np.asarray(ds['Irmp'])
            prmp = np.asarray(ds['Prmp'])
            phi = np.asarray(ds['phi'])
            phid = np.asarray(ds['phi_dot'])
            omega0 = np.asarray(ds['omega0'])
            omega = np.asarray(ds['omega'])
            w = np.asarray(ds['W'])
            wv = np.asarray(ds['W_vac'])
            psim = np.asarray(ds['Psi_-'])
            psip = np.asarray(ds['Psi_+'])
            psivm = np.asarray(ds['Psi_v-'])
            psivp = np.asarray(ds['Psi_v+'])
            deltap = np.asarray(ds['DeltaP'])
            deltap0 = np.asarray(ds['DeltaP0'])
            lock = 0
        else:
            time = np.append(time, np.asarray(ds['time']), 0)
            irmp = np.append(irmp, np.asarray(ds['Irmp']), 0)
            prmp = np.append(prmp, np.asarray(ds['Prmp']), 0)
            phi = np.append(phi, np.asarray(ds['phi']), 0)
            phid = np.append(phid, np.asarray(ds['phi_dot']), 0)
            omega0 = np.append(omega0, np.asarray(ds['omega0']), 0)
            omega = np.append(omega, np.asarray(ds['omega']), 0)
            w = np.append(w, np.asarray(ds['W']), 0)
            wv = np.append(wv, np.asarray(ds['W_vac']), 0)
            psim = np.append(psim, np.asarray(ds['Psi_-']), 0)
            psip = np.append(psip, np.asarray(ds['Psi_+']), 0)
            psivm = np.append(psivm, np.asarray(ds['Psi_v-']), 0)
            psivp = np.append(psivp, np.asarray(ds['Psi_v+']), 0)
            deltap = np.append(deltap, np.asarray(ds['DeltaP']), 0)
            deltap0 = np.append(deltap0, np.asarray(ds['DeltaP0']), 0)
            
    ncfilename = 'Outputs/ncFiles1/' + str(k) + '.nc'    
    ncfile = nc.Dataset(ncfilename, 'w')
    xires = ncfile.createDimension('i_res', mpol.shape[0])
    xtime = ncfile.createDimension('time', time.shape[0])

    ympol = ncfile.createVariable('m_pol', int, ('i_res',))
    ypsin = ncfile.createVariable('PsiN_res', np.double, ('i_res',))
    ytime = ncfile.createVariable('time', np.double, ('time',))
    yirmp = ncfile.createVariable('Irmp', np.double, ('time',))
    yprmp = ncfile.createVariable('Prmp', np.double, ('time',))
    yphi = ncfile.createVariable('phi', np.double, ('time', 'i_res'))
    yphid = ncfile.createVariable('phi_dot', np.double, ('time', 'i_res'))
    yomega0 = ncfile.createVariable('omega0', np.double, ('time', 'i_res'))
    yomega = ncfile.createVariable('omega', np.double, ('time', 'i_res'))
    yw = ncfile.createVariable('W', np.double, ('time', 'i_res'))
    ywv = ncfile.createVariable('W_vac', np.double, ('time', 'i_res'))
    ypsim = ncfile.createVariable('Psi_-', np.double, ('time', 'i_res'))
    ypsip = ncfile.createVariable('Psi_+', np.double, ('time', 'i_res'))
    ypsivm = ncfile.createVariable('Psi_v-', np.double, ('time', 'i_res'))
    ypsivp = ncfile.createVariable('Psi_v+', np.double, ('time', 'i_res'))
    ydeltap = ncfile.createVariable('DeltaP', np.double, ('time', 'i_res'))
    ydeltap0 = ncfile.createVariable('DeltaP0', np.double, ('time', 'i_res'))

    ympol[:] = mpol
    ypsin[:] = psin
    ytime[:] = time
    yirmp[:] = irmp
    yprmp[:] = prmp
    yphi[:,:] = phi
    yphid[:,:] = phid
    yomega0[:,:] = omega0
    yomega[:,:] = omega
    yw[:,:] = w
    ywv[:,:] = wv
    ypsim[:,:] = psim
    ypsip[:,:] = psip
    ypsivm[:,:] = psivm
    ypsivp[:,:] = psivp
    ydeltap[:,:] = deltap
    ydeltap0[:,:] = deltap0

    ncfile.close()

    os.system('touch Outputs/ncFiles1/Index')
            
    
