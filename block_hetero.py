#!/usr/bin/env python
import numpy as np
import scipy.interpolate
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema


def read_PREM():
    '''
    '''
    PREM = np.loadtxt('PREM_1s.csv',delimiter=',')

    rad =  PREM[:,0]
    rho =  PREM[:,2]
    vs  =  PREM[:,5]
    vp  =  PREM[:,3]

    fine_rad = np.linspace(rad.min(),rad.max(),1000)

    f_rho = scipy.interpolate.interp1d(rad,rho)
    f_vs = scipy.interpolate.interp1d(rad,vs)
    f_vp = scipy.interpolate.interp1d(rad,vp)

    int_rho = f_rho(fine_rad)
    int_vp = f_vp(fine_rad)
    int_vs = f_vs(fine_rad)

    return fine_rad,int_rho,int_vs,int_vp

def make_bounds(radmin,radmax,colatmin,colatmax):
    '''
    Makes bounds for cylindrical symmetry
    '''
    rad_domain = np.arange(radmin,radmax,2.)
    colat_domain = np.arange(colatmin,colatmax,1)

    colat_array,rad_array = np.meshgrid(colat_domain,rad_domain)
    return colat_array,rad_array

def make_slant_bounds(r1,alpha,L,w,phi,dvp,dvs,drho):
    '''
    cartesian slab with tilt angle phi from horizontal
    '''
    phi = np.radians(-1*phi)
    c, s = np.cos(phi), np.sin(phi)
    #R = np.matrix('{} {}; {} {}'.format(c, -1*s, s, c))

    x1,y1 = r1*np.sin(np.radians(alpha)), r1*np.cos(np.radians(alpha))
    x2,y2 = x1+L*np.cos(phi), y1-L*np.sin(phi)

    xlim = np.linspace(x1,x2,num=1000)
    print xlim.min(),xlim.max()
    ylim = np.linspace(y1,y2,num=1000)
    print ylim.min(),ylim.max()
    xtick = abs(xlim[0]-xlim[1])
    ytick = abs(ylim[0]-ylim[1])
    off_int = int(w/(2*np.hypot(xtick,ytick)))

    xx,yy = np.meshgrid(xlim,ylim)
    triu = np.triu_indices(1000,off_int)
    tril = np.tril_indices(1000,-1*off_int)
    vs_anom = np.ones(xx.shape)
    vp_anom = np.ones(xx.shape)
    rho_anom = np.ones(xx.shape)
    vs_anom[triu] = 0.
    vs_anom[tril] = 0.
    vp_anom[triu] = 0.
    vp_anom[tril] = 0.
    rho_anom[triu] = 0.
    rho_anom[tril] = 0.
    colat_array = np.degrees(np.arctan(xx/yy))
    radius_array = np.sqrt(xx**2+yy**2)

    colat_array = np.reshape(colat_array,(colat_array.size,1))
    radius_array = np.reshape(radius_array,(colat_array.size,1))
    vp_array = np.flipud(np.reshape(vp_anom,(colat_array.size,1)))
    vs_array = np.flipud(np.reshape(vs_anom,(colat_array.size,1)))
    rho_array = np.flipud(np.reshape(rho_anom,(colat_array.size,1)))

    rad,rho,vs,vp = read_PREM()

    for idx,ii in enumerate(radius_array):
        vel_ind = np.abs(rad-ii).argmin()
        if vp_array[idx]:
           vp_array[idx] = vp[vel_ind]*(1+dvp/100.)
           vs_array[idx] = vs[vel_ind]*(1+dvs/100.)
           rho_array[idx] = rho[vel_ind]*(1+drho/100.)
        elif vp_array[idx] == False:
           vp_array[idx] = vp[vel_ind]
           vs_array[idx] = vs[vel_ind]
           rho_array[idx] = rho[vel_ind]

    colat_array = np.reshape(colat_array,xx.shape)
    radius_array = np.reshape(radius_array,xx.shape)
    vp_array = np.reshape(vp_array,xx.shape)*1000
    vs_array = np.reshape(vs_array,xx.shape)*1000
    rho_array = np.reshape(rho_array,xx.shape)*1000

    return colat_array, radius_array, vp_array, vs_array, rho_array

    #return colat_array,rad_array,vp_array,vs_array,rho_array

def slant_boxes(r1,alpha,w,phi,dvp,dvs,drho):
    '''
    make slanted slab by stacking boxes next together
    '''
    phi = np.radians(phi)
    x1,y1 = int(r1*np.sin(np.radians(alpha))), int(r1*np.cos(np.radians(alpha)))

    p_rad,p_rho,p_vs,p_vp = read_PREM()
    x_list = []
    y_list = []
    for ii in range(15):
       x2,y2 = x1+w,y1-w
       #xlim = np.linspace(x1,x2)
       xlim = np.round(np.arange(x1,x2,2))
       ylim = np.round(np.linspace(y2,y1,2))
       xx,yy = np.meshgrid(xlim,ylim)
       #c,r = np.degrees(np.arctan(xx,yy)),np.sqrt(xx**2+yy**2)
       xx = np.reshape(xx,(xx.size,1))
       yy = np.reshape(yy,(yy.size,1))
       x_list.append(xx)
       y_list.append(yy)
       x1,y1 = x1+w,y1-w*np.sin(phi)

    '''
    c,r = np.degrees(np.arctan(xx,yy)),np.sqrt(xx**2+yy**2)
    c,r = np.reshape(c,(c.size,1)),np.reshape(r,(r.size,1))
    idx = np.abs(rad-rad_array[ii][0]).argmin()
    vp,vs,rho = assign_velocity(c,r,dvp,dvs,drho)
    c_a.append(c[:,0])
    r_a.append(r[:,0])
    vp_a.append(vp[:,0])
    vs_a.append(vs[:,0])
    rho_a.append(rho[:,0])
    '''
    x = np.array(x_list)
    y = np.array(y_list)
    x = np.reshape(x,(x.size,1))
    y = np.reshape(y,(y.size,1))
    c,r = np.degrees(np.arctan(x/y)), np.sqrt(x**2+y**2)

    vp,vs,rho = assign_velocity(c,r,dvp,dvs,drho)
    return np.round(c),np.round(r),vp,vs,rho

    #return np.array(c_a),np.array(r_a),np.array(vp_a),np.array(vs_a),np.array(rho_a)

def assign_velocity(colat_array,rad_array,dvp,dvs,drho,**kwargs):
    '''
    Enter perturbations in percent values.

    dvp = 4 means 4% vp increase with respect to PREM

    '''

    rad,rho,vs,vp = read_PREM()

    vp_array = np.zeros(rad_array.shape)
    vs_array = np.zeros(rad_array.shape)
    rho_array = np.zeros(rad_array.shape)

    for ii in range(rad_array.shape[0]):
        idx = np.abs(rad-rad_array[ii][0]).argmin()
        vp_array[ii,:] = vp[idx]*(1+dvp/100.)*1000
        vs_array[ii,:] = vs[idx]*(1+dvs/100.)*1000
        rho_array[ii,:] = rho[idx]*(1+drho/100.)*1000

    return vp_array,vs_array,rho_array

def write_output(colat_array,rad_array,vp_array,vs_array,rho_array,**kwargs):
    '''
    write files for AxiSEM to read.
    '''

    name = kwargs.get('name','slab.sph')
    slant = kwargs.get('slant',True)

    rad,rho,vs,vp = read_PREM()

    b_idx = np.abs(rad-3480.).argmin()
    t_idx = np.abs(rad-6371.).argmin()


    print 'Beginning output_datafiles'

    rad = np.reshape(rad_array,(rad_array.size,1))
    theta = np.reshape(colat_array,(rad_array.size,1))
    rho = np.reshape(rho_array,(rad_array.size,1))
    vp = np.reshape(vp_array,(rad_array.size,1))
    vs = np.reshape(vs_array,(rad_array.size,1))

    with open(name,'w') as f_handle:
        f_handle.write(str(rad.size)+'\n')
        np.savetxt(f_handle,np.hstack((rad,theta,vp,vs,rho)),delimiter='   ',fmt='%1.3f')

dvp = 0
dvs = -2
drho = 2

#rholist = np.arange(1,7,2)
#vslist = np.arange(-1,-7,-2)
thmin = np.arange(3,9,1)
thmax = np.arange(9,15,1)
#rhogrid,vsgrid = np.meshgrid(rholist,vslist)

intheta = zip(thmin,thmax)
#inpair = zip(vsgrid.reshape(vsgrid.size),rhogrid.reshape(rhogrid.size))

radmax = 4571
radmin = 4561
diff= str(radmax-radmin)

#thmin = 4
#thmax = 9

'''
colat_array, rad_array = make_bounds(radmin,radmax,thmin,thmax)
vp_array, vs_array, rho_array = assign_velocity(colat_array,rad_array,
                                            dvp,dvs,drho)
write_output(colat_array,rad_array,vp_array,vs_array,rho_array,
             name='dvs_{}_drho_{}_thmin_{}_thmax_{}_diff_{}.sph'.format(
                   str(dvs),str(drho),str(thmin),str(thmax),diff))


for ii in intheta:
    colat_array, rad_array = make_bounds(radmin,radmax,ii[0],ii[1])
    vp_array, vs_array, rho_array = assign_velocity(colat_array,rad_array,
                                                dvp,dvs,drho)
    write_output(colat_array,rad_array,vp_array,vs_array,rho_array,
             name='dvs_{}_drho_{}_dvp_{}_thmin_{}_thmax_{}_diff_{}.sph'.format(
                   str(dvs),str(drho),str(dvp),str(ii[0]),str(ii[1]),diff))

for ii in inpair:
    colat_array, rad_array = make_bounds(radmin,radmax,4,9)
    vp_array, vs_array, rho_array = assign_velocity(colat_array,rad_array,
                                                dvp,ii[0],ii[1])
    write_output(colat_array,rad_array,vp_array,vs_array,rho_array,
             name='dvs_{}_drho_{}_dvp_{}_thmin_{}_thmax_{}_diff_{}.sph'.format(
                   str(ii[0]),str(ii[1]),str(1),str(thmin),str(thmax),diff))

'''

'''
ang_list = [(4,9),
            (3,8),
            (3,9),
            (4,10)]

for ii in ang_list:
    colat_array, rad_array = make_bounds(4561,4571,ii[0],ii[1])
    vp_array, vs_array, rho_array = assign_velocity(colat_array,rad_array,
                                            dvp,dvs,drho)
    write_output(colat_array,rad_array,vp_array,vs_array,rho_array,
            name='min_'+str(ii[0])+'_max_'+str(ii[1])+'.sph')
'''


#Make rectangular slab at angle
phi = 30
r1=4471
alpha=4
L=150
w=10
dvp=0
dvs=-2
drho=2

#c_box,r_box,vp_box,vs_box,rho_box = slant_boxes(r1,alpha,w,phi,dvp,dvs,drho)
c_box,r_box,vp_box,vs_box,rho_box = make_slant_bounds(r1,alpha,L,w,phi,dvp,dvs,drho)
c_box = np.reshape(c_box,(c_box.size,1))
r_box = np.reshape(r_box,(r_box.size,1))
vp_box = np.reshape(vp_box,(vp_box.size,1))
vs_box = np.reshape(vs_box,(vs_box.size,1))
rho_box = np.reshape(rho_box,(rho_box.size,1))

prem_rad,prem_rho,prem_vs,prem_vp = read_PREM()
rad = np.arange(prem_rad.min(),prem_rad.max(),1.0)
f_rho = interp1d(prem_rad,prem_rho)
f_vs = interp1d(prem_rad,prem_vs)
f_vp = interp1d(prem_rad,prem_vp)

vs = f_vs(rad)
vp = f_vp(rad)
rho = f_vp(rad)
colat = np.linspace(0,15,num=1000)
prem_c, prem_r = np.meshgrid(colat,rad)
prem_vs = np.flipud(np.transpose(np.tile(vs,(len(colat),1)))*1000)
prem_vp = np.flipud(np.transpose(np.tile(vp,(len(colat),1)))*1000)
prem_rho = np.flipud(np.transpose(np.tile(rho,(len(colat),1)))*1000)
prem_r = np.flipud(prem_r)

#prem_c = np.reshape(prem_c,(prem_c.size,1))
#prem_r = np.reshape(prem_r,(prem_r.size,1))
#prem_vs = np.reshape(prem_vs,(prem_vs.size,1))
#prem_vp = np.reshape(prem_vp,(prem_vp.size,1))
#prem_rho = np.reshape(prem_rho,(prem_rho.size,1))
c_ind = prem_c[0,:]
r_ind = prem_r[:,0]
#PREM = np.hstack((prem_c,prem_r,prem_vs,prem_vp,prem_rho))

#c,r,vp,vs,rho = slant_boxes(r1,alpha,w,phi,dvp,dvs,drho)
#c,r,vp,vs,rho  make_slant_bounds(r1,alpha,L,w,phi,dvp,dvs,drho):
for idx,ii in enumerate(c_box):
    vert = np.argmin(np.abs(r_ind-r_box[idx]))
    hori = np.argmin(np.abs(c_ind-c_box[idx]))

    prem_vp[vert,hori] = vp_box[idx]
    prem_vs[vert,hori] = vs_box[idx]
    prem_rho[vert,hori] = rho_box[idx]

out_c = np.reshape(prem_c,(prem_c.size,1))
out_r = np.reshape(prem_r,(prem_r.size,1))
out_vp = np.reshape(prem_vp,(prem_vp.size,1))
out_vs = np.reshape(prem_vs,(prem_vs.size,1))
out_rho = np.reshape(prem_rho,(prem_rho.size,1))


write_output(out_c,out_r,out_vp,out_vs,out_rho,name=str(phi)+'degree_{}vs_{}rho.sph'.
             format(str(dvs),str(drho)))
#write_output(c,r,vp,vs,rho,name='header.sph')

'''
for idx,ii in enumerate(c_box):

    try:
        ind = int(s.pop())
    except KeyError:
        continue
    PREM[ind,2] = vs_box[idx]
    PREM[ind,3] = vp_box[idx]
    PREM[ind,4] = rho_box[idx]

for phi in phi_list:
   c,r,vp,vs,rho = slant_boxes(r1,alpha,w,phi,dvp,dvs,drho)
   for idx, ii in
   if phi < 0:
      phi = 'neg_'+str(abs(phi))
   write_output(c,r,vp,vs,rho,name=str(phi)+'_degree.sph')
'''
