#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : takeoff_angle_amplitude.py
Purpose : Find P,SV,SH amplitudes for different takeoff angles/azimuths.
          page 108-109. Aki/richards. equation (4.88). chapter : Elastic waves
          from a point dislocation source
Creation Date : 11-04-2017
Last Modified : Tue 11 Apr 2017 11:53:39 AM EDT
Created By : Samuel M. Haugland

==============================================================================
'''

import numpy as np
from matplotlib import pyplot as plt
from subprocess import call
from os import listdir
import h5py
import obspy
import seispy
from sys import stdout
from obspy.taup import TauPyModel
import mplstereonet
import itertools

def main():
    '''
    homedir = '/home/samhaug/work1/SP_brazil_data/2007-07-21-mw60-western-brazil-4/'
    Mxyz = cmt2mxyz(homedir+'CMTSOLUTION')
    st = obspy.read(homedir+'sparse_R.pk')
    stat_list = []
    rat_list = []

    for tr in st:
        azS,toaS = taup_angles(tr,['S'])
        azS1800P,toaS1800P = taup_angles(tr,['S1800P'])
        amp_S1800P = find_sv_amp(Mxyz,azS1800P,toaS1800P)
        amp_S = find_sv_amp(Mxyz,azS,toaS)

        stat_list.append(tr.stats.station)

        rat_list.append(amp_S/amp_S1800P)
        print amp_S/amp_S1800P
    write_ratio('ratio.dat',zip(stat_list,rat_list))
    '''
    beachball('/home/samhaug/work1/SP_brazil_data/2007-07-21-mw60-western-brazil-4/')

def beachball(homedir):

    def main(homedir):
        st = obspy.read(homedir+'sparse_R.pk')
        Mxyz = cmt2mxyz(homedir+'CMTSOLUTION')
        fig = plt.figure(figsize=(4,8))
        ax1 = fig.add_subplot(311, projection='stereonet')
        ax1.set_title('P')
        ax2 = fig.add_subplot(312, projection='stereonet')
        ax2.set_title('SV')
        ax3 = fig.add_subplot(313, projection='stereonet')
        ax3.set_title('SH')
        p_coords = stereonet_coords(Mxyz,find_p_amp)
        sv_coords = stereonet_coords(Mxyz,find_sv_amp)
        sh_coords = stereonet_coords(Mxyz,find_sh_amp)
        plot_coords(p_coords,ax1)
        plot_coords(sv_coords,ax2)
        plot_coords(sh_coords,ax3)
        ray_coords = get_ray_coordinates(st)
        for ii in ray_coords:
            ax1.pole(90+ii[0],ii[1],markersize=4.0,color='limegreen',mew=0.)
            ax2.pole(90+ii[0],ii[1],markersize=4.0,color='limegreen',mew=0.)
            ax3.pole(90+ii[0],ii[1],markersize=4.0,color='limegreen',mew=0.)

        plt.show()

    def stereonet_coords(Mxyz,func):
        r = np.linspace(0,90,num=190)
        t = np.linspace(0,360,num=460)
        coords = list(itertools.product(r,t))
        mag = []
        for ii in range(0,len(coords)):
            mag.append(func(Mxyz,np.radians(coords[ii][1]),np.radians(coords[ii][0])))
        for ii in range(len(mag)):
            if mag[ii] <= 0:
                coords[ii] = 0
        coords = [i for i in coords if i != 0]
        coords = np.array(coords)
        return coords

    def plot_coords(coords,ax):
        ax.pole(90+coords[:,1],coords[:,0],color='k',marker='o',rasterized=True)
        ax.set_azimuth_ticklabels([])

    def get_ray_coordinates(st):
        model = TauPyModel(model='prem50')
        ray_coord = []
        for tr in st:
            evdp = tr.stats.sac['evdp']
            gcarc = tr.stats.sac['gcarc']
            arrivals = model.get_travel_times(source_depth_in_km=evdp,
                                              distance_in_degree=gcarc,
                                              phase_list=['S1800P'])
            ang = arrivals[0].takeoff_angle
            ray_coord.append([tr.stats.sac['az'],ang])
        return ray_coord

    main(homedir)

def write_ratio(fname,ratio_list):
    with open(fname,'w') as f:
        for ii in ratio_list:
            f.write("{} {}\n".format(ii[0],ii[1]))

def taup_angles(tr,phase_list):
    model = TauPyModel(model='prem50')
    gcarc = tr.stats.sac['gcarc']
    evdp = tr.stats.sac['evdp']
    arr = model.get_travel_times(distance_in_degree=gcarc,source_depth_in_km=evdp,
                           phase_list=phase_list)
    toa = arr[0].takeoff_angle
    az = tr.stats.sac['az']
    return np.radians(toa),np.radians(az)

def orig_cmt2mxyz(file):
	li = open(file,'r').readlines()
	M = []
	M.append([])
	M[0].append(float(li[ 7].split(':')[1]))
	M[0].append(float(li[10].split(':')[1]))
	M[0].append(float(li[11].split(':')[1]))
	M.append([])
	M[1].append(float(li[10].split(':')[1]))
	M[1].append(float(li[ 8].split(':')[1]))
	M[1].append(float(li[12].split(':')[1]))
	M.append([])
	M[2].append(float(li[11].split(':')[1]))
	M[2].append(float(li[12].split(':')[1]))
	M[2].append(float(li[ 9].split(':')[1]))

	N = []
	N.append([])
	N[-1].append( M[1][1])
	N[-1].append(-M[1][2])
	N[-1].append( M[0][1])
	N.append([])
	N[-1].append(-M[1][2])
	N[-1].append( M[2][2])
	N[-1].append(-M[0][2])
	N.append([])
	N[-1].append( M[0][1])
	N[-1].append(-M[0][2])
	N[-1].append( M[0][0])
	return N

def cmt2mxyz(file):
    '''modification on what jeroen did in orig_cmt2mxyz'''
    '''See box 4.4, aki richards page 113'''
    a = np.genfromtxt(file,skip_header=7)[:,1]
    Mrr = a[0]
    Mtt = a[1]
    Mpp = a[2]
    Mrt = a[3]
    Mrp = a[4]
    Mtp = a[5]
    Mxyz = np.array([[Mrr,Mrt,Mrp],[Mrt,Mtt,Mtp],[Mrp,Mtp,Mpp]])
    new_Mxyz = Mxyz.copy()
    new_Mxyz[0,0] = Mxyz[1,1]
    new_Mxyz[0,1] = -Mxyz[1,2]
    new_Mxyz[0,2] = Mxyz[0,1]

    new_Mxyz[1,0] = -Mxyz[1,2]
    new_Mxyz[1,1] = Mxyz[2,2]
    new_Mxyz[1,2] = -Mxyz[0,2]

    new_Mxyz[2,0] = Mxyz[0,1]
    new_Mxyz[2,1] = -Mxyz[0,2]
    new_Mxyz[2,2] = Mxyz[0,0]
    return new_Mxyz

def find_p_amp(Mxyz,az,toa):
    '''az:azimuthal direction, toa:takeoff angle'''
    Ap = 0
    si,ci = np.sin(toa),np.cos(toa)
    sf,cf = np.sin(az),np.cos(az)
    g = [si*cf, si*sf, ci]
    v = [ci*cf, ci*sf,-si]
    h = [-sf, cf, 0.]
    #for k1 in xrange(3):
    #    for k2 in xrange(3):
    #        Ap += Mxyz[k1][k2]*g[k1]*g[k2]
    Ap = np.dot(np.array(g),np.dot(np.array(Mxyz),np.array(g)))
    return Ap

def find_sv_amp(Mxyz,az,toa):
    Asv = 0
    si,ci = np.sin(toa),np.cos(toa)
    sf,cf = np.sin(az),np.cos(az)
    g = [si*cf, si*sf, ci]
    v = [ci*cf, ci*sf,-si]
    h = [-sf, cf, 0.]
    #for k1 in xrange(3):
    #    for k2 in xrange(3):
    #        Asv += Mxyz[k1][k2]*g[k1]*v[k2]
    Asv = np.dot(np.array(v),np.dot(np.array(Mxyz),np.array(g)))
    return Asv

def find_sh_amp(Mxyz,az,toa):
    Ash = 0
    si,ci = np.sin(toa),np.cos(toa)
    sf,cf = np.sin(az),np.cos(az)
    g = [si*cf, si*sf, ci]
    v = [ci*cf, ci*sf,-si]
    h = [-sf, cf, 0.]
    #for k1 in xrange(3):
    #    for k2 in xrange(3):
    #        Ash += Mxyz[k1][k2]*g[k1]*h[k2]
    Ash = np.dot(np.array(h),np.dot(np.array(Mxyz),np.array(g)))
    return Ash

main()
