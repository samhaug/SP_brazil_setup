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

def main():
    homedir = '/home/samhaug/work1/SP_brazil_data/2007-07-21-mw60-western-brazil-4/'
    Mxyz = cmt2mxyz(homedir+'CMTSOLUTION')
    st = obspy.read(homedir+'sparse_R.pk')
    stat_list = []
    rat_list = []

    for tr in st:
        azS,toaS = taup_angles(tr,['S'])
        azS1800P,toaS1800P = taup_angles(tr,['S1800P'])
        #print toaS,toaS1800P

        S_Asv = find_sv_amp(Mxyz,azS,toaS)
        #print S_Asv
        S1800P_Asv = find_sv_amp(Mxyz,azS1800P,toaS1800P)
        #print S1800P_Asv
        stat_list.append(tr.stats.station)
        rat_list.append(S_Asv/S1800P_Asv)

    write_ratio('ratio.dat',zip(stat_list,rat_list))

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

def cmt2mxyz(file):
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

def find_p_amp(Mxyz,az,toa):
    '''az:azimuthal direction, toa:takeoff angle'''
    Ap = 0
    si,ci = np.sin(toa),np.cos(toa)
    sf,cf = np.sin(az),np.cos(az)
    g = [si*cf, si*sf, ci]
    v = [ci*cf, ci*sf,-si]
    h = [-sf, cf, 0.]
    for k1 in xrange(3):
        for k2 in xrange(3):
            Ap += Mxyz[k1][k2]*g[k1]*g[k2]
    return Ap

def find_sv_amp(Mxyz,az,toa):
    Asv = 0
    si,ci = np.sin(toa),np.cos(toa)
    sf,cf = np.sin(az),np.cos(az)
    g = [si*cf, si*sf, ci]
    v = [ci*cf, ci*sf,-si]
    h = [-sf, cf, 0.]
    for k1 in xrange(3):
        for k2 in xrange(3):
            Asv += Mxyz[k1][k2]*g[k1]*v[k2]
    return Asv

def find_sh_amp(Mxyz,az,toa):
    Ash = 0
    si,ci = np.sin(toa),np.cos(toa)
    sf,cf = np.sin(az),np.cos(az)
    g = [si*cf, si*sf, ci]
    v = [ci*cf, ci*sf,-si]
    h = [-sf, cf, 0.]
    for k1 in xrange(3):
        for k2 in xrange(3):
            Ash += Mxyz[k1][k2]*g[k1]*h[k2]
    return Ash

main()
