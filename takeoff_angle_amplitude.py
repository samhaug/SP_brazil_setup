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
    Mxyz = cmt2mxyz('CMTSOLUTION')
    st = obspy.read('')

    Asv = find_sv_amp(Mxyz,az,toa)

def taup_angles(tr,phase_list):
    model = TauPyModel(model='prem50')
    gcarc = tr.stats.sac['gcarc']
    evdp = tr.stats.sac['evdp']
    arr = model.get_travel_times(distance_in_degree=gcarc,source_depth_in_km=,evdp,
                           phase_list=phase_list)
    toa = arr[0].takeoff_angle
    az = tr.stats.sac['az']
    return toa,az

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

def find_p_amp(M,az,toa):
    '''az:azimuthal direction, toa:takeoff angle'''
    si,ci = np.sin(toa),np.cos(toa)
    sf,cf = np.sin(az),np.cos(az)
    g = [si*cf, si*sf, ci]
    v = [ci*cf, ci*sf,-si]
    h = [-sf, cf, 0.]
    for k1 in xrange(3):
        for k2 in xrange(3):
            Ap  += Mxyz[k1][k2]*g[k1]*g[k2]
    return Ap

def find_sv_amp(M,az,toa):
    si,ci = np.sin(toa),np.cos(toa)
    sf,cf = np.sin(az),np.cos(az)
    g = [si*cf, si*sf, ci]
    v = [ci*cf, ci*sf,-si]
    h = [-sf, cf, 0.]
    for k1 in xrange(3):
        for k2 in xrange(3):
            Asv += Mxyz[k1][k2]*g[k1]*v[k2]
    return Asv

def find_sh_amp(M,az,toa):
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
