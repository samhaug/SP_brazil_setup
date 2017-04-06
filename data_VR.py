#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : data_VR.py
Purpose : get data into V and R
Creation Date : 05-04-2017
Last Modified : Wed 05 Apr 2017 02:11:47 PM EDT
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

def main():
    stz,ste,stn = read_data()
    str,stt = rotate(stn,ste)

    stz.write('sparse_Z.pk',format='PICKLE')
    str.write('sparse_R.pk',format='PICKLE')

def read_data():
    ste = obspy.read('BHE.pk')
    stn = obspy.read('BHN.pk')
    stz = obspy.read('BHZ.pk')
    ste.filter('highpass',freq=1./20)
    stn.filter('highpass',freq=1./20)
    stz.filter('highpass',freq=1./20)
    stn.interpolate(50)
    stz.interpolate(50)
    ste.interpolate(50)
    for idx,tr in enumerate(ste):
        ste[idx] = ste[idx].slice(ste[idx].stats.starttime,ste[idx].stats.starttime+650)
        stn[idx] = stn[idx].slice(stn[idx].stats.starttime,stn[idx].stats.starttime+650)
        stz[idx] = stz[idx].slice(stz[idx].stats.starttime,stz[idx].stats.starttime+650)
    return stz,ste,stn

def rotate(stn,ste):
    str,stt = seispy.rotate.rotate_ne_rt(stn,ste)
    return str,stt


main()
