#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : rotate_svaxi.py
Purpose : rotate all streams to lq and save. run in svaxi directory
Creation Date : 28-03-2017
Last Modified : Tue 04 Apr 2017 04:30:34 PM EDT
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
from obspy.taup import TauPyModel
model = TauPyModel(model='prem50')

def main():
    str = obspy.read('*0.R')
    rotate_stream(str)

def rotate_stream(str):
    sts = str.copy()
    stp = str.copy()
    stp_norm = str.copy()
    h = str[0].stats.sac['evdp']
    for idx,tr in enumerate(str):
        stp[idx].stats.sac['o'] += -202.
        stp_norm[idx].stats.sac['o'] += -202.
        sts[idx].stats.sac['o'] += -202.
        str[idx].stats.sac['o'] += -202.
        arrivals = model.get_travel_times(source_depth_in_km=h,
                                          distance_in_degree=tr.stats.sac['gcarc'],
                                          phase_list=['S','S1800P'])
        i_s = np.radians(arrivals[0].incident_angle)
        i_p = np.radians(arrivals[1].incident_angle)
        stp[idx].data = tr.data/np.sin(i_p)
        sts[idx].data = tr.data/np.sin(i_s)
        a = seispy.data.phase_window(sts[idx],phase=['S'])
        stp_norm[idx].data = stp[idx].data/np.max(np.abs(a.data))
    stp.write('stp.pk',format='PICKLE')
    stp_norm.write('stp_norm.pk',format='PICKLE')
    sts.write('sts.pk',format='PICKLE')

main()
