#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : svaxi_setup.py
Purpose : setup svaxi simulations such that S1800P has vertical component amp
Creation Date : 06-04-2017
Last Modified : Thu 06 Apr 2017 12:31:40 PM EDT
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

def read_synth():
    print('setup')
    st = obspy.read('*0.R')
    st = seispy.filter.range_filter(st,(56,72))
    for tr in st:
        tr.stats.sac['o'] += -202
    st = seispy.data.normalize_on_phase(st,phase=['S'],min=True)
    for tr in st:
        tr.data *= 1./np.tan(np.radians(17))
    st.write('stv.pk',format='PICKLE')

read_synth()
