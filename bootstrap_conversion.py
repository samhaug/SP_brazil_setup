#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : bootstrap_conversion.py
Purpose : Use to bootstrap the aligned h5py phases. Need to use align_on_conversion.py
Creation Date : 09-04-2017
Last Modified : Sun 16 Apr 2017 08:23:42 PM EDT
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
    conv = read_h5()
    #a = np.genfromtxt('ratio.dat')
    #rat = a[:,1]
    #conv *= rat[:,None]
    bootstrap(conv,1000)

def read_h5():
    f = h5py.File('align_conversion.h5','r')
    conv = f['wave'][...]
    return conv

def bootstrap(conv,iterations):
    print 'bootstrap'
    boot_list = []
    for ii in range(iterations):
        single = []
        for jj in range(0,conv.shape[0]):
            idx = np.random.randint(0,conv.shape[0])
            single.append(conv[idx])
        single = np.array(single).mean(axis=0)
        boot_list.append(single)
    f = h5py.File('align_bootstrap.h5','w')
    f.create_dataset('boot',data = np.array(boot_list))
    f.close()

main()
