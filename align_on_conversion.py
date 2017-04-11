#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : align_on_conversion.py
Purpose : align traces on S1800P by drag-and-drop
Creation Date : 07-04-2017
Last Modified : Sun 09 Apr 2017 06:49:35 PM EDT
Created By : Samuel M. Haugland

==============================================================================
'''

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import pylab as p
from subprocess import call
from os import listdir
import h5py
import obspy
import seispy
from obspy.taup import TauPyModel
model = TauPyModel(model='prem50')


def main():
    st = read_stream()
    data_dict,gcarc_list = window_phase(st)
    fig,ax = plt.subplots(figsize=(10,15))
    for idx,keys in enumerate(data_dict):
        #e = tr.stats.npts/tr.stats.sampling_rate
        #l = tr.stats.station
        #t = np.linspace(0,850,num=tr.stats.npts)
        ax.plot(gcarc_list[idx]+data_dict[keys],color='k',label=keys,pickradius=50,picker=True)

    ax.axvline(len(data_dict[keys])/2.,color='b')
    ax.axvline(len(data_dict[keys])*(1./6),color='b')
    ax.axvline(len(data_dict[keys])*(2./6),color='b')
    ax.axvline(len(data_dict[keys])*(4./6),color='b')
    ax.axvline(len(data_dict[keys])*(5./6),color='b')
    dragh = DragHandler(figure=fig)
    plt.show()
    roll_list = zip(dragh.station,dragh.shift)
    #print roll_list
    data_dict = roll_data(roll_list,data_dict)
    write_h5(data_dict)

def write_h5(data_dict):
    array = np.array([data_dict[keys] for keys in data_dict])
    f = h5py.File('align_conversion.h5','w')
    f.create_dataset('wave',data=array)
    f.close()

def roll_data(roll_list,data_dict):
    for ii in roll_list:
        data_dict[ii[0]] = np.roll(data_dict[ii[0]],ii[1])
    return data_dict

def window_phase(st):
    data_dict = {}
    gcarc_list = []
    for idx,tr in enumerate(st):
        d = seispy.data.phase_window(st[idx],phase=['S1800P'],window=(-30,30))
        data_dict[tr.stats.station] = d.data
        gcarc_list.append(tr.stats.sac['gcarc'])
    return data_dict,gcarc_list

def read_stream():
    stz = obspy.read('/home/samhaug/work1/SP_brazil_data/2007-07-21-mw60-western-brazil-4/sparse_Z.pk')
    str = obspy.read('/home/samhaug/work1/SP_brazil_data/2007-07-21-mw60-western-brazil-4/sparse_R.pk')
    str.filter('highpass',freq=1./20)
    stz.filter('highpass',freq=1./20)
    stz = seispy.data.align_on_phase(stz,phase=['P'])
    for idx,tr in enumerate(str):
        d = seispy.data.phase_window(str[idx],['S'],window=(-30,30))
        e = np.abs(d.data).max()
        stz[idx].data *= 1./e
        tr.stats.location = tr.stats.sac['gcarc']
    stz.interpolate(50)
    stz.sort(['location'])
    stz = seispy.data.align_on_phase(stz)
    return stz

class DragHandler(object):
    """ A simple class to handle Drag n Drop.

    This is a simple example, which works for Text objects only.
    """
    def __init__(self, figure=None) :
        if figure is None : figure = p.gcf()
        # simple attibute to store the dragged text object
        self.dragged = None
        self.station = []
        self.shift = []

        # Connect events and callbacks
        figure.canvas.mpl_connect("pick_event", self.on_pick_event)
        figure.canvas.mpl_connect("button_release_event", self.on_release_event)

    def on_pick_event(self, event):

        self.dragged = event.artist
        self.xdata = event.artist.get_data()[0]
        self.ydata = event.artist.get_data()[1]
        self.pick_pos = event.mouseevent.xdata
        self.station.append(event.artist.get_label())
        return True

    def on_release_event(self, event):

        newx = event.xdata
        newy = np.roll(self.ydata,int(newx-self.pick_pos))
        self.dragged.set_data(self.xdata,newy)
        self.dragged = None
        p.draw()
        self.shift.append(int(newx-self.pick_pos))
        return True

main()




