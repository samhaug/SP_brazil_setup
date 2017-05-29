#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : align_on_conversion.py
Purpose : align traces on S1800P by drag-and-drop
Creation Date : 07-04-2017
Last Modified : Mon 29 May 2017 10:15:33 AM EDT
Created By : Samuel M. Haugland

==============================================================================
'''

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import pylab as p
import h5py
import obspy
import seispy
from obspy.taup import TauPyModel
model = TauPyModel(model='prem50')

def main():
    stz,str = read_stream()
    stz = window_phase(stz)
    fig,ax = plt.subplots(figsize=(8,10))
    for idx,tr in enumerate(stz):
        t = np.linspace(-30,30,num=tr.stats.npts)
        ax.plot(t,tr.stats.sac['gcarc']+(tr.data/tr.data.max()),color='k',
                pickradius=20,label=tr.stats.station,picker=True,lw=3)
    ax.axvline(0)
    ax.axvline(10)
    ax.axvline(-10)
    ax.set_xlim(-30,30)
    dragh = DragHandler(stz,figure=fig)
    plt.show()
    stz = roll_data(stz,dragh.station,dragh.shift)
    stz.write('Z_align.pk',format='PICKLE')
    str.write('R.pk',format='PICKLE')

def roll_data(st,stat_list,shift_list):
    for idx,tr in enumerate(st):
        if tr.stats.station in stat_list:
            i = stat_list.index(tr.stats.station)
            st[idx].data = np.roll(tr.data,shift_list[i])
    return st

def window_phase(st):
    data_dict = {}
    gcarc_list = []
    for idx,tr in enumerate(st):
        st[idx] = seispy.data.phase_window(tr,phase=['S1800P'],window=(-30,30))
    return st

def read_stream():
    stz = obspy.read('/home/samhaug/work1/SP_brazil_data/2007-07-21-mw60-western-brazil-4/sparse_Z.pk')
    str = obspy.read('/home/samhaug/work1/SP_brazil_data/2007-07-21-mw60-western-brazil-4/sparse_R.pk')
    str.interpolate(50)
    stz.interpolate(50)
    str.filter('highpass',freq=1./50)
    stz.filter('highpass',freq=1./50)
    stz = seispy.data.align_on_phase(stz,phase=['P'])
    for idx,tr in enumerate(str):
    #    d = seispy.data.phase_window(str[idx],['S'],window=(-30,-30))
    #    e = np.abs(d.data).max()
    #    stz[idx].data *= 1./e
        tr.stats.location = tr.stats.sac['gcarc']
        stz[idx].stats.location = stz[idx].stats.sac['gcarc']
    stz.sort(['location'])
    stz = seispy.data.align_on_phase(stz,phase=['P'])
    return stz,str

class DragHandler(object):
    """ A simple class to handle Drag n Drop.

    This is a simple example, which works for Text objects only.
    """
    def __init__(self,st,figure=None) :
        if figure is None : figure = p.gcf()
        # simple attibute to store the dragged text object
        self.dragged = None
        self.station = []
        self.shift = []
        self.srate = st[0].stats.sampling_rate

        # Connect events and callbacks
        figure.canvas.mpl_connect("pick_event", self.on_pick_event)
        figure.canvas.mpl_connect("button_release_event", self.on_release_event)

    def on_pick_event(self, event):
        #print type(event.artist.get_label())

        self.dragged = event.artist
        self.xdata = event.artist.get_data()[0]
        self.ydata = event.artist.get_data()[1]
        self.pick_pos = event.mouseevent.xdata
        self.station.append(event.artist.get_label())
        return True

    def on_release_event(self, event):

        newx = event.xdata
        try:
            newy = np.roll(self.ydata,int((newx-self.pick_pos)*self.srate))
            self.dragged.set_data(self.xdata,newy)
            self.dragged = None
            p.draw()
            self.shift.append(int((newx-self.pick_pos)*self.srate))
            return True
        except (TypeError,AttributeError) as e:
            return True

main()

