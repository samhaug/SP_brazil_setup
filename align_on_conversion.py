#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : align_on_conversion.py
Purpose : align traces on S1800P by drag-and-drop
Creation Date : 07-04-2017
Last Modified : Sat 08 Apr 2017 05:39:15 PM EDT
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
    fig,ax = plt.subplots(figsize=(10,15))
    for idx,tr in enumerate(st):
        e = tr.stats.npts/tr.stats.sampling_rate
        l = tr.stats.station
        #t = np.linspace(0,850,num=tr.stats.npts)
        gcarc = tr.stats.sac['gcarc']
        ax.plot(gcarc+tr.data,color='k',label=l,pickradius=50,picker=True)

    dragh = DragHandler(figure=fig)
    plt.show()
    print zip(dragh.station,dragh.shift)

def read_stream():
    st = obspy.read('/home/samhaug/work1/SP_brazil_data/2007-07-21-mw60-western-brazil-4/sparse_Z.pk')
    st.normalize()
    st = seispy.data.align_on_phase(st,phase=['P'])
    for idx,tr in enumerate(st):
        st[idx] = seispy.data.phase_window(st[idx],['P'],window=(-50,800))
        tr.stats.location = tr.stats.sac['gcarc']
    st.interpolate(50)
    st.sort(['location'])
    return st

class DragHandler(object):
    """ A simple class to handle Drag n Drop.

    This is a simple example, which works for Text objects only.
    """
    def __init__(self, figure=None) :
        """ Create a new drag handler and connect it to the figure's event system.
        If the figure handler is not given, the current figure is used instead
        """
        if figure is None : figure = p.gcf()
        # simple attibute to store the dragged text object
        self.dragged = None
        self.station = []
        self.shift = []

        # Connect events and callbacks
        figure.canvas.mpl_connect("pick_event", self.on_pick_event)
        figure.canvas.mpl_connect("button_release_event", self.on_release_event)

    def on_pick_event(self, event):
        " Store which text object was picked and were the pick event occurs."

        self.dragged = event.artist
        self.xdata = event.artist.get_data()[0]
        self.ydata = event.artist.get_data()[1]
        self.pick_pos = event.mouseevent.xdata
        self.station.append(event.artist.get_label())
        return True

    def on_release_event(self, event):
        " Update text position and redraw"

        #if self.dragged is not None :
        newx = event.xdata
        newy = np.roll(self.ydata,int(newx-self.pick_pos))
        self.dragged.set_data(self.xdata,newy)
        self.dragged = None
        p.draw()
        self.shift.append(int(newx-self.pick_pos))
        return True


main()




