#!/home/samhaug/anaconda2/bin/python

import numpy as np
import obspy
import seispy

st = obspy.read('st_Z.pk')
newst = obspy.core.Stream()

for idx in range(0,6):
    a = st[idx::6].copy()
    a = seispy.data.slant(a,0.5)
    seispy.plot.section(a,picker=True)
    for tr in a:
        newst.append(tr)
    print len(newst)
newst.write('sparse.pk',format='PICKLE')

