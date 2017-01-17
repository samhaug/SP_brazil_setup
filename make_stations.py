#!/usr/bin/env python


from numpy import arange
stat_file = open('STATIONS','w')

lat = arange(0,90,1)

for ii in lat:
    stat_file.write('STAT_0_{} II {} 0.0 0.0 0.0\n'.format(str(ii),str(ii)))
    stat_file.write('STAT_90_{} II {} 90.0 0.0 0.0\n'.format(str(ii),str(ii)))
    stat_file.write('STAT_45_{} II {} 45.0 0.0 0.0\n'.format(str(ii),str(ii)))
stat_file.close()
