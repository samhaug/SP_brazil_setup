from numpy import pi ,sin, cos
import obspy
import numpy as np
from geopy.distance import great_circle
import ses3drotation
import seispy

def make_R(theta, u):
    '''
    Make rotation matrix
    '''
    return [[cos(theta) + u[0]**2 * (1-cos(theta)),
             u[0] * u[1] * (1-cos(theta)) - u[2] * sin(theta),
             u[0] * u[2] * (1 - cos(theta)) + u[1] * sin(theta)],
            [u[0] * u[1] * (1-cos(theta)) - u[2] * sin(theta),
             cos(theta) + u[1]**2 * (1-cos(theta)),
             u[1] * u[2] * (1 - cos(theta)) + u[0] * sin(theta)],
            [u[0] * u[2] * (1-cos(theta)) - u[1] * sin(theta),
             u[1] * u[2] * (1-cos(theta)) - u[0] * sin(theta),
             cos(theta) + u[2]**2 * (1-cos(theta))]]

def rotate(pointToRotate, euler_pole, theta):
    '''
    Rotate point euler pole by angle theta. Right hand rule applies

    pointToRotate = [x,y,z] of point to rotate

    euler_pole = [x,y,z] of euler pole
    '''
    u= []
    squaredSum = 0

    for i,f in zip(euler_pole, [0,0,1]):
        u.append(f-i)
        squaredSum += (f-i) **2
    u = [i/squaredSum for i in u]
    r = make_R(theta, u)
    rotated = []

    for i in range(3):
        rotated.append((np.sum([r[j][i]*pointToRotate[j] for j in range(3)])))

    return rotated

def find_euler_pole(station_coord):
    '''
    station_coord = [x,y,z] of station
    '''
    z = [0,0,1]
    euler_pole = np.cross(station_coord,z)

    return euler_pole

def convert_cart(lat,lon):
    '''
    convert lat/lon coordinates to normalized [x,y,z]
    '''
    x = cos(np.radians(lat))*cos(np.radians(lon))
    y = cos(np.radians(lat))*sin(np.radians(lon))
    z = sin(np.radians(lat))

    return [x,y,z]

def stream_2_coord(st):
    '''
    Make list of coordinates for each station
    '''
    stat_coord = []
    for tr in st:
        stat_coord.append(convert_cart(tr.stats.sac['stla'],
                         tr.stats.sac['stlo']))
    event_coord = convert_cart(st[0].stats.sac['evla'],
                               st[0].stats.sac['evlo'])
    return stat_coord, event_coord

def rotate_stream(st):
    '''
    '''
    stat_coord, event_coord = stream_2_coord(st)

    ep = find_euler_pole(event_coord)
    gc = great_circle((90,0),(st[0].stats.sac['evla'],
                  st[0].stats.sac['evlo'])).km/111.195
    rotate_list = []
    for ii in stat_coord:
        rotate_list.append(rotate(ii,ep,gc))
    return rotate_list,stat_coord,event_coord,ep

def cart_to_lonlat(coord):
    '''
    convert xyz to lat
    '''
    lat = np.degrees(np.arcsin(coord[-1]))
    lon = np.degrees(np.arccos(coord[0]/np.arccos(np.radians(lat))))
    return lat,lon



def full_rotate_stream(st):
    '''
    all in one
    '''

    strot = st.copy()
    rotphi = 90-st[0].stats.sac['evla']
    evvec = convert_cart(st[0].stats.sac['evla'],st[0].stats.sac['evlo'])
    ep = find_euler_pole(evvec)

    for tr in strot:
        colat = 90-tr.stats.sac['stla']
        lon = tr.stats.sac['stlo']
        colatnew, lonnew = ses3drotation.rotate_coordinates(ep,rotphi,colat,lon)
        tr.stats.sac['stla'] = 90-colatnew
        tr.stats.sac['stlo'] = lonnew
        colat = 90-tr.stats.sac['evla']
        lon = tr.stats.sac['evlo']
        colatnew, lonnew = ses3drotation.rotate_coordinates(ep,rotphi,colat,lon)
        tr.stats.sac['evla'] = 90-colatnew
        tr.stats.sac['evlo'] = lonnew
    return strot

'''
rotate_list, stat_coord,event_coord,ep = rotate_stream(st)
latlon = []
for ii in stat_coord:
    latlon.append(cart_to_lonlat(ii))
print latlon
seispy.convert.axisem_stations(st)
'''
st = obspy.read('super_sparse.pk')
seispy.mapplot.source_reciever_plot(st)
new = full_rotate_stream(st)
seispy.convert.axisem_stations(new)

#stat_coords = convert_cart(-8.08,-71)
#ep = find_euler_pole(stat_coords)
#theta = np.radians(45.)

#out = rotate(point,ep,theta)






