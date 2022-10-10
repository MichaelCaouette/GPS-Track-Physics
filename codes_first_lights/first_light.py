# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 15:59:16 2022

goal:
    Load and do some basic stuff with the gps data
    
A lot of inspiration from 
https://www.youtube.com/watch?v=9Q8nEA_0ccg
Florian Wilhelm - Handling GPS Data with Python

And his github repo:
    https://github.com/FlorianWilhelm/gps_data_with_python

@author: client
"""


import gpxpy
import numpy as np
import matplotlib.pyplot as plt

# path = "C://Users/client//Desktop//GPS track physics"
path = "C://Users//mcaoue2\Desktop\GPS track physics 2022_09_06"
filename = "TRK4.GPX"
#filename = "WPTS.GPX" 

filepath = path + '/' + filename

# Get the t, x, y, z (AKA motion!)
with open(filepath) as fh:
    gpx_file = gpxpy.parse(fh)
    
print("File has {} track(s).".format(len(gpx_file.tracks)))   
print("Track has {} segment(s).".format(len(gpx_file.tracks[0].segments))) 
# So we take the only track and segment. 
segment = gpx_file.tracks[0].segments[0]
print('The segment has {} points'.format(len(segment.points)))

# The structure is that each point has the t, x, y, z information.
# If we want an array for each degree of freedom, we have to restructure it

t_init_dt = segment.points[0].time
t_init_sec = t_init_dt.timestamp()

## Pre-allocation
#ele, lat, lon, t = np.zeros((4,len(segment.points) ))
#for i,p, in enumerate(segment.points):
#    ele[i] = p.elevation
#    lat[i] = p.latitude
#    lon[i] = p.longitude
#    t[i]   = p.time.timestamp() - t_init_sec

# Method to remove outliers
ele_max = 5000 # Maximum elevation
ele, lat, lon, t = [], [], [], []
for p in segment.points:
    if p.elevation <= ele_max:
        ele.append(p.elevation)
        lat.append(p.latitude)
        lon.append(p.longitude)
        t.append(p.time.timestamp() - t_init_sec)
ele = np.array(ele)
lat = np.array(lat)
lon = np.array(lon)
t   = np.array(t  )
        
# Check this simply
plt.figure()
plt.subplot(211)
plt.title(filename + '  Starting at ' + str(t_init_dt))
plt.plot(lat, lon)
plt.plot(lat[0], lon[0], 'or')
plt.subplot(212)
plt.plot(t/60 ,ele)
plt.xlabel('Time elapsed (min)')
    
# 3D trajectory
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(lat, lon, ele)
plt.title(filename + '  Starting at ' + str(t_init_dt))















    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

















