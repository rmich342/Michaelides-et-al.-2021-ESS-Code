#!/usr/bin/env python
# coding: utf-8

# In[ ]:




import h5py
import numpy as np
import glob
import pandas as pd
import matplotlib.pyplot as plt
from pyproj import Proj, transform, Transformer

#simple function for reading in files
def read_h5(fname, vnames=[]):
        """ Simple HDF5 reader. """
        with h5py.File(fname, 'r') as f:
            return [f[v][:] for v in vnames]

#step1: get segment of track and interpolate over different distances
        
#files corresponding to example crossovers
f1 = 'processed_ATL06_20190514145800_07080303_003_01_gt1l_spot1_A.h5'
f2 = 'processed_ATL06_20190609022711_10970305_003_01_gt3r_spot6_D.h5'
xovers_40 = pd.read_csv('xovers_40_full_period_low_slopes_filtered.csv') #read in source file to get a sample crossover
i = xovers_40.index[(xovers_40['source_1'] == f1) & (xovers_40['source_2'] == f2)][0] #find crossover in dataframe
#load in crossover point
x_c = xovers_40['x'][i]
y_c = xovers_40['y'][i]

#read in data from source file of interest
lon1, lat1, h1, sigma_1 = read_h5(f1, ['lon', 'lat', 'h_elv', 's_elv'])

#project into UTM:
transformer = Transformer.from_crs("epsg:4326", "epsg:32606")
x1, y1 = transformer.transform(lat1, lon1)

#get distance from ntersection point:
distance_1 = np.sqrt((x1 - x_c)**2 + (y1 - y_c)**2)

#get points within 20m of intersection
lon_20_1f = x1[distance_1 <= 20]
lat_20_1f = y1[distance_1 <= 20]
h_20_1f = h1[distance_1 <= 20]
sigma_20_1f = sigma_1[distance_1 <= 20]
#get along-track distance from furthest point
along_track_20_1 = np.sqrt((lon_20_1f - lon_20_1f[0])**2 + (lat_20_1f - lat_20_1f[0])**2)
#get along-track distance of crossing point
c_at_20_1 = np.sqrt((x_c - lon_20_1f[0])**2 + (y_c - lat_20_1f[0])**2)
#calculate linear interpolation of along-track points
p1_20, v1_20 = np.polyfit(along_track_20_1, h_20_1f, 1, w = sigma_20_1f**(-1), cov= 'unscaled')
h1_20 = np.polyval(p1_20, along_track_20_1) #evaluate for all points
h1_20_c = np.polyval(p1_20, c_at_20_1) #evaluate at crossing location


#repeat for 40m
lon_40_1f = x1[distance_1 <= 40]
lat_40_1f = y1[distance_1 <= 40]
h_40_1f = h1[distance_1 <= 40]
sigma_40_1f = sigma_1[distance_1 <= 40]
#get along-track distance from furthest point
along_track_40_1 = np.sqrt((lon_40_1f - lon_40_1f[0])**2 + (lat_40_1f - lat_40_1f[0])**2)
#get along-track distance of crossing point
c_at_40_1 = np.sqrt((x_c - lon_40_1f[0])**2 + (y_c - lat_40_1f[0])**2)
#interpolate to crossing point
p1_40, v1_40 = np.polyfit(along_track_40_1, h_40_1f, 1, w = sigma_40_1f**(-1), cov= 'unscaled')
h1_40 = np.polyval(p1_40, along_track_40_1)
h1_40_c = np.polyval(p1_40, c_at_40_1)


#repeat for 60m
lon_60_1f = x1[distance_1 <= 60]
lat_60_1f = y1[distance_1 <= 60]
h_60_1f = h1[distance_1 <= 60]
sigma_60_1f = sigma_1[distance_1 <= 60]
#get along-track distance from furthest point
along_track_60_1 = np.sqrt((lon_60_1f - lon_60_1f[0])**2 + (lat_60_1f - lat_60_1f[0])**2)
#get along-track distance of crossing point
c_at_60_1 = np.sqrt((x_c - lon_60_1f[0])**2 + (y_c - lat_60_1f[0])**2)
#interpolate to crossing point
p1_60, v1_60 = np.polyfit(along_track_60_1, h_60_1f, 1, w = sigma_60_1f**(-1), cov= 'unscaled')
h1_60 = np.polyval(p1_60, along_track_60_1)
h1_60_c = np.polyval(p1_60, c_at_60_1)


#repeat for 80m interpolation
lon_80_1f = x1[distance_1 <= 80]
lat_80_1f = y1[distance_1 <= 80]
along_track_80_1 = np.sqrt((lon_80_1f - lon_80_1f[0])**2 + (lat_80_1f - lat_80_1f[0])**2)
c_at_80_1 = np.sqrt((x_c - lon_80_1f[0])**2 + (y_c - lat_80_1f[0])**2)
h_80_1f = h1[distance_1 <= 80]
sigma_80_1f = sigma_1[distance_1 <= 80]
p1_80, v1_80 = np.polyfit(along_track_80_1, h_80_1f, 1, w = sigma_80_1f**(-1), cov= 'unscaled')
h1_80 = np.polyval(p1_80, along_track_80_1)
h1_80_c = np.polyval(p1_80, c_at_80_1)

#repeat for 100m interpolation
lon_100_1f = x1[distance_1 <= 100]
lat_100_1f = y1[distance_1 <= 100]
along_track_100_1 = np.sqrt((lon_100_1f - lon_100_1f[0])**2 + (lat_100_1f - lat_100_1f[0])**2)
c_at_100_1 = np.sqrt((x_c - lon_100_1f[0])**2 + (y_c - lat_100_1f[0])**2)
h_100_1f = h1[distance_1 <= 100]
sigma_100_1f = sigma_1[distance_1 <= 100]
p1_100, v1_100 = np.polyfit(along_track_100_1, h_100_1f, 1, w = sigma_100_1f**(-1), cov= 'unscaled')
h1_100 = np.polyval(p1_100, along_track_100_1)
h1_100_c = np.polyval(p1_100, c_at_100_1)
usetex = True
fig = plt.figure(figsize= (15, 10))
ax1 = fig.add_subplot(2, 1, 1)

#step 2: plotting

#load in crossover files and filter for short-period crossovers
xovers_20= pd.read_csv('xovers_full_period_low_slopes_new_filtered.csv')
xovers_20_s = xovers_20[xovers_20['dt (days)'] <= 14]

xovers_40= pd.read_csv('xovers_40_full_period_low_slopes_filtered.csv')
xovers_40_s = xovers_40[xovers_40['dt (days)'] <= 14]

xovers_60= pd.read_csv('xovers_60_full_period_low_slopes_filtered.csv')
xovers_60_s = xovers_60[xovers_60['dt (days)'] <= 14]

xovers_80= pd.read_csv('xovers_80_full_period_low_slopes_filtered.csv')
xovers_80_s = xovers_80[xovers_80['dt (days)'] <= 14]

xovers_100= pd.read_csv('xovers_100_full_period_low_slopes_filtered.csv')
xovers_100_s = xovers_100[xovers_100['dt (days)'] <= 14]

#get array of standard deviation of short-period crossovers
sds_s = [np.std(xovers_20_s['dh'].values*100), np.std(xovers_40_s['dh'].values*100), np.std(xovers_60_s['dh'].values*100),                                            np.std(xovers_80_s['dh'].values*100), np.std(xovers_100_s['dh'].values*100)]


r = [20, 40, 60, 80, 100] #interpolation radii
#top plot: short-period standard deviation vs radius
plt.scatter(r, sds_s, c='blue', label = None)
plt.ylim([0, 25])
plt.xlabel('radius (m)', fontsize = 20)
plt.ylabel('short-period $\sigma$ (cm)',fontsize = 20)

ax1.tick_params(axis='both', which='major', labelsize=18)
ax1.tick_params(axis='y', colors='blue')
ax1.yaxis.label.set_color('blue')
ax1.spines['left'].set_color('blue')
#calculated median propogated uncertainty
rms = [np.median(xovers_20['sigma_dh'].values*100), np.median(xovers_40['sigma_dh'].values*100),               np.median(xovers_60['sigma_dh'].values*100), np.median(xovers_80['sigma_dh'].values*100),               np.median(xovers_100['sigma_dh'].values*100)]

plt.annotate('a)', (.01, .9), xycoords = 'axes fraction', fontsize = 18)
#add median propogated uncertainty plot
ax2 = ax1.twinx()
plt.scatter(r, rms, c = 'r', label = None)
plt.xlabel('radius (m)', fontsize = 20)
plt.ylabel('median uncertainty (cm)', fontsize = 20 )

ax2.tick_params(axis='both', which='major', labelsize=18)
ax2.tick_params(axis = 'y', colors = 'red')
plt.ylim([0, 25])
ax2.spines['right'].set_color('red')
ax2.yaxis.label.set_color('red')

ax3 = fig.add_subplot(2, 1, 2)
#plot each interpolation with x-axis centered around the crossover point
plt.plot(along_track_100_1 - c_at_100_1, h_100_1f,c = 'blue', marker = 'o', label = 'ATL06') #plot ATL06 tracks
plt.plot(along_track_40_1 - c_at_40_1, h1_40, c = 'orange', marker = 'o',label = 'linear interpolation, r = 40') #80m inteprolation
plt.plot(along_track_60_1 - c_at_60_1,h1_60, c = 'red', marker = 'o',label = 'linear interpolation, r = 60')
plt.plot(along_track_80_1 - c_at_80_1,h1_80, c = 'green', marker = 'o',label = 'linear interpolation, r = 80')
plt.plot(along_track_100_1 - c_at_100_1,h1_100, c = 'purple', marker = 'o',label = 'linear interpolation, r = 100')

#plot interpolated height at crossover location for each radius
plt.plot(0, h1_20_c,  marker = '*', c = 'blue', markersize =10, label = None)
plt.plot(0, h1_40_c,  marker = '*', c = 'orange',markersize =10, label = None)
plt.plot(0, h1_60_c, marker = '*', c = 'r',markersize =10,label = None)
plt.plot(0, h1_80_c,  marker = '*', c = 'green',markersize =10, label = None)
plt.plot(0, h1_100_c, marker = '*', c = 'purple',markersize =10,label = None)#interpolated crossing point(40)
plt.legend(prop={'size': 18})
plt.xlabel('along-track distance (m)', fontsize = 20)
plt.ylabel('elevation (m)', fontsize = 20)
plt.tick_params(axis='both', which='major', labelsize=18)
plt.xlim([-100, 100])
plt.xticks(np.arange(-100, 101, 20))
plt.annotate('b)', (.01, .9), xycoords = 'axes fraction', fontsize = 18)
plt.savefig('ess_fig3.pdf')
#fig.tight_layout(pad=3.0)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




