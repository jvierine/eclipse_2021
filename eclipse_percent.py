#!/usr/bin/env python
#
# Calculate the totality of the eclipse.
#
import ephem
import math
from operator import itemgetter
import numpy as n
import matplotlib.pyplot as plt
import scipy.io as sio

from mpl_toolkits.basemap import Basemap, cm


def check_non_zero(x):
    return x > 0

def intersection(r0,r1,d,n_s=100):
    A1=n.zeros([n_s,n_s])
    A2=n.zeros([n_s,n_s])
    I=n.zeros([n_s,n_s])
    x=n.linspace(-2.0*r0,2.0*r0,num=n_s)
    y=n.linspace(-2.0*r0,2.0*r0,num=n_s)
    xx,yy=n.meshgrid(x,y)
    A1[n.sqrt((xx+d)**2.0+yy**2.0) < r0]=1.0
    n_sun=n.sum(A1)
    A2[n.sqrt(xx**2.0+yy**2.0) < r1]=1.0
    S=A1+A2
    I[S>1]=1.0
    eclipse=n.sum(I)/n_sun
    return(eclipse)

def plot_map(p,lons,lats,t0,alt=0,lat_0=69,lon_0=16):
    fig = plt.figure(figsize=(8,8))
    plt.style.use('dark_background')
    # create polar stereographic Basemap instance.
    m = Basemap(projection='gnom',lon_0=lon_0,lat_0=lat_0, resolution='l',width=15.e6,height=15.e6)
    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    m.drawcountries()
    
    lon,lat=n.meshgrid(lons,lats)
    x, y = m(lon, lat) # compute map proj coordinates.
    # draw filled contours.
    cs = m.pcolormesh(x,y,1.0-p[0,0,:,:],vmin=0,vmax=1.0,cmap="inferno")
    # add colorbar.
    cbar = m.colorbar(cs,location='bottom',pad="5%")
    plt.title("Fraction of sunlight at %1.2f km\n%s (UTC)"%(alt/1e3,t0))
    fname="eclipse-%f.png"%(float(t0))
    plt.text(1e5,1e5,"University of Tromso",color="black")
    print("saving %s"%(fname))
    plt.savefig(fname)

# calculate the fraction that sun is eclipsed at given altitudes, latitude, and longitude
#
# returns eclipse fraction (0..1) and time (seconds after t0) 
def get_eclipse(t0,n_t=60*10,dt=60.0,alts=n.linspace(0,600e3,num=600),lats=[69.0],lons=[16.02]):
    # Location
    obs = ephem.Observer()
    n_alts=len(alts)
    n_lats=len(lats)
    n_lons=len(lons)
    
    p=n.zeros([n_t,n_alts,n_lats,n_lons])
    times=n.arange(n_t)*dt
    dts=[]
    for ti,t in enumerate(times):
#        print("Time %1.2f (s)"%(t))
        for ai,alt in enumerate(alts):
            for lai,lat in enumerate(lats):
                for loi,lon in enumerate(lons):
                    #obs.lon, obs.lat = '-1.268304', '51.753101'#'16.02', '78.15' # ESR
                    obs.lon, obs.lat = '%1.2f'%(lon), '%1.2f'%(lat) # ESR
                    obs.elevation=alt
                    obs.date= t0#(ephem.date(ephem.date(t0)+t*ephem.second))
                    sun, moon = ephem.Sun(), ephem.Moon()
                    
                    # Output list
                    results=[]
                    seps=[]
                    sun.compute(obs)
                    moon.compute(obs)
                    r_sun=(sun.size/2.0)/3600.0
                    r_moon=(moon.size/2.0)/3600.0
                    s=n.degrees(ephem.separation((sun.az, sun.alt), (moon.az, moon.alt)))
                    percent_eclipse=0.0

                    if s < (r_moon+r_sun):
#                        print("eclipsed")
                        if s < 1e-3:
                            percent_eclipse=1.0
                        else:
                            percent_eclipse=intersection(r_moon,r_sun,s,n_s=100)

                    if n.degrees(sun.alt) <= r_sun:
                        if n.degrees(sun.alt) <= -r_sun:
                            percent_eclipse=1.0
                        else:
                            percent_eclipse=1.0-((n.degrees(sun.alt)+r_sun)/(2.0*r_sun))*(1.0-percent_eclipse)

#                    print(percent_eclipse)
                    p[ti,ai,lai,loi]=percent_eclipse
        dts.append(obs.date)
    return(p,times,dts)

#(ephem.date(+t*ephem.second))
def create_map(t0=ephem.date((2015, 3, 20, 10, 10, 0)),lat_0=69,lon_0=16.,alt=150e3,nlats=200,nlons=500):
    n_t=1
    dt=10*60.0    
    lats=n.linspace(-70,90,num=nlats)
    lons=n.linspace(0,360,num=nlons)
    alts=[alt]
    
    p,t,dts=get_eclipse(t0,dt=dt,n_t=n_t,alts=alts,lats=lats,lons=lons)
    plot_map(p,lons,lats,t0=dts[0],alt=alts[0],lat_0=lat_0,lon_0=lon_0)

    
#create_map(t0=ephem.date((2015,3,20,10,0,0)),lat_0=69.58,lon_0=19.23)
#for i in range(20):
##    create_map(t0=ephem.date((2015,3,20,7,0,0))+ephem.second*900*i,lat_0=69.58,lon_0=19.23)
#    create_map(t0=ephem.date((2017,8,21,16,15,0))+ephem.second*900*i,lat_0=38.496077,lon_0=360-90.152751,nlats=300,nlons=600)

def create_matlab_file(t0,lat,lon,alts):
    dt=60.0
    n_t=60*24.0
    p,t,dts=get_eclipse(t0,dt=dt,n_t=n_t,alts=alts,lats=[lat],lons=[lon])
    plt.pcolormesh(t,alts/1e3,p[:,:,0,0].T,vmin=0,vmax=1.0)
    plt.xlabel("Time (s since midnight %s)"%(dts[0]))
    plt.ylabel("Range (km)")
    plt.colorbar()
    plt.show()
    fname="eclipse-%d-%1.2f-%1.2f.mat"%(int(t0),lat,lon)
    print("saving %s"%(fname))
    sio.savemat(fname,{"eclipse":p[:,:,0,0],"t":t,"t0":"%s"%(dts[0]),"dates":dts})

#create_matlab_file((2015, 3, 20, 0, 0, 0),lat=69.58,lon=19.23,alts=n.arange(600)*1e3)

if __name__ == "__main__":
#    p,t,dts=get_eclipse(t0,dt=dt,n_t=n_t,alts=alts,lats=lats,lons=lons)

    p,t,dts=get_eclipse(t0=ephem.date((2021, 6, 10, 10, 10, 0)),n_t=1,dt=60.0,alts=[300e3],lons=n.linspace(0,360,num=100),lats=n.linspace(0,90,num=90))

    print(p.shape)
    plt.pcolormesh(p[0,0,:,:])
    plt.colorbar()
    plt.show()
