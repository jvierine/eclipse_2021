import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import eclipse_percent as ep
import ephem
import numpy as n
import stuffr

def plot_eclipse(t0=ephem.date((2021, 6, 10, 10, 40, 0))):
    #
    lons=n.linspace(0,360,num=100)
    lats=n.linspace(-90,90,num=90)
    p,t,dts=ep.get_eclipse(t0=t0,n_t=1,dt=60.0,alts=[300e3],lons=lons,lats=lats)
    
    fig=plt.figure(figsize=(10,4))
    #prj=ccrs.Orthographic(central_longitude=0,central_latitude=90)
    #prj=ccrs.Orthographic(central_longitude=0,central_latitude=90)
    #prj=ccrs.PlateCarree()#Geodetic()#RotatedPole()
    prj=ccrs.RotatedPole(pole_latitude=0)
#    fig=plt.figure(figsize=(10,6))
#    ax = plt.axes(projection=prj)
    ax = fig.add_subplot(1,1,1,projection=prj)
#    ax.set_extent([0, 360,0, 90])#,ccrs.Geodetic())
    ax.set_global()
    plt.pcolormesh(lons,lats,1.0-p[0,0,:,:],transform=ccrs.PlateCarree(),cmap="nipy_spectral")
    cb=plt.colorbar()
#    cb.set_label("Fraction of solar disc blocked")
    cb.set_label("Fraction of solar disc visible")    
    
    ax.coastlines(color="gray")
    ax.gridlines()

    print(ax.get_xlim())
    print(ax.get_ylim())    
    #ax.set_ylim([0,90])
    
    haarp_p=[62.391667, -145.15]
    ski_p=[69.34844744991165, 20.36342820500161]
    
    print("plotting red line")
    plt.plot([haarp_p[1],ski_p[1]],[haarp_p[0],ski_p[0]],color="red",transform=ccrs.Geodetic())
    
    unix_t0=ephem.date((1970, 1, 1, 0, 0, 0))
    unix_t=(t0-unix_t0)/ephem.second
    utc_str=stuffr.unix2datestr(unix_t)
    plt.title("%s UTC"%(utc_str))
    


if __name__ == "__main__":
    unix_t0=ephem.date((1970, 1, 1, 0, 0, 0))
    dt=900.0
    for i in range(24*2*2):
        t0=ephem.date((2021, 6, 10, 0, 0, 0))+ephem.second*i*dt
        unix_t=(t0-unix_t0)/ephem.second
        print(unix_t)
        print(stuffr.unix2datestr(unix_t))
        plot_eclipse(t0)
        plt.savefig("eclipse_%02d.png"%(i),dpi=200)
        plt.clf()
        plt.close()
