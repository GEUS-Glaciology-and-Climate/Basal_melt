# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 15:54:11 2021

@author: nbk
"""

#Make figures 1D, 1E, 1F and 2A from Karlsson et al., 2021, Nature Comms.

#Import necessary python packages

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as mcols
import netCDF4 as nc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#load file of heat sources - can be downloaded from  https://doi.org/10.22008/FK2/PLNUEO as heatsources_Karlssonetal2021.nc
filename='heatsources.nc'
ds=nc.Dataset(filename)

#Coordinate system information in meta-data. Projection polar_stereographic, epsg:3413
x = ds['x'][:]
y = ds['y'][:]

#Geothermal flux - mean of Fox Maule et al. (2009), Martos et al. (2018) and Shapiro and Ritzwoller (2004)
gf_mean = ds['gf_mean'][:]

#Friction heat from Elmer/Ice, see Gillet-Chaulet et al., 2012
fric_heat = ds['fric_heat'][:]

#Viscous heat dissipation from surface melt water, see Mankoff and Tulaczyk 2017
vhd_heat = ds['srfwater_heat'][:]

#
#****Calculate basal melt rates in m/year*****
#
rho = 918 # density of ice
L = 334*1e3 #latent heat of ice J/kg

basalmelt_gf = (gf_mean)/(rho*L)*(365.25*24*60*60)
basalmelt_fric = fric_heat/(rho*L)*(60*60*24*365.35)
basalmelt_vhd = vhd_heat/(rho*L)*(60*60*24*365.35)

#Apply mask of basal conditions from MacGregor et al., 2016
basalmask = ds['basalstate'][:]

basalmelt_gf=basalmelt_gf*basalmask
basalmelt_fric=basalmelt_fric*basalmask
basalmelt_vhd=basalmelt_vhd*basalmask

#
#****Calculate total basal melt rates in m/year*****
#

totalbasalmelt=basalmelt_gf+basalmelt_fric+basalmelt_vhd
totalbasalmelt[totalbasalmelt==0]=np.nan

#****Basal melt terms for each region, see Table 1 in manuscript
gf_reg= [0.51, 0.71, 1.27, 0.44, 0.56, 0.66, 1.17] #geothermal flux melt per region
fmelt_reg=[1.21, 2.4 , 1.02, 0.64, 2.05, 2.22, 1.33] #friction melt per region
vhd_reg=[0.54, 0.74, 0.47, 0.4 , 0.84, 0.83, 1.36] #surface water melt melt per region
tot_reg=fmelt_reg+vhd_reg+gf_reg # total melt per region
tottm=np.sum(tot_reg)


########## MAKE PLOTS ################

fig = plt.figure(1)
fig.set_size_inches(12.5, 8.5)

# color scale limits
bmmin=1e-3;
bmmax=1e0;
#colormap
cmap2 = plt.get_cmap('PuBu')
lincol=[cmap2(0.66)]

ax1 = plt.subplot(1,3,1)   
CS1=ax1.pcolormesh(x,y,basalmelt_gf,norm=mcols.LogNorm(bmmin,bmmax),cmap=cmap2,shading='auto')
CS1b=ax1.contour(x,y,basalmelt_gf,levels=[0],colors=lincol,linewidths=1)
ax1.axis('square')
ax1.set_title('Geothermal')
h1,_ = CS1b.legend_elements()
ax1.legend([h1[0]], ['0 m per year'],loc='upper right',fontsize=13)
ax1.set_xlim(-6.5e5,8.65e5)
ax1.set_xlabel('X Easting (m)')
ax1.set_ylabel('Y Northing (m)')

ax2 = plt.subplot(1,3,2)   
CS2 = ax2.pcolormesh(x,y,basalmelt_fric,norm=mcols.LogNorm(bmmin,bmmax),cmap=cmap2,shading='auto')
CS2b = ax2.contour(x,y,basalmelt_fric,levels=[0.01],colors='C0',linewidths=1)
ax2.axis('square')
ax2.set_title('Friction')
lines = [ CS2b.collections[0]]
labels = ['0.01 m per year']
ax2.legend(lines,labels,loc='upper right',fontsize=13)
ax2.set_xlim(-6.5e5,8.65e5)
ax2.set_xlabel('X Easting (m)')

ax3 = plt.subplot(1,3,3)   
CS3=ax3.pcolormesh(x,y,basalmelt_vhd,norm=mcols.LogNorm(bmmin,bmmax),cmap=cmap2,shading='auto')
ax3.axis('square')
ax3.set_title('Surface melt water')
ax3.set_xlim(-6.5e5,8.65e5)
ax3.set_xlabel('X Easting (m)')

pos1 = ax1.get_position()
pos3 = ax3.get_position()

colorbar_axes = plt.gcf().add_axes([pos1.x0, 0.13, pos3.x1-pos1.x0, 0.04]) #x0, y0, width, height
colorbar = plt.colorbar(CS2, colorbar_axes, orientation='horizontal')
colorbar.set_label(r'Basal melt rate (m per year)',fontsize=20)
colorbar.ax.tick_params(labelsize=16)

#*****total basal melt

fig = plt.figure(2)
fig.set_size_inches(12.5, 8.5)

cmap3 = plt.get_cmap('Reds')

CSt=plt.pcolormesh(x,y,totalbasalmelt,cmap=cmap3,vmin=0, vmax=0.25,shading='auto')
plt.axis('square')
plt.xlim(-6.5e5,8.65e5)
plt.title('Total basal melt')
plt.xlabel('X Easting (m)')
plt.ylabel('Y Northing (m)')

colorbar = plt.colorbar(CSt)
colorbar.set_label(r'Total basal melt (m per year)',fontsize=20)
colorbar.ax.tick_params(labelsize=16)

#define pie drawing
#see https://stackoverflow.com/questions/45266955/adding-pie-chart-at-given-coordinates-to-cartopy-projection
colors = ['C0','black','grey','yellow','magenta','purple']
def plot_pie_inset(data,ilon,ilat,ax,width):
    ax_sub= inset_axes(ax, width=width, height=width, loc=10, 
                       bbox_to_anchor=(ilon, ilat),
                       bbox_transform=ax.transData,
                       borderpad=0)
    wedges,texts= ax_sub.pie(data,colors=colors)
    ax_sub.set_aspect("equal")   

resizer=2

#Coordinates for pies
xpsout=[669662,11408,566135,136736,-147229,154543,-158850]
ypsout=[-2190297,-2178831,-1514127,-1258608,-1587010,-2948794,-3030985]
#draw pies for each region
for k in range(np.shape(xpsout)[0]):
    piesize=np.sqrt((fmelt_reg[k]+gf_reg[k]+vhd_reg[k])/np.pi)/resizer
    plot_pie_inset([fmelt_reg[k]/tot_reg[k],gf_reg[k]/tot_reg[k],vhd_reg[k]/tot_reg[k]],xpsout[k],ypsout[k],fig.axes[0],piesize)

#draw total pie chart
piesize=np.sqrt(tottm/np.pi)/resizer
plot_pie_inset([np.sum(fmelt_reg)/tottm,np.sum(gf_reg)/tottm,np.sum(vhd_reg)/tottm],574306,-2856985,fig.axes[0],piesize)

