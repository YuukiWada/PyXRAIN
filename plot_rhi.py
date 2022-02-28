#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PyXRAIN import xrain
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import math
import sys
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


if __name__ == "__main__":
    # parameter settings
    azimuth = float(sys.argv[1])
    input_list = sys.argv[2]
    
    if (azimuth<0) or (azimuth>=360):
        print("Error: Azimuth angle should be between 0 and 360 degree.")
        exit()

    az_index = round((azimuth-0.6)/1.2)
    radar = list()
    elv_original = list()
    f = open(input_list, "r")
    input_files = f.readlines()

    if len(input_files)!=12:
        print("Error: The number of input files should be 12.")
        exit()

    for i in range(len(input_files)):
        radar.append(xrain(input_files[i].replace("\n",""),False))
        elv_original.append(radar[i].elv)
    f.close()

    if radar[0].par[1]!="12":
        print("Error: this input file is not RZH0 (horizontal reflectivity).")
        exit()

    elv_order = sorted(list(set(elv_original)))
    
    if len(elv_order)!=12:
        print("Error: Provided file sets have a duplicated zimuth angle.")
        exit()

    elv_border=list()
    elv_border.append(round(elv_order[0]-(elv_order[1]-elv_order[0])/2.0,3))
    for i in range(11):
        elv_border.append(round((elv_order[i+1]+elv_order[i])/2.0,3))
    elv_border.append(round(elv_order[11]+(elv_order[11]-elv_order[10])/2.0,3))
    
    rng_num = radar[0].range_num
    rng_min = radar[0].range_min
    rng_step = radar[0].range_step

    elv = np.tile(np.array(elv_border),(rng_num+1,1)).transpose()
    r = np.tile(rng_min+np.arange(rng_num+1)*rng_step/1000.0,(len(elv_border),1))
    z = r*np.tan(elv*math.pi/180.0) 
    rhi = np.zeros((12, rng_num))
    for i in range(12):
        index = elv_original.index(elv_order[i])
        rhi[i,:] = radar[index].ppi[az_index,:]
    
    # color settings
    colors = ['#FFFFFF', '#A0D2FF', '#218CFF', '#0041FF', '#FAF500', '#FF9900', '#FF2800', '#B40068']
    values = range(len(colors))
    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append((v/vmax,c))
    custom_color = LinearSegmentedColormap.from_list('custom_cmap', color_list)

    # plotting
    fig = plt.figure()
    ax = fig.add_subplot(111)
    p = ax.pcolormesh(r, z, rhi, vmin=0.0, vmax=50.0, cmap=custom_color)
    pp = fig.colorbar(p, ax=ax, orientation="vertical")
    pp.set_clim(0.0,50.0)
    pp.set_label("dBZ", fontname="Arial", fontsize=10)
    plt.xlabel("Radius [km]")
    plt.ylabel("Height [km]")
    ax.set_xlim(0.0, 60.0)
    ax.set_ylim(0.0, 6.0)
    #ax.set_aspect('equal')
    plt.show()

