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
    input_file = sys.argv[1]
    radar = xrain(input_file)
    if radar.par[1]!="12":
        print("Error: this input file is not RZH0 (horizontal reflectivity).")
        exit()
    x = radar.x
    y = radar.y
    zh = radar.ppi

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
    p = ax.pcolormesh(x, y, zh, vmin=0.0, vmax=50.0, cmap=custom_color)
    pp = fig.colorbar(p, ax=ax, orientation="vertical")
    pp.set_clim(0.0,50.0)
    pp.set_label("dBZ", fontname="Arial", fontsize=10)
    plt.xlabel("x [km]")
    plt.ylabel("y [km]")
    ax.set_xlim(-60.0, 60.0)
    ax.set_ylim(-60.0, 60.0)
    ax.set_aspect('equal', 'datalim')
    plt.show()

    
