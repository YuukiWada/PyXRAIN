#!/usr/bin/env python
from PyXRAIN import composite
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import geopandas as gpd
import numpy as np
import math
import sys
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

def color():
    cmap = mpl.colors.ListedColormap(["#f5f5ff", "#b9ddff", "#64a8ff", "#4971ff", "#fff866", "#ffb35b", "#ff6653"])
    cmap.set_under("#9e9e9e")
    cmap.set_over("#c9518f")
    bounds = [0.1, 1.0, 5.0, 10.0, 20.0, 30.0, 50.0, 80.0]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    return cmap, norm

if __name__ == "__main__":
    input_file = sys.argv[1]
    lon_range = [134.0,144.0]
    lat_range = [31.0,37.0] 

    radar = composite(input_file)
    x,y = radar.mesh
    rain = radar.comp
    edge = radar.edge
    date = radar.date
    time = radar.time
    cmap, norm = color()

    x_index = [radar.lon_index(lon_range[0]), radar.lon_index(lon_range[1])]
    y_index = [radar.lat_index(lat_range[0]), radar.lat_index(lat_range[1])]    
    x = x[y_index[0]:y_index[1]+1,x_index[0]:x_index[1]+1]
    y = y[y_index[0]:y_index[1]+1,x_index[0]:x_index[1]+1]
    rain = rain[y_index[0]:y_index[1],x_index[0]:x_index[1]]
    
    # plotting
    fig = plt.figure(figsize=(8,5))
    df = gpd.read_file('./japan.geojson')
    ax = fig.add_subplot(111)
    p = ax.pcolormesh(x, y, rain, cmap=cmap, norm=norm)
    pp = fig.colorbar(p, ax=ax, orientation="vertical", extend="both")
    pp.set_label("mm/h", fontname="Arial", fontsize=10)
    plt.xlabel("Longitude [deg]")
    plt.ylabel("Latitude [deg]")
    plt.ylim(lat_range[0], lat_range[1])
    plt.xlim(lon_range[0], lon_range[1])
    plt.title(date+" "+time)
    df.plot(ax=ax, color="none", edgecolor="black", facecolor="white", linewidth=0.5)
    plt.show()
