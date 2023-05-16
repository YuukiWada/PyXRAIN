#!/usr/bin/env python
from PyXRAIN import composite
import numpy as np
import math
import sys
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_folder = sys.argv[2]
    output_file = output_folder+"/"+input_file.split("/")[-1]+".csv"

    radar = composite(input_file)
    x,y = radar.mesh
    rain = radar.comp
    edge = radar.edge
    
    print("Data range")
    print("N"+str(np.round(radar.edge[0][0],decimals=3))+"-"+str(np.round(radar.edge[0][1],decimals=3))+"deg ("+str(rain.shape[1])+")")
    print("E"+str(np.round(radar.edge[1][0],decimals=3))+"-"+str(np.round(radar.edge[1][1],decimals=3))+"deg ("+str(rain.shape[0])+")")
    
    np.savetxt(output_file, rain, delimiter=",", fmt="%.2f")

