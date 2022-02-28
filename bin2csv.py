#!/usr/bin/env python
# -*- coding: utf-8 -*- 
from PyXRAIN import xrain
import numpy as np
import sys

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_folder = sys.argv[2]
    output_file = output_folder+"/"+input_file.split("/")[-1]+".csv"
    radar = xrain(input_file)
    data = radar.ppi
    np.savetxt(output_file, data, delimiter=",", fmt="%.2f")
