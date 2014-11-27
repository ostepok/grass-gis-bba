#!/usr/bin/env python

#%module
#% description: Bundle block adjustment prototype.
#%end

#%option G_OPT_M_DIR
#% key: bingo_data_dir
#% description: Path to the file of Bingo-F protocol.
#% required: yes
#%end

#%option G_OPT_M_DIR
#% key: protocol_dir
#% description: Directory where the protocols will be generated.
#% required: yes
#%end

#%flag
#% key: f
#% description: Adjust as free network.
#%end


import os
import sys
import numpy as np

import grass.script as grass


from plots import PlotScene
from readpar import InputData
from parsers import BingoParser
from helmert import HelmertTransform, Test
from relative_or import RelativeOrientation
from bba import RunBundleBlockAdjutment, RotMatBLUH, GetBLUHAngles


def main():
    
    # read data from Bingo-F protocols
    in_dt = InputData(BingoParser(options["bingo_data_dir"]))

    phs = in_dt.GetPhotos()
    pts = in_dt.GetPoints()
    cam_m, distor = in_dt.GetCameras().values()[0].GetParams()
    

    # compute relative orientations and merged them into common system
    ros = RelativeOrientation(in_dt)

    np.set_printoptions(precision=3, suppress=True)

    #PlotScene(pts, phs)
    #PlotRelOrs(pts, phs)

    # performs helmert transformation into world coordinate system 
    HelmertTransform(pts, phs)

    # run bundle block adjustment
    ph_errs, pt_errs = RunBundleBlockAdjutment(in_dt, options['protocol_dir'], flags['f'])   

    PlotScene(pts, phs)
    return 

if __name__ == "__main__":
    options, flags = grass.parser()
    sys.exit(main())
