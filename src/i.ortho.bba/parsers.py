import os
from math import pi
import numpy as np

from bba import RotMatBLUH
from readpar import Point, Camera, Photo

from relative_or import undistortPoints

# boring functions parsing input for bingo

class BingoParser:
    def __init__(self, path):   
        self.path = path
    
    def Parse(self):

        self.cams, self.gcps = self._readCameraFile()
        self.phs = self._readPhotosFile()
        self.pts = self.gcps.copy()
        self._readPhotoPointsFile(self.pts, self.phs)

        self.res_eo = self._readResult(self.phs, self.pts)

        self.cams.values()[0].SetChipSize((11,8))

    def GetData(self):
        return self.cams, self.phs, self.pts, self.gcps

    def GetResult(self):
        return self.res_eo

    def _readResult(self, phs, pts):
        f = open(self.path + '/old-itera.dat')

        while 1:
            line = f.readline()
            if not line:
                break
            if line[0] == "*" or not line.strip():
                continue

            lin=line.rstrip("")
            l = lin.split()
            if not l:
                continue
            if "ORIA" in l[0]:

                ph_id = int(l[1])

                eo = map(float, l[2:-1])
                eo[3:] = map(lambda d : d / 200 * pi, eo[3:])

                r =  RotMatBLUH(*eo[3:])
                t = np.array([eo[:3]])


                eo = np.hstack((r, t.T))

                phs[ph_id].SetResultExtOr(eo=eo)

            elif "CORD" in l[0]:
                pt_id = int(l[1])
                c = map(float, l[2:5])


                pts[pt_id].SetResultCoords(c)

        return phs

    def _readCameraFile(self):
        f = open(self.path + '/geoin.dat')

        cams = {}
        gcps = {}

        while 1:
            line = f.readline()
            if not line:
                break
           
            lin=line.strip()
            
            l = lin.split()
            if not l:
                continue
           
            if l[0] == '*' or l[0] == 'C':
                continue
            #TODO implemented for just one camera 
            if 'CAPA' in l[0]:
                fl = float(l[3])
                x0 = float(l[4])
                y0 = float(l[5])

            if 'ECCA' in l[0]:
                lge = float(l[2])
                lgn = float(l[3])
                lgh = float(l[4])

                distor = np.array([0, 0, 0, 0], dtype=float)
                cam_m = np.array([[fl, 0, x0],
                                  [0, fl, y0],
                                  [0, 0, 1]])

                self.cam = cams[l[1]] = Camera(l[0], cam_m, distor)

            if 'CONT' in l[0] or 'CHCK' in l[0]:
                control = False
                if 'CHCK' in l[0]:
                    control = True

                gcp_id = int(l[1])
                e = float(l[2])
                n = float(l[3])
                h = float(l[4])

                c = np.array([e, n, h])

                gcps[gcp_id] = Point(gcp_id)
                gcps[gcp_id].SetGcp(c, control)
                

        return cams, gcps

    def _readPhotosFile(self):
        f = open(self.path + '/gps.dat')

        phs = {}
        count = 0

        while 1:
            line = f.readline()
            if not line:
                break
            if line[0] == "#":
                continue

            lin=line.rstrip("\r\n")
            l = lin.split()
            if not l:
                continue
            if "LINE" in l[0]:
                continue

            ph_id = int(l[0])
            #TOOD just one cam
            cam = self.cam

            eo = map(float, l)
            phs[ph_id] = Photo(ph_id, self.cam) 
                               #np.array(eo), TODO
                               #None
                     
        return phs

    def _readPhotoPointsFile(self, pts, phs):
        f = open(self.path + '/image.dat')

        count = 0

        new_ph = True

        #TODO works with just one camere
        while 1:
            line = f.readline()
            if not line:
                break
            if line[0] == "#":
                continue

            lin=line.strip()
            l = lin.split()
            if not l:
                continue

            if '-99' == l[0]:
                new_ph = True
                continue

            if new_ph:
                #if len(l) > 1:
                #TODO camera
                ph_id = int(l[0])
                new_ph = False
                continue

            try:
                pt_id = int(l[0])
            except:
                pt_id = int(l[0][1:])

            #TODO 
            c = np.array([float(l[1]), float(l[2])])

            if not pts.has_key(pt_id):
                pts[pt_id] = Point(pt_id)

            pts[pt_id].AddPhotoPoint(phs[ph_id], c)

