import os
import sys
import numpy as np

import grass.script as grass
from relative_or import compute3Dpoints

try:
    from grass.lib.imagery import *
except ImportError, e:
    sys.stderr.write(_("Loading ctypes libs failed"))
from ctypes import *


# this file gathers functions dealing with transformation of common relative coordinate system 
# into world coordinate system (for in depth description see chapters 1.6.7, 3.2.3, 4))

def HelmertTransform(pts, phs):
    # computes 
    ComputeCoords(pts)

    gcp_pts = []
    ro_pts = []

    # get GCPS for computation of Helemert transformation coefficients 
    for pt in pts.itervalues():

        # only non control GCPS are used (control GCPS are treated as tie points)
        gcp_c, control = pt.GetGcp()
        if gcp_c is not None and not control and pt.GetCoords() is not None:
            gcp_pts.append(gcp_c)
            ro_pts.append(pt.GetCoords())
    
    ro_pts = np.array(ro_pts)
    gcp_pts = np.array(gcp_pts)
    
    # computes helmert transformation coefficients 
    hm_p, hm_s, r = ComputeHelmertTransf(ro_pts, gcp_pts)

    # transforms all object points by Helmert transformation coefficients
    for pt in pts.itervalues():
        if pt.GetCoords() is not None:
            p = GeoRefOr(pt.GetCoords(), r)
            pt.SetCoords(p)

    # transforms all exterior orientations of photos
    for ph in phs.itervalues():
        if ph.GetEO() is not None:
            eo = ph.GetEO()
            t1 = eo[:, 3]
            r1 = eo[:, :3]
            ph.HelmertTransformEO(hm_p, hm_s)

def ComputeCoords(pts):
    # median is used for merging corresponding object points coordinates  
    # determined by different relative orientation 
    # 
    for pt in pts.itervalues():
        ro_pts = pt.GetRelOrsCoords()

        #for ro_pt in ro_pts.itervalues():
        #    plot_pts.append(ro_pt)

        tmp = []
        for ro_pt in ro_pts.itervalues():
            tmp.append(ro_pt)

        if tmp:
            m = np.median(np.array(tmp), axis=0)
            pt.SetCoords(m)
            #plot_pts.append(m)

def Test(ros, in_dt, cam_m, distor, rel_or_phs):
    ph_ids_to_ro_id = {}


    phs = in_dt.GetPhotos()
    pts = in_dt.GetPoints()

    for i_ro_id, ro_phs in enumerate(rel_or_phs):


        #tie_point_ids, tie_pts = in_dt.GetTiePoints(ro_phs)
        ro = ros[0][i_ro_id]
        tie_pts = ro.tie_pts

        eo_ph = phs[ro_phs[0]].GetP()
        eo_trans = phs[ro_phs[1]].GetP()
        out_pts = compute3Dpoints(eo_ph, eo_trans, tie_pts[:,:2], tie_pts[:,2:])

        for i_pt, pt_add_id in enumerate(ro.pt_ids):
            pts[pt_add_id].AddRelOrCoords(i_ro_id, out_pts[i_pt])

    ComputeCoords(pts)

def ComputeHelmertTransf(pts1, pts2):
    # calls GRASS imagery library function I_compute_georef_equations_or 
    # through ctypes interface to calculate transf. coefficients   
    def RotMatCToRotMat(rm_c):
        
        rm = np.zeros((3, 4))
        for i in range(9):
            rm[i/3, i%3] = rm_c[i]

        rm[0,3] = rm_c[9]
        rm[1,3] = rm_c[10]
        rm[2,3] = rm_c[11]

        s = rm_c[12]
        return rm, s

    def CToRotMatRotMat(rm, scale, shift):
        c_double_p = POINTER(c_double)
        DArray = c_double * 15

        rm_c = DArray()
        for i in range(9):
            rm_c[i] = rm[i/3, i%3]
        #shift = [0,0,0]
        #scale = 0

        rm_c[9] =  shift[0]
        rm_c[10] = shift[1]
        rm_c[11] = shift[2]
        rm_c[12] = scale
        rm_c[13] = scale
        rm_c[14] = scale

        return rm_c

    r1 = c_double * 15
    r2 = c_double * 15 
    r1 = r1()
    r2 = r2()
            
    pts_c = Control_Points_3D()
    pts_c.count = len(pts1)

    c_double_p = POINTER(c_double)
    c_int_p = POINTER(c_int)

    DArray = c_double * len(pts1)
    pts_c.e1 = DArray()
    pts_c.n1 = DArray()
    pts_c.z1 = DArray()

    pts_c.e2 = DArray()
    pts_c.n2 = DArray()
    pts_c.z2 = DArray()

    IArray = c_int * len(pts1)
    pts_c.status = IArray()

    for i, pt2 in enumerate(pts2):
        pt1 = pts1[i]

        pts_c.e1[i] = c_double(pt1[0]) 
        pts_c.n1[i] = c_double(pt1[1]) 
        pts_c.z1[i] = c_double(pt1[2])

        pts_c.e2[i] = c_double(pt2[0])
        pts_c.n2[i] = c_double(pt2[1])
        pts_c.z2[i] = c_double(pt2[2])
        s = c_int(1)
        pts_c.status[i] = s

    I_compute_georef_equations_or(byref(pts_c),
                                         r1, r2)

    P, s = RotMatCToRotMat(r1)
    #P1, s1 = RotMatCToRotMat(r2)
    r2 = CToRotMatRotMat(P[:,:3], s, P[:,3])

    return P, s, r2


def GeoRefOr(ph_ob, rm_c):
    # transforms point by Helmert transformation coefficients
    e = c_double(0) 
    n = c_double(0)
    z = c_double(0) 
                
    I_georef_or(ph_ob[0], ph_ob[1], ph_ob[2], 
                byref(e), byref(n), byref(z), rm_c)

    #print "%.2f %.2f %.2f" % (e.value, n.value, z.value)
    return np.array([e.value, n.value, z.value])

