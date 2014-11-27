import os
import sys
import numpy as np
import scipy as sc


from plots import PlotRelOrs, PlotAllRelOrs
import grass.script as grass
from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE
from math import sqrt
import cv2
from copy import deepcopy

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


# this file comprise functions  performing relative orientation by five point algorithm 
# and than merges individual relative coordinates systems into the common one (for in depth description  (1.6.5, 1.6.6, 3.2.3, 4))

def solvePnP(in_dt):
    # computes exterior orientation from photo points and corresponding object points (not used)
    phs = in_dt.GetPhotos()
    pts = in_dt.GetPoints()
    cam_m, distor = in_dt.GetCameras().values()[0].GetParams()
    import cv2


    for ph in phs.itervalues():
        img_pts = []
        obj_pts = []
        pts = ph.GetPoints()
        i = 0
        for pt in pts.itervalues():
            gcp = pt.GetGcp()[0]
            phpts = pt.GetPhotoPoints()
            phpt = phpts[ph]
            img_pts.append([phpt[0], phpt[1]])
            obj_pts.append([gcp[0], gcp[1], gcp[2]])
            i += 1
            #if i == 4:
            #    break
        img_pts = np.array([img_pts])
        obj_pts = np.array([obj_pts])
    
        ret, r, t = cv2.solvePnP(obj_pts, img_pts, cam_m, distor, flags=cv2.EPNP)
        ret, r, t = cv2.solvePnP(obj_pts, img_pts, cam_m, distor, rvec=r, tvec=t, useExtrinsicGuess=True)

        r = cv2.Rodrigues(r)[0]
        t = t[:,0]
        c = - r.T.dot(t)
        eo = ph.GetEO()
        c_eo = eo[:,3]

        p = np.zeros((3, 4))
        p[:, :3] = r
        p[:, 3] = t

        ph.SetP(p)
    return

    #cv2.solvePnP(objectPoints, imagePoints, cameraMatrix, distCoeffs[, rvec[, tvec[, useExtrinsicGuess[, flags]]]])  retval, rvec, tvec
class RelOr:
    # stores needed information about relative orientation 
    def __init__(self, id, ph1_id, ph2_id, pm, pt_ids, sps_pts, tie_pts):
        self.id = id
        self.ph1_id = ph1_id
        self.ph2_id = ph2_id
        self.pm = pm
        self.pt_ids = pt_ids
        self.sps_pts = sps_pts
        self.tie_pts = tie_pts
        self.scale = None

def AddRelOrToIndex(ph1_id, ph2_id, ro_id, ph_ids_to_ro_id):
    #  the index are dictionaries in dictionary, it is 
    # used for search of all corresponding relative 
    # orientation where the photo belongs  
    # e. g.  ph_ids_to_ro_id[551] returns 
    #  {552 : RelOr class object, 560 : RelOr class object}
    # it means that  the photo 551 is included in two computed relative orientations 

    # the function adds new relative orientation into the index

    if ph1_id not in ph_ids_to_ro_id:
        ph_ids_to_ro_id[ph1_id] = {ph2_id : ro_id}
    else:
        ph_ids_to_ro_id[ph1_id][ph2_id] = ro_id

    if ph2_id  not in ph_ids_to_ro_id:
        ph_ids_to_ro_id[ph2_id] = {ph1_id : ro_id}
    else:
        ph_ids_to_ro_id[ph1_id][ph2_id] = ro_id


def RelativeOrientation(in_dt):
    ph_ids_to_ro_id = {}
    ros = []

    print "computing relative orientations"

    cam_m, distor = in_dt.GetCameras().values()[0].GetParams()

    phs = in_dt.GetPhotos()
    pts = in_dt.GetPoints()

    # first pair with the highest number of tie points is returned 
    ro_phs, phs_outer_envelope, phs_oriented =  GetFirstPair(in_dt)

    i_ro_id = 0
    # it goes through the graph of relative orientations until all 
    # available photos are  relatively oriented and transformed in
    # the common system  
    while True:
        # newly oriented photos
        ph1, ph2 = ro_phs

        ph1_id = ph1.GetId()
        ph2_id = ph2.GetId()

        print "ph1"
        print ph1_id


        print "ph2"
        print ph2_id

        # their common tie points
        tie_point_ids, tie_pts = in_dt.GetTiePoints(ro_phs)

        # relative orientation is computed and object 
        # point coordinates are determined by triangulation  
        sps_pts, P, p1, p2 = computeRelativeOrientation(tie_pts, 
                                                        cam_m, 
                                                        distor)

        tie_pts = np.hstack((p1[0, ...] , p2[0, ...]))

        #  RelOr  structure is created
        ro = RelOr(i_ro_id, ro_phs[0], ro_phs[1], P, tie_point_ids, sps_pts, tie_pts)
        ros.append(ro)


        # first relative orientation determines common relative orientation system 
        if i_ro_id == 0:
            pm1 = np.eye(3, 4)
            phs[ro_phs[0]].SetP(pm1)
            phs[ro_phs[1]].SetP(ro.pm)
            out_pts = ro.pt_ids
        
            ro.scale = 1

            for i_pt, pt_id in enumerate(ro.pt_ids):
                pts[pt_id].AddRelOrCoords(i_ro_id, sps_pts[i_pt][:3])
       
        # another orientations are transformed into the common system
        else:
            if ph_ids_to_ro_id.has_key(ro_phs[0]):
                conn_ph = phs[ro_phs[0]]
            elif ph_ids_to_ro_id.has_key(ro_phs[1]):
                conn_ph = phs[ro_phs[1]]
    
            if ph_ids_to_ro_id.has_key(ro.ph1_id):
                ro_conn = ros[ph_ids_to_ro_id[ro.ph1_id].values()[0]]
            else:
                ro_conn = ros[ph_ids_to_ro_id[ro.ph2_id].values()[0]]

            eo_ph, eo_trans = ConnectRelOr(phs, ro, ro_conn)
            out_pts = compute3Dpoints(eo_ph, eo_trans, tie_pts[:,:2], tie_pts[:,2:])

            for i_pt, pt_add_id in enumerate(ro.pt_ids):
                pts[pt_add_id].AddRelOrCoords(i_ro_id, out_pts[i_pt])

        
        #if (ph1_id, ph2_id) == (6985, 6986):


        apts =  {}
        for pid in ro.pt_ids:
            apts[pts[pid]] = pts[pid]

        aphs =  {}
        for phid in ro_phs:
            aphs[phs[phid]] = phs[phid]

        PlotRelOrs(apts, aphs, i_ro_id)


        PlotAllRelOrs(pts, phs)

        # new relative orientation is added into the index
        AddRelOrToIndex(ro_phs[0], ro_phs[1], i_ro_id, ph_ids_to_ro_id)


        # next pair to be oriented is returned from all available pairs 
        # where one photo is already transformed into common system 
        # and another is not. It is chosen the  pair with highest number of shared tie points. 
        ret = NextPair(ros, ph_ids_to_ro_id, phs_outer_envelope, phs_oriented)

        # all pairs are already merged, terminate the loop
        if not ret:
            break
        ro_phs = ret
 
        i_ro_id += 1


    #createRelOrProtocol(protocol_path, ros, cam_m, distor, in_dt)
        

    # computes points, whose object coordinates have not been determined 
    # by relative orientation process 
    computeMissingPoints(in_dt)
    return ros, ph_ids_to_ro_id

def computeMissingPoints(in_dt):

    phs = in_dt.GetPhotos()
    pts = in_dt.GetPoints()

    for pt in pts.itervalues():
        if not pt.GetRelOrsCoords():

            ph_pts = pt.GetPhotoPoints()
            ph_pts_list = ph_pts.keys()
            found = False
            for i, ph_pt1 in enumerate(ph_pts_list[:-1]):
                if ph_pt1.GetP() is not None:
                    for ph_pt2 in ph_pts_list[i + 1:]:
                        if ph_pt2.GetP() is not None:

                            p1 = ph_pt1.GetP()
                            p2 = ph_pt2.GetP()

                            tie_pts1 = np.array([ph_pts[ph_pt1]])
                            tie_pts2 = np.array([ph_pts[ph_pt2]])

                            cam_m1, distor1 = ph_pt1.GetCamera().GetParams()
                            cam_m2, distor2 = ph_pt2.GetCamera().GetParams()

                            tp1_u = undistortPoints(tie_pts1, cam_m1, distor1)
                            tp2_u = undistortPoints(tie_pts2, cam_m1, distor1)

                            out_pts = compute3Dpoints(p1, p2, tp1_u, tp2_u)
                            pt.AddRelOrCoords(-1, out_pts[0])
                        found = True
                        break
                if found:
                    break

def GetFirstPair(in_dt):
    # returns pair with highest number of tie points in the set

    sorted_pairs = in_dt.GetPhotosConnectivity()
    phs_pair = sorted_pairs[0][0]

    ph1, ph2 = phs_pair
    phs_outer_envelope = [ph1, ph2]
    phs_oriented  = [ph1, ph2]
    

    # envelope holds information about 
    # candidate photos for next orientation 
    # phs_oriented is list of already  oriented  photos
    return phs_pair, phs_outer_envelope, phs_oriented

def NextPair(ros, ph_ids_to_ro_id, phs_outer_envelope, phs_oriented):
    # returns pairs, for next relative orientation 
    conn_ph_tie_pts_n = 0
    conn_ph = None

    for i_out_ph, out_ph in enumerate(phs_outer_envelope[:]):
        neigh_phs = out_ph.GetNeighboroughPhotos()
        
        ph_in_envelope = False
        
        for ph, tie_pts_n in dict(neigh_phs).iteritems():
            if ph in phs_oriented:
                continue

            ph_in_envelope = True
            if tie_pts_n > conn_ph_tie_pts_n:
                conn_ph = ph
                env_ph = out_ph
                conn_ph_tie_pts_n = tie_pts_n 


        if not ph_in_envelope:
            phs_outer_envelope.remove(out_ph)

    if not conn_ph:
        return None

    phs_outer_envelope.append(conn_ph)
    phs_oriented.append(conn_ph)

    return (env_ph, conn_ph)


def ConnectRelOr(phs, add_rel_or, conn_rel_or):
    # transforms relative orientation add_rel_or into common system

    def GetPivotPhoto(rel_or1, rel_or2):
        # returns photo which is present in both relative orientations 

        if rel_or1.ph1_id == rel_or2.ph1_id:
            pivot_ph_id = rel_or1.ph1_id
        elif rel_or1.ph1_id == rel_or2.ph2_id:
            pivot_ph_id = rel_or1.ph1_id
        else:
            pivot_ph_id = rel_or1.ph2_id

        return pivot_ph_id

    def computeScale(add_rel_or, conn_rel_or, phs):
        #  computes scale for transformation of add_rel_or relative orientation 
        #  into conn_rel_or
        pivot_ph_id = GetPivotPhoto(add_rel_or, conn_rel_or)

        if add_rel_or.ph1_id == pivot_ph_id:
            eo_add_ro = np.eye(3, 4)
        else:
            eo_add_ro = add_rel_or.pm

        if conn_rel_or.ph1_id == pivot_ph_id:
            eo_conn_ro = np.eye(3, 4)
        else:
            eo_conn_ro = conn_rel_or.pm
        scales = []

        t_conn = eo_conn_ro[:, 3]
        R_conn = eo_conn_ro[:, :3] 

        t_add = eo_add_ro[:, 3]
        R_add = eo_add_ro[:, :3] 

        for i_pt_add, pt_add_id in enumerate(add_rel_or.pt_ids):
            if pt_add_id in conn_rel_or.pt_ids:
                i_pt_conn = conn_rel_or.pt_ids.index(pt_add_id)
                pt_conn = conn_rel_or.sps_pts[i_pt_conn, :3]
                pt_conn = R_conn.dot(pt_conn) + t_conn

                pt_add = add_rel_or.sps_pts[i_pt_add, :3]
                pt_add = R_add.dot(pt_add) + t_add

                scales.append(np.linalg.norm(pt_conn) / np.linalg.norm(pt_add))

        scales = np.array(scales)
        med = np.median(scales)

        return med

    pivot_ph_id = GetPivotPhoto(add_rel_or, conn_rel_or)

    pivot_ph = phs[pivot_ph_id]

    pm_ro = add_rel_or.pm
    ro_pt_ids = add_rel_or.pt_ids
    t_ro = pm_ro[:, 3]
    R_ro = pm_ro[:, :3]

    #ro_pts = add_rel_or.sps_pts[:,:3]
    if pivot_ph == add_rel_or.ph2_id:
        add_ph = phs[add_rel_or.ph1_id]
    else:
        add_ph = phs[add_rel_or.ph2_id]

    eo_ph = pivot_ph.GetP()
    t_ph = eo_ph[:, 3]
    R_ph = eo_ph[:, :3]

    scale = computeScale(add_rel_or, conn_rel_or, phs)

    add_rel_or.scale = scale * conn_rel_or.scale

    t_trans = R_ro.dot(t_ph) + add_rel_or.scale * t_ro
    R_trans = R_ro.dot(R_ph)
    
    eo_trans = np.hstack((R_trans, np.expand_dims(t_trans, axis=1)))
    add_ph.SetP(eo_trans)

    return eo_ph, eo_trans 

def undistortPoints(pts, cam_m, distor):
    # returns normalized coordinates of photo points    
    tp = np.expand_dims(pts, axis=0)
    tp_u = cv2.undistortPoints(tp, cam_m, distor)
    return tp_u

def UndistorTiePoints(tie_pts, cam_m, distor):
    # returns coordinates in photo system without effect of distortion 
    tp1_u = undistortPoints(tie_pts[:, :2], cam_m, distor)
    tp2_u = undistortPoints(tie_pts[:, 2:4], cam_m, distor)

    shape = (1, len(tie_pts), 2)
    tp1 = np.empty(shape)
    tp2 = np.empty(shape)
    tp1[0,:,0] = cam_m[0, 0] * tp1_u[0,:,0] + cam_m[0, 2]
    tp1[0,:,1] = cam_m[1, 1] * tp1_u[0,:,1] + cam_m[1, 2]

    tp2[0,:,0] = cam_m[0, 0] * tp2_u[0,:,0] + cam_m[0, 2]
    tp2[0,:,1] = cam_m[1, 1] * tp2_u[0,:,1] + cam_m[1, 2]

    return tp1, tp2, tp1_u, tp2_u

def computeRelativeOrientation(tie_pts, cam_m, distor):
    #  computes relative orientation  using OpenCV 5 point algorithms

    tp1, tp2, tp1_u, tp2_u  = UndistorTiePoints(tie_pts, cam_m, distor)

    # eight point algorithm 
    #F, mask = cv2.findFundamentalMat(tp1[0, :], tp2[0, :], param1=0.1, param2=0.95, method = cv2.FM_RANSAC)
    #em = cam_m.T.dot(F).dot(cam_m)
    
    em, mask = cv2.findEssentialMat(tp1[0, :], tp2[0, :], 
                                    threshold=0.05, 
                                    prob=0.95, 
                                    focal=cam_m[0, 0], 
                                    pp=(cam_m[0, 2], cam_m[1, 2]))


    F = np.linalg.inv(cam_m.T).dot(em).dot(np.linalg.inv(cam_m))

    # optimal solution for triangulation of  object points
    p1, p2 = cv2.correctMatches(em, tp1_u, tp2_u)

    pts, R, t, mask = cv2.recoverPose(em, p1, p2)
    P = np.hstack((R, t))
    pm1 = np.eye(3, 4)

    pts_sps = compute3Dpoints(pm1, P, p1, p2)

    return pts_sps, P, p1, p2

def computeError(phs_pts, sps_pts, R, t, cam_m, distor):

    rvec, J = cv2.Rodrigues(R)
    #cam_m = np.eye(3)

    phs_pts_rep, J = cv2.projectPoints(sps_pts, rvec, t, cam_m, distor)
    sh = phs_pts_rep.shape
    phs_pts_rep = phs_pts_rep.reshape((sh[1], sh[0], sh[2]))

    s = phs_pts[0,...] - phs_pts_rep[0,...]
    errs = []
    for i in s:
        errs.append(np.linalg.norm(i))

    return np.array(errs)

def compute3Dpoints(P1, P2, npts1, npts2):
    # computes object point coordinates from photo points correspondence using DLT method 
    ptsh3d = cv2.triangulatePoints(P1, P2, npts1.T, npts2.T).T

    pts_sps = deepcopy(ptsh3d[:,:3])
    for i, pt in enumerate(ptsh3d):
        pt = (pt / pt[3])[:3]
        pts_sps[i, :] = pt[:3]

    return pts_sps

