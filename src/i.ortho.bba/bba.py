import os
import sys
import numpy as np
from math import atan2, sin, cos, sqrt, pi
from bba_derivs import x_derivs, y_derivs

import scipy.misc as misc
import scipy.sparse as sparse
from  scipy.sparse.linalg import splu

from protocols import CreateIterationProtocolMeasurements, CreateParametersIterationProtocol, reprojection_errors, \
 GenerateLatexPhotosIterationTable, GenerateLatexPhotosAnglesIterationTable, GenerateLatexPointsIterationTable, \
 CreateParametersLatex, CreatePhootoParametersLatex

from collections import namedtuple
from sympy.matrices import *
import cv2

# bundle block adjustment (see chapters 1.7, 2.1.6, and 3.2.4)
# computes distances 
#  between projective center coordinates  and object points  form Bingo-F results and current results 
# (e.g. after some adjustment iteration) 
# it also computes differences of angles 
def ComputeDiffs(pts, phs):
    gcp_pts = []
    ro_pts = []

    pt_errs = []
    no_coords = []
    pt_errs = [] #np.ma.zeros((len(pts), 2))
    #pt_errs.mask = np.ma.zeros((len(pts), 2))
    for i_pt, pt in enumerate(pts.itervalues()):
        o_pt, control = pt.GetGcp()
        if  (o_pt is not None and not control):
            continue

        if pt.GetCoords() is not None and pt.GetResultCoords() is not None: 
            pt_errs.append([pt.GetId(), abs(np.linalg.norm(pt.GetCoords() - pt.GetResultCoords()))])
        #else:
            #pt_errs.mask[i_pt, 1] = True
         
    pt_errs = np.array(pt_errs)

    skip_pts = []
    for err in pt_errs:
        if err[1] > 0.3:
            skip_pts.append(int(err[0]))

    ph_errs = []
    ph_angle_errs = []

    for ph in phs.itervalues():
        if ph.GetEO() is not None:
            eo = ph.GetEO()
            res_eo = ph.GetResultExtOr()
            
            a = np.array([[1,0,0],[0,1,0],[0,0,-1]])
            angles = GetBLUHAngles(a.dot(eo[:,:3]))

            res_coords = res_eo[:,3]
            res_angles = GetBLUHAngles(np.array(res_eo))

            ph_errs.append([ph.GetId(), 
                            abs(np.linalg.norm(eo[:,3] - res_coords))] )
            ph_angle_errs.append(np.append(ph.GetId(), angles - res_angles))
        else:
            no_coords.append(ph.GetId())

    ph_errs = np.array(ph_errs)

    ph_angle_errs = np.array(ph_angle_errs)

    return ph_errs, pt_errs, ph_angle_errs
# rotation matrices / euler angles  functions 


def RotMatBLUH(ph, om, ka):
    rotm = np.zeros((3,3))
    _rotMatBLUH(ph, om, ka, rotm, sin, cos)

    return rotm

def GetBLUHAngles(rm):

    ka = atan2(rm[0, 1], rm[1, 1])
    ph = atan2(rm[2, 0], rm[2, 2])
    om = atan2(-rm[2, 1], sqrt(rm[1, 1] * rm[1, 1] + rm[0, 1] * rm[0, 1]))

    return np.array([ph, om, ka])

def _rotMatBLUH(ph, om, ka, rotm, sin, cos):

    so = sin(om)
    co = cos(om)
    sp = sin(ph)
    cp = cos(ph)
    sk = sin(ka)
    ck = cos(ka)

    rotm[0][0] = ck * cp + sk * so * sp
    rotm[1][0] = - sk * cp + ck * so * sp
    rotm[2][0] = co * sp

    rotm[0][1] = sk * co
    rotm[1][1] = ck * co
    rotm[2][1] = - so
    
    rotm[0][2] =  - ck * sp + sk * so * cp  
    rotm[1][2] =  sk * sp + ck * so * cp
    rotm[2][2] = co * cp


# collinearity equations
def xp_col_eq(xi_c, yi_c, f, xo_c, yo_c, zo_c ,xo_pt, yo_pt, zo_pt, r, r1_dis, **kwargs):
    x_und = _xp_col_eq(xi_c, f, xo_c, yo_c, zo_c ,xo_pt, yo_pt, zo_pt, r)
    y_und = _yp_col_eq(yi_c, f, xo_c, yo_c, zo_c ,xo_pt, yo_pt, zo_pt, r)
    return x_und
    #dist_sq = x_und ** 2 +  y_und ** 2

    #return x_und + x_und * dist_sq * r1_dis / 1000.0

def yp_col_eq(xi_c, yi_c, f, xo_c, yo_c, zo_c, xo_pt, yo_pt, zo_pt, r, r1_dis, **kwargs):
    x_und = _xp_col_eq(xi_c, f, xo_c, yo_c, zo_c ,xo_pt, yo_pt, zo_pt, r)
    y_und = _yp_col_eq(yi_c, f, xo_c, yo_c, zo_c ,xo_pt, yo_pt, zo_pt, r)
    return y_und
    #dist_sq = x_und ** 2 +  y_und ** 2

    #return y_und *(1 + dist_sq * r1_dis)

def _xp_col_eq(xi_c, f, xo_c, yo_c, zo_c ,xo_pt, yo_pt, zo_pt, r, **kwargs):
    return xi_c - f * (r[0][0]*(xo_pt-xo_c) + r[0][1]*(yo_pt-yo_c) + r[0][2]*(zo_pt-zo_c)) / \
                      (r[2][0]*(xo_pt-xo_c) + r[2][1]*(yo_pt-yo_c) + r[2][2]*(zo_pt-zo_c)) 

def _yp_col_eq(yi_c, f, xo_c, yo_c, zo_c, xo_pt, yo_pt, zo_pt, r, **kwargs):
    return yi_c - f * (r[1][0]*(xo_pt-xo_c) + r[1][1]*(yo_pt-yo_c) + r[1][2]*(zo_pt-zo_c)) / \
                      (r[2][0]*(xo_pt-xo_c) + r[2][1]*(yo_pt-yo_c) + r[2][2]*(zo_pt-zo_c))

# returns all parameters of collinearity equations defined 
# by photo point and measured photo point coordinates. The photo point can be given as combination of photo and point objects (see redpar.py) 
# because Point object represents all photo points and given Photo object specifies single photo point
# The function returns photo point coordinates, object point coordinates, exterior and interior orientation of photo.
 
def GetParamsVals(ph, pt, free_net):

    params = {}

    cam = ph.GetCamera()

    cam_m, distor = cam.GetParams()

    params['r1_dis'] = distor[0] 
    params['f'] = cam_m[0, 0]
    params['xi_c'] =  cam_m[0, 2]
    params['yi_c'] = cam_m[1, 2]

    ph_eo = ph.GetEO()

    #ph_p = ph.GetP()
    #r = ph_p[:,:3]
    #params['rvect'] = cv2.Rodrigues(r)[0]
    #params['tvect'] = ph_p[:,3]
    #params['distor'] = distor
    #params['cam_m'] = cam_m

    o_c = ph_eo[:,3]
    eo_r = ph_eo[:,:3]

    # reflection caused by different handness of Opencv and BLUH systems 
    a = np.array([[1,0,0],[0,1,0],[0,0,-1]])
    angles = GetBLUHAngles(a.dot(eo_r))
    for i, ang in enumerate(['ph', 'om', 'ka']):
        params[ang] = angles[i]

    for i, ax in enumerate(['xo_c', 'yo_c', 'zo_c']):
        params[ax] = o_c[i]

    o_pt, control = pt.GetGcp()
    if not free_net and (o_pt is not None and not control):
        params["gcp"] = True 
    else:
        params["gcp"] = False
        o_pt = pt.GetCoords()

    if o_pt is None:
        return None

    for i, ax in enumerate(['xo_pt', 'yo_pt', 'zo_pt']):
        params[ax] = o_pt[i]

    phpts = pt.GetPhotoPoints()

    i = 0
    for phpt in phpts.iterkeys():
        if phpt.GetP() is not None:
            i += 1

    if i < 2:
        return None


    i_pt = phpts[ph]
    params['xi_pt'] = i_pt[0]
    params['yi_pt'] = i_pt[1]

    return params
    #.evalf(<args>).

BbaIndexes = namedtuple('BbaIndexes', 'Lrows_idxs,' 
                                      'pts_idxs,  idxs_pts, '
                                      'phs_idxs, idxs_phs, '
                                      'cams_idxs, idxs_cams, '
                                      'cams_n, photos_n, tie_pts_n, '
                                      'image_points_n')


# builds L and L_XO vectors for least squares iteration  (see chapter 1.7.)
def GetL(in_dt, free_net):
    cams = in_dt.GetCameras()
        
    L = []
    L_XO = []


    phs = in_dt.GetPhotos()

    #TODO for testing
    use_photos = phs.keys()
    #use_photos = [598, 597]

    cams_idxs = {}
    idxs_cams = []

    phs_idxs = {}
    idxs_phs = []
    
    pts_idxs = {}
    idxs_pts = []
    
    Lrows_idxs = []

    photos_n = 0
    image_points_n = 0
    tie_pts_n = 0
    cams_n = 0

    for cam in cams.itervalues():
        phs = cam.GetPhotos()

        for ph in phs.itervalues():

            ph_pts = ph.GetPoints()

            photo_points_n = 0
            if ph.GetId() not in use_photos: 
                continue

            for pt in ph_pts.itervalues():

                #if not i use_photos: 
                pt_phs = pt.GetPhotoPoints().keys()
                if len((set(pt_phs) & set(use_photos))) < 2:
                    continue

                v = GetParamsVals(ph, pt, free_net)

                if v is not None:
                    if cam not in idxs_cams:
                        idxs_cams.append(cam) 
                        cams_idxs[cam] = cams_n
                        cams_n += 1

                    if ph not in idxs_phs:
                        idxs_phs.append(ph) 
                        phs_idxs[ph] = photos_n
                        photos_n += 1

                    if pt not in pts_idxs and not v['gcp']:
                        pts_idxs[pt] = tie_pts_n
                        idxs_pts.append(pt) 
                        tie_pts_n +=  1

                    Lrows_idxs.append((cam, ph, pt, v))

                    r = RotMatBLUH(v['ph'], v['om'], v['ka'])

                    #pt = np.array([[v['xo_pt'], v['yo_pt'], v['zo_pt']]])
                    #ip = cv2.projectPoints(pt, v['rvect'], v['tvect'], v['cam_m'], v['distor'])[0]
                    #L_XO.append(ip[0, 0,0])
                    #L_XO.append(ip[0,0,1])
                    L_XO.append(xp_col_eq(r=r, **v))
                    L_XO.append(yp_col_eq(r=r, **v))

                    L.append(v['xi_pt'])
                    L.append(v['yi_pt'])

                    photo_points_n += 1
                else:
                    pass

        image_points_n += photo_points_n
        bba_idxs = BbaIndexes(Lrows_idxs, 
                              pts_idxs, idxs_pts, phs_idxs, idxs_phs, cams_idxs, idxs_cams, 
                              cams_n, photos_n, tie_pts_n, image_points_n)
    return L, L_XO, bba_idxs


# updates Camera, Photo and  Point objects with adjusted parameters after BBA iteration
def AdjustData(Xa, bba_idxs, unk, apo):
    bi = bba_idxs

    cam_step = len(unk.cam)
    ph_step = len(unk.ph)
    pt_step = len(unk.pt)

    for cam_idx in range(bi.cams_n):
        cam = bi.idxs_cams[cam_idx]
        cam_x_col = apo.cam_0col + cam_idx * cam_step

        cam_io = {}
        for i, un in enumerate(unk.cam):
            cam_io[un] = Xa[cam_x_col + i]

        cam.SetIO(cam_io)

    for pt_idx in range(bi.tie_pts_n):

        pt = bi.idxs_pts[pt_idx]
        pt_x_col = apo.pt_0col + pt_idx * pt_step

        pt_coords = np.zeros((3))

        for i in range(pt_step):
            pt_coords[i] = Xa[pt_x_col + i]
        pt.SetCoords(pt_coords)

    for ph_idx in range(bi.photos_n):
        ph = bi.idxs_phs[ph_idx]
        ph_x_col = apo.ph_0col + ph_idx * ph_step

        ph_eo = np.zeros((ph_step))
        for i in range(ph_step):
            ph_eo[i] = Xa[ph_x_col + i]

        # reflection back TODO error prone, merge it
        a = np.array([[1,0,0],[0,1,0],[0,0,-1]])
        eo_r = a.T.dot(RotMatBLUH(*ph_eo[:3]))
        eo = np.hstack((eo_r, np.array([ph_eo[3:6]]).T))
        ph.SetEO(eo)


Unknowns = namedtuple('Unknowns', 'cam, ph, pt')


# runs BBA
def RunBundleBlockAdjutment(in_dt, protocol_dir, free_net):

    prot_fd = open(os.path.join(protocol_dir, "parameters_prot.txt") , 'w')
    prot_mesurements_fd = open(os.path.join(protocol_dir, "measur_prot.txt"), 'w')
    
    phs = in_dt.GetPhotos()
    pts = in_dt.GetPoints()

    cam_unk = []
    ph_unk = ['ph', 'om', 'ka', 'xo_c', 'yo_c', 'zo_c']
    pt_unk = ['xo_pt', 'yo_pt', 'zo_pt']

    unk =  Unknowns(cam_unk, ph_unk, pt_unk)


    ph_errs, pt_errs, ph_angles_errs = ComputeDiffs(in_dt.pts, in_dt.phs)

    i = 0
    while True:

         #   cam_unk = ['f', 'xi_c', 'yi_c']
         #   unk =  Unknowns(cam_unk, ph_unk, pt_unk)

        #if i == 10:
        #    free_net = True
        #    ph_unk =  []
        #    pt_unk = ['xo_pt', 'yo_pt', 'zo_pt']
        #    cam_unk = ['f', 'xi_c', 'yi_c']
        #    unk =  Unknowns(cam_unk, ph_unk, pt_unk)

        print ("iteration: %d" % i)
        Xa, dx, e_taylor, e_repro, bi, apo, l, cov_x, cov_l = RunBBAIteration(in_dt, unk, free_net)
        CreateParametersIterationProtocol(i, cov_x, Xa, dx, bi, apo, prot_fd, unk, in_dt.GetGcps(), free_net, l)
        CreateIterationProtocolMeasurements(i, cov_l, e_repro, bi, prot_mesurements_fd)

        it_ph_errs, it_pt_errs, it_ph_angles_errs = ComputeDiffs(pts, phs)
        it_ph_errs = np.expand_dims(it_ph_errs[:,1], axis=0)
        it_pt_errs = np.expand_dims(it_pt_errs[:,1], axis=0)
        it_ph_angles_errs = it_ph_angles_errs[:,1:]

        ph_errs = np.ma.hstack((ph_errs, it_ph_errs.T)) 
        pt_errs = np.ma.hstack((pt_errs, it_pt_errs.T))
        ph_angles_errs = np.ma.hstack((ph_angles_errs, it_ph_angles_errs)) 

        #iteration termination criterion 
        if np.max(np.absolute(dx))  < 0.0001:
            break
        i += 1
    GenerateLatexPhotosIterationTable(ph_errs, os.path.join(protocol_dir, "itera_photos.txt"))
    GenerateLatexPhotosAnglesIterationTable(ph_angles_errs, os.path.join(protocol_dir, "itera_angles.txt"))
    GenerateLatexPointsIterationTable(pt_errs, os.path.join(protocol_dir, "itera_points.txt"))
    CreateParametersLatex(i, cov_x, Xa, dx, bi, apo, os.path.join(protocol_dir, "points_res.txt"), unk, in_dt.GetGcps(), free_net, l)
    CreatePhootoParametersLatex(i, cov_x, Xa, dx, bi, apo, os.path.join(protocol_dir, "photo_res.txt"), unk, in_dt.GetGcps(), free_net, l)
    reprojection_errors(e_repro, "/home/ostepok/Desktop/itera_repro.txt")
    prot_fd.close()
    prot_mesurements_fd.close()
    return ph_errs, pt_errs

AdjustedParametersOrder = namedtuple('AdjustedParameters', 'cam_0col, ph_0col, pt_0col, cols_n')


def GetAdjParamsOrder(unk, bi):
    cam_zero_col = 0
    ph_zero_col = bi.cams_n * len(unk.cam)
    pt_zero_col = bi.photos_n * len(unk.ph) + ph_zero_col

    col_n = ph_zero_col + pt_zero_col + bi.tie_pts_n * len(unk.pt)

    return cam_zero_col, ph_zero_col, pt_zero_col, col_n


# represents one least square iteration 
def RunBBAIteration(in_dt, unk, free_net):

    # build vectors L (measured photo points), L_XO (photo points computed from parameters)
    L, L_XO, bi = GetL(in_dt, free_net)    
    L_XO = np.array(L_XO)
    L = np.array(L)

    cam_step = len(unk.cam)
    ph_step = len(unk.ph)
    pt_step = len(unk.pt)

    apo = AdjustedParametersOrder(*GetAdjParamsOrder(unk, bi))

    A = np.zeros((len(L), apo.cols_n))
    X0 = np.zeros((apo.cols_n))


    # builds design matrix A and parameters vector X0 according 
    prev_cam = None
    prev_ph = None
    for i, d in enumerate(bi.Lrows_idxs):
        cam, ph, pt, v = d
        if prev_ph is None:
            i_ph_idx = 0
            prev_ph = ph
        elif prev_ph != ph:
            prev_ph = ph
            i_ph_idx += 1

        if prev_cam is None:
            i_cam_idx = 0
            prev_cam = cam
        elif prev_ph != ph:
            prev_cam = cam
            i_cam_idx += 1

        row_x = i * 2
        row_y = i * 2 + 1


        # fist columns of design matrix represent partial derivatives respect to
        # interior orientation parameters of cameras (if they are adjusted)
        cam_idx = apo.cam_0col + bi.cams_idxs[cam] * cam_step
        for i_un, un in enumerate(unk.cam):
            c = cam_idx + i_un
            A[row_x, c] = x_derivs[un](**v)
            A[row_y, c] = y_derivs[un](**v)

            # put values of parameter in corresponding position to design matrix column
            X0[c] = v[un]

        # in the middle, there are derivatives respect to  exterior orientations of photos   
        ph_idx = apo.ph_0col + bi.phs_idxs[ph] * ph_step
        for i_un, un in enumerate(unk.ph):
            c = ph_idx + i_un
            A[row_x, c] = x_derivs[un](**v)
            A[row_y, c] = y_derivs[un](**v)
            X0[c] = v[un]

        # last columns comprise derivatives respect to object coordinates of object points
        if not v['gcp']:
            pt_idx = apo.pt_0col + bi.pts_idxs[pt] * pt_step
            for i_un, un in enumerate(unk.pt):
                c = pt_idx + i_un
                A[row_x, c] = x_derivs[un](**v)
                A[row_y, c] = y_derivs[un](**v)
                X0[c] = v[un]



    l = L - L_XO
    
    # identity matrix is used for weights
    P = np.eye(A.shape[0])

    N = A.T.dot(P).dot(A)
    n = A.T.dot(P).dot(l)
    #Ni = np.linalg.inv(N)

    # computes pseudoinverse, which allows to deal with datum deficiency 
    # allowing to adjust free network (see ch. 1.4.4.) 
    Ni = np.linalg.pinv(N)
    dx = Ni.dot(n)

    Xa = X0 + dx
    e = A.dot(dx) - l

    # the data structures are updated with adjusted data
    AdjustData(Xa, bi, unk, apo)

    # compute reprojection error
    L1, L_XO1, bi = GetL(in_dt, free_net)
    e_repro = L_XO1 - L

    print "test"
    print np.amax(e - e_repro)

    # computes covariance matrices 
    r = A.shape[0] - A.shape[1]
    sigma0_2 = (e.T.dot(P).dot(e)) / r
    cov_x =  sigma0_2 * Ni
    cov_l= A.dot(cov_x).dot(A.T)

    L_XO1 = np.array(L_XO)
    L1 = np.array(L)    
    l = L1 - L_XO1
    return Xa, dx, e, e_repro, bi, apo, l, cov_x, cov_l

