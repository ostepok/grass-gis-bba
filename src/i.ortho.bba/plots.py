import os
import sys
import numpy as np
from copy import deepcopy

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


# helper functions for plotting graphs for the thesis

def plotPts(pts_list, lines_list):
    fig = plt.figure()

    i_col = 1
    cols = len(pts_list)
    for i, pts in  enumerate(pts_list):
        if pts is not None:
            ax = fig.add_subplot(1, cols, i_col, projection='3d')
            s = 40
            if pts.shape[1] < 6:
                ax.scatter(pts[:,0], pts[:,1], pts[:,2], c='g', s=s, marker='^')
            else:
                ax.scatter(pts[:,0], pts[:,1], pts[:,2], c=pts[:,3:6], s=s, marker='^')

            i_col += 1

        if lines_list[i] is not None:
            for lns in lines_list[i]:
                for ln in lns:

                    ln = np.reshape(ln, (-1, 3))
                    ax.plot(ln[:,0], ln[:,1], ln[:,2], lw=2)

    plt.show()


def PlotResults(pts, phs):
    plot_pts = []
    for pt in pts.itervalues():
        if pt.GetResultCoords() is not None:
            plot_pts.append(pt.GetResultCoords())

    lines = []
    for ph in phs.itervalues():
         
        if  ph.GetResultExtOr() is not None:
            eo = ph.GetResultExtOr()
            lines.append(GetPhotoPyramids(ph, 18))

    plot_pts = [np.array(plot_pts)]


    plotPts(plot_pts, lines)


def PlotRelOrs(pts, phs, ro):

    plot_pts_result = []
    plot_pts_ro = []

    inv_coord = 0

    # quick hack - when no result, put prev point in order to fit to numpy array and prev point is used in roder to be close to other points
    prev_pt = None

    for pt in pts.itervalues():
        if pt.GetRelOrsCoords().has_key(ro):
            
            if pt.GetResultCoords() is None:
                plot_pts_result.append(prev_pt)
            else:
                plot_pts_result.append(pt.GetResultCoords())
                prev_pt =  pt.GetResultCoords()
            
            pt_coords = deepcopy(pt.GetRelOrsCoords()[ro])
            pt_coords[inv_coord] = pt_coords[inv_coord] * -1
            plot_pts_ro.append(pt_coords)
            

    eo_result = []
    eo_ro = []

    for ph in phs.itervalues():
         #ph.GetResultExtOr() is not None

            #res_eo = ph.GetResultExtOr()
            #res_pts.append(res_eo[:, 3])
            pyr = GetPhotoPyramids(ph, 18)
            for ln in pyr:
                ln = np.reshape(ln, (-1, 3))
                ln[:,inv_coord] = -1*ln[:,inv_coord]

            eo_ro.append(pyr)

            #r = res_eo[:, :3]
            #eo_res = np.hstack((r, np.array([ph.GetEO()[:,3]]).T))
            #eo_res = np.array([ph.GetEO()[:,3]]).T
            #lines.append(GetCameraLines(eo_res, inv=True))
            
            eo_result.append(GetPhotoPyramids(ph, 7, ph.GetResultExtOr(), True))

    plot_pts_result = np.array(plot_pts_result)

    print plot_pts_result.shape

    plot_pts_ro = np.array(plot_pts_ro)
    plot_list = [plot_pts_result, plot_pts_ro]

    lines = [eo_result, eo_ro]


    #if ro:
    #    plot_list.append(ro)


    plotPts(plot_list, lines)


def PlotAllRelOrs(pts, phs):

    plot_pts_result = []
    plot_pts_ro = []

    inv_coord = 0
    for pt in pts.itervalues():
        if pt.GetRelOrsCoords():
            c = np.random.rand(3,1)
            plot_pts_result.append(pt.GetResultCoords())
            for ro, coords in deepcopy(pt.GetRelOrsCoords()).iteritems():
                pt_coords = np.append(coords, c)
                pt_coords[inv_coord] = pt_coords[inv_coord] * -1
                plot_pts_ro.append(pt_coords)



    eo_result = []
    eo_ro = []

    for ph in phs.itervalues():
        if ph.GetEO() is not None:
         #ph.GetResultExtOr() is not None

            #res_eo = ph.GetResultExtOr()
            #res_pts.append(res_eo[:, 3])            
            pyr = GetPhotoPyramids(ph, 18)
            for ln in pyr:
                ln = np.reshape(ln, (-1, 3))
                ln[:,inv_coord] = -1*ln[:,inv_coord]

            eo_ro.append(pyr)

            #r = res_eo[:, :3]
            #eo_res = np.hstack((r, np.array([ph.GetEO()[:,3]]).T))
            #eo_res = np.array([ph.GetEO()[:,3]]).T
            #lines.append(GetCameraLines(eo_res, inv=True))
            
            eo_result.append(GetPhotoPyramids(ph, 7, ph.GetResultExtOr(), True))



    plot_pts_result = np.array(plot_pts_result)
    plot_pts_ro = np.array(plot_pts_ro)
    plot_list = [plot_pts_ro]#,  plot_pts_result]
    lines = [eo_ro]#, eo_result]

    #if ro:
    #    plot_list.append(ro)
    plotPts(plot_list, lines)


def PlotScene(pts, phs):
    #PlotResults(pts, phs)
    #return

    #PlotANgles(pts, phs)
    #return
    plot_pts = []
    for pt in pts.itervalues():
         #pt.GetResultCoords() is not None and 
        if pt.GetGcp()[0] is not None:
            #res_pts.append(pt.GetResultCoords())
            plot_pts.append(pt.GetGcp()[0])

    res_pts = []
    for pt in pts.itervalues():
         #pt.GetResultCoords() is not None and 
        if pt.GetCoords() is not None:
            #res_pts.append(pt.GetResultCoords())
            plot_pts.append(pt.GetCoords())

    lines = []
    for ph in phs.itervalues():
         
        if ph.GetEO() is not None and ph.GetResultExtOr() is not None:
            res_eo = ph.GetResultExtOr()
            res_pts.append(res_eo[:, 3])
            plot_pts.append(ph.GetEO()[:,3])
            #lines.append(GetCameraLines(ph.GetEO()))

            lines.append(GetPhotoPyramids(ph, 18))

            r = res_eo[:, :3]
            eo_res = np.hstack((r, np.array([ph.GetEO()[:,3]]).T))
            #lines.append(GetCameraLines(eo_res, inv=True))
            
            #lines.append(GetCameraLines(ph.GetEO()))
    plot_pts = np.array(plot_pts)
    res_pts = np.array(res_pts)

    plot_list = [plot_pts]
    #if ro:
    #    plot_list.append(ro)

    plotPts(plot_list, [lines])

def PlotANgles(pts, phs):

    plot_pts = []
    for pt in pts.itervalues():
         #pt.GetResultCoords() is not None and 
        if pt.GetGcp()[0] is not None:
            #res_pts.append(pt.GetResultCoords())
            plot_pts.append(pt.GetGcp()[0])

    res_pts = []
    for pt in pts.itervalues():
         #pt.GetResultCoords() is not None and 
        if pt.GetCoords() is not None:
            #res_pts.append(pt.GetResultCoords())
            plot_pts.append(pt.GetCoords())

    lines = []
    for ph in phs.itervalues():
         
        if ph.GetEO() is not None and ph.GetResultExtOr() is not None:
            res_eo = ph.GetResultExtOr()
            res_pts.append(res_eo[:, 3])
            plot_pts.append(ph.GetEO()[:,3])
            #lines.append(GetCameraLines(ph.GetEO()))

            #lines.append(GetPhotoPyramids(ph, 18))

            r = res_eo[:, :3]
            eo_res = np.hstack((r, np.array([ph.GetEO()[:,3]]).T))
            lines.append(GetCameraLines(eo_res, inv=True))
            
            lines.append(GetCameraLines(ph.GetEO()))
    plot_pts = np.array(plot_pts)
    res_pts = np.array(res_pts)

    plot_list = [plot_pts]
    #if ro:
    #    plot_list.append(ro)

    plotPts(plot_list, [lines])

def GetCameraLines(eo, inv=False):

    size = 1
    apts = np.array([
                  [size, 0, 0],
                  [0, size, 0],
                  [0, 0, -size]])
    origin_mat = np.tile(eo[:,3], (apts.shape[0] , 1))
    axis_pts = (eo[:,:3].T.dot(apts.T)).T + origin_mat
    

    if not inv:
        apts = np.array([
                      [size, 0, 0],
                      [0, size, 0],
                      [0, 0, size]])
    
        origin_mat = np.tile(eo[:,3], (apts.shape[0] , 1))
        axis_pts = (eo[:,:3].T.dot(apts.T)).T + origin_mat
        
    """
    if not inv:
        dz = axis_pts[:, 2] - eo[2,3]
        axis_pts[:, 2]  = eo[2,3] - dz
    """

    return  np.hstack((origin_mat, axis_pts))

def GetPhotoPyramids(ph, scale, eo=None, inv=False):
    cam = ph.GetCamera()
    if eo is None:
        eo = ph.GetEO()

    chs_x, chs_y = cam.GetChipSize() 

    chs_x /=  (2.0 * scale)
    chs_y /=  (2.0 * scale)

    f = cam.GetParams()[0][0, 0]
    f /= (2.0 * scale)
    if inv:
        f = f * -1

    apts = [np.array(
                    [[chs_x, chs_y, -f], 
                     [chs_x, -chs_y, -f], 
                     [-chs_x, -chs_y, -f],
                     [-chs_x, chs_y, -f],  
                     [chs_x, chs_y, -f]]
                    ),
            np.array([[chs_x, chs_y, -f], [0, 0, 0]]),
            np.array([[-chs_x, chs_y, -f], [0, 0, 0]]),
            np.array([[-chs_x, -chs_y, -f], [0, 0, 0]]),
            np.array([[chs_x, -chs_y, -f], [0, 0, 0]])
            ]


    for i, ln in enumerate(apts):
        origin_mat = np.tile(eo[:,3], (ln.shape[0] , 1))
        ln = eo[:,:3].T.dot(ln.T).T + origin_mat
        apts[i] = ln.reshape((1, -1))

    """
    if not inv:
        dz = axis_pts[:, 2] - eo[2,3]
        axis_pts[:, 2]  = eo[2,3] - dz
    """

    return  apts