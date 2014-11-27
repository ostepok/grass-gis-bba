import numpy as np
from math import pi

# File contains helper functions for priniting various protocols

# helper functions for creation of latex tables for the thesis

def GenerateLatexPointsIterationTable(diff, out_file):
    _genDiff(diff, out_file, "Pt. no.")

def GenerateLatexPhotosIterationTable(diff, out_file):
    _genDiff(diff, out_file, "Ph. no.")


def reprojection_errors(e_repro, prot_fd):


    fd = open(prot_fd, "w")


    e_repro = np.absolute(e_repro)
    e_repro = e_repro * (4592 / 25.1) # transformation to pixels  
    amax = np.amax(e_repro)
    amin = np.amin(e_repro)
    avg = np.average(e_repro)

    fd.write("max. & %8.6f \\\\ \hline \n" % amax)
    fd.write("min. & %8.6f \\\\ \hline \n" % amin)
    fd.write("avg.  & %8.6f \\\\ \hline \n" % avg)

    perc = [90, 75, 50, 25, 10] 
    for p in perc:
        pr = np.percentile(e_repro, p, axis=0)
        fd.write("perc. %d  & %8.6f \\\\ \hline \n" % (p, pr))

    fd.close()
    return p

def GetStats(diff):
    adiff = np.absolute(diff)
    avg = list(np.average(adiff, axis = 0))
    amax = list(np.amax(adiff, axis = 0))

    avg[0] = "avg."
    amax[0] = "max."

    stats = [avg, amax]

    return stats

def GenerateLatexPhotosAnglesIterationTable(diff, out_file):
    fd = open(out_file, "w")
    r, c = diff.shape 
    diff = diff[diff[:,0].argsort(axis=0)]

    fd.write("\\begin{center} \n")
    fd.write("\\rowcolors{1}{Gray}{white} \n")
    fd.write("\\begin{longtable}{| c || \n")

    ncols = c 
    for i in range(3): 
        for j in range( (ncols - 1) / 3):     
            if j == 0: 
                fd.write("c || ")
            else:
                fd.write("c | ")

    fd.write("} \n")

    fd.write("\hline \n")

    fd.write("Photo no. & ") 

    for i in range(3): 
        for j in range( (ncols - 1) / 3):     
            if j == 0: 
                fd.write("\\textbf{Bef. BBA}}")
            else:
                fd.write("\\textbf{It. %d}}" % (j))
            
            if i * (j-1) < 3 * (ncols - 1): 
                fd.write(" & ")

    fd.write("\\\\ \hline \n")

    new_diff = np.array(diff)

    width = (diff.shape[1] - 1)  / 3
    for i, col in enumerate(diff.T):
        if  i == 0:
            continue
        mod_col = (i + 2) % 3
        step = (i - 1) / 3 
        shift = mod_col * width 
        new_diff[:, shift + step + 1] = col 

    for i_row, row in enumerate(new_diff):
        mod_row = (i_row) % 3

        for i, col in enumerate(row):
            if i == 0: 
                fd.write("%d & " % col)
                continue
            
            fd.write("%6.3f " % (col * 180 / pi))
            
            if i < (len(row) - 1): 
                fd.write(" & ")

        fd.write("\\\\ \hline \n")

    stats = GetStats(new_diff)
    for i_row, row in enumerate(stats):
        mod_row = (i_row) % 3

        for i, col in enumerate(row):
            if i == 0: 
                fd.write("%s & " % col)
                continue
            
            fd.write("%6.3f " % (col * 180 / pi))
            
            if i < (len(row) - 1): 
                fd.write(" & ")

        fd.write("\\\\ \hline \n")

    fd.write("\hline \n")
    fd.write("\end{longtable}\n")
    fd.write("\end{center}\n")

    fd.close()

def _genDiff(diff, out_file, name):
    fd = open(out_file, "w")
    r, c = diff.shape 

    diff = diff[diff[:,0].argsort(axis=0)]

    fd.write("\\begin{longtable}{| c || \n")
 
    for i in range(c - 1):
        fd.write("c | ")

    fd.write("} \n")

    fd.write("\hline \n")
    for i in range(c):     
        if i == 0: 
            fd.write(name)        
        elif i == 1: 
            fd.write("Bef. I.")
        else:
            fd.write("I. %d " % (i - 1))
        if i != c - 1: 
                fd.write("&")

    fd.write("\\\\ \hline \n")

    for row in diff:
        for i, col in enumerate(row):
            if i == 0: 
                fd.write("%d " % col)
            else: 
                fd.write("%6.3f " % col)

            if i != c - 1: 
                fd.write(" & ")

        fd.write("\\\\ \hline \n")
    
    stats = GetStats(diff)
    for row in stats:
        for i, col in enumerate(row):
            if i == 0: 
                fd.write("%s " % col)
            else: 
                fd.write("%6.3f " % col)

            if i != len(row) - 1: 
                fd.write("&")

        fd.write("\\\\ \hline \n")


    fd.write("\hline \n")
    fd.write("\end{longtable}\n")
    fd.write("\end{center}\n")

    fd.close()

def CreateParametersLatex(i_num, cov, Xa, dx, bi, apo, prot_fd, unk, gcps, free_net, l):

    fd = open(prot_fd, "w")
    sigmas = np.sqrt(np.diagonal(cov))


    def _getDt(Xa, dx, sigmas, unk_l, x_col, pt):
        
        dt = [pt.GetId()]
        for i, un in enumerate(unk_l):
            val = Xa[x_col + i]
            s = sigmas[x_col + i]

            if un in ['ph', 'om', 'ka']:
                val = val * 180 / pi
                s = s  * 180 / pi

            dt.append(val)
            dt.append(s)

        return dt

    #diff = diff[numpy.lexsort(diff[:,0])]

    dts = []

    pt_step = len(unk.pt)
    for pt_idx in range(bi.tie_pts_n):

        pt = bi.idxs_pts[pt_idx]
        pt_x_col = apo.pt_0col + pt_idx * pt_step


        dts.append(_getDt(Xa, dx, sigmas, unk.pt, pt_x_col, pt))
    
    dts = np.array(dts)
    dts = dts[dts[:,0].argsort(axis=0)]

    for d in dts:
        for i, c in enumerate(d):
            if i == 0: 
                fd.write("%d " % (c))
            else:
                fd.write("%10.3f " % (c))

            if i < len(d) - 1:
                fd.write(" & ")

        fd.write(" \\\\ \hline \n")

    fd.close()

def CreatePhootoParametersLatex(i_num, cov, Xa, dx, bi, apo, prot_fd, unk, gcps, free_net, l):

    fd = open(prot_fd, "w")
    sigmas = np.sqrt(np.diagonal(cov))


    def _printFeature(Xa, dx, sigmas, unk_l, x_col, ph):
        fd.write("%d &" % (ph.GetId()))
        for i, un in enumerate(unk_l):
            val = Xa[x_col + i]
            s = sigmas[x_col + i]

            if un in ['ph', 'om', 'ka']:
                val = val * 180 / pi
                s = s  * 180 / pi
            fd.write("%10.3f & %10.3f " % (val, s))

            if i < len(unk_l) - 1:
                fd.write(" & ")

        fd.write(" \\\\ \hline \n")
    ph_step = len(unk.ph)


    for ph_idx in range(bi.photos_n):
        ph = bi.idxs_phs[ph_idx]
        ph_x_col = apo.ph_0col + ph_idx * ph_step

        _printFeature(Xa, dx, sigmas, unk.ph, ph_x_col, ph)


    fd.close()


# functions for creation of adjustment protocols

def CreateIterationProtocolMeasurements(i_num, cov_l, e_repro, bi, prot_fd):
    
    prot_fd.write(_('\n\n\n\n*********************************************************************'))
    prot_fd.write(_('\n\nIteration: %d\n\n' % (i_num + 1)))

    sigmas = np.sqrt(np.diagonal(cov_l))


    prev_ph = None
    for i, d in enumerate(bi.Lrows_idxs):
        cam, ph, pt, v = d
        i1 = i * 2
        i2 = i * 2 + 1
        if prev_ph is None or prev_ph != ph:
            prev_ph = ph
            prot_fd.write("\n\nPhoto %10d:\n" % (ph.GetId()))


        phpt = pt.GetPhotoPoints()[ph]
        prot_fd.write("%10d: %10.4f %10.4f " % (pt.GetId(), phpt[0], phpt[1], ))

        prot_fd.write("%10.4f %10.4f" % (sigmas[i1], sigmas[i2]))
        prot_fd.write("%10.4f %10.4f" % (e_repro[i1], e_repro[i2]))

        prot_fd.write("\n")



def CreateParametersIterationProtocol(i_num, cov, Xa, dx, bi, apo, prot_fd, unk, gcps, free_net, l):

    prot_fd.write(_('\n\n\n\n*********************************************************************'))
    prot_fd.write(_('\n\nIteration: %d\n\n' % (i_num + 1)))

    prot_fd.write(_('Ground control point differences:\n'))

    sigmas = np.sqrt(np.diagonal(cov))
    #sigmas = None
    prot_fd.write("%10s, %10s\n" % ("gcp id", "diff"))
    for gcp in gcps.itervalues():
        if gcp.GetCoords() is None or ( not gcp.GetGcp()[1] and not free_net):
            continue

        dist = np.linalg.norm(gcp.GetGcp()[0] - gcp.GetCoords())
        bingo_dist = np.linalg.norm(gcp.GetGcp()[0] - gcp.GetResultCoords())

        prot_fd.write("%10d, %10.4f" % (gcp.GetId(), dist))
        prot_fd.write("%10d, %10.4f, %10.4f" % (gcp.GetId(), dist, bingo_dist))
        prot_fd.write("\n")

    prot_fd.write("\n\n")

    cam_step = len(unk.cam)
    ph_step = len(unk.ph)
    pt_step = len(unk.pt)

    def _printParamsLine(x, unk_l, x_col):

        for i, un in enumerate(unk_l):
            val = x[x_col + i]
            if un in ['ph', 'om', 'ka']:
                val = val * 200 / pi

            prot_fd.write("%15.4f" % (val))

    def _printParamCaptions(unk_l):
            for i, un in enumerate(unk_l):
                prot_fd.write("%15s" % (un))
            prot_fd.write("\n\n")

    def _printFeature(Xa, dx, sigmas, unk_l, x_col):
        _printParamsLine(Xa, unk_l, x_col)
        prot_fd.write("\n")
        _printParamsLine(dx, unk_l, x_col)
        prot_fd.write("\n")
        _printParamsLine(sigmas, unk_l, x_col)
        prot_fd.write("\n\n")

    prot_fd.write(_('Interior orientaitons adjustment results\n'))

    for cam_idx in range(bi.cams_n):
        if cam_idx == 0:
            _printParamCaptions(unk.cam)

        cam = bi.idxs_phs[cam_idx]
        prot_fd.write(_("Camera %d\n" % cam.GetId()))
        cam_x_col = apo.cam_0col + cam_idx * cam_step

        _printFeature(Xa, dx, sigmas, unk.cam, cam_x_col)

    prot_fd.write(_('Exterior orientaitons adjustment results\n'))

    for ph_idx in range(bi.photos_n):
        if ph_idx == 0:
            _printParamCaptions(unk.ph)

        ph = bi.idxs_phs[ph_idx]
        prot_fd.write(_("Photo %d\n" % ph.GetId()))
        ph_x_col = apo.ph_0col + ph_idx * ph_step

        _printFeature(Xa, dx, sigmas, unk.ph, ph_x_col)

    prot_fd.write(_('Object points orientaitons adjustment results\n'))

    for pt_idx in range(bi.tie_pts_n):
        if pt_idx == 0:
            _printParamCaptions(unk.pt)

        pt = bi.idxs_pts[pt_idx]
        prot_fd.write(_("Point %d\n" % pt.GetId()))
        pt_x_col = apo.pt_0col + pt_idx * pt_step

        _printFeature(Xa, dx, sigmas, unk.pt, pt_x_col)