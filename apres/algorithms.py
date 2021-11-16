import numpy as np
import time

from apres.math import tetra_volume_sign

def compute_intersections(q, p, step=200):
    """subdivide and compute intersection shell"""
    dim_rays, dim_faces = q.shape[1:], q.shape[-2:]
    r = np.empty(dim_rays)

    ind = np.indices(dim_faces)
    I, J = ind[0].ravel(), ind[1].ravel()

    n = len(I)

    for i in range(0, n, step):
        print('Completed {:0.2f} %\n'.format(i/n*100))
        i0 = I[i:i+step]
        j0 = J[i:i+step]
        q_sub = q[:, :, :, i0, j0]
        p_sub = p
        r_tmp, idx = intersect_line_triangle(q_sub, p_sub, tails=False, dtype=np.float32)
        r[:, :, i0, j0] = r_tmp

    return r

def intersect_line_triangle(q, p, tails=True, plot=False, dtype=None, step=200, fast=True):
    """calculate intersection between n rays (q) and m faces (p) using Moller-Trumbore algorithm
    inputs:
        q: [2 x 3 x n1 x n2 x ...] rays where nr = prod(n1, n2, ...)
        p: [3 x 3 x nt] faces

    outputs:
        r: [3 x nt x n1 x n2 x ...] intersections
        ind: [nt x n1 x n2 x ...] collision indices
    optional:
        tails: exclude additional intersections of a given ray
        dtype: define data dtype to optimize memory consumption
    """
    ti = time.time()

    # Reshape to Broadcast
    dim_p = p.shape
    dim_q = q.shape
    nt = p.shape[2]
    nr = np.prod(q.shape[2:])
    print('{:} rays upon {:} faces'.format(nr, nt))

    q = np.reshape(q, (2, 3, 1, nr)).astype(dtype)
    p = np.reshape(p, (3, 3, nt, 1)).astype(dtype)

    print('q: {:}\t{:0.2f} MB'.format(q.shape, q.nbytes/1e6))
    print('p: {:}\t{:0.2f} MB'.format(p.shape, p.nbytes/1e6))

    print('Computing tetrahedra volumes...')
    s1 = tetra_volume_sign(q[0], p[0], p[1], p[2])
    s2 = tetra_volume_sign(q[1], p[0], p[1], p[2])

    # Identify Intersection possibilities
    ind = (s1 != s2)
    print('Identified {:} possibilities from {:} combinations!'.format(
        ind.sum(), ind.size))

    # Verify Intersection Subspace
    if fast:
        idx = np.argwhere(ind)
        I, J = idx[:, 0], idx[:, 1]
        s3 = tetra_volume_sign(
            q[0][:, 0, J], q[1][:, 0, J], p[0][:, I, 0], p[1][:, I, 0])
        s4 = tetra_volume_sign(
            q[0][:, 0, J], q[1][:, 0, J], p[1][:, I, 0], p[2][:, I, 0])
        s5 = tetra_volume_sign(
            q[0][:, 0, J], q[1][:, 0, J], p[2][:, I, 0], p[0][:, I, 0])
        ind[I, J] = (s3 == s4) & (s4 == s5)
    else:
        s3 = tetra_volume_sign(q[0], q[1], p[0], p[1])
        s4 = tetra_volume_sign(q[0], q[1], p[1], p[2])
        s5 = tetra_volume_sign(q[0], q[1], p[2], p[0])
        ind *= (s3 == s4) & (s4 == s5)
    print('Verified {:} intersections!'.format(ind.sum()))

    # Compute Intersection Properties
    if fast:
        idx = np.argwhere(ind)
        I, J = idx[:, 0], idx[:, 1]
        n = np.cross(p[1][:, I, 0]-p[0][:, I, 0], p[2][:, I, 0] -
                     p[0][:, I, 0], axisa=0, axisb=0, axis=0)
        t = np.sum((p[0][:, I, 0]-q[0][:, 0, J])*n, axis=0) / \
            np.sum((q[1][:, 0, J]-q[0][:, 0, J])*n, axis=0)
        r = q[0][:, 0, J]+t*(q[1][:, 0, J]-q[0][:, 0, J])
    else:
        with np.errstate(divide='ignore', invalid='ignore'):
            n = np.cross(p[1]-p[0], p[2]-p[0], axisa=0, axisb=0, axis=0)
            t = np.sum((p[0]-q[0])*n, axis=0)/np.sum((q[1]-q[0])*n, axis=0)
            r = q[0] + t * (q[1]-q[0])
    print('r: {:}\t{:0.2f} MB'.format(r.shape, r.nbytes/1e6))

    # Cluster intersections by rays
    cluster = cluster_intersections(ind)

    R = np.empty(dim_q[1:])
    R.fill(np.nan)
    if not tails:
        for ray in cluster:
            i, faces = cluster[ray]
            # i = cluster[ray]
            ind_t = np.argsort(t[i])  # mark for removal by parametric extent

            faces_rem = np.array(faces)[ind_t][1:]
            idx_keep = np.array(i)[ind_t][0]

            # r0 = r[:,idx_keep]

            ind[faces_rem, ray] = False
            # r[:, idx_rem] = np.nan

            # convert linear ray index to ray matrix index
            idx = np.unravel_index(ray, dim_q[2:])
            for j in range(3):
                R[j][idx] = r[j, idx_keep]
        print('Refined to {:} intersections!'.format(ind.sum()))

    # Exclude invalid intersections
    idx = np.argwhere(ind)
    print('Refined r: {:}\t{:0.2f} MB'.format(r.shape, r.nbytes/1e6))

    tf = time.time()
    print('Finished in {:0.2f} seconds!'.format(tf-ti))

    return R, idx  # return intersections and face indices


def cluster_intersections(ind):
    """cluster (face, ray) intersections by hashable ray indices

    inputs:
        ind: face x ray intersection indices

    ouputs:
        cluster[ray] = list of ([idx], [face])"""
    cluster = {}
    for i, (face, ray) in enumerate(np.argwhere(ind)):
        try:
            cluster[ray][0].append(i)
            cluster[ray][1].append(face)
        except:
            cluster[ray] = ([i], [face])
    print('Clustered {:} rays'.format(len(cluster)))
    return cluster
