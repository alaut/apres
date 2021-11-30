import numpy as np

def define_rays(p, rmax=1000, n_r=1, n_th=4, n_z=100):
    """given faces, define rays to test for intersections"""
    xdata, ydata, zdata = p[:, 0], p[:, 1], p[:, 2]

    r = np.linspace(0, rmax, n_r+1)
    th = np.linspace(0, 2*np.pi, n_th, endpoint=False)
    zq = np.linspace(zdata.min()-1, zdata.max()+1, n_z)
    R, TH, Z = np.meshgrid(r, th, zq, indexing='ij')

    X = R*np.cos(TH)
    Y = R*np.sin(TH)

    q1 = np.stack((X[:-1, :, :], Y[:-1, :, :], Z[:-1, :, :]))
    q2 = np.stack((X[1::, :, :], Y[1::, :, :], Z[1::, :, :]))
    q = np.stack((q1, q2))
    print('q: {:}'.format(q.shape))

    return q