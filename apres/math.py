import numpy as np


def curvilinear2cartesian(data, R_0):
    """convert toroidal coordaintes to curvilinear projections"""
    data['s'] = R_0(data['phi'])*data['phi']
    data['x'] = data['r']*np.cos(data['theta'])
    data['y'] = data['r']*np.sin(data['theta'])


def cartesian2curvilinear(x, y, z, R_0=0, ccw=1, dphi=0):
    """return toroidal/curvilinear coordinates from cartesian coordinates
    https://en.wikipedia.org/wiki/Toroidal_and_poloidal_coordinates"""
    R = np.sqrt(x**2+y**2)
    phi = np.arctan2(y, x)

    if callable(R_0):
        R_0 = R_0(phi)

    r = np.sqrt(z**2+(R-R_0)**2)
    theta = np.arctan2(z, R-R_0)

    th = ccw*(phi-dphi)
    th[th < 0] += 2*np.pi

    s = R_0*th
    h = r*np.cos(theta)
    v = r*np.sin(theta)

    return {'r': r, 'theta': theta, 'phi': phi, 's': s, 'h': h, 'v': v}

def scalar_triple_product(a, b, c, verbose=False):
    """return the scalar triple product of vectors a, b, and c"""
    if verbose:
        print(f'stp ~ {a.shape} * ( {b.shape} x {c.shape} )')

    # classic approach
    base = np.cross(b, c, axisa=0, axisb=0, axis=0)
    stp = np.sum(a*base, axis=0)
    if verbose:
        print(f'stp: {stp.nbytes/1e6} MB')

    return stp


def tetra_volume_sign(a, b, c, d, verbose=False):
    """return positivity of tetrahedra volume defined by vertices a, b, c, d
    https://en.wikipedia.org/wiki/Tetrahedron"""
    sign = np.sign(scalar_triple_product(a-d, b-d, c-d)) == 1
    if verbose:
        print(f'sign: {sign.nbytes/1e6} MB')
    return sign
