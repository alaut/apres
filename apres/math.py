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
