from dataclasses import dataclass
import numpy as np


@dataclass
class IrregularRoundedPolygon():

    N = 20              # major sectors
    L = 4402.5          # major length (mm)

    n = (1, 4)          # minor sectors
    l = (3000, 1600)    # minor lengths (mm)

    R = 100e3           # major radius (mm)
    rho = 70.079e3      # minor radius (mm)

    def __post_init__(self):

        k = self.n[0]*(self.l[0]+self.L)/(self.n[1] *
                                          (self.l[1]+self.L) + self.n[0]*(self.l[0]+self.L))

        a = 2*np.pi/self.N
        b = a/self.n[0]*k
        c = a/self.n[1]*(1-k)

        dphi = np.array([0] + self.N*([b]*self.n[0] + [c]*self.n[1]))

        # magnet centers
        self.phi_m0 = b - np.cumsum(dphi)

        # drift centers
        self.phi_d = (self.phi_m0[0:-1] + self.phi_m0[1:])/2

        self.phi_m = self.phi_m0[1:]
        # dphi = dphi[1:]

    def radius(self, phi):
        """return reference radius given toroidal angle"""
        # print(R.shape)
        l = self.L
        R = self.R
        l1, l2 = self.l

        d = l/2/R
        d1 = l1/R
        # d2 = l2/R

        # validate generalized coordaintes to defined polar domain
        a = d + d1 - 2*np.pi
        b = d + d1
        phi = ((phi - a) % (b-a)) + a

        r = np.empty(phi.shape)
        r.fill(np.nan)
        ind = np.zeros(phi.shape, dtype=bool)

        phi_max = self.phi_m0[0:-1] - d
        phi_min = self.phi_m0[1:] + d

        # iterate over sectors
        for i in range(self.phi_m.size):
            print('sector: {:}'.format(i))

            theta = phi - self.phi_m[i]
            psi = phi - self.phi_d[i]

            # magnet boundaries
            ind_m = np.abs(theta) < d

            # drift boundaries
            ind_d = (phi <= phi_max[i]) & (phi >= phi_min[i])

            # drift section
            dn = phi_max[i]-phi_min[i]
            tmp = R*np.cos(dn/2)/np.cos(psi)
            r[ind_d] = tmp[ind_d]

            # bend section
            h = R*np.cos(d)-np.sqrt(R**2*np.cos(d)**2-R**2+self.rho**2)
            tmp = h*np.cos(theta)+np.sqrt(h**2*np.cos(theta)
                                          ** 2-h**2+self.rho**2)
            r[ind_m] = tmp[ind_m]

            ind += ind_d

        return {'R': r, 'drifts': ind}

@dataclass
class RegularRoundedPolygon:
    R: 100
    rho: 70
    n: 5
    dphi: float = 0
    verbose: bool = False
    # def rounded_polygon(R, rho, n, verbose=False, dphi=0):
    """return r(phi) for rounded polygon

    inputs:
        R: characteristic radius
        rho: corner radius
        n: number of sides

    optional:
        dphi: polar offset
        verbose: print statements

    outputs:
        r: function handle returning radius as a function of angle
        ind: drift section indices"""

    def __post__init__(self):

        alpha = 2*np.pi/n  # sector angle

        L = 2*np.pi*(R-rho)/n  # straight length
        l = 2*np.pi*rho/n  # bend length

        s = np.sqrt(2*rho**2*(1-np.cos(alpha)))  # bend cord length
        # cropped corner edge length
        b = np.abs(s*np.sin(alpha/2)/np.sin(alpha))

        h = (b+L/2)/np.tan(alpha/2)

        gamma = 2*np.arctan(L/(2*h))    # drift angle
        beta = alpha-gamma  # bend angle

        delta = (alpha-beta)/2
        k = rho*np.sin(delta)/np.sin(beta/2)

        if verbose:
            print('L: {:0.2f}'.format(L))
            print('l: {:0.2f}'.format(l))
            print('alpha: {:0.2f} deg'.format(np.degrees(alpha)))
            print('beta:  {:0.2f} deg'.format(np.degrees(beta)))
            print('gamma: {:0.2f} deg'.format(np.degrees(gamma)))

    def radius(phi):
        psi = (phi + dphi) % alpha  # subsector angle

        theta = psi-gamma/2
        vphi = psi-beta/2-gamma

        r_drift = h/np.cos(theta)
        r_bend = k*np.cos(vphi) + np.sqrt(k**2*np.cos(vphi)**2 - k**2+rho**2)

        r_drift[np.isnan(r_drift)] = 0
        r_bend[np.isnan(r_bend)] = 0
        r = r_drift*(psi <= gamma) + r_bend*(psi > gamma)

        return r, psi <= gamma
