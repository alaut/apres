from dataclasses import dataclass
import numpy as np
from stl import mesh

import matplotlib.pyplot as plt

from apres.math import cartesian2curvilinear
from apres.algorithms import compute_intersections
from apres.plotting import plot_projections


@dataclass
class Aperture:

    files: list

    reference: float = 0
    origin: tuple = (0, 0, 0)

    nr: int = 1
    nz: int = 100
    nth: int = 4
    rmax: float = 1000

    intersections = {}
    tri = {}

    def __post_init__(self):
        self.load_stl()

    def load_stl(self):
        """load stl files"""
        for file in self.files:
            print('Loading {:} ...'.format(file))
            self.tri[file] = mesh.Mesh.from_file(file)

    def show(self, KEYS=[('x', 'y')]):

        if self.reference is not None:
            phi = np.linspace(0, 2*np.pi, 100_000)

            r = {
                'x': self.reference(phi)*np.cos(phi),
                'y': self.reference(phi)*np.sin(phi),
                'z': np.zeros(phi.shape),
            }

        for k1, k2 in KEYS:
            fig, ax = plt.subplots(constrained_layout=True)
            ax.set_aspect('equal')
            ax.set_xlabel(k1)
            ax.set_ylabel(k2)
            for name, tri in self.tri.items():
                x1 = getattr(tri, k1).ravel()
                x2 = getattr(tri, k2).ravel()
                ax.plot(x1, x2, ',', label=name)
                # ax.annotate(text=name, xy=(np.mean(x1), np.mean(x2)))

            if self.reference is not None:
                ax.plot(r[k1], r[k2], 'm-')

    def linearize(self):

        for file, tri in self.tri.items():

            print('linearizing {:} ...'.format(file))
            tri.vectors -= np.array(self.origin)
            data = cartesian2curvilinear(
                x=tri.x.T,
                y=tri.y.T,
                z=tri.z.T,
                R_0=self.reference,
                ccw=-1,
                dphi=0,
            )

            tri.x, tri.y, tri.z = data['h'].T, data['v'].T, data['s'].T
            self.tri[file] = tri

        self.reference = None

    def inflate(self):

        for file, tri in self.tri.items():

            p = np.stack([tri.x.T, tri.y.T, tri.z.T], axis=1)
            print('p: {:}\t{:}'.format(p.shape, file))

            q = self.define_rays(p)

            r = compute_intersections(q, p, step=1000)

            np.save(f"{file}", r)
            self.intersections[file] = r

    def plot_aperture(self):
        keys = [('z', 'x'), ('z', 'y')]

        fig, axes = plt.subplots(
            len(keys), 1, constrained_layout=True, num='aperture', sharex=True)

        for file, tri in self.tri.items():
            r = self.intersections[file]

            hits = {'x': r[0, 0].T, 'y': r[1, 0].T, 'z': r[2, 0].T}
            data = {'x': tri.x.T, 'y': tri.y.T, 'z': tri.z.T}

            plot_projections(data, keys=keys, axes=axes,
                             style='mesh', equal=True)
            plot_projections(hits, keys=keys, axes=axes,
                             style='line', spec='-')
            # axes[0].set_title(file)

    def define_rays(self, p):
        """given faces, define rays to test for intersections"""
        xdata, ydata, zdata = p[:, 0], p[:, 1], p[:, 2]

        r = np.linspace(0, self.rmax, self.nr+1)
        th = np.linspace(0, 2*np.pi, self.nth, endpoint=False)
        zq = np.linspace(zdata.min()-1, zdata.max()+1, self.nz)
        R, TH, Z = np.meshgrid(r, th, zq, indexing='ij')

        X = R*np.cos(TH)
        Y = R*np.sin(TH)

        q1 = np.stack((X[:-1, :, :], Y[:-1, :, :], Z[:-1, :, :]))
        q2 = np.stack((X[1::, :, :], Y[1::, :, :], Z[1::, :, :]))
        q = np.stack((q1, q2))
        print('q: {:}'.format(q.shape))

        return q
