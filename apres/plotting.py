import matplotlib.pyplot as plt
import numpy as np

from matplotlib.tri import Triangulation

def plot_projections(data, keys, num=None, axes=None, equal=False, style='mesh', spec='k,'):
    """construct subplots from data dictionary and list of key tuples

    inputs:
        data: dictionary of plotable mesh vertices
        keys: list of tuples describing axes bases

    optional:
        num: None
        axes: None
        equal: False
        style: 'mesh

    ouputs:
        axes:"""

    if axes is None:
        if True:
            fig, axes = plt.subplots(
                len(keys), 1, num=num, constrained_layout=True)
            if len(keys) == 1:
                axes = [axes]
        else:
            axes = []
            for (k1, k2) in keys:
                fig, ax = plt.subplots(
                    num='{:}-{:}-{:}'.format(k1, k2, num), constrained_layout=True)
                axes.append(ax)

        for ax, (k1, k2) in zip(axes, keys):
            ax.set_xlabel(k1)
            ax.set_ylabel(k2)
            if equal:
                ax.set_aspect('equal')

    for ax, (k1, k2) in zip(axes, keys):
        if style == 'points':
            ax.plot(data[k1].ravel(), data[k2].ravel(), spec)
        elif style == 'mesh':
            import matplotlib as mpl
            mpl.rcParams['agg.path.chunksize'] = 10000

            con = np.reshape(range(data[k1].size), (data[k1].shape[1], 3))
            x = data[k1].ravel(order='F')
            y = data[k2].ravel(order='F')
            tri = Triangulation(x, y, con)
            ax.triplot(tri, spec, lw=0.25)
        elif style == "line":
            ax.plot(data[k1], data[k2], spec)
        else:
            raise NameError('unrecognized plotting style')

    return axes