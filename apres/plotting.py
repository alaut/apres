from matplotlib.tri import Triangulation
import numpy as np
import matplotlib.pyplot as plt


def plot_projections(data, keys, axes=None, equal=False, style='mesh', spec='k,'):
    """plot dataset projections given data dict and list of key pairs"""

    if axes is None:
        axes = []
        for (k1, k2) in keys:
            fig, ax = plt.subplots(
                num=f'projection: {k1}-{k2}', constrained_layout=True)
            ax.set_xlabel(k1)
            ax.set_ylabel(k2)
            axes.append(ax)

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
