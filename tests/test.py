import matplotlib.pyplot as plt
from glob import glob
import numpy as np

from apres.apres import Aperture
from apres.shapes import IrregularRoundedPolygon

ps = Aperture(
    files=glob('./tests/models/*.stl'),
    # reference=100e3,
    reference=IrregularRoundedPolygon().radius,
    origin=(0, 0, 1260),
    nz=500,
    nth=16,
)

ps.show(
    keys=[('y', 'x')],
    phi=np.radians(np.linspace(3, -9, 100)),
)

ps.dphi = np.radians(3)
ps.linearize()

ps.show(
    keys=[('x', 'y'), ('z', 'x'), ('z', 'y')],
    equal=False,
)

ps.inflate()
ps.plot_aperture(
    keys=[('z', 'x'), ('z', 'y'), ('x', 'y')],
)

plt.show()
print('done!')
