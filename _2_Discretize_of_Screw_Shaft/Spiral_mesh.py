# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import sys
import pathlib
sys.path.append(str(pathlib.Path(__file__).resolve().parent.parent))
from _9_Mechanics.Spiral import Spiral

R = 0.05
SP0 = Spiral(alp=0e-4, l=0.6*np.pi*2, r=0.8, eta=np.array([0.00141421356, 0.00141421356]), R=0.052)

th = np.deg2rad(np.linspace(-0.05, 0.05, int(1e1)))
alp = np.deg2rad(np.linspace(44, 46, int(1e1)))

xyz = SP0.get_mesh(th, alp)
x = xyz[:,:,0]
y = xyz[:,:,1]

fig = plt.figure()
ax = fig.add_subplot(111)
plt.scatter(x, y)
ax.set_aspect('equal', adjustable='box')

segs1 = np.stack((x,y), axis=2)
segs2 = segs1.transpose(1,0,2)
plt.gca().add_collection(LineCollection(segs1))
plt.gca().add_collection(LineCollection(segs2))
plt.show()

