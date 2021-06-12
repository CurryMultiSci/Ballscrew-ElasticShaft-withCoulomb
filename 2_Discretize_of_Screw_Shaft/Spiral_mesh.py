# -*- coding: utf-8 -*-

# from . import Spiral

import numpy as np
from Spiral import Spiral
import matplotlib.pyplot as plt

R = 0.05
SP0 = Spiral(alp=0e-4, l=0.6*np.pi*2, r=0.8, eta=np.array([0.00141421356, 0.00141421356]), R=0.052)


th = np.deg2rad(np.linspace(-0.01, 0.01, int(1e1)))
alp = np.deg2rad(np.linspace(44, 46, int(1e1)))

xyz = SP0.get_mesh(th, alp)
