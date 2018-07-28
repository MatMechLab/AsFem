#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: FlyFox
@email: walkandthinker@gmail.com
@discuss: QQ group: 797998860
"""

from Mesh import *

mesh1d=Mesh1D(10,0.0,1.0,1)
mesh1d.CreateMesh()
# mesh1d.PlotMesh()

test=np.zeros((6,90))
