#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: FlyFox
@email: walkandthinker@gmail.com
@discuss: QQ group: 797998860
"""

from Mesh import *

mesh2d=Mesh2D(20,20,0.0,2.0,0.0,2.0,'quad4')
mesh2d.CreateMesh()
mesh2d.SplitBCMesh()
mesh2d.PlotMesh()