[mesh]
  type=gmsh
  file=grains50.msh
[end]

[dofs]
name=c d ux uy
[end]


[nonlinearsolver]
  type=nr
  maxiters=25
  r_rel_tol=1.0e-10
  r_abs_tol=5.0e-7
[end]

[ics]
  [constd]
    type=const
    dof=c
    params=0.0
  [end]
[end]

[output]
  type=vtu
  interval=2
[end]

[timestepping]
  type=be
  dt=5.0e-4
  time=1.0e2
  adaptive=true
  optiters=4
  growthfactor=1.2
  cutfactor=0.85
  dtmax=1.0e-1
[end]

[projection]
scalarmate=vonMises
[end]


[bcs]
  [fixux]
    type=dirichlet
    dofs=ux uy
    value=0.0
    boundary=left
  [end]
  [fixuy]
    type=dirichlet
    dofs=ux
    value=0.0
    boundary=right
  [end]
  [flux]
    type=neumann
    dofs=c
    value=-0.2
    boundary=bottom top
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]

//************************************************
[elmts]
  [myelmt1]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate1
    domain=1
  [end]
  [myelmt2]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate2
    domain=2
  [end]
  [myelmt3]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate3
    domain=3
  [end]
  [myelmt4]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate4
    domain=4
  [end]
  [myelmt5]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate5
    domain=5
  [end]
  [myelmt6]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate6
    domain=6
  [end]
  [myelmt7]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate7
    domain=7
  [end]
  [myelmt8]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate8
    domain=8
  [end]
  [myelmt9]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate9
    domain=9
  [end]
  [myelmt10]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate10
    domain=10
  [end]
  [myelmt11]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate11
    domain=11
  [end]
  [myelmt12]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate12
    domain=12
  [end]
  [myelmt13]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate13
    domain=13
  [end]
  [myelmt14]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate14
    domain=14
  [end]
  [myelmt15]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate15
    domain=15
  [end]
  [myelmt16]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate16
    domain=16
  [end]
  [myelmt17]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate17
    domain=17
  [end]
  [myelmt18]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate18
    domain=18
  [end]
  [myelmt19]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate19
    domain=19
  [end]
  [myelmt20]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate20
    domain=20
  [end]
  [myelmt21]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate21
    domain=21
  [end]
  [myelmt22]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate22
    domain=22
  [end]
  [myelmt23]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate23
    domain=23
  [end]
  [myelmt24]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate24
    domain=24
  [end]
  [myelmt25]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate25
    domain=25
  [end]
  [myelmt26]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate26
    domain=26
  [end]
  [myelmt27]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate27
    domain=27
  [end]
  [myelmt28]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate28
    domain=28
  [end]
  [myelmt29]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate29
    domain=29
  [end]
  [myelmt30]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate30
    domain=30
  [end]
  [myelmt31]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate31
    domain=31
  [end]
  [myelmt32]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate32
    domain=32
  [end]
  [myelmt33]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate33
    domain=33
  [end]
  [myelmt34]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate34
    domain=34
  [end]
  [myelmt35]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate35
    domain=35
  [end]
  [myelmt36]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate36
    domain=36
  [end]
  [myelmt37]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate37
    domain=37
  [end]
  [myelmt38]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate38
    domain=38
  [end]
  [myelmt39]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate39
    domain=39
  [end]
  [myelmt40]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate40
    domain=40
  [end]
  [myelmt41]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate41
    domain=41
  [end]
  [myelmt42]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate42
    domain=42
  [end]
  [myelmt43]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate43
    domain=43
  [end]
  [myelmt44]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate44
    domain=44
  [end]
  [myelmt45]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate45
    domain=45
  [end]
  [myelmt46]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate46
    domain=46
  [end]
  [myelmt47]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate47
    domain=47
  [end]
  [myelmt48]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate48
    domain=48
  [end]
  [myelmt49]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate49
    domain=49
  [end]
  [myelmt50]
    type=diffusionfracture
    dofs=c d ux uy
    mate=mymate50
    domain=50
  [end]
[end]


[mates]
  [mymate1]
    type=diffusionfracturemate
    params= 1.5000e+01  3.8538e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate2]
    type=diffusionfracturemate
    params= 1.5000e+01  2.3588e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate3]
    type=diffusionfracturemate
    params= 1.5000e+01  2.6718e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate4]
    type=diffusionfracturemate
    params= 1.5000e+01  3.5439e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate5]
    type=diffusionfracturemate
    params= 1.5000e+01  4.6147e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate6]
    type=diffusionfracturemate
    params= 1.5000e+01  4.5724e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate7]
    type=diffusionfracturemate
    params= 1.5000e+01  4.4332e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate8]
    type=diffusionfracturemate
    params= 1.5000e+01  5.0123e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate9]
    type=diffusionfracturemate
    params= 1.5000e+01  3.9020e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate10]
    type=diffusionfracturemate
    params= 1.5000e+01  4.6325e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate11]
    type=diffusionfracturemate
    params= 1.5000e+01  4.2351e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate12]
    type=diffusionfracturemate
    params= 1.5000e+01  3.5028e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate13]
    type=diffusionfracturemate
    params= 1.5000e+01  4.5522e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate14]
    type=diffusionfracturemate
    params= 1.5000e+01  4.2147e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate15]
    type=diffusionfracturemate
    params= 1.5000e+01  4.5536e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate16]
    type=diffusionfracturemate
    params= 1.5000e+01  3.9821e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate17]
    type=diffusionfracturemate
    params= 1.5000e+01  5.6622e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate18]
    type=diffusionfracturemate
    params= 1.5000e+01  4.3349e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate19]
    type=diffusionfracturemate
    params= 1.5000e+01  3.9489e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate20]
    type=diffusionfracturemate
    params= 1.5000e+01  4.6540e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate21]
    type=diffusionfracturemate
    params= 1.5000e+01  2.7760e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate22]
    type=diffusionfracturemate
    params= 1.5000e+01  3.8087e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate23]
    type=diffusionfracturemate
    params= 1.5000e+01  3.1807e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate24]
    type=diffusionfracturemate
    params= 1.5000e+01  3.2900e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate25]
    type=diffusionfracturemate
    params= 1.5000e+01  5.1266e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate26]
    type=diffusionfracturemate
    params= 1.5000e+01  3.8458e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate27]
    type=diffusionfracturemate
    params= 1.5000e+01  4.1310e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate28]
    type=diffusionfracturemate
    params= 1.5000e+01  4.1958e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate29]
    type=diffusionfracturemate
    params= 1.5000e+01  4.6036e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate30]
    type=diffusionfracturemate
    params= 1.5000e+01  4.0516e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate31]
    type=diffusionfracturemate
    params= 1.5000e+01  4.7021e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate32]
    type=diffusionfracturemate
    params= 1.5000e+01  4.0598e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate33]
    type=diffusionfracturemate
    params= 1.5000e+01  6.2156e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate34]
    type=diffusionfracturemate
    params= 1.5000e+01  5.2309e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate35]
    type=diffusionfracturemate
    params= 1.5000e+01  4.0867e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate36]
    type=diffusionfracturemate
    params= 1.5000e+01  3.8019e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate37]
    type=diffusionfracturemate
    params= 1.5000e+01  3.7469e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate38]
    type=diffusionfracturemate
    params= 1.5000e+01  3.6825e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate39]
    type=diffusionfracturemate
    params= 1.5000e+01  5.2089e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate40]
    type=diffusionfracturemate
    params= 1.5000e+01  3.4652e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate41]
    type=diffusionfracturemate
    params= 1.5000e+01  3.3110e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate42]
    type=diffusionfracturemate
    params= 1.5000e+01  5.5391e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate43]
    type=diffusionfracturemate
    params= 1.5000e+01  3.2730e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate44]
    type=diffusionfracturemate
    params= 1.5000e+01  2.9201e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate45]
    type=diffusionfracturemate
    params= 1.5000e+01  4.7679e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate46]
    type=diffusionfracturemate
    params= 1.5000e+01  3.4930e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate47]
    type=diffusionfracturemate
    params= 1.5000e+01  3.9983e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate48]
    type=diffusionfracturemate
    params= 1.5000e+01  4.1094e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate49]
    type=diffusionfracturemate
    params= 1.5000e+01  2.5553e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
  [mymate50]
    type=diffusionfracturemate
    params= 1.5000e+01  4.3878e-02  2.7000e-03  1.0000e-02  1.0000e-06  2.1000e+02  2.5000e-01  0.0000e+00
  [end]
[end]
