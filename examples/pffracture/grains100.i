[mesh]
  type=gmsh
  file=grains100.msh
[end]


[dofs]
name=d ux uy
[end]

[nonlinearsolver]
  type=nr
  maxiters=20
  r_rel_tol=5.0e-10
  r_abs_tol=5.5e-7
  //solver=superlu
[end]

[ics]
  [constd]
    type=const
    dof=d
    params=0.0
  [end]
[end]

[output]
  type=vtu
  interval=10
[end]

[timestepping]
  type=be
  dt=2.0e-5
  time=5.0e-1
  adaptive=true
  optiters=5
  growthfactor=1.2
  cutfactor=0.85
  dtmin=1.0e-12
  dtmax=1.5e-3
[end]

[bcs]
  [fixux]
    type=dirichlet
    dofs=ux
    value=0.0
    boundary=left
  [end]
  [fixuy]
    type=dirichlet
    dofs=uy
    value=0.0
    boundary=bottom
  [end]
  [load]
    type=dirichlet
    dofs=uy
    value=1.0*t
    boundary=top
  [end]
[end]

[projection]
scalarmate=vonMises
[end]

[job]
  type=transient
  debug=dep
[end]

//***************************************************
[elmts]
  [myelmt1]
    type=miehefrac
    dofs=d ux uy
    mate=mymate1
    domain=1
  [end]
  [myelmt2]
    type=miehefrac
    dofs=d ux uy
    mate=mymate2
    domain=2
  [end]
  [myelmt3]
    type=miehefrac
    dofs=d ux uy
    mate=mymate3
    domain=3
  [end]
  [myelmt4]
    type=miehefrac
    dofs=d ux uy
    mate=mymate4
    domain=4
  [end]
  [myelmt5]
    type=miehefrac
    dofs=d ux uy
    mate=mymate5
    domain=5
  [end]
  [myelmt6]
    type=miehefrac
    dofs=d ux uy
    mate=mymate6
    domain=6
  [end]
  [myelmt7]
    type=miehefrac
    dofs=d ux uy
    mate=mymate7
    domain=7
  [end]
  [myelmt8]
    type=miehefrac
    dofs=d ux uy
    mate=mymate8
    domain=8
  [end]
  [myelmt9]
    type=miehefrac
    dofs=d ux uy
    mate=mymate9
    domain=9
  [end]
  [myelmt10]
    type=miehefrac
    dofs=d ux uy
    mate=mymate10
    domain=10
  [end]
  [myelmt11]
    type=miehefrac
    dofs=d ux uy
    mate=mymate11
    domain=11
  [end]
  [myelmt12]
    type=miehefrac
    dofs=d ux uy
    mate=mymate12
    domain=12
  [end]
  [myelmt13]
    type=miehefrac
    dofs=d ux uy
    mate=mymate13
    domain=13
  [end]
  [myelmt14]
    type=miehefrac
    dofs=d ux uy
    mate=mymate14
    domain=14
  [end]
  [myelmt15]
    type=miehefrac
    dofs=d ux uy
    mate=mymate15
    domain=15
  [end]
  [myelmt16]
    type=miehefrac
    dofs=d ux uy
    mate=mymate16
    domain=16
  [end]
  [myelmt17]
    type=miehefrac
    dofs=d ux uy
    mate=mymate17
    domain=17
  [end]
  [myelmt18]
    type=miehefrac
    dofs=d ux uy
    mate=mymate18
    domain=18
  [end]
  [myelmt19]
    type=miehefrac
    dofs=d ux uy
    mate=mymate19
    domain=19
  [end]
  [myelmt20]
    type=miehefrac
    dofs=d ux uy
    mate=mymate20
    domain=20
  [end]
  [myelmt21]
    type=miehefrac
    dofs=d ux uy
    mate=mymate21
    domain=21
  [end]
  [myelmt22]
    type=miehefrac
    dofs=d ux uy
    mate=mymate22
    domain=22
  [end]
  [myelmt23]
    type=miehefrac
    dofs=d ux uy
    mate=mymate23
    domain=23
  [end]
  [myelmt24]
    type=miehefrac
    dofs=d ux uy
    mate=mymate24
    domain=24
  [end]
  [myelmt25]
    type=miehefrac
    dofs=d ux uy
    mate=mymate25
    domain=25
  [end]
  [myelmt26]
    type=miehefrac
    dofs=d ux uy
    mate=mymate26
    domain=26
  [end]
  [myelmt27]
    type=miehefrac
    dofs=d ux uy
    mate=mymate27
    domain=27
  [end]
  [myelmt28]
    type=miehefrac
    dofs=d ux uy
    mate=mymate28
    domain=28
  [end]
  [myelmt29]
    type=miehefrac
    dofs=d ux uy
    mate=mymate29
    domain=29
  [end]
  [myelmt30]
    type=miehefrac
    dofs=d ux uy
    mate=mymate30
    domain=30
  [end]
  [myelmt31]
    type=miehefrac
    dofs=d ux uy
    mate=mymate31
    domain=31
  [end]
  [myelmt32]
    type=miehefrac
    dofs=d ux uy
    mate=mymate32
    domain=32
  [end]
  [myelmt33]
    type=miehefrac
    dofs=d ux uy
    mate=mymate33
    domain=33
  [end]
  [myelmt34]
    type=miehefrac
    dofs=d ux uy
    mate=mymate34
    domain=34
  [end]
  [myelmt35]
    type=miehefrac
    dofs=d ux uy
    mate=mymate35
    domain=35
  [end]
  [myelmt36]
    type=miehefrac
    dofs=d ux uy
    mate=mymate36
    domain=36
  [end]
  [myelmt37]
    type=miehefrac
    dofs=d ux uy
    mate=mymate37
    domain=37
  [end]
  [myelmt38]
    type=miehefrac
    dofs=d ux uy
    mate=mymate38
    domain=38
  [end]
  [myelmt39]
    type=miehefrac
    dofs=d ux uy
    mate=mymate39
    domain=39
  [end]
  [myelmt40]
    type=miehefrac
    dofs=d ux uy
    mate=mymate40
    domain=40
  [end]
  [myelmt41]
    type=miehefrac
    dofs=d ux uy
    mate=mymate41
    domain=41
  [end]
  [myelmt42]
    type=miehefrac
    dofs=d ux uy
    mate=mymate42
    domain=42
  [end]
  [myelmt43]
    type=miehefrac
    dofs=d ux uy
    mate=mymate43
    domain=43
  [end]
  [myelmt44]
    type=miehefrac
    dofs=d ux uy
    mate=mymate44
    domain=44
  [end]
  [myelmt45]
    type=miehefrac
    dofs=d ux uy
    mate=mymate45
    domain=45
  [end]
  [myelmt46]
    type=miehefrac
    dofs=d ux uy
    mate=mymate46
    domain=46
  [end]
  [myelmt47]
    type=miehefrac
    dofs=d ux uy
    mate=mymate47
    domain=47
  [end]
  [myelmt48]
    type=miehefrac
    dofs=d ux uy
    mate=mymate48
    domain=48
  [end]
  [myelmt49]
    type=miehefrac
    dofs=d ux uy
    mate=mymate49
    domain=49
  [end]
  [myelmt50]
    type=miehefrac
    dofs=d ux uy
    mate=mymate50
    domain=50
  [end]
  [myelmt51]
    type=miehefrac
    dofs=d ux uy
    mate=mymate51
    domain=51
  [end]
  [myelmt52]
    type=miehefrac
    dofs=d ux uy
    mate=mymate52
    domain=52
  [end]
  [myelmt53]
    type=miehefrac
    dofs=d ux uy
    mate=mymate53
    domain=53
  [end]
  [myelmt54]
    type=miehefrac
    dofs=d ux uy
    mate=mymate54
    domain=54
  [end]
  [myelmt55]
    type=miehefrac
    dofs=d ux uy
    mate=mymate55
    domain=55
  [end]
  [myelmt56]
    type=miehefrac
    dofs=d ux uy
    mate=mymate56
    domain=56
  [end]
  [myelmt57]
    type=miehefrac
    dofs=d ux uy
    mate=mymate57
    domain=57
  [end]
  [myelmt58]
    type=miehefrac
    dofs=d ux uy
    mate=mymate58
    domain=58
  [end]
  [myelmt59]
    type=miehefrac
    dofs=d ux uy
    mate=mymate59
    domain=59
  [end]
  [myelmt60]
    type=miehefrac
    dofs=d ux uy
    mate=mymate60
    domain=60
  [end]
  [myelmt61]
    type=miehefrac
    dofs=d ux uy
    mate=mymate61
    domain=61
  [end]
  [myelmt62]
    type=miehefrac
    dofs=d ux uy
    mate=mymate62
    domain=62
  [end]
  [myelmt63]
    type=miehefrac
    dofs=d ux uy
    mate=mymate63
    domain=63
  [end]
  [myelmt64]
    type=miehefrac
    dofs=d ux uy
    mate=mymate64
    domain=64
  [end]
  [myelmt65]
    type=miehefrac
    dofs=d ux uy
    mate=mymate65
    domain=65
  [end]
  [myelmt66]
    type=miehefrac
    dofs=d ux uy
    mate=mymate66
    domain=66
  [end]
  [myelmt67]
    type=miehefrac
    dofs=d ux uy
    mate=mymate67
    domain=67
  [end]
  [myelmt68]
    type=miehefrac
    dofs=d ux uy
    mate=mymate68
    domain=68
  [end]
  [myelmt69]
    type=miehefrac
    dofs=d ux uy
    mate=mymate69
    domain=69
  [end]
  [myelmt70]
    type=miehefrac
    dofs=d ux uy
    mate=mymate70
    domain=70
  [end]
  [myelmt71]
    type=miehefrac
    dofs=d ux uy
    mate=mymate71
    domain=71
  [end]
  [myelmt72]
    type=miehefrac
    dofs=d ux uy
    mate=mymate72
    domain=72
  [end]
  [myelmt73]
    type=miehefrac
    dofs=d ux uy
    mate=mymate73
    domain=73
  [end]
  [myelmt74]
    type=miehefrac
    dofs=d ux uy
    mate=mymate74
    domain=74
  [end]
  [myelmt75]
    type=miehefrac
    dofs=d ux uy
    mate=mymate75
    domain=75
  [end]
  [myelmt76]
    type=miehefrac
    dofs=d ux uy
    mate=mymate76
    domain=76
  [end]
  [myelmt77]
    type=miehefrac
    dofs=d ux uy
    mate=mymate77
    domain=77
  [end]
  [myelmt78]
    type=miehefrac
    dofs=d ux uy
    mate=mymate78
    domain=78
  [end]
  [myelmt79]
    type=miehefrac
    dofs=d ux uy
    mate=mymate79
    domain=79
  [end]
  [myelmt80]
    type=miehefrac
    dofs=d ux uy
    mate=mymate80
    domain=80
  [end]
  [myelmt81]
    type=miehefrac
    dofs=d ux uy
    mate=mymate81
    domain=81
  [end]
  [myelmt82]
    type=miehefrac
    dofs=d ux uy
    mate=mymate82
    domain=82
  [end]
  [myelmt83]
    type=miehefrac
    dofs=d ux uy
    mate=mymate83
    domain=83
  [end]
  [myelmt84]
    type=miehefrac
    dofs=d ux uy
    mate=mymate84
    domain=84
  [end]
  [myelmt85]
    type=miehefrac
    dofs=d ux uy
    mate=mymate85
    domain=85
  [end]
  [myelmt86]
    type=miehefrac
    dofs=d ux uy
    mate=mymate86
    domain=86
  [end]
  [myelmt87]
    type=miehefrac
    dofs=d ux uy
    mate=mymate87
    domain=87
  [end]
  [myelmt88]
    type=miehefrac
    dofs=d ux uy
    mate=mymate88
    domain=88
  [end]
  [myelmt89]
    type=miehefrac
    dofs=d ux uy
    mate=mymate89
    domain=89
  [end]
  [myelmt90]
    type=miehefrac
    dofs=d ux uy
    mate=mymate90
    domain=90
  [end]
  [myelmt91]
    type=miehefrac
    dofs=d ux uy
    mate=mymate91
    domain=91
  [end]
  [myelmt92]
    type=miehefrac
    dofs=d ux uy
    mate=mymate92
    domain=92
  [end]
  [myelmt93]
    type=miehefrac
    dofs=d ux uy
    mate=mymate93
    domain=93
  [end]
  [myelmt94]
    type=miehefrac
    dofs=d ux uy
    mate=mymate94
    domain=94
  [end]
  [myelmt95]
    type=miehefrac
    dofs=d ux uy
    mate=mymate95
    domain=95
  [end]
  [myelmt96]
    type=miehefrac
    dofs=d ux uy
    mate=mymate96
    domain=96
  [end]
  [myelmt97]
    type=miehefrac
    dofs=d ux uy
    mate=mymate97
    domain=97
  [end]
  [myelmt98]
    type=miehefrac
    dofs=d ux uy
    mate=mymate98
    domain=98
  [end]
  [myelmt99]
    type=miehefrac
    dofs=d ux uy
    mate=mymate99
    domain=99
  [end]
  [myelmt100]
    type=miehefrac
    dofs=d ux uy
    mate=mymate100
    domain=100
  [end]
[end]


[mates]
  [mymate1]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  3.1011e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate2]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  4.4367e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate3]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.7297e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate4]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.6858e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate5]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  5.9324e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate6]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  3.4924e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate7]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  4.1837e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate8]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  5.0700e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate9]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.1565e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate10]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  5.7606e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate11]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.4381e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate12]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.0442e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate13]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  9.3933e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate14]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.1852e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate15]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.6262e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate16]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  2.8022e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate17]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  9.8475e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate18]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.4214e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate19]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.9002e+00  0.0000e+00  0.0000e+00
  [end]
  [mymate20]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  6.1131e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate21]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.1512e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate22]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.1764e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate23]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.1095e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate24]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.3116e+00  0.0000e+00  0.0000e+00
  [end]
  [mymate25]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.5362e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate26]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.0077e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate27]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  4.5168e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate28]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.9168e+00  0.0000e+00  0.0000e+00
  [end]
  [mymate29]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  6.8999e+00  0.0000e+00  0.0000e+00
  [end]
  [mymate30]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  9.2069e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate31]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.3373e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate32]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.7659e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate33]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.1850e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate34]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.5187e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate35]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  6.3566e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate36]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  2.3824e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate37]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.4971e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate38]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.7067e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate39]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  4.1485e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate40]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  9.4432e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate41]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  3.6465e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate42]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.7228e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate43]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.2074e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate44]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  3.6011e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate45]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  6.9048e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate46]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  9.6256e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate47]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.6433e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate48]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.7657e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate49]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.0302e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate50]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  2.0627e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate51]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  9.7051e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate52]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.4453e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate53]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  9.3438e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate54]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  9.9087e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate55]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  3.5163e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate56]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.6197e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate57]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  5.1217e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate58]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.0648e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate59]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  8.6263e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate60]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.0666e+00  0.0000e+00  0.0000e+00
  [end]
  [mymate61]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  8.8185e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate62]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.2498e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate63]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.4538e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate64]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  6.9887e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate65]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.7438e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate66]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.2720e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate67]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  2.5430e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate68]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  2.3035e+00  0.0000e+00  0.0000e+00
  [end]
  [mymate69]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.3735e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate70]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.0754e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate71]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.7837e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate72]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.1686e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate73]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.7238e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate74]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  3.5303e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate75]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  6.1089e+00  0.0000e+00  0.0000e+00
  [end]
  [mymate76]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.3263e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate77]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.9050e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate78]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  2.8570e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate79]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  6.8064e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate80]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  9.9181e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate81]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  4.2778e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate82]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.6750e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate83]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  3.3200e+00  0.0000e+00  0.0000e+00
  [end]
  [mymate84]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.3916e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate85]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.3416e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate86]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.0537e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate87]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  5.2265e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate88]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.4679e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate89]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  8.2594e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate90]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.7846e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate91]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.6900e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate92]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.1479e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate93]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  2.0052e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate94]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  9.3230e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate95]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  2.5974e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate96]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.4064e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate97]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.3883e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate98]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  8.1535e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate99]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.5731e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate100]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.7815e+02  0.0000e+00  0.0000e+00
  [end]
[end]
