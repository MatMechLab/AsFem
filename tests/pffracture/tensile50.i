[mesh]
	type=gmsh
	file=rect50.msh
[end]

[dofs]
	name=ux uy d
[end]

[projection]
	name=vonMises Hydro sigxx sigyy sigxy
[end]

[elmts]
	[block1]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac1
		domain=1
	[end]
	[block2]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac2
		domain=2
	[end]
	[block3]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac3
		domain=3
	[end]
	[block4]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac4
		domain=4
	[end]
	[block5]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac5
		domain=5
	[end]
	[block6]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac6
		domain=6
	[end]
	[block7]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac7
		domain=7
	[end]
	[block8]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac8
		domain=8
	[end]
	[block9]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac9
		domain=9
	[end]
	[block10]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac10
		domain=10
	[end]
	[block11]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac11
		domain=11
	[end]
	[block12]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac12
		domain=12
	[end]
	[block13]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac13
		domain=13
	[end]
	[block14]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac14
		domain=14
	[end]
	[block15]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac15
		domain=15
	[end]
	[block16]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac16
		domain=16
	[end]
	[block17]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac17
		domain=17
	[end]
	[block18]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac18
		domain=18
	[end]
	[block19]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac19
		domain=19
	[end]
	[block20]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac20
		domain=20
	[end]
	[block21]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac21
		domain=21
	[end]
	[block22]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac22
		domain=22
	[end]
	[block23]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac23
		domain=23
	[end]
	[block24]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac24
		domain=24
	[end]
	[block25]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac25
		domain=25
	[end]
	[block26]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac26
		domain=26
	[end]
	[block27]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac27
		domain=27
	[end]
	[block28]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac28
		domain=28
	[end]
	[block29]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac29
		domain=29
	[end]
	[block30]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac30
		domain=30
	[end]
	[block31]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac31
		domain=31
	[end]
	[block32]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac32
		domain=32
	[end]
	[block33]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac33
		domain=33
	[end]
	[block34]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac34
		domain=34
	[end]
	[block35]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac35
		domain=35
	[end]
	[block36]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac36
		domain=36
	[end]
	[block37]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac37
		domain=37
	[end]
	[block38]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac38
		domain=38
	[end]
	[block39]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac39
		domain=39
	[end]
	[block40]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac40
		domain=40
	[end]
	[block41]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac41
		domain=41
	[end]
	[block42]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac42
		domain=42
	[end]
	[block43]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac43
		domain=43
	[end]
	[block44]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac44
		domain=44
	[end]
	[block45]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac45
		domain=45
	[end]
	[block46]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac46
		domain=46
	[end]
	[block47]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac47
		domain=47
	[end]
	[block48]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac48
		domain=48
	[end]
	[block49]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac49
		domain=49
	[end]
	[block50]
		type=miehefrac
		dofs=ux uy d
		mate=anisofrac50
		domain=50
	[end]
[end]

[mates]
	[anisofrac1]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 70.6307 69.513 89.0213
	[end]
	[anisofrac2]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 102.789 79.3597 90.2422
	[end]
	[anisofrac3]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 60.2582 82.2597 86.6815
	[end]
	[anisofrac4]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 25.7617 74.1068 78.3665
	[end]
	[anisofrac5]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 109.174 52.8944 87.3305
	[end]
	[anisofrac6]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 15.0278 88.9393 79.8049
	[end]
	[anisofrac7]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 67.15 86.1478 86.6108
	[end]
	[anisofrac8]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 43.8886 79.8927 91.9686
	[end]
	[anisofrac9]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 90.1691 87.5166 90.9584
	[end]
	[anisofrac10]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 63.7935 76.8016 91.0318
	[end]
	[anisofrac11]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 98.0412 68.6126 86.2793
	[end]
	[anisofrac12]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 74.5166 78.5469 94.3578
	[end]
	[anisofrac13]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 81.2352 92.5299 93.8862
	[end]
	[anisofrac14]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 101.077 50.5222 87.9804
	[end]
	[anisofrac15]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 113.313 66.2884 88.8689
	[end]
	[anisofrac16]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 89.6553 61.0419 79.6044
	[end]
	[anisofrac17]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 113.378 71.7552 94.8064
	[end]
	[anisofrac18]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 34.8805 92.9461 87.5508
	[end]
	[anisofrac19]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 69.5792 82.4483 91.4032
	[end]
	[anisofrac20]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 82.4822 87.1484 87.1982
	[end]
	[anisofrac21]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 65.9923 65.9311 98.3645
	[end]
	[anisofrac22]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 32.1263 87.858 90.656
	[end]
	[anisofrac23]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 20.7254 100.907 87.1601
	[end]
	[anisofrac24]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 73.5442 71.1666 80.8558
	[end]
	[anisofrac25]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 69.2283 88.1798 96.1777
	[end]
	[anisofrac26]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 56.6284 52.8001 84.9201
	[end]
	[anisofrac27]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 37.6265 64.2583 89.8539
	[end]
	[anisofrac28]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 81.931 43.6826 83.9687
	[end]
	[anisofrac29]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 104.449 73.6778 91.0602
	[end]
	[anisofrac30]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 69.2348 62.5783 92.3833
	[end]
	[anisofrac31]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 49.3836 73.9274 79.3551
	[end]
	[anisofrac32]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 80.4433 73.6108 80.8321
	[end]
	[anisofrac33]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 57.7375 50.3341 87.0999
	[end]
	[anisofrac34]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 65.1828 61.2264 84.393
	[end]
	[anisofrac35]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 86.9656 56.7425 85.4594
	[end]
	[anisofrac36]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 93.9858 69.0521 89.4861
	[end]
	[anisofrac37]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 109.161 67.7694 79.5168
	[end]
	[anisofrac38]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 83.4841 100.784 87.9079
	[end]
	[anisofrac39]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 68.3246 61.6974 91.0974
	[end]
	[anisofrac40]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 117.712 79.9803 87.771
	[end]
	[anisofrac41]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 97.0538 94.1388 91.1305
	[end]
	[anisofrac42]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 88.4038 91.8819 84.8044
	[end]
	[anisofrac43]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 58.0559 67.5473 89.2036
	[end]
	[anisofrac44]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 72.0202 100.153 99.9895
	[end]
	[anisofrac45]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 51.9394 79.6886 92.8106
	[end]
	[anisofrac46]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 73.2603 65.336 84.2291
	[end]
	[anisofrac47]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 69.6811 54.7638 97.5388
	[end]
	[anisofrac48]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 84.8913 58.0527 105.614
	[end]
	[anisofrac49]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 41.4414 73.0687 86.8858
	[end]
	[anisofrac50]
		type=anisopffrac
		params=210 90 0.3 0.2 0.0027 0.004 1e-05 31.4383 70.6422 91.4489
	[end]
[end]

[bcs]
	[fixUx]
		type=dirichlet
		dof=ux
		boundary=left
		value=0.0
	[end]
	[fixUy]
		type=dirichlet
		dof=uy
		boundary=bottom
		value=0.0
	[end]
	[loadUy]
		type=dirichlet
		dof=uy
		boundary=top
		value=1.0*t
	[end]
[end]

[timestepping]
	type=be
	dt=1.0e-5
	dtmax=4.0e-4
	dtmin=5.0e-9
	endtime=1.0e2
	opts=4
	adaptive=true
[end]

[nonlinearsolver]
	type=newtonls
	r_abs_tol=5.0e-7
	r_rel_tol=1.5e-8
	maxiters=25
[end]

[job]
	type=transient
	debug=dep
	projection=true
[end]
