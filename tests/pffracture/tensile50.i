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
		params=210 0.3 0.0027 0.004 1e-05 76.6975 93.1389 84.9454
	[end]
	[anisofrac2]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 96.3187 97.9709 93.3117
	[end]
	[anisofrac3]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 85.1263 85.1316 90.1527
	[end]
	[anisofrac4]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 102.216 94.6241 111.044
	[end]
	[anisofrac5]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 69.9686 97.7003 78.2755
	[end]
	[anisofrac6]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 84.6326 89.65 97.6428
	[end]
	[anisofrac7]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 85.6275 83.1001 74.6402
	[end]
	[anisofrac8]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 100.568 94.8636 81.9489
	[end]
	[anisofrac9]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 106.832 89.2857 97.897
	[end]
	[anisofrac10]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 83.9148 91.0929 75.4978
	[end]
	[anisofrac11]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 115.398 83.2321 109.496
	[end]
	[anisofrac12]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 92.4453 66.5475 108.367
	[end]
	[anisofrac13]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 80.6359 95.9231 99.1303
	[end]
	[anisofrac14]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 85.2934 97.673 81.0156
	[end]
	[anisofrac15]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 97.5807 102.749 86.9126
	[end]
	[anisofrac16]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 98.7526 98.8432 97.4714
	[end]
	[anisofrac17]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 86.3579 79.4671 81.3739
	[end]
	[anisofrac18]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 101.784 75.9145 95.6981
	[end]
	[anisofrac19]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 95.7734 91.3447 87.8987
	[end]
	[anisofrac20]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 83.8471 75.8511 74.1238
	[end]
	[anisofrac21]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 86.8876 75.3361 80.7084
	[end]
	[anisofrac22]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 91.0596 78.9156 96.8479
	[end]
	[anisofrac23]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 90.3676 95.7127 91.4149
	[end]
	[anisofrac24]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 84.8891 91.1558 80.3176
	[end]
	[anisofrac25]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 90.1254 95.2175 97.6514
	[end]
	[anisofrac26]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 86.6582 86.7502 118.426
	[end]
	[anisofrac27]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 102.958 110.469 97.3493
	[end]
	[anisofrac28]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 96.1886 88.2264 77.9534
	[end]
	[anisofrac29]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 102.41 83.1243 78.1723
	[end]
	[anisofrac30]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 108.039 92.3597 89.6973
	[end]
	[anisofrac31]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 73.6351 86.3908 83.8451
	[end]
	[anisofrac32]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 80.5891 78.3667 93.8978
	[end]
	[anisofrac33]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 87.9151 96.1757 87.0381
	[end]
	[anisofrac34]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 90.9724 90.1396 85.8023
	[end]
	[anisofrac35]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 79.2475 102.696 90.1446
	[end]
	[anisofrac36]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 82.3906 73.658 89.5184
	[end]
	[anisofrac37]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 94.8149 77.0486 96.4735
	[end]
	[anisofrac38]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 105.727 77.3053 109.33
	[end]
	[anisofrac39]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 90.2339 101.649 71.9136
	[end]
	[anisofrac40]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 95.3894 88.2516 96.5774
	[end]
	[anisofrac41]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 109.06 80.1763 94.3119
	[end]
	[anisofrac42]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 85.1165 87.4295 103.403
	[end]
	[anisofrac43]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 95.7218 82.3792 89.3299
	[end]
	[anisofrac44]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 72.9248 86.8734 88.6534
	[end]
	[anisofrac45]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 83.5275 88.3733 98.0869
	[end]
	[anisofrac46]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 109.82 103.227 105.622
	[end]
	[anisofrac47]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 86.1197 82.9479 81.6163
	[end]
	[anisofrac48]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 88.1777 89.8504 90.4295
	[end]
	[anisofrac49]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 93.1657 84.8125 87.3945
	[end]
	[anisofrac50]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 75.159 82.1861 79.9291
	[end]
[end]

[bcs]
	[fixUx]
		type=dirichlet
		dof=ux
		boundary=bottom
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
	dtmax=1.0e-4
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
