[mesh]
	type=gmsh
	file=rect5.msh
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
[end]

[mates]
	[anisofrac1]
		type=anisopffrac
		params=210 50 0.3 0.2 0.0027 0.004 1e-05 49.6071 84.6307 100.967
	[end]
	[anisofrac2]
		type=anisopffrac
		params=210 50 0.3 0.2 0.0027 0.004 1e-05 61.7844 68.9158 85.559
	[end]
	[anisofrac3]
		type=anisopffrac
		params=210 50 0.3 0.2 0.0027 0.004 1e-05 23.3104 65.2141 82.5777
	[end]
	[anisofrac4]
		type=anisopffrac
		params=210 50 0.3 0.2 0.0027 0.004 1e-05 70.8536 76.8472 91.5993
	[end]
	[anisofrac5]
		type=anisopffrac
		params=210 50 0.3 0.2 0.0027 0.004 1e-05 54.2964 82.8409 88.7877
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
