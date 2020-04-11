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
		params=210 0.3 0.0027 0.004 1e-05 69.6709 91.9911 76.7659
	[end]
	[anisofrac2]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 108.75 86.7641 74.0175
	[end]
	[anisofrac3]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 78.3269 99.7136 79.3489
	[end]
	[anisofrac4]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 82.0308 81.6665 108.364
	[end]
	[anisofrac5]
		type=anisopffrac
		params=210 0.3 0.0027 0.004 1e-05 90.7335 87.1857 91.866
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
