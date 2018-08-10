This fold contains the scripts for poisson equation.

# Poisson1D.py
Implement the 1d poisson equation which looks like:
$$A\nabla\cdot(\nabla\phi)=F$$
the related weak form for this problem is:
$$R^{I}_{\phi}=-\int_{\Omega}FN^{I}dV-\int_{\Omega}A\nabla\phi\nabla N^{I}dV$$

the comparison between analytical solution and FEM results looks like:
![poisson1d](compare.png)