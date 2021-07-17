# Model Predictive Control of a double pendulum on a cart

I, Vinicius (reach me at vinicius.fontes@aluno.puc-rio.br), have headed the efforts to build and publish this work, which also belongs to the following authors:
- Bruno
- Nicholas Casaprima
- Walisson Pinto

This work has been comissioned by dr. Helon Ayala for the class "Control of Mechanical Systems" taught at the Pontifical Catholic University of Rio de Janeiro (1st semester of 2021). For further information, check https://github.com/helonayala/

This repo consists in the following set of Matlab® Livescripts with the implementation of an MPC of a double pendulum on a cart:
 
 - Model.mlx: Derivation of the dynamics of the system using Euler-Lagrange;
 - MPC_singleshooting.mlx: Implementation of the MPC usign single shooting (rather slow);
 - MPC_multipleshooting.mlx: Same as above, but using multiple shooting (significantly faster).

The in-depth explanation of the code is meant to aid the reader in every step of the implementation, for teaching purposes, hence some choices have been adopted to favor readability rather than perfomance. The parametrization of the controller was set to values that show the application of the MPC method, but could be further optimized (try it yourself by changing the weight matrices).

Derivation and parameters are based on [1]. Matlab implementation was based on [2], where significant changes were made to make it easier to read. An optimization tool [3] has been used to improve performance.

References:

[1] K. Furuta, T. Okutani, H. Sone, Computer control of a double inverted pendulum, Computers & Electrical Engineering, Volume 5, Issue 1, 1978, Pages 67-84, ISSN 0045-7906, https://doi.org/10.1016/0045-7906(78)90018-6.

[2] Mohamed W. Mehrez. Workshop: MPC and MHE implementation in MATLAB using Casadi. Further information available at:
https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi

[3] Andersson, J. A. E.; Gillis, J.; Horn, G.; Rawlings, J. B. & Diehl, M.
CasADi -- A software framework for nonlinear optimization and optimal control 
Mathematical Programming Computation, In Press, 2018