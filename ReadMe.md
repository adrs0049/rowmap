Python interface for the ROWMAP integrator
===================================

Numerical methods
-----------------------

### ODE integrator

The ROWMAP integrator is an integrator suited for large stiff systems of ODEs.  It's code can be found in the folder rowmap. The rowmap solver is implemented in fortran. A python wrapper using the old scipy integrator interface is provided.

**TODO:**

1. Move the python rowmap interface to the new scipy integrator interface.


### MOL (Method of lines) 
Integrates a system of PDEs using a method of lines approach. A functor doing the discretization must be provided. 

