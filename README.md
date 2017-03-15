# Fins

## MATLAB

### Quasi Orange free
 
Adaptive IMEX time stepping.

Fourth order _BWENO_.

_Centered FD_ is unstable.

---

class::      ins.m

setup(fS)

getDt(fS)
		  
regression:: caseConvergence.m

main::       main(fS)

% form the LHS matrix u with implicit time stepping 

formLHS(fS)

% form the LHS matrix for pressure Poisson eq

formLHSP(fS)


---
BC:

  1. no-slip
  2. periodic ;
  3. dirichlet ;
  4. inflow; 
  5. outflow;

