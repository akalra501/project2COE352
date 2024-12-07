## Project Objective

This project aims to solve the 1D heat transfer problem using the Galerkin Finite Element Method (FEM). The problem is defined by the partial differential equation:

\[ 
u_t - u_{xx} = f(x, t), \quad (x, t) \in (0, 1) \times (0, 1),
\]

with the following conditions:

- **Initial Condition**:  
  \[
  u(x, 0) = \sin(\pi x),
  \]

- **Dirichlet Boundary Conditions**:  
  \[
  u(0, t) = u(1, t) = 0,
  \]

- **Source Term**:  
  \[
  f(x, t) = (\pi^2 - 1)e^{-t} \sin(\pi x),
  \]

- **Analytical Solution**:  
  \[
  u(x, t) = e^{-t} \sin(\pi x).
  \]

### Key Features:
1. **Generalized Framework**:  
   The code is designed to allow different source terms \( f(x, t) \) and boundary conditions to be easily implemented.

2. **1D Lagrange Basis Functions**:  
   Second-order Gaussian quadrature is used to numerically integrate the weak form of the source term.

3. **Elemental Matrices**:  
   The mass and stiffness matrices are built by mapping from the physical space (\(x\)) to the reference space (\(\xi\)), and integration is performed in the parent space.

4. **Node Flexibility**:  
   While the default setup uses \( N = 11 \) nodes, the code is generalized to handle any number of spatial nodes.

This project demonstrates the application of the Galerkin FEM approach for time-dependent partial differential equations (PDEs), emphasizing modularity, flexibility, and numerical accuracy.
