## Project Objective

This project aims to solve the 1D heat transfer problem using the Galerkin Finite Element Method (FEM). The problem is defined by the partial differential equation:
u_t - u_{xx} = f(x, t)

with the following conditions:

- **Initial Condition**:  
  u(x, 0) = sin(pi*x),

- **Dirichlet Boundary Conditions**:  
  u(0, t) = u(1, t) = 0,

- **Source Term**:  
  f(x, t) = (pi^2 - 1)e^{-t} sin(pi*x),

- **Analytical Solution**:  
  u(x, t) = e^{-t} * sin(pi*x).

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


### Implementation:
The main Python script can be run after cloning the repository. The number of nodes, timesteps, and the choice of the Euler Method can be modified by the user.

## Problem 1:
The weak form of the equation was derived simply by multiplying the equation by a weak function and then integrating it by parts. The boundary conditions were also applied.

## Problem 2: 
The problem was solved  using the Forward Euler time derivative discretization. The initial time steps entered was 551(as mentioned in the problem), dt was 1/timeSteps = 1/551. This graph was plotted and is saved in the repository. The next step was to increase the time-step(dt) to observe instability. To do this the timeSteps variable had to be reduced as it is in the denominator. The timeSteps variable was decreased and instability was observed at each time-step(dt). First decreased to 500, then to 400, then to 300, then to 200. At timeSteps = 200 instability was present. This meant that the instability occurred at a timeSteps value between 200 and 300. Next, the timeSteps value was set to 250, and instability was still present, meaning the instability occurred between a timeSteps value of 250 and 300. The plots were then observed at timeSteps 260, 270, and 280. Instability was observed between timeSteps values of 270 and 280. Then the timeSteps value was decreased from 280 by 1, and instability was found to occur at timeSteps = 278 or the time-step value of 1/278. The plot of the Forward Euler solution at 278 is included in the repo. 

The effect of lower nodes was also observed. A plot with numNodes(N) = 6 and timeSteps(t) = 551 was observed. As shown in the plot, the Forward Euler solution undershot the actual solution when the number of spatial nodes was decreased.

## Problem 3
Problem 3 was solved using the Backward Euler time derivative discretization. The plots were constructed with time-steps of 1/551, and 1/278(same as forward Euler). However, the problem also asked us to solve the problem where the time-step was equal to 11. The plots with time-steps of 1/551, and 1/278 were both stable indicating the unconditional stability of the Backward Euler Method. However, when the time-step value was equal to or greater than the spatial time-step value instability was observed. This is because the approximation of the derivative becomes less accurate, leading to lesser accuracy overall.

_The plots for all the Forward and Backward Euler solutions are saved in the repository._



