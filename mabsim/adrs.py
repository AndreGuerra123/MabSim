from fenics import *
import numpy as np;

set_log_active(False)

T = 1.0            # final time
num_steps = 10    # number of time steps
dt = T / num_steps # time step size
K = 0.001        # reaction rate

#Unit in mesh
mesh =  UnitIntervalMesh(50)

# Define function space for system of concentrations
P1 = FiniteElement('P', interval, 1)
element = MixedElement([P1, P1, P1])
V = FunctionSpace(mesh, element)

# Define test functions
v_1, v_2, v_3 = TestFunctions(V)

# Define functions for velocity and concentrations
u = Function(V)
u_n = Function(V)

# Split system functions to access components
u_1, u_2, u_3 = split(u)
u_n1, u_n2, u_n3 = split(u_n)

# Define source terms
f_1 = Constant(1)
f_2 = Constant(1)
f_3 = Constant(0)

# Define expressions used in variational forms
k = Constant(dt)
K = Constant(K)
w = Constant(0.005)

# Define variational problem
F = ((u_1 - u_n1) / k)*v_1*dx + w*u_1.dx(0)*v_1*dx + K*u_1*u_2*v_1*dx  \
  + ((u_2 - u_n2) / k)*v_2*dx + w*u_2.dx(0)*v_2*dx + K*u_1*u_2*v_2*dx  \
  + ((u_3 - u_n3) / k)*v_3*dx + w*u_3.dx(0)*v_3*dx - K*u_1*u_2*v_3*dx \
  - f_1*v_1*dx - f_2*v_2*dx - f_3*v_3*dx


t = 0

for n in range(num_steps):

    # Update current time
    t += dt

    # Solve variational problem for time step
    solve(F == 0, u);

    # Update previous solution
    u_n.assign(u)

#Return final state array

print(u.vector().get_local())

