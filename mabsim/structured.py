import os;
os.environ["OMPI_MCA_btl_base_warn_component_unused"] = "0"
from fenics import *;
set_log_active(False)
import sys,json;
args =  json.loads(sys.argv[1]);
import numpy as np;
import pandas as pd;

try:
    #Defining initial conditions and inputs

    iter = args['iter']
    input = np.array(args['input'])

    #Defining parameters

    w = Constant(args['vo'])
    # kf1 = Constant(args['kf1'])
    E1max = Constant(args['E1max'])
    z1max = Constant(args['z1max'])
    w1 = Constant(args['w1'])
    Kd1 = Constant(args['Kd1'])

    #This system is from a continuous process evaluated in real time, so the retrieved dt is the actual time elapsed

    T = args['dt']     # final time
    num_steps = 10    # number of time steps
    dt = Constant(T / num_steps) # time step size

    #Defining Mesh
    m = 10;
    mesh = UnitIntervalMesh(m)

    #Defining function space for system of concentrations

    P1  = FiniteElement('P', interval, 1)
    element = MixedElement([P1,P1,P1])
    V = FunctionSpace(mesh,element)

    #Defining Test Functions

    v_1,v_2,v_3 = TestFunctions(V)

    #Defining  Velocity and Concentration Functions

    u = Function(V)
    u_n = Function(V)

    if iter > 1:
        u_input = np.load('u.npy')
        u.vector().set_local(u_input)

    #Split for Access
    s = 3
    u_1, u_2, u_3 = split(u)
    u_n1, u_n2, u_n3 = split(u_n)

    #Defining Source Terms

    f_1,f_2,f_3 = input

    f_1 = Constant(f_1)
    f_2 = Constant(f_2)
    f_3 = Constant(f_3)

    #Defining variational problem
    K= Kd1;  #jUST FOR NOW

    F = ((u_1 - u_n1) / dt)*v_1*dx + w*u_1.dx(0)*v_1*dx + K*u_1*u_2*v_1*dx  \
        + ((u_2 - u_n2) / dt)*v_2*dx + w*u_2.dx(0)*v_2*dx + K*u_1*u_2*v_2*dx  \
        + ((u_3 - u_n3) / dt)*v_3*dx + w*u_3.dx(0)*v_3*dx - K*u_1*u_2*v_3*dx \
        - f_1*v_1*dx - f_2*v_2*dx - f_3*v_3*dx

    t = 0
    for n in range(num_steps):

        # Update current time
        t += dt

        # Solve variational problem for time step
        solve(F == 0, u)

        # Update previous solution
        u_n.assign(u)

    # Save state to file
    u_input = u.vector().get_local()
    os.remove('u.npy')
    np.save('u',u_input)

    #  Output variables values
    out_1 = u_input[-s:]
    out_2 = np.nan_to_num(out_1)
    out_3 = np.around(out_2,6)
    out_4 = json.dumps(out_3.tolist())
    print(out_4)


except Exception as e:

    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)





