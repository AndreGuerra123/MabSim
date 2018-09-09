import os;
os.environ["OMPI_MCA_btl_base_warn_component_unused"] = "0"
from fenics import *;
set_log_active(False)
import sys,json;
import numpy as np;

args =  json.loads(sys.argv[1]);

try:
    n_species = 77;
    n_steps_t = 100;
    n_steps_z = 100;

    #Defining initial conditions and parameters
    #{iter,dt,w,mab,E1max,z1max,w1,Kd1,vo,kf1}

    iter = args['iter']
    T = args['dt']
    dt = Constant(T/n_steps_t) 
    mab0 = args['mab0']
    w = Constant(args['vo'])
    kf1 = Constant(args['kf1'])
    E1max = Constant(args['E1max'])
    z1max = Constant(args['z1max'])
    w1 = Constant(args['w1'])
    Kd1 = Constant(args['Kd1'])

    #Defining Mesh
    mesh = UnitIntervalMesh(n_steps_z)

    #Defining function space for system of concentrations

    P1  = FiniteElement('P', interval, 1)
    elements = MixedElement([P1]*n_species)
    V = FunctionSpace(mesh,elements)

    #Defining Test Functions

    v = TestFunctions(V)

    #Defining  Velocity and Concentration Functions

    u = Function(V)
    u_n = Function(V)

    if iter > 1:
        u_input = np.load('u.npy')
        u.vector().set_local(u_input)

    #Defining variational problem
    K= Kd1;  #jUST FOR NOW

    F = ((u_1 - u_n1) / dt)*v_1*dx + w*u_1.dx(0)*v_1*dx + K*u_1*u_2*v_1*dx  
        + ((u_2 - u_n2) / dt)*v_2*dx + w*u_2.dx(0)*v_2*dx + K*u_1*u_2*v_2*dx  
        + ((u_3 - u_n3) / dt)*v_3*dx + w*u_3.dx(0)*v_3*dx - K*u_1*u_2*v_3*dx 
        

        - Constant(mab0)*v[0]*dx 

    t = 0
    for n in range(n_steps_t):

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


except Exception as e:

    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)





