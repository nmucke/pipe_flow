from __future__ import print_function
from dolfin import *
from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
import pdb

#def pipe_flow():
# Right hand side f, zero-vector
c = 1227
rho_ref = 870
p_ref = 5e5
k = c*c*rho_ref-p_ref
w_l = 1e6
x_l = 500.
D = 20.
A = 3.1415*0.5*D*0.5*D



# Initial
u0 = Constant(0.)
p0 = Constant(0.)

N = 2048
a = 0.
b = 2000.

# Geometry setup and meshing
# Create domain
mesh = IntervalMesh(N,a,b)

# Specify boundary
def inlet(x, on_boundary):
    return on_boundary and near(x[0], a)

def outlet(x, on_boundary):
    return on_boundary and near(x[0], b)

# Set up discrete setting
# Finite Element degree
deg_velocity = 2
deg_pressure = 2
# Define Element
V_el = FiniteElement("P", mesh.ufl_cell(), deg_velocity)
Q_el = FiniteElement("P", mesh.ufl_cell(), deg_pressure)
# Mixed Element
W_el = MixedElement([V_el, Q_el])
# Define function space
V = FunctionSpace(mesh, "P", deg_velocity)
Q = FunctionSpace(mesh, "P", deg_pressure)
W = FunctionSpace(mesh, W_el)
# Define trial and test functions
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
# Define solution
w_sol = Function(W)
u_sol = Function(V)
p_sol = Function(Q)
# Interpolate initial datum into Finite Element space
u_old = interpolate(u0, V)
p_old = interpolate(p0, Q)

heaviside = interpolate(Expression('x[0] >= 0.0 ? 1.0 : 0.0',degree=1),V)
dirac = interpolate(Expression('x[0]-x_l == 0.0 ? 1.0 : 0.0',x_l=x_l,degree=2),V)

# Time discretization
# Theta method
theta = 1
# Final Time
T_end = 1.5
# number of time dofs
time_dofs = 7500
# Time stepsize
dt = T_end / time_dofs

# Set up boundary data
# Define boundary values
bv_inlet = Expression(('.1'), degree=2)
bv_outlet = Constant(0.)
# Define boundary conditions
bc_inlet = DirichletBC(W.sub(0), bv_inlet, inlet)
bc_outlet = DirichletBC(W.sub(1), bv_outlet, outlet)
bcs = [bc_inlet, bc_outlet]


# Define variational problem
if theta == 1:
    F = 1/dt*(p-p_old)*q*dx + k*u.dx(0)*q*dx + c*c/A*w_l*dirac*q*dx + 1/dt*(u-u_old)*v*dx + c*c/k*p.dx(0)*v*dx - c*c/(k*A)*u*w_l*dirac*v*dx


#elif theta == 0:

#else:

L = lhs(F)
R = rhs(F)


# Time loop
for t_ind in range(time_dofs):
    # Get current time
    time = t_ind * dt

    # Solve linear system
    solve(L==R, w_sol, bcs)
    #pdb.set_trace()
    (u_sol, p_sol) = w_sol.split()
    # Plot solution
    # plot(u_sol, key='u', title='Numerical solution')
    # u_sol.rename('Temperature', 'Temperature')
    # Update old solution
    assign(u_old, u_sol)
    assign(p_old, p_sol)

plot(u_sol)
plt.title('Velocity')
plt.show()

plot(p_sol)
plt.title('Pressure')
plt.show()
'''
u = ufun.sub(0).vector().array_view()
u = np.reshape(u,(N+1,N+1))

v = ufun.sub(1).vector().array_view()
v = np.reshape(v,(N+1,N+1))

VV = dolfin.VectorFunctionSpace(mesh, 'CG', 1)

U = dolfin.project(ufun,VV)
u = U.sub(0).compute_vertex_values()
u = np.reshape(u,(N+1,N+1))

v = U.sub(1).compute_vertex_values()
v = np.reshape(v,(N+1,N+1))


#plt.subplot(1, 3, 3)
plt.figure()
plt.imshow(np.power(u,2)+np.power(v,2),origin='lower',interpolation='bilinear')
plt.colorbar()
plt.show()


t0 = time.time()
sol = burgers(nu=1e-3,u_vel=0.5,v_vel=0.5,a=2,N=99,t_end=1.,Nt=100)
t1 = time.time()
print(t1-t0)

x = np.linspace(-1.0,1.0,99+1)
y = np.linspace(-1.0,1.0,99+1)
X,Y = np.meshgrid(x, y) # grid of point
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, np.power(sol[-1,:,:,0],2)+np.power(sol[-1,:,:,1],2), rstride=1, cstride=1,
                      cmap='viridis',linewidth=0, antialiased=False)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
fig.colorbar(surf, shrink=0.5, aspect=15)
plt.show()

plt.figure()
plt.imshow(np.power(sol[-1,:,:,0],2)+np.power(sol[-1,:,:,1],2),origin='lower',interpolation='bilinear')
plt.colorbar()
plt.show()
'''

