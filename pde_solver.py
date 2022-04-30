import numpy as np
from math import pi
import scipy.sparse.linalg
import matplotlib.pyplot as plt

def finite_difference_method(mx,mt,L,T,k,method ,boundary_condition , initial_condition_func ,l_conditional_value = 0,r_conditional_value = 0 ):
# the boundary condition need to be provided
# ms: the number of gridpoints the space
# mx; The mesh point in space
# initial_condition_func : func to get the initial condition
    if boundary_condition == 'dirichlet':
    #no need change for this type
    #(n-1)(n-1) matrix
        matrix_size = mx -1

    elif boundary_condition == 'periodic':
    # (n)(n) matrix
        matrix_size = mx
    elif boundary_condition == 'neumann':
        matrix_size = mx + 1

    #As calling lambda will be a function. Therefore , use lamda as the name of lambda
    x = np.linspace(0, L, mx + 1)  # mesh points in space
    t = np.linspace(0, T, mt + 1)
    delta_x = x[1] - x[0]
    delta_t = t[1] - t[0]
    lamda = (k*delta_t)/delta_x**2
    # Set up the solution variables
    u_j = np.zeros(x.size)        # u at current time step
    u_jp1 = np.zeros(matrix_size)      # u at next time step
    vec = np.zeros(matrix_size)
    # Set initial condition
    for i in range(0, mx+1):
        u_j[i] = initial_condition_func(x[i])

    if method == 'forward':
        matrix_A = np.identity(matrix_size)
        matrix_B = tridiagonal_matrix(matrix_size,1 - 2*lamda ,lamda , lamda)
    elif method == 'backward':
        matrix_A = tridiagonal_matrix(matrix_size,1 + 2*lamda ,-lamda , -lamda)
        matrix_B = np.identity(matrix_size)
    elif method == 'cranker':
        matrix_A = tridiagonal_matrix(matrix_size, 1 + lamda,-lamda/2 , -lamda/2 )
        matrix_B = tridiagonal_matrix(matrix_size, 1 - lamda, lamda/2, lamda / 2)


    #no need change for dirichlet
    #(n-1)(n-1) matrix

    if boundary_condition == 'periodic':
    # (n)(n) matrix
        matrix_A[0,-1] = matrix_A[0,1]
        matrix_A[-1, 0] = matrix_A[0,1]
        matrix_B[0, -1] = matrix_B[0, 1]
        matrix_B[-1, 0] = matrix_B[0, 1]
        #new_u_j = u_j[:-2] #not sure]
        #new_u_j = np.append(new_u_j, new_u_j[-1])

        #u_j remain the same
    elif boundary_condition == 'neumann':
        matrix_A[0,1] = 2*matrix_A[0,1]
        matrix_A[-1, -2] = 2 * matrix_A[-1, -2]
        matrix_B[0, 1] = 2 * matrix_B[0, 1]
        matrix_B[-1, -2] = 2 * matrix_B[-1, -2]
        print(matrix_B)
        #new_u_j = u_j
    def get_u_j(u_j ,boundary_condition):
        if boundary_condition == 'dirichlet':
            new_u_j = u_j[1:-1]
        elif boundary_condition == 'periodic':
            new_u_j = u_j[:-2]  # not sure]
            new_u_j = np.append(new_u_j, new_u_j[-1])
        elif boundary_condition == 'neumann':
            new_u_j = u_j
        return new_u_j
    def solve_pde_step(u_j,matrix_A,matrix_B,method,boundary_condition):


        if boundary_condition == 'dirichlet':
            vec[0] = lamda* l_conditional_value
            vec[-1] = lamda* r_conditional_value
            new_u_j = u_j[1:-1]
            l_n = 1
            r_n = -1
        elif boundary_condition == 'periodic':
            new_u_j = u_j[:-1]  # not sure]

            l_n = None
            r_n = -1
        elif boundary_condition == 'neumann':
            vec[0] = 2*lamda *delta_x * -l_conditional_value
            vec[-1] = 2*lamda* delta_x * r_conditional_value
            new_u_j = u_j

            l_n = None
            r_n = None

        u_jp1 = scipy.sparse.linalg.spsolve(matrix_A, np.matmul(matrix_B,new_u_j)+vec)
        if boundary_condition == 'neumann':
            u_jp1[0] = l_conditional_value
            u_jp1[-1] = r_conditional_value
        elif boundary_condition == 'periodic':
            u_jp1[0] = l_conditional_value
            u_j[-1] = r_conditional_value
        u_j[l_n:r_n] = u_jp1
        return u_j



    u_j[0] = l_conditional_value
    u_j[-1] = r_conditional_value
    for j in range(0, mt):
        # Forward Euler timestep at inner mesh points
        # PDE discretised at position x[i], time t[j]
        # Boundary conditions
        u_j = solve_pde_step(u_j, matrix_A, matrix_B, method, boundary_condition)

    print(u_j)
    return x , u_j
def tridiagonal_matrix(n,diagonal,upper_diagonal,lower_diagonal ):
    tridiagonal_matrix =np.eye(n, n, k=-1)*lower_diagonal + \
                        np.eye(n, n) * diagonal + np.eye(n, n, k=1) * upper_diagonal
    return tridiagonal_matrix


def u_I(x):
    # initial temperature distribution
    y = np.sin(pi*x/L)
    return y

def u_exact(x,t ,kappa , L):
    # the exact solution
    y = np.exp(-kappa*(pi**2/L**2)*t)*np.sin(pi*x/L)
    return y
kappa = 1.0   # diffusion constant
L=2.0         # length of spatial domain
T=0.5         # total time to solve for
mx = 50     # number of gridpoints in space
mt = 1000   # number of gridpoints in time

# Solve for each method
c_x, c_u = finite_difference_method(mx, mt, L, T, kappa, "backward",'neumann', u_I )

"""
plotting exact solution
"""
xx = np.linspace(0, L, 250)
exact = u_exact(xx, T, kappa, L)
plt.plot(xx, exact, label='exact')

plt.plot(c_x, c_u, label='crank')

plt.legend()
plt.xlabel('x')
plt.ylabel('u(x,0.5)')
plt.show()