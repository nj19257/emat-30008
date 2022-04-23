import numpy as np
from scipy import optimize
from solve_ode import solve_ode
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math
from new_numerial_shooting import get_most_repeated_value_index



def predator_prey(r ,t,*args):
    """
    The predator-prey equation function
        Parameter:
            x = x is the number of prey (for example, rabbits).
            y : y is the number of some predator (for example, foxes).
            a ,b, d : are positive real parameters describing the interaction of the two species.

    :return: [dxdt, dydt ]  which represent the instantaneous growth rates of the two populations;
    """
    a = args[0]
    b = args[1]
    d = args[2]
    dxdt=r[0]*(1-r[0])- (a*r[0]*r[1])/(d+r[0])
    dydt=b*r[1]*(1-r[1]/r[0])
    return np.array([dxdt, dydt])

def get_x_y_with_dxdt_equal_zero(a,b,d):
    dxdt=x*(1-x)- (a*x*y)/(d+x)

#Assume b = 0.5 as 0.5 > 0.26
# predator_prey(0.4 ,y ,1 , 0.5,0.1 )

def f(x ,y ,a ,b ,d ):
    #dxdt=r[0]*(1-r[0])- (1*r[0]*r[1])/(0.1+r[0])
    dxdt = x * (1 - x) - (a * x * y) / (d + x)
    #dydt = b * y(1 - y / x)
    return dxdt
def f_der(x ,y ,a ,b ,d ):
    dydt = b * y*(1 - y / x)
    return dydt


def solve(u0, p):
    # solve the combined problem
    u, info, ier, msg = fsolve(predator_prey , u0, args=(p,), full_output=True)
    if ier == 1:
        print("Root finder found the solution u={} after {} function calls; the norm of the final residual is {}".format(u, info["nfev"], np.linalg.norm(info["fvec"])))
        return u
    else:
        print("Root finder failed with error message: {}".format(msg))
        return None
def solve_dxdt(x):
    change=x[:-1]-x[1:]
    asign = np.sign(change)
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
    signchange[0]=0
    index_change= np.where(signchange==1)
    return np.array(index_change)

def get_dxdt_0(x,index_change):
    #dxdt_0=
    n=np.array(index_change)
    print(x[n])
    print(x[n+1])

def get_phase_condition(f,initial_x,target_time ,row,tol,*args , n=200):
    print(args)
    x ,t= solve_ode( f, initial_x , 0 , target_time , 0.001 , 'rk4'  ,*args, n=n ,plot= 'plot_x_y')
    print(x)
    index = solve_dxdt(x[:, row])
    x=x[index,row]
    print(x)
    i=np.size(x)
    print(i)
    for k in range(i):
        diff[k]=x[k+2]-x[k]
        print(diff)







'''f1 = lambda x, y : x * (1 - x) - (1 * x * y) / (0.1 + x)
fder = lambda x, y : 0.5 * y*(1 - y / x)
rng = np.random.default_rng()
x = rng.standard_normal(100)
y= np.arange(-50, 50)
print(optimize.newton(f1,x , fprime=fder, args=(y, ) , maxiter=200))'''
solution,t = solve_ode( predator_prey,[1,1] , 0 , 200 , 0.001 , 'rk4'  ,1, 0.15,0.1, n=2000,plot= 'plot_x_y')
print(solution)
#plt.plot(solution[:,0], solution[:,1])
#plt.ylabel('y')
#plt.xlabel('x')
#plt.show()
#u0=[10,10]
#solve(u0,[1,0.4,0.1])
#print(predator_prey([0.27015621 ,0.27015621] ,1 ,1, 0.4,0.1))
#find dxdt=0= change of x over time = 0
#dxdt = solve_dxdt(solution[:,0])
#get_dxdt_0(solution,dxdt)
#get_phase_condition(predator_prey, [1,1] , 500 , 0, 1e-10, 1, 0.15,0.1, n=1000)