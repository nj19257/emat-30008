import numpy as np
from scipy import optimize
from solve_ode import solve_ode
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import random


def predator_prey(r , t , *args):
    """
    The predator-prey equation function
        Parameter:
            x = x is the number of prey (for example, rabbits).
            y : y is the number of some predator (for example, foxes).
            a ,b, d : are positive real parameters describing the interaction of the two species.

    :return: [dxdt, dydt ]  which represent the instantaneous growth rates of the two populations;
    """
    #print(args)
    a = args[0]
    b = args[1]
    d = args[2]
    dxdt=r[0]*(1-r[0])- (a*r[0]*r[1])/(d+r[0])
    dydt=b*r[1]*(1-r[1]/r[0])
    return np.array([dxdt, dydt])

def solve(u0, *args):
    """
    Used to find saddle point
    :param u0:
    :param p:
    :return:
    """
    t=(1,) #let t=1
    args= t +args
    # solve the combined problem
    u, info, ier, msg = fsolve(predator_prey , u0, args=args, full_output=True)
    if ier == 1:
        print("Root finder found the solution u={} after {} function calls; the norm of the final residual is {}".format(u, info["nfev"], np.linalg.norm(info["fvec"])))
        return u
    else:
        print("Root finder failed with error message: {}".format(msg))

def single_preiodic_orbit(t):
    output, i=get_most_repeated_value_index(np.diff(t))
    print(output)
    return  output

def get_most_repeated_value_index(input,dp=5):
    input = np.around(input, dp)
    unique, counts = np.unique(input, return_counts=True )
    #print(counts)
    output = unique[np.where(counts == np.max(counts))]
    #print(output)
    index=np.where(input == output)
    #print(index)
    return output,np.array(index)

def solve_dxdt(x,ignore_first=True):
    change=x[:-1]-x[1:]
    asign = np.sign(change)
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
    if ignore_first == True :
        signchange[0]=0
    index_change= np.where(signchange==1)
    print(index_change)
    return np.array(index_change)

def get_phase_condition(f,initial_x,target_time ,row,tol,*args , n=200):
    print(args)
    x ,t= solve_ode( f, initial_x , 0 , target_time , 0.001 , 'rk4'  ,*args, n=n )
    t_interval=t[3]-t[0]
    index = solve_dxdt(x[:, row])
    #correct the index as the currnt index is after it change direction
    index = index -1
    x_odd = x[index[0, 1::2]]
    x_drop, i_odd = get_most_repeated_value_index(x_odd[:,0], dp=4)
    print(i_odd)
    #print(np.random.choice(i_odd[0]))
    x_odd, t = solve_ode(f, x_odd[np.random.choice(i_odd[0])], 0, t_interval, 0.001, 'rk4', *args, n=n)
    x_even = x[index[0, ::2]]
    x_drop, i_even = get_most_repeated_value_index(x_even[:,0], dp=4)
    x_even, t = solve_ode(f, x_even[np.random.choice(i_even[0])], 0, t_interval, 0.001, 'rk4', *args, n=n)
    #x, t = solve_ode(f, x_odd[random.choice(i_odd)], 0, target_time, 0.001, 'rk4', *args, n=n)
    finial_odd_i = solve_dxdt(x_odd[:, row],ignore_first=False)
    finial_even_i = solve_dxdt(x_even[:, row],ignore_first=False)
    print(finial_even_i)
    print(x_even)
    phase_condition1=x_odd[finial_odd_i[0,-1]]
    phase_condition2=x_even[finial_even_i[0,-1]]

    return phase_condition1, phase_condition2

def get_test_most_repeated_value(f ,x0 , target_time ):
    solution, t = solve_ode(predator_prey, [1, 1], 0, 500, 0.001, 'rk4', a, b, d, n=500)

def initial_guesss(f, x0 , target_time ):
    n = target_time*10
    solution, t = solve_ode(predator_prey, x0, 0, target_time, 0.001, 'rk4', 1, 0.15, 0.1, n=n)
    #get most repeated x values
    output, i = get_most_repeated_value_index(solution[:,0])
    print(output)
    print(i)
    #get time interval between most repeated x values
    cycle_time = single_preiodic_orbit(t[i])
    values_guess = np.append( solution[np.random.choice(i[0])] , cycle_time)
    return values_guess

"""

u0=[1,1]
saddle_point = solve(u0,1,0.15,0.1)
print(saddle_point)
#plot x_y diagram when a=1 , d = 0.1 and b=0.45 (Thus, b >0.26)
a=1
d=0.1
b=0.15
#plot figure of x againist y with b= 0.15
solution1,t1 = solve_ode( predator_prey,[1,1] , 0 , 500 , 0.001 , 'rk4'  ,a, b, d, n=500 , plot='plot_x_y')
def unique_rows(data):
    uniq = np.unique(data.view(data.dtype.descr * data.shape[1]))
    return uniq.view(data.dtype).reshape(-1, data.shape[1])
#plot figure of x againist y with b= 0.45"""
initial_guesss(predator_prey,[1,1],500)

"""

a=1
d=0.1
b=0.15
#plot figure of x againist y with b= 0.15
solution1,t1 = solve_ode( predator_prey,[1,1] , 0 , 500 , 0.001 , 'rk4'  ,a, b, d, n=500 , plot='plot_x_y')
o,i=get_most_repeated_value_index(solution1[:,0])
#print(solution1)
#need to rewrite this part in word
print('y')
print(o)
print(np.diff(t1[i]))
print(single_preiodic_orbit(t1[i]))
plt.plot(solution1[58,0], solution1[58,1], 'go', label="Manually found orbit")
plt.plot(solution1[:,0], solution1[:,1])
plt.ylabel('y')
plt.xlabel('x')
plt.legend()
plt.show()"""
"""
#start Q2 get phase_condition
i=solve_dxdt(solution1[:,0])

print(i[0,::2])
x=solution1[i[0,::2]]
print(get_most_repeated_value_index(x[:,0]))
print(x)
result1,result2 = get_phase_condition(predator_prey, [1,1] , 1000 , 0, 1e-10, 1, 0.15,0.1, n=5000)
print(result1)
print(result2)
print(solve_ode( predator_prey,[0.06864949 ,0.15707303] , 0 , 25.0501002 , 0.001 , 'rk4'  ,a, 0.15, d, n=100 , plot='plot_x_y')) """

def shooting(f,U,*args):
    g = lambda U, *args: [
        *(U[:-1] - solve_ode(f, U[:-1], 0, U[-1], 0.001, 'rk4', *args)[0][-1]),
        f(U[:-1], 1, *args)[0],  # dx/dt(0) = 0
    ]
    fsolve_roots = fsolve(g, initial_guess, parameter)
    print(fsolve_roots)

a=1
d=0.1
b=0.15
parameter=(a,b,d)
initial_guess = [0.06861238 , 0.15582958,25.05]
#print(solve_ode( predator_prey,[1,1] , 0 , 100 , 0.001 , 'rk4'  ,a, b, d, n=100 , plot='plot_x_y')[0][-1])
g = lambda U,*args : [
    *(U[:-1] - solve_ode(predator_prey,U[:-1], 0 , U[-1] , 0.001 , 'rk4' ,*args )[0][-1]),
    predator_prey(U[:-1],1,*args)[0] ,# dx/dt(0) = 0
]

G = lambda U : [
    *(U[:-1] - solve_ode(predator_prey_test,U[:-1], 0 , U[-1] , 0.001 , 'rk4'  )[0][-1]),
    U[0]*(1-U[0])- (1*U[0]*U[1])/(0.1+U[0]) ,# dx/dt(0) = 0
]


def predator_prey_test(r , t ):
    """
    The predator-prey equation function
        Parameter:
            x = x is the number of prey (for example, rabbits).
            y : y is the number of some predator (for example, foxes).
            a ,b, d : are positive real parameters describing the interaction of the two species.

    :return: [dxdt, dydt ]  which represent the instantaneous growth rates of the two populations;
    """
    #print(args)
    a = 1
    d = 0.1
    b = 0.15
    dxdt=r[0]*(1-r[0])- (a*r[0]*r[1])/(d+r[0])
    dydt=b*r[1]*(1-r[1]/r[0])
    return np.array([dxdt, dydt])



