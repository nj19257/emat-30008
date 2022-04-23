import numpy as np
from scipy import optimize
from solve_ode_finial import solve_ode
from scipy.optimize import fsolve
from arbitary_differential_equation import *
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
    a = args[0]
    b = args[1]
    d = args[2]
    dxdt=r[0]*(1-r[0])- (a*r[0]*r[1])/(d+r[0])
    dydt=b*r[1]*(1-r[1]/r[0])
    return np.array([dxdt, dydt])

def get_most_repeated_value_index(input,dp=5):
    counts=0
    #to ensure it will only give up will a
    while np.max(counts) < 1 :
        print(dp)
        rounded_input = np.around(input, dp)
        unique, counts = np.unique(rounded_input, return_counts=True)
        print(counts)
        output = unique[np.where(counts == np.max(counts))]
        print(output)
        print(np.max(counts))
        index=np.where(np.in1d(rounded_input , output))
        print(index)
        dp = dp - 1
    return output, np.array(index)


def single_preiodic_orbit(f, initial_guess, target_time, *args, n=500):
    print(args)
    solution, t = solve_ode(f, initial_guess, target_time, 'rk4', *args, n=n)
    print(solution)
    o, i = get_most_repeated_value_index(solution[:, 0])
    print(i)
    output, i = get_most_repeated_value_index(np.diff(t[i]))
    return output

def solve_dxdt(x,ignore_first=True):
    change=x[:-1]-x[1:]
    asign = np.sign(change)
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
    if ignore_first == True :
        signchange[0]=0
    index_change= np.where(signchange==1)
    print(index_change)
    return np.array(index_change)

def get_phase_condition(f,initial_x,target_time ,row, tol ,*args , n=200):
    print(args)
    dp=5
    i_odd= []
    i_even=[]
    x ,t= solve_ode( f, initial_x , target_time ,  'rk4'  ,*args, n=n )
    t_interval=t[3]-t[0]
    index = solve_dxdt(x[:, row])
    #correct the index as the currnt index is after it change direction
    index = index -1
    x_odd = x[index[0, 1::2]]
    print(x_odd)
    x_drop, i_odd = get_most_repeated_value_index(x_odd[:,0], dp=5)
    print(i_odd)
    #print(np.random.choice(i_odd[0]))
    x_odd, t = solve_ode(f, x_odd[np.random.choice(i_odd[0])],  t_interval,  'rk4', *args, n=n)
    x_even = x[index[0, ::2]]
    x_drop, i_even = get_most_repeated_value_index(x_even[:,0], dp=7)
    print(i_even)
    x_even, t = solve_ode(f, x_even[np.random.choice(i_even[0])], t_interval,  'rk4', *args, n=n)
    #x, t = solve_ode(f, x_odd[random.choice(i_odd)], 0, target_time, 0.001, 'rk4', *args, n=n)
    finial_odd_i = solve_dxdt(x_odd[:, row],ignore_first=False)
    finial_even_i = solve_dxdt(x_even[:, row],ignore_first=False)
    print(finial_even_i)
    print(x_even)
    phase_condition1=x_odd[finial_odd_i[0,-1]]
    phase_condition2=x_even[finial_even_i[0,-1]]
    print('The phase condition is {} and {}'.format(phase_condition1, phase_condition2))

    return phase_condition1, phase_condition2


def shooting(ode, u0,*args):
    """
    A function that uses numerical shooting to find limit cycles of
    a specified ODE.

    Parameters
    ----------
    ode : function
        The ODE to apply shooting to. The ode function should take
        a single parameter (the state vector) and return the
        right-hand side of the ODE as a numpy.array.
    u0 : numpy.array
        An initial guess at the initial values for the limit cycle.
    args : Any extra arguments .

    Returns
    -------
    Returns a numpy.array containing the corrected initial values
    for the limit cycle. If the numerical root finder failed, the
    returned array is empty.
    """
    # Here is the code that does the shooting
    g = lambda u0, *args: [
        *(u0[:-1] - solve_ode(ode, u0[:-1],  u0[-1],  'rk4', *args)[0][-1]),
        ode(u0[:-1], 1, *args)[0] ,  # dx/dt(0) = 0
    ]
    fsolve_roots = fsolve(g, u0, *args)
    print(fsolve_roots)
def main():
    #plot x_y diagram when a=1 , d = 0.1 and b=0.45 (Thus, b >0.26)
    a=1
    d=0.1
    b=0.15
    #plot figure of x againist y with b= 0.15
    solution1,t1 = solve_ode( predator_prey,[1,1] ,  500 ,  'rk4'  ,a, b, d, n=500 , plot='plot_x_y')

    #plot figure of x againist y with b= 0.45
    b=0.45
    solution,t = solve_ode( predator_prey,[1,1] ,  100 ,  'rk4'  ,a, b, d, n=100 , plot='plot_x_y')

    #get starting condition
    print(single_preiodic_orbit(predator_prey, [1,1], 200, 1,0.15,0.1, n=400))


    plt.plot(solution1[58,0], solution1[58,1], 'go', label="Manually found orbit")
    plt.plot(solution1[:,0], solution1[:,1])
    plt.ylabel('y')
    plt.xlabel('x')
    plt.legend()
    plt.show()
    print(solution1[58])
    starting_condition=solution1[58]
    print(starting_condition.tolist())
    #get phase condition

    #print(i[0,::2])
    #x=solution1[i[0,::2]]
    #print(get_most_repeated_value_index(x[:,0]))
    result1,result2 = get_phase_condition(predator_prey, starting_condition.tolist() ,1000 , 0, 1e-10, 1, 0.15,0.1, n=5000)

    #root finding
    parameter=(a,0.15,d)
    initial_guess =np.append(result1,25.05)
    print(initial_guess.tolist())
    shooting(predator_prey, initial_guess.tolist(),parameter)

    #testing arbitrary equation
    shooting(d2xdt2_equals_minus_x, [1,1,10])

if __name__ == "__main__":
    main()