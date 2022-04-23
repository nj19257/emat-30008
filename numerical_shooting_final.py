import numpy as np
from scipy import optimize
from solve_ode_finial import solve_ode
from scipy.optimize import fsolve
from arbitary_differential_equation import *
from root_finding_final import newton
import matplotlib.pyplot as plt
import sys
import warnings

def single_preiodic_orbit(t):
    output, i=get_most_repeated_value_index(np.diff(t))
    print(output)
    return  output

def initial_guesss(f, x0 , target_time ,*args):
    n = target_time*20
    solution, t = solve_ode(f, x0, target_time, 'rk4', *args , n=n)
    try:
        #get most repeated x values
        output, i = get_most_repeated_value_index(solution[:,0])
        #get time interval between most repeated x values
        cycle_time = np.random.choice(single_preiodic_orbit(t[i]))
    except ValueError:
        print('The values of input is not enough to find a good initial guess , try to input a larger target time')
        sys.exit([])
    values_guess = np.append( solution[np.random.choice(i[0])] , cycle_time)
    print(solution[np.random.choice(i[0])])
    print(values_guess)
    return values_guess

def get_most_repeated_value_index(input,dp=5):
    input = np.around(input, dp)
    unique, counts = np.unique(input, return_counts=True )
    #print(counts)
    output = unique[np.where(counts == np.max(counts))]
    #print(output)
    index=np.where(1122222222222222222222222222222200000000000000000000000000000 == output)
    #print(index)
    return output,np.array(index)


def shooting(ode, u0,*args , method='fsolve' ,plot=False):
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


    g = lambda u0, *args: np.array([
        *(u0[:-1] - solve_ode(ode, u0[:-1],  u0[-1],  'rk4', *args)[0][-1]),
        ode(u0[:-1], 1, *args)[0] ,  # dx/dt(0) = 0
    ])
    warnings.simplefilter("error")
    try :
        if method == 'fsolve':
            root = fsolve(g, u0, args)
        elif method == 'newton':
            root = newton(g,u0,*args)
        else:
            print('not valid input for method')
            root = np.array([])
        if plot == True:
            print(root[:-1])
            print(root[-1])
            solution, t = solve_ode(ode, root[:-1], root[-1], 'rk4', *args ,n=1000)
            print(solution)
            plt.plot(root[0], root[1], 'go', label="Manually found orbit" )
            plt.plot(solution[:, 0], solution[:, 1])
            plt.ylabel('y')
            plt.xlabel('x')
            plt.legend()
            plt.show()
    except ValueError:
            print('test')
            root = np.array([])

    return root



def predator_prey(r, t, *args):
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

    dxdt = r[0] * (1 - r[0]) - (a * r[0] * r[1]) / (d + r[0])
    dydt = b * r[1] * (1 - r[1] / r[0])
    return np.array([dxdt, dydt])


def main():
    # plot x_y diagram when a=1 , d = 0.1 and b=0.45 (Thus, b >0.26)
    a = 1
    d = 0.1
    b = 0.15
    #plot figure of x againist y with b= 0.15
    #solution1, t1 = solve_ode(predator_prey, [1, 1], 500, 'rk4', a, b, d, n=500, plot='plot_x_y')

    # plot figure of x againist y with b= 0.45
    b = 0.45
    #solution, t = solve_ode(predator_prey, [1, 1], 100, 'rk4', a, b, d, n=100, plot='plot_x_y')
    a = 1
    d = 0.1
    parameter = (a, 0.15, d)
    #shooting(hopf_bifurcation_normal, [0.01, 0.01, 3, 4], *(2, -1))
    solution, t = solve_ode(hopf_bifurcation_normal, [1, 1], 100, 'rk4', 2,-1, n=1000, plot='plot_x_y')
    initial_guess = initial_guesss(hopf_bifurcation_normal1 , [1.4, 0] ,1500,  2      )
    # initial_guess =np.append(result1,25.05)
    print(initial_guess)
    print(initial_guess)
    print(shooting(hopf_bifurcation_normal, initial_guess, *( 0.02020202 ,-1.  ) ,method = 'fsolve', plot=True))
    print(shooting(predator_prey, [0.7,0.3,100], *parameter ,method = 'fsolve', plot=True))


if __name__ == "__main__":
    main()