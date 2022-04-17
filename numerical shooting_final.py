import numpy as np
from scipy import optimize
from solve_ode_finial import solve_ode
from scipy.optimize import fsolve
from arbitary_differential_equation import *
import matplotlib.pyplot as plt

def single_preiodic_orbit(f, initial_guess, target_time, *args, n=500, plot=False ):
    print(args)
    #Numerically find root so that x(0)-x(T) = 0
    #get list of x(t)
    solution, t = solve_ode(f, initial_guess, target_time, 'rk4', *args, n=n)
    n=len(solution)
    print(sum(solution[0]))
    #root finding
    timecycle= 0
    #find x(T)
    for i in range(1,n):
        if sum(abs(solution[i] - initial_guess)) < 1e-7:
            print(solution[i])
            timecycle = t[i]
            print(timecycle)
    if plot == True:
        plt.plot(initial_guess[0], initial_guess[1], 'go', label="Manually found orbit")
        plt.plot(solution[:,0], solution[:,1])
        plt.ylabel('y')
        plt.xlabel('x')
        plt.legend()
        plt.show()
    return timecycle


def shooting(ode, u0,*args ,plot=False):
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
    print(args[0])
    tf = fsolve_roots[-1]
    print(fsolve_roots)
    timecycle = single_preiodic_orbit(predator_prey, fsolve_roots[:-1], tf, *args[0], n= int(tf*10) , plot = plot)

    return [fsolve_roots[0] , fsolve_roots[1] , timecycle]



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
    # plot figure of x againist y with b= 0.15
    solution1, t1 = solve_ode(predator_prey, [1, 1], 500, 'rk4', a, b, d, n=500, plot='plot_x_y')

    # plot figure of x againist y with b= 0.45
    b = 0.45
    solution, t = solve_ode(predator_prey, [1, 1], 100, 'rk4', a, b, d, n=100, plot='plot_x_y')
    a = 1
    d = 0.1
    parameter = (a, 0.15, d)
    # initial_guess =np.append(result1,25.05)
    # print(initial_guess.tolist())
    shooting(predator_prey, [0.71072924, 0.23452027, 23], parameter , plot=True)


if __name__ == "__main__":
    main()