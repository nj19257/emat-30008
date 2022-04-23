import numpy as np
import scipy
import matplotlib.pyplot as plt
from arbitary_differential_equation import *
from numerical_shooting_final import shooting
from scipy.optimize import fsolve ,root

def continuation(ode , x0, pars, vary_par_index, vary_par_range, vary_par_number, discretisation, solver=fsolve ):
    """
    Takes a function and applies a specified continuation method over a range of parameter values for a selected
    parameter. Returns the solution of the function for each parameter value.

    :param method: Numerical continuation method to use. 'natural' for natural parameter continuation and 'pseudo' for
    pseudo-arclength continuation.
    :param function: Function to vary the parameter of
    :param x0: Array containing guess for initial state of the function
    :param pars: Array containing the function parameters
    :param vary_par: Index of the parameter to vary in pars
    :param vary_par_range: Array specifying the inclusive range to vary the parameter in the form [a, b] where a is the
    lower bound and b is the upper bound.
    :param vary_par_number: Number of equally spaced values to split the range into
    :param discretisation: Discretisation to be applied to function. Use shootingG for shooting problems
    :param solver: Solver to use for root finding. fsolve is suggested as it has the best performance
    :return: Returns par_list, x where par_list is a list of parameters and x is a list of the respective discretised
    function solutions.
    """
    print(str(discretisation))
    def solution_solver(ode ,x0 ,*args):
        if discretisation == shooting :
            output = shooting(ode, x0, *args )
        else:
            output = discretisation(solver(ode,x0, args))
        return output

    vary_par_list , pars_list = get_parameters_list(pars ,vary_par_index,vary_par_range,vary_par_number)

    solution, test = pseudo_arclength_function(ode ,x0,pars_list,vary_par_index, solution_solver ,discretisation )
    #solution = nature_contunuation_function(ode ,x0,pars_list, solution_solver  )

    print(test)
    #if plot ==True :
    norm_x_nat = scipy.linalg.norm(solution, axis=1, keepdims=True)
    #norm_x_ps = scipy.linalg.norm(x_ps[:, :-1], axis=1, keepdims=True)
    print(norm_x_nat)
    plt.plot(pars_list, norm_x_nat, label='Nat. Param.')
    #plt.plot(par_list_ps, norm_x_ps[:, 0], label='Pseudo-arclength')
    plt.xlabel('c')
    plt.ylabel('||x||')
    plt.legend()
    plt.show()

def get_parameters_list(pars ,vary_par_index,vary_par_range,vary_par_number):
    vary_par_list = np.linspace(vary_par_range[0], vary_par_range[1], vary_par_number)
    n = len(vary_par_list)
    w = np.size(pars)
    pars_list = np.zeros((n ,w))
    for i in range(n):
        pars[vary_par_index] = vary_par_list[i]
        pars_list[i] = pars
    return vary_par_list , pars_list
def nature_contunuation_function(ode ,x0, pars_list , solution_solver):
    n=len(pars_list)
    solutions=[]
    for i in range(n):
        x0 = solution_solver(ode, x0, *pars_list[i])
        print(x0)
        solutions.append(x0)
        x0 = np.round(x0, 2)
    return np.array(solutions)

#def pseudo_arclength_function(ode ,x0, pars_list ,vary_par_index, solution_solver):
def pseudo_arclength_function(ode ,x0, pars_list ,vary_par_index, solution_solver ,discretisation):
    par_0 = pars_list[0]
    par_1 = pars_list[1]
    range1 = pars_list[0]
    range2 = pars_list[-1]
    vary_par_0=par_0[vary_par_index]
    vary_par_1 = par_1[vary_par_index]
    #par_1[vary_par_index] = vary_par_1
    print(par_0)
    pars_list = []
    solutions =[]
    statement = False
    n = len(par_1)
    #for shoot only
    # get the current true value
    print(x0)
    x0 = solution_solver(ode, x0, *par_0)
    x0 = np.round(x0, 2)
    v_0=np.concatenate((x0,vary_par_0), axis=None)

    pars_list.append(par_0)
    solutions.append(x0)
    # let predicted x1 = x0 and get true x1
    x1 = solution_solver(ode, x0, *par_1)
    print(x1)
    v_1 = np.concatenate((x1, vary_par_1), axis=None)
    print(v_0)
    print(v_1)
    pars_list.append(par_1)
    solutions.append(x1)
    find_root = False
    while statement == False:
        #get secant
        secant = v_1-v_0
        print(secant)
        print('y')
        PAeq = lambda v1: np.dot((v1 - predicted_v1), secant)
        def get_par(vary_par_1):
            print(vary_par_1)
            par_1[vary_par_index] = vary_par_1
            print(par_1)
            return par_1


        v_0= v_1
        par_0=par_1

        predicted_v1 = sum_with_roll(v_1, secant)
        print(predicted_v1)
        print(PAeq(predicted_v1))
        #for example v1[u0,vary_par]
        print(solution_solver(ode, predicted_v1[:-1], *get_par(vary_par_1)))
        print('end')

        #predicted_x1 = np.array(x1 + secant_x)
        solve = lambda v1: np.append(
            ode(v1[n:], *v1[:n]),
            PAeq(v1)
        )
        #solve = lambda v1: np.append(
        #    solution_solver(ode, v1[:-1], *get_par(v1[-1])),
        #    PAeq(v1)
        #)
        print(predicted_v1)
        v_1 = fsolve(solve, predicted_v1)
        solutions.append(v_1[n:])
        pars_list.append(v_1[:n])
        print('start')
        print(v_1)
        check_value_x = v_1 - predicted_v1
        if predicted_v1[n:] > range2 or predicted_v1[n:]< range1:
            statement = True
    solutions = np.array(solutions)
    pars_list=np.array(pars_list)
    print(pars_list[:, vary_par_index])
    return solutions, pars_list[:, vary_par_index]
def pseudo_arclength_method(ode ,predicted_x1 ):
    return np.append(solution_solver(ode, predicted_x1, *pars1), dot_product(x0, predicted_x1, secant))

test = lambda u0, *args: np.array([
    *(solution_solver(ode, x1, *pars1)),
    dot_product(x0, x1, secant) ,
    ])   # dx/dt(0) = 0
 #   pseudo_arclength_equation
def sum_with_roll(a, b, offset=0):
    e = np.zeros(a.shape)
    e[tuple(map(slice, b.shape))] = b
    return a + np.roll(e, offset)

def pseudo(x, delta_x, p, delta_p):

    """
    Pseudo-arclength equation. Takes state vector (x), secant of the state vector (delta_x), parameter value (p),
    and secant of the parameter value (delta_p), returns the value of the pseudo-arclength equation for these values

    param x: Function solution (state vector)
    :param delta_x: Secant of the function solution (state vector)
    :param p: Parameter value
    :param delta_p: Secant of the parameter value
    :return: Returns the value of the pseudo-arclength equation for the given values
    """

    # calculates predicted values of x and p by adding the secant
    x_pred = x + delta_x
    p_pred = p + delta_p

    # pseudo-arclength equation
    arc = np.dot(x - x_pred, delta_x) + np.dot(p - p_pred, delta_p)
    return arc

def cubic(x, *args):
    #print(args)
    c = args[0]
    return x ** 3 - x + c

def dot_product(X1, X2, x):
    return np.vdot(secant(X1, X2), (x - x_tilde(X1, X2)))

#print(pseudo_arclength_function(hopf_bifurcation_normal ,[1.4, 0, 6.3], get_parameters_list([2, -1] ,0,[2 , -1],40)[1] ,shooting))
#[ 1.41421356e+00 -7.80704240e-11  2.51327412e+01] hopf normal
continuation(cubic, np.array([1]), [2], 0, [-2, 2], 200,
                                       discretisation=lambda x: x, solver=fsolve)

#continuation(hopf_bifurcation_normal, [1.4, 0, 6.3], [2,-1], 0, [2, 0], 200,
#                                       discretisation=shooting , solver=fsolve)