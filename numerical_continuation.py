import numpy as np
import scipy
import matplotlib.pyplot as plt
from arbitary_differential_equation import *
from numerical_shooting import shooting
from scipy.optimize import fsolve ,root
from solve_ode import solve_ode
import warnings

def numerical_continuation(funs , x0, pars, vary_par_index, vary_par_range, vary_par_number, discretisation, solver=fsolve ,method = 'natural' , plot = False):
    """
    Takes a function and applies a specified continuation method over a range of parameter values for a selected
    parameter. Returns the solution of the function for each parameter value.
    :param funs: Function to vary the parameter of
    :param x0: Array containing guess for initial state of the function
    :param pars: Array containing the function parameters
    :param vary_par_index: Index of the parameter to vary in pars
    :param vary_par_range: Array specifying the inclusive range to vary the parameter in the form [a, b]
    :param vary_par_number: Number of equally spaced values to split the range into
    :param discretisation: Discretisation to be applied to function. Use shooting for shooting problems
    :param solver: Solver to use for root finding. fsolve is suggested as it has the best performance
    :param method: Numerical continuation method to use. 'natural' for natural parameter continuation and 'pseudo'
    :param plot: The statement of plot figure or not
    pseudo-arclength continuation.



     :return: Returns par_list, x where par_list is a list of parameters and x is a list of the respective discretised
    function solutions.
    """

    def solution_solver(ode ,x0 ,*args):
        if discretisation == shooting :
            output = shooting(ode, x0, *args )
            #output = shootingG(funs)
        else:
            output = solver(discretisation(ode),x0, args)
        return output
    def get_parameters_list(pars ,vary_par_index,vary_par_range,vary_par_number):
        vary_par_list = np.linspace(vary_par_range[0], vary_par_range[1], vary_par_number)
        n = len(vary_par_list)
        w = np.size(pars)
        pars_list = np.zeros((n ,w))
        for i in range(n):
            pars[vary_par_index] = vary_par_list[i]
            pars_list[i] = pars
        return vary_par_list , pars_list
    def nature_contunuation_function(funs ,x0, pars_list , solution_solver ):
        n=len(pars_list)
        solutions=[]
        warnings.simplefilter("ignore")
        for i in range(n):
            x0 = solution_solver(funs, x0, *pars_list[i])
            solutions.append(x0)
            x0 = np.round(x0, 2)
        return np.array(solutions)

    def pseudo_arclength_function(funs ,x0, pars_list ,vary_par_index, solution_solver ,discretisation):
        #initially used fsolve twice
        par_0 = pars_list[0]
        par_1 = pars_list[1]
        range_min = min([pars_list[0][vary_par_index] , pars_list[-1][vary_par_index]])
        range_max = max([pars_list[0][vary_par_index] , pars_list[-1][vary_par_index]])
        vary_par_0=par_0[vary_par_index]
        vary_par_1 = par_1[vary_par_index]
        pars_list = []
        solutions =[]
        statement = False
        n = len(par_0)
        # get the current true value
        x0 = np.round(x0, 2)
        x0 = solution_solver(funs, x0, *par_0)
        v_0=np.concatenate((x0,vary_par_0), axis=None)
        pars_list.append(vary_par_0)
        solutions.append(x0)
        # let predicted x1 = x0 and get true x1
        x1 = solution_solver(funs, np.round(x0, 2), *par_1)
        v_1 = np.concatenate((x1, vary_par_1), axis=None)
        pars_list.append(vary_par_1)
        solutions.append(x1)
        find_root = False
        g = shootingG(funs)
        def get_par(vary_par_1):
            par_1[vary_par_index] = vary_par_1
            return par_1
        if discretisation == shooting:
            solve = lambda v1: np.append(
                g(v1[:-1], *get_par(v1[-1])), #np.round added if error arise this might be problem
                PAeq(v1)
            )
        else:
            solve = lambda v1: np.append(
                funs(v1[:-1], *get_par(v1[-1])),
                PAeq(v1)
            )

        while statement == False:
            secant = v_1 - v_0
            PAeq = lambda v1: np.dot((v1 - predicted_v1), secant)
            v_0 = v_1
            par_0 = par_1
            predicted_v1 = v_1+ secant
            predicted_v1[:-1] = np.round(predicted_v1[:-1],2) #this might be the problem
            v_1 = solver(solve, predicted_v1)
            solutions.append(v_1[:-1])
            pars_list.append(v_1[-1])
            check_value_x = v_1 - predicted_v1
            if v_1[-1] > range_max or v_1[-1] < range_min:
                statement = True

        solutions = np.array(solutions)
        pars_list=np.array(pars_list)
        return solutions, pars_list
    #get the list of time in the time range
    vary_par_list , pars_list = get_parameters_list(pars ,vary_par_index,vary_par_range,vary_par_number)

    if method == 'natural' :
        solution = nature_contunuation_function(funs, x0, pars_list, solution_solver )
        label = 'natural parameter continuation'
    elif method == 'pseudo-arclength':
        #to replace the original vary_par_list
        solution, vary_par_list = pseudo_arclength_function(funs, x0, pars_list, vary_par_index, solution_solver, discretisation)
        label = 'pseudo-arclength continuation'

    if plot ==True :
        #this mean that it is a ode and the solution will be included will a timecycle
        if discretisation == shooting:
            norm_x = scipy.linalg.norm(solution[:, :-1], axis=1, keepdims=True)
        else:
            norm_x = scipy.linalg.norm(solution, axis=1, keepdims=True)
        plt.plot(vary_par_list, norm_x, label=label)
        plt.xlabel('c')
        plt.ylabel('||x||')
        plt.legend()
        plt.show()
    return solution, vary_par_list




def shootingG(ode):
    """
    Constructs the shooting root finding problem for given ODE
    :param ode: Function defining ODE(s) to solve in the form f(t,x,*args) which returns derivative value at (t,x)
    :return: Returns the function, G,  whose root solves the shooting problem.
    """
    #g = lambda u0, *args: np.array([
    #    *(u0[:-1] - solve_ode(ode, u0[:-1],  u0[-1],  'rk4', *args)[0][-1]),
    #    ode(u0[:-1], 1, *args)[0] ,  # dx/dt(0) = 0
    #])

    def G(u0,*args , phase_condition = lambda u0 ,*args: ode(u0[:-1], 1, *args)[0] ):
        return np.array([
        *(u0[:-1] - solve_ode(ode, u0[:-1],  u0[-1],  'rk4', *args)[0][-1]),
        phase_condition(u0,*args) ,  # dx/dt(0) = 0
    ])


    return G
def cubic(x, *args):
    #print(args)
    c = args[0]
    return x ** 3 - x + c

def main():
    #tested all three equation and work with no problems
    #print(pseudo_arclength_function(hopf_bifurcation_normal ,[1.4, 0, 6.3], get_parameters_list([2, -1] ,0,[2 , -1],40)[1] ,shooting))
    #[ 1.41421356e+00 -7.80704240e-11  2.51327412e+01] hopf normal
    #Using cubic equation
    numerical_continuation(cubic, np.array([1,1]), [2], 0, [-2, 2], 200,
                                          discretisation=lambda x: x, solver=fsolve ,method = 'natural' , plot=True)
    numerical_continuation(cubic, np.array([1,1]), [2], 0, [-2, 2], 200,
                                           discretisation=lambda x: x, solver=fsolve ,method = 'pseudo-arclength', plot=True)
    #Using Hopf bifurcation nor form equation
    numerical_continuation(hopf_bifurcation_normal, [1.4, 0, 6.3], [2, -1], 0, [2, 0], 50,
                                           discretisation=shooting , solver=fsolve ,method = 'natural', plot=True)
    numerical_continuation(hopf_bifurcation_normal, [1.4, 0, 6.3], [2, -1], 0, [2, 0], 50,
                                           discretisation=shooting , solver=fsolve ,method = 'pseudo-arclength', plot=True)
    numerical_continuation(hopf_bifurcation_normal, [1.4, 0, 6.3], [2,-1], 0, [2, 0], 50,
                                           discretisation=shooting , solver=fsolve ,method = 'pseudo-arclength')
    #Using modified_hopf_bifurcation_normal
    numerical_continuation(modified_hopf_bifurcation_normal, [1.4, 0 ,19], [2], 0, [2, -1], 50,
                                           discretisation=shooting , solver=fsolve ,method = 'natural', plot=True)
    numerical_continuation(modified_hopf_bifurcation_normal, [1.4, 0, 6.3], [2], 0, [2, -1], 50,
                                           discretisation=shooting , solver=fsolve ,method = 'pseudo-arclength', plot=True)
if __name__ == "__main__":
    main()