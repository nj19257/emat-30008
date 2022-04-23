# importing packages
import numpy as np
import math
import matplotlib.pyplot as plt
from arbitary_differential_equation import *
import time

def euler_step(f, x, t, h,*args): #work
    """
     single Euler step function
        Parameter:
            f(function) = function used ( need to define more detail)
            x_current : the current value of x
            t_current : the current timestep
            h : interval between step (step size) (assumed to be constant for the sake of simplicity) is then given by t_{n+1} = h + t_n
            args : Any extra arguments .


    :return:
    [ x_{n+1} , t_{n+1} ] or [x_new , t_new] , which it the x and t for the next step
    """
    x_new = x + h * f(x, t,*args)
    t_new = t + h
    return [x_new, t_new]

def RK4_step(f, x, t, h,*args): #work
    """
     single 4th-order Runge-Kutta step function
        Parameter:
            f(function) = function used ( need to define more detail)
            x_current : the current value of x
            t_current : the current timestep
            h : interval between step (step size) (assumed to be constant for the sake of simplicity) is then given by t_{n+1} = h + t_n
            args : Any extra arguments .

    :return:
    [ x_{n+1} , t_{n+1} ] or [x_new , t_new] , which x_new and t_new are the value of x and t for the next timestep
    """

    k1 = f(x, t,*args)
    k2 = f(x + h * (k1/2) , t + h / 2, *args)
    k3 = f(x + h * (k2/2), t + h / 2, *args)
    k4 = f(x + h * k3, t + h, *args)
    x_new = x + (1 / 6) * h * (k1 + 2 * k2 + 2 * k3 + k4)
    t_new = t + h
    return [x_new, t_new]

def solve_to( f ,x , t , t1 , deltat_max , Method,*args): # work
    """
     Function which solves from x_n , t_n to x_{n+1} , t_{n+1} in steps no bigger than deltat_max
        Parameter:
            f(function) = function used ( need to define more detail)
            x_current(x) : the current value of x
            t_current(t) : the current timestep
            t1 : the value of next timestep based on the current timestep
            deltat_max : The maximum timestep allowed
            Method : Which method to solve the ODE
            args : Any extra arguments .
    :return:
    x : the x value for the target timestep
    """
    Method=str.lower(Method) #to ignore error cause by upper_case input
    #Determine which method to be used
    if Method == 'euler':
        method=euler_step
    if Method == 'rk4':
        method = RK4_step

    while (t + deltat_max) < t1:  # steps while t value is < t2
        x, t = method(f, x, t, deltat_max, *args)
    else:
        x, t = method(f, x, t, t1 - t, *args)  # bridges gap between last step and t2 using step size t2-t_prev_step
    return x

def solve_ode( f , x ,  t_end ,  Method  , *args, deltat_max =1e-3 ,n=5 ,t_0=0 , plot=False , timer=False):
    """
    The function that use solve_to by a first-order numerical procedure for ODE (such as forward euler and 4th-order Runge-Kutta method)
    to generate a series of numerical solution estimates x1,x2,x3,...
        Parameter:
            f(function) = function used ( need to define more detail)
            x_current : the current value of x
            t_0 : the initial timestep (default start time is 0 )
            t_end : the target timestep
            Method : Which method to use
            deltat_max : The maximum timestep allowed (default value is 1e-3 )
            n : number of step to reach the target timestep(the default timesteps is 5)
            args : Any extra arguments .

    :return:
    x_list : the values of x corresponding to t(time)
    t : the time of x in x_list
    """
    if timer == True:
        start = time.time()
    t = np.linspace(t_0, t_end, n )
    t_n = len(t)
    x_n = np.size(x)
    x_list = np.zeros((t_n,x_n))
    x_list[0]= x

    for i in range(t_n-1):
        x_list[i+1]=solve_to(f, x_list[i], t[i], t[i+1], deltat_max, Method,*args)

    if timer == True:
        end = time.time()
        print("Time taken for %s = %f sec" % (Method, end - start))
    #plot figure of x against y
    if plot == 'plot_x_y':
        plot_x_y(x_list, Method)
    # plot figure of x against t
    elif plot == 'plot_x_t':
        plot_x_t(x_list[:,0],t, Method, 'x')
    # plot figure of y against t
    elif plot == 'plot_y_t':
        plot_y_t(x_list[:,1],t, Method, 'y')
    elif plot == False:
        # No changes are made
        pass
    else:
        print('The input string for plot parameter is incorrect')
    return x_list ,t





def error_depend_deltat_t( f ,x  , t_end ,h, Method , *args , plot=False , timer=False):
    """
    Function that find out the error depending on a range of delta t
        Parameter:
            f(function) = function used ( need to define more detail)
            x_current : the current value of x
            t_end : the target timestep
            h : The list of values for maximum timesteps allowed for testing
            Method : Which method to use (Forward Euler or RK4 Method)
            n : number of step to reach the target timestep(the default timesteps is 5)


    :return:
    list_error =
    """
    if timer == True:
        start = time.time()

    list_error = np.zeros(len(h))
    for i in range(len(h)):
        x_list ,t =solve_ode(f, x , t_end , Method,*args , deltat_max=h[i])
        if f == dxdt_equal_x:
            abs_error=find_error_dxdt_result(x_list[-1], t_end)

        list_error[i]=abs_error
        #print(x_list[-1])
    if timer == True:
        end = time.time()
        print("Time taken for %s = %f sec" % (Method, end - start))
    #Therefore , the plot function wont take into account of the method
    if plot == True:
        plot_loglog(h,list_error,Method)
    return list_error

def find_error_dxdt_result( x ,t ):
    abs_error = abs(x-np.exp(t))
    return abs_error



def plot_loglog(list_h,list_error , Method):
    Method=str.lower(Method) #to ignore error cause by upper_case input
    if Method == 'euler':
        plt.loglog(list_h, list_error, label ="Euler's Method")
    elif Method == 'rk4':
        plt.loglog(list_h, list_error, label ="RK4's Method")

    plt.legend()
    plt.xlabel('h')
    plt.ylabel('Absolute Error')

def plot_x_y( x_list ,Method):

    plt.plot(x_list[:,0], x_list[:,1], label=Method)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

def plot_x_t( x_list ,t,Method ,xlabel):

    plt.plot(t, x_list, label=Method)
    plt.xlabel('t')
    plt.ylabel(xlabel)
    plt.show()

def same_error(list_error_1,list_error_2 , h , tol=1e-5):
    error=abs(list_error_1-list_error_2)
    error_min_i = np.where(error == np.amin(error))
    #error_stat = (error < tol)
    output = None

    print(error_min_i)
    if error[error_min_i] < tol :  #At least one value smaller than tol
        index = np.where(error < tol)
        output=h[index]
    else:
        print('no step-sizes that give same error in the two method')
    return output


"""
This section is used to answer Q3
"""
def main():
    list_h=np.logspace(-1 , -6 , 100 ) #Therefore, h will be always smaller than the default timestep
    list_error_rk4 = error_depend_deltat_t(dxdt_equal_x,1 ,1 ,list_h , 'rk4', plot=True , timer = True)
    list_error_euler = error_depend_deltat_t(dxdt_equal_x,1  ,1 ,list_h ,  'euler', plot=True , timer = True)
    plt.show()
    print(same_error(list_error_rk4,list_error_euler , list_h))
    n = 500
    euler_output,t=solve_ode(d2xdt2_equals_minus_x,[1,1]  ,100 ,   'euler' , n=n , timer=True)
    rk4_output,t=solve_ode(d2xdt2_equals_minus_x,[1,1]  ,100 ,   'rk4' , n=n  ,timer=True )
    true_answer=d2xdt2_equals_minus_x_true(t)
    #plot_xlabel_ylabel('t', x, x , y ,Method)
    plt.subplot(2, 1, 1)
    plt.plot(t, true_answer[0], alpha=0.15 , lw=5, label='True')
    plt.plot(t, euler_output[:,0],'--', label='Euler')
    plt.plot(t, rk4_output[:,0],':', label='RK4')

    plt.legend()
    plt.xlabel('t')
    plt.ylabel('x')

    plt.subplot(2, 1, 2)
    plt.plot(t, true_answer[1], alpha=0.15 , lw=5, label='True')
    plt.plot(t, euler_output[:,1],'--', label='Euler')
    plt.plot(t, rk4_output[:,1],':', label='RK4')


    plt.legend()
    plt.xlabel('t')
    plt.ylabel('y (dx/dt)')

    plt.show()

    plt.plot(true_answer[1], true_answer[0], alpha=0.15 , lw=5, label='True')
    plt.plot(euler_output[:,1], euler_output[:,0],'--', label='Euler')
    plt.plot(rk4_output[:,1], rk4_output[:,0],':', label='RK4')
    plt.legend()
    plt.xlabel('y (dx/dt)')
    plt.ylabel('x')
    plt.show()

if __name__ == "__main__":
    main()