#importing packages
import numpy as np
import math


# definition of parameters
# f = the func
def euler_step( f ,x, t , h ):
    """
     single Euler step function 
        Parameter:
            f(function) = function used ( need to define more detail)
            x_current : the current value of x
            t_current : the current timestep
            h : interval between step (step size) (assumed to be constant for the sake of simplicity) is then given by t_{n+1} = h + t_n
            
            
    :return: 
    [ x_{n+1} , t_{n+1} ] or [x_new , t_new] , which it the x and t for the next step
    """
    x_new = x + h * f(x,t)
    t_new = t + h
    return [x_new , t_new]

def RK4_step( f ,x, t , h ):
    k1 = f(x,t)
    k2 = f( x+h*k1/2 , t+h/2)
    k3 = f( x+h*k2/2 , t+h/2)
    k4 = f( x+h*k3, t+h)
    x_new = x + (1/6)*h*(k1+2*k2+2*k3+k4)
    t_new = t + h
    return [x_new , t_new]

def solve_to( f ,x , t , h , deltat_max , Method):
    Method=str.lower(Method)
    if deltat_max < h :
        if Method == 'euler':
            while deltat_max < h :
                [x,t]=euler_step( f ,x , t , deltat_max )
                h -= deltat_max  #in order to get out of loop
                #print(t)
                #print(x)
        elif Method == 'rk4':
            while deltat_max < h :
                [x, t] = RK4_step(f, x, t, deltat_max)
                h -= deltat_max  # in order to get out of loop

    else:
        if Method == 'euler':
            [x, t] = euler_step(f, x, t, h)
        elif Method == 'rk4':
            [x, t] = RK4_step(f, x, t, h)


    return x
def solve_ode( f ,x , t_0 , t_end , n , deltat_max , Method):
    """
    single Euler step function
        Parameter:
            f(function) = function used ( need to define more detail)
            x_current : the current value of x
            t_0 : the initial timestep
            t_end : the target timestep
            n = number of step to reach the target timestep

    :return:
    """

    t = np.linspace(t_0, t_end, n)

    x_list = np.zeros(n)
    x_list[0]=x

    for i in range(n-1):
        h=t[i+1]-t[i]
        x_list[i+1]=solve_to(f, x_list[i], t[i], h, deltat_max, Method)

    return x_list








def dxdt(x,t):
    dxdt= x
    return dxdt

#print(solve_to(dxdt,1 ,0 ,1, 0.001 , 'rk4' ))
x=solve_ode(dxdt,1 ,0 ,1 ,100, 0.0001 , 'euler')
#[x1,t1]=euler_step( dxdt ,1 , 0.5 , 0.1 )
#print(RK4_step( dxdt ,1 , 0.5 , 0.1 ))
print(x)