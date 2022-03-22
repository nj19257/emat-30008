#importing packages
import numpy as np
import math
import matplotlib.pyplot as plt


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

    return x_list ,t
def find_error_dxdt( x_list ,t ):
    n=len(x_list)
    abs_error= np.zeros(n)
    for i in range(n):
        abs_error[i] = abs(x_list[i]-np.exp(t[i]))


    return abs_error

def find_error_dxdt_result( x ,t ):
    abs_error = abs(x-np.exp(t))


    return abs_error

def error_depend_deltat_t( f ,x , t_0 , t_end , n ,h, Method):
    list_error = np.zeros(len(h))
    for i in range(len(h)):
         x_list,t =solve_ode(f, x, t_0, t_end, n, h[i], Method)
         abs_error=find_error_dxdt_result(x_list[-1], t_end)
         list_error[i]=abs_error
           #print(abs_error)

    return list_error


#def test( f ,x , t_0 , t_end , n , deltat_max , Method):
   # n= math.ceil((t_end-t_0)/deltat_max)
   # x_list = np.zeros(n)
  # x_list[0]=x
  #  t=t_0
   # for i in range(n-1):
  #      x_list[i+1] = euler_step(f, x, t, deltat_max)
   # return x_list ,t
def dxdt(x,t):
    dxdt= x
    return dxdt

#print(solve_to(dxdt,1 ,0 ,1, 0.001 , 'rk4' ))
#x,t=solve_ode(dxdt,1 ,0 ,1 ,100, 0.0001 , 'euler')
#error ,h=find_error_dxdt(x,t)
#print(x)
list_h=np.logspace(-1 , -6 , 1000 )
list_error = error_depend_deltat_t(dxdt,1 ,0 ,1 ,5,list_h ,  'euler')
print(list_error)
#[x1,t1]=euler_step( dxdt ,1 , 0.5 , 0.1 )
#print(RK4_step( dxdt ,1 , 0.5 , 0.1 ))
plt.loglog(list_h,list_error,label="Euler's Method")
plt.xlabel('h')
plt.ylabel('Absolute Error')
plt.show()
print(len(list_error))
print(len(list_h))