import numpy as np
def dxdt_equal_x(x,t):
    dxdt= x
    return dxdt

def d2xdt2_equals_minus_x( r , t):
    """
    Function defining system of  ODEs dxdt = y, dy/dt = -x
    :param t: t value
    :param r: array [x, y]
    :return: returns value of dxdt and dy/dt at (t,u)
    """
    dxdt = r[1]
    dydt = -r[0]

    return np.array([dxdt, dydt])
def d2xdt2_equals_minus_x_true(t):
    """
    The general solution for 2nd order ODE d2xdt2_equals_minus_x()
    :return: Return the true value of x and y at a given time
    """
    x = np.sin(t) + np.cos(t)
    y = np.cos(t) - np.sin(t)
    return np.array([x, y])

def hopf_bifurcation_normal(u,t,*args):
    #print(u)
    #print(t)
    beta = args[0]
    sigma = args[1]
    du1dt = beta*u[0] - u[1] + sigma*u[0]*(u[0]**2+u[1]**2)
    du2dt = u[0] + beta*u[1] + sigma*u[1]*(u[0]**2+u[1]**2)
    return np.array([du1dt , du2dt])

def hopf_bifurcation_normal1(u,t,*args):
    #print(u)
    #print(t)
    beta = args[0]
    du1dt = beta*u[0] - u[1] - u[0]*(u[0]**2+u[1]**2)
    du2dt = u[0] + beta*u[1] - u[1]*(u[0]**2+u[1]**2)
    return np.array([du1dt , du2dt])

#when sigma is -1
def hopf_bifurcation_general(u,*args):
    beta = args[0]
    theta = 1
    u1_t= math.sqrt(beta)*math.cos(t+theta)
    u2_t = math.sqrt(beta) * math.sin(t + theta)
    return  np.array([u1_t,u2_t])

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

def modified_hopf_bifurcation_normal(u,t,beta):
    du1dt = beta*u[0] - u[1] + u[0]*(u[0]**2+u[1]**2) - u[0]*(u[0]**2+u[1]**2)**2
    du2dt = u[0] + beta*u[1] + u[1]*(u[0]**2+u[1]**2) - u[1]*(u[0]**2+u[1]**2)**2
    return np.array([du1dt , du2dt])
