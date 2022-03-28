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