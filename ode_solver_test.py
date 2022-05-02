import math
import unittest
import numpy as np
from solve_ode import solve_ode
from numerical_shooting import shooting
from arbitary_differential_equation import *

class test_shooting(unittest.TestCase):
    def test_wrong_dimensions_initial_values(self):
        try:
            solve_ode(d2xdt2_equals_minus_x,[1,1,1]  ,1000 ,   'euler' , n=100 )
            print('solve_ode function unable to detect error of initial values with wrong dimensions: test failed')
        except ValueError:
            print('solve_ode function able to detect error of initial values with wrong dimensions: test passed')


    def test_wrong_input_method(self):
        try:
            solve_ode(d2xdt2_equals_minus_x, [1, 1], 1000, 'euler1', n=100)
            print('solve_ode function unable to detect error of input method which is not valid: test failed')
        except UnboundLocalError:
            print('solve_ode function able to detect error of input method which is not valid: test passed')

    def test_undefine_ODE_input(self):
        try:
            solve_ode(abc,[1,1,1]  ,1000 ,   'euler' , n=100 )
            print('solve_ode function unable to detect error of undefined ODE as input: test failed')
        except NameError:
            print('solve_ode function able to detect error of undefined ODE as input: test passed')

    def test_rk4_method_accaracy(self):
        rk4_true =[ 0.70259117,  1.22734088]
        rk4 = solve_ode(d2xdt2_equals_minus_x,[1,1]  ,50 ,   'rk4' , n=100)[0][-1]
        error = abs(rk4 - rk4_true)
        if np.all(error < 1e-8) == True:
            print('RK4 method in solve_ode is resulting in reliable value : test passed')
        else:
            print('RK4 method in solve_ode is resulting in unreliable value : test failed')
    def test_euler_method_accaracy(self):
        euler_true = [0.72035466, 1.25842017]
        euler = solve_ode(d2xdt2_equals_minus_x,[1,1]  ,50 ,   'euler' , n=100)[0][-1]
        error = abs(euler - euler_true)
        if np.all(error < 1e-8) == True:
            print('Forward-Euler method in solve_ode is resulting in reliable value : test passed')
        else:
            print('Forward-Euler method in solve_ode is resulting in unreliable value : test failed')

if __name__ == "__main__":
    unittest.main()