import math
import unittest
import numpy as np
from solve_ode import solve_ode
from numerical_shooting import shooting
from arbitary_differential_equation import *


class test_shooting(unittest.TestCase):
    def test_wrong_dimensions_initial_values(self):
        try:
            shooting(hopf_bifurcation_normal, [0.01, 0.01, 3 ,4],*(2, -1))
            print('shooting function unable to detect error of initial values with wrong dimensions: test failed')
        except ValueError:
            print('shooting function able to detect error of initial values with wrong dimensions: test passed')

    def test_root_finding_not_converge(self):
        result = shooting(predator_prey, [0,0,10], *(1, 0.15, 0.1), method='newton')
        true_result = np.array([])
        if true_result.size == result.size :
            print('shooting function able to return a empty array when the numerical root finder failed: test passed')

        else:
            print('shooting function able to return a empty array when the numerical root finder failed: test failed')

    def test_wrong_input_method(self):
        try:
            result = shooting(predator_prey, [0.06861238, 0.15582958, 25.05], *(1, 0.15, 0.1), method='1234')
            print('shooting function able to return a empty array when the input method is not valid: test failed')
        except NameError:
            print('shooting function able to return a empty array when the input method is not valid: test passed')


    def test_undefine_ODE_input(self):
        try:
            shooting(abc, [0.06861238 , 0.15582958,25.05], *(1, 0.15, 0.1), method='newton')
            print('shooting function unable to detect error of undefined ODE as input: test failed')
        except NameError:
            print('shooting function able to detect error of undefined ODE as input: test passed')

    def test_newton_method_accaracy(self):
        newton = shooting(predator_prey, [0.06861238, 0.15582958, 25.05], *(1, 0.15, 0.1), method='newton')
        true_result = shooting(predator_prey, [0.06861238, 0.15582958, 25.05], *(1, 0.15, 0.1), method='fsolve')
        error = abs(newton - true_result)
        if np.all(error < 1e-8) == True:
            print('Newton root-finding method in shooting function is resulting in reliable value : test passed')
        else:
            print('Newton root-finding method in shooting function is resulting in unreliable value : test failed')

if __name__ == "__main__":
    unittest.main()
