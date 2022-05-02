import math
import unittest
import numpy as np
from numerical_continuation import numerical_continuation
from numerical_shooting import shooting
from scipy.optimize import fsolve
from arbitary_differential_equation import *
import warnings

class test_numerical_continuation(unittest.TestCase):
    def test_wrong_dimensions_initial_values(self):
        try:
            numerical_continuation(hopf_bifurcation_normal, [1.4, 0, 6.3,6], [2,-1], 0, [2, 0], 50,discretisation=shooting , solver=fsolve ,method = 'natural')
            print('numerical_continuation function unable to detect error of initial values with wrong dimensions: test failed')
        except ValueError:
            print('numerical_continuation function able to detect error of initial values with wrong dimensions: test passed')

    def test_wrong_input_method(self):
        try:
            numerical_continuation(hopf_bifurcation_normal, [1.4, 0, 6.3], [2,-1], 0, [2, 0], 50,discretisation=shooting , solver=fsolve ,method = 'wrong')
            print('numerical_continuation function able to return a empty array when the input method is not valid: test failed')
        except UnboundLocalError:
            print('numerical_continuation function able to return a empty array when the input method is not valid: test passed')


    def test_undefine_ODE_input(self):
        try:
            numerical_continuation(abc, [1.4, 0, 6.3], [2,-1], 0, [2, 0], 50,discretisation=shooting , solver=fsolve ,method = 'natural')
            print('numerical_continuation function unable to detect error of undefined ODE as input: test failed')
        except NameError:
            print('numerical_continuation function able to detect error of undefined ODE as input: test passed')

    def test_natural_method_accaracy(self):
        natural_true = [ 8.59272445e-08, -9.56612151e-23,  6.28318531e+00]
        natural = numerical_continuation(hopf_bifurcation_normal, [1.4, 0, 6.3], [2,-1], 0, [2, 0], 20,discretisation=shooting , solver=fsolve ,method = 'natural')[0][-1]
        error = abs(natural - natural_true)
        if np.all(error < 1e-8) == True:
            print('Natural method with discretisation=shooting in numerical_continuation function is resulting in reliable value : test passed')
        else:
            print('Nature method with discretisation=shooting in numerical_continuation function is resulting in unreliable value : test failed')

    def test_pseudo_arclength_method_accaracy(self):
        pseudo_arclength_true = [-1.41595553e+00, -2.97167863e-12,  6.28318531e+00]
        pseudo_arclength = numerical_continuation(hopf_bifurcation_normal, [1.4, 0, 6.3], [2,-1], 0, [2, 0], 50,discretisation=shooting , solver=fsolve ,method = 'pseudo-arclength')[0][-1]
        error = abs(pseudo_arclength - pseudo_arclength_true)
        if np.all(error < 1e-8) == True:
            print('pseudo arclength method with discretisation=shooting in numerical_continuation function is resulting in reliable value : test passed')
        else:
            print('pseudo arclength method with discretisation=shooting in numerical_continuation function is resulting in unreliable value : test failed')
    def test_natural_method_accaracy_normal(self):
        warnings.simplefilter("ignore")
        natural_true = [0.57777267, 0.57777265]
        natural = numerical_continuation(cubic, np.array([1,1]), [2], 0, [-2, 2], 200, discretisation=lambda x: x, solver=fsolve ,method = 'natural')[0][-1]
        error = abs(natural - natural_true)
        if np.all(error < 1e-8) == True:
            print('Natural method with discretisation!=shooting in numerical_continuation function is resulting in reliable value : test passed')
        else:
            print('Nature method with discretisation!=shooting in numerical_continuation function is resulting in unreliable value : test failed')

    def test_pseudo_arclength_method_accaracy_normal(self):
        pseudo_arclength_true = [-1.52273608, -1.52273608]
        pseudo_arclength = numerical_continuation(cubic, np.array([1,1]), [2], 0, [-2, 2], 200, discretisation=lambda x: x, solver=fsolve ,method = 'pseudo-arclength')[0][-1]
        error = abs(pseudo_arclength - pseudo_arclength_true)
        if np.all(error < 1e-8) == True:
            print('pseudo arclength method with discretisation!=shooting in numerical_continuation function is resulting in reliable value : test passed')
        else:
            print('pseudo arclength method with discretisation!=shooting in numerical_continuation function is resulting in unreliable value : test failed')
if __name__ == "__main__":
    unittest.main()