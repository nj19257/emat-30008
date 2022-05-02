import math
import unittest
import numpy as np
from pde_solver import finite_difference_method ,u_I
from numerical_shooting import shooting
from arbitary_differential_equation import *


class test_shooting(unittest.TestCase):
    def test_wrong_type_parameters(self):
        try:
            finite_difference_method(50.5, 1000.1, 2, 0.5, 1, "forward",'neumann', u_I )
            print('finite_difference_method function unable to detect error of parameters mx and mt not in the form of integer: test failed')
        except ValueError:
            print('finite_difference_method function able to detect error of parameters mx and mt not in the form of integer: test passed')

    def test_finite_difference_method_not_converge(self):
        try:
            finite_difference_method(500, 1000, 2, 0.5, 1, "forward", 'neumann', u_I)
            print('finite_difference_method function unable to raise error when it dont converge: test failed')
        except ValueError:
            print('finite_difference_method function able to raise error when it dont converge: test passed')
    def test_wrong_input_method(self):
        try:
            finite_difference_method(50, 1000, 2, 0.5, 1, "wrong",'neumann', u_I )
            print('finite_difference_method function unable to detect error of input method which is not valid: test failed')
        except UnboundLocalError:
            print('finite_difference_method function able to detect error of input method which is not valid: test passed')

    def test_wrong_input_boundary_condition(self):
        try:
            finite_difference_method(50, 1000, 2, 0.5, 1, "forward",'wrong', u_I )
            print('finite_difference_method function unable to detect error of input boundary_condition which is not valid: test failed')
        except UnboundLocalError:
            print('finite_difference_method function able to detect error of input boundary_condition which is not valid: test passed')

    def test_undefine_ODE_input(self):
        try:
            finite_difference_method(50, 1000, 2, 0.5, 1, "forward",'wrong', wrong )
            print('finite_difference_method function unable to detect error of undefined Initial condition function as input: test failed')
        except NameError:
            print('finite_difference_method function able to detect error of undefined Initial condition function as input: test passed')

    def test_accaracy_forward_with_neumann(self):
        forward_with_neumann_true =0.01827892
        forward_with_neumann = finite_difference_method(50, 1000, 2, 0.5, 1, "forward",'neumann', u_I )[1][-2]
        error = abs(forward_with_neumann - forward_with_neumann_true)
        if np.all(error < 1e-8) == True:
            print('Forward method in finite_difference_method with neumann as boundary condition is resulting in reliable value : test passed')
        else:
            print('RK4 method in solve_ode is resulting with neumann as boundary condition in unreliable value : test failed')
    def test_accaracy_backward_with_neumann(self):
        backward_with_neumann_true =0.03035395
        backward_with_neumann = finite_difference_method(50, 1000, 2, 0.5, 1, "backward",'neumann', u_I )[1][-2]
        error = abs(backward_with_neumann - backward_with_neumann_true)
        if np.all(error < 1e-8) == True:
            print('Backward method in finite_difference_method with neumann as boundary condition is resulting in reliable value : test passed')
        else:
            print('Backward method in finite_difference_method with neumann as boundary condition is resulting in unreliable value : test failed')
    def test_accaracy_crank_with_neumann(self):
        backward_with_neumann_true =0.02426942
        backward_with_neumann = finite_difference_method(50, 1000, 2, 0.5, 1, "crank",'neumann', u_I )[1][-2]
        error = abs(backward_with_neumann - backward_with_neumann_true)
        if np.all(error < 1e-8) == True:
            print('Crank method in finite_difference_method with neumann as boundary condition is resulting in reliable value : test passed')
        else:
            print('Crank method in finite_difference_method with neumann as resulting in unreliable value : test failed')


if __name__ == "__main__":
    unittest.main()