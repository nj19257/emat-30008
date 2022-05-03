import numpy as np
from scipy.linalg import solve


def jacobian_matrix(f, x, tol=1e-8, *args):
    """
    Numerical approximation of the Jacobian of the given function at the point x.
    :param f: Function to have its Jacobian approximated.
    :param x: Values to approximate Jacobian at
    :param tol: Tolerance
    :param args: Any additional args to be passed to function
    :return: Returns numerical approximation of the Jacobian for a given function at x.
    """
    n = len(x)
    function = f(x, *args)
    jacobian_matrix = np.zeros((n, n))
    for i in range(n):
        Dxj = (abs(x[i]) * tol) if x[i] != 0 else tol
        x_new = [(xi if k != i else xi + Dxj) for k, xi in enumerate(x)]
        jacobian_matrix[:, i] = (f(x_new, *args) - function) / Dxj
    return jacobian_matrix


def newton(g, x0, *args):
    J =jacobian_matrix(g, x0, 1e-8, *args)
    x1_minus_x0 = solve(J, -g(x0, *args))
    x1 = x1_minus_x0 + x0
    while not np.allclose(x1, x0, atol=1e-12):
        x0=x1
        J = jacobian_matrix(g, x0, 1e-8,  *args)
        x1_minus_x0 = solve(J, -g(x0, *args))
        x1 = x1_minus_x0 + x0
    return x1