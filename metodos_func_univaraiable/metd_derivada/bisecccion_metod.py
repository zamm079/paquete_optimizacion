import numpy as np


def central_difference_1(f, x, delta_x):
    return (f(x + delta_x) - f(x - delta_x)) / (2 * delta_x)

def bisection_method(funcion, a, b, epsilon, delta_x):
    x1 = a
    x2 = b
    max_iter=10000
    if (central_difference_1(funcion, a, delta_x) < 0) and (central_difference_1(funcion, b, delta_x) > 0):
        epsilon = epsilon
    else:
        raise ValueError("La función no cumple con la condición")
    
    iteraciones = 0

    while abs(x1 - x2) > epsilon and iteraciones < max_iter:
        z = (x1 + x2) / 2
        f_z = central_difference_1(funcion, z, delta_x)

        if abs(f_z) <= epsilon:
            return z, z 

        if f_z < 0:
            x1 = z
        else:
            x2 = z

        iteraciones += 1

    return [x1, x2]