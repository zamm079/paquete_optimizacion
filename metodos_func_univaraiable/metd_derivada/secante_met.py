import numpy as np


def delta_x_func(x):
    return 0.01 * abs(x) if abs(x) > 0.01 else 0.0001

def central_difference_1(f, x, delta_x):
    return (f(x + delta_x) - f(x - delta_x)) / (2 * delta_x)

def secant_method(funcion, a, b, epsilon, delta_x):
    x1 = a
    x2 = b
    max_iter=10000
    if (central_difference_1(funcion, a, delta_x) < 0) and (central_difference_1(funcion, b, delta_x) > 0):
        epsilon = epsilon
    else:
        raise ValueError("La función no cumple con la condición")
    
    i = 0

    while abs(x1 - x2) > epsilon and i < max_iter:
        z = x2 - (central_difference_1(funcion, x2, delta_x)/((central_difference_1(funcion, x2, delta_x)-central_difference_1(funcion, x1, delta_x))/(x2-x1)))
        f_z = central_difference_1(funcion, z, delta_x)

        if abs(f_z) <= epsilon:
            return z, z 

        if f_z < 0:
            x1 = z
        else:
            x2 = z

        i += 1

    return [x1, x2] 


