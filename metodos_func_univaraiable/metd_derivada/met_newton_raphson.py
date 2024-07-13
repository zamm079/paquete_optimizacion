import numpy as np


def central_difference_1(f, x, delta_x):
    return (f(x + delta_x) - f(x - delta_x)) / (2 * delta_x)

def central_difference_2(f, x, delta_x):
    return (f(x + delta_x) - (2 * f(x)) + f(x - delta_x)) / (delta_x ** 2)

def newton_raphson_method(funcion, i_guess, delta_fun, epsilon):
    x = i_guess
    k = 1
    max_iter=10000
    while k < max_iter:
        #step1
        delta_x = delta_fun(x)
        f_derivada1 = central_difference_1(funcion, x, delta_x)
        #step2
        f_derivada2= central_difference_2(funcion, x, delta_x)
        
        if abs(f_derivada1) < epsilon:
            return x
        #step 3
        x_k1 = x - f_derivada1 / f_derivada2
        #step 4
        if abs(x_k1 - x) < epsilon:
            return x_k1
        
        x = x_k1
        k += 1
    
    return x


