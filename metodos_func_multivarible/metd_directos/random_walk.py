import numpy as np


def met_random_walk(funcion,x0,epsilon,max_iter):

    def gen_aleatorio(xk):
        return xk + np.random.uniform(-epsilon,epsilon,size=xk.shape)

    x_mejor = x0
    xk = x0
    iteraciones = 0
    while iteraciones < max_iter:
        xk1 = gen_aleatorio(xk)
        if funcion(xk1) < funcion(x_mejor):
            x_mejor = xk1
        xk = xk1
        iteraciones += 1

    return x_mejor

