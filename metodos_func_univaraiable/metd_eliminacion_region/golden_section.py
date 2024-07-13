import numpy as np


def golden_section_search(funcion, a, b, epsilon):
    
    golden = 0.618
    golden2 = 1 - golden

    w1 = a + golden2 * (b - a)
    w2 = a + golden * (b - a)

    f_x1 = funcion(w1)
    f_x2 = funcion(w2)

    while b - a > epsilon:
        
        if f_x1 < f_x2:
            b = w2
            w2 = w1
            w1 = a + golden2 * (b - a)
            f_x2 = f_x1
            f_x1 = funcion(w1)
        else:
            a = w1
            w1 = w2
            w2 = a + golden * (b - a)
            f_x1 = f_x2
            f_x2 = funcion(w2)

    return [a, b]