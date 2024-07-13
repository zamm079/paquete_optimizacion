import numpy as np
import matplotlib.pyplot as plt


def fibonacci(n):
    if n == 0:
        return 1
    elif n == 1:
        return 1
    else:
        a, b = 1, 1
        for _ in range(2, n + 1):
            a, b = b, a + b
        return b

def fibonacci_search(a,b,n,funcion):
    #step 1
    L = b-a
    k=2
    while k != n:
        #step 2
        Lk = fibonacci(n-k+1)/fibonacci(n+1)
        x1 = a + Lk
        x2 = b - Lk
        #step 3
        fx1 = funcion(x1)
        fx2 = funcion(x2)
        if fx1 > fx2:
            a = x1
        if fx1 < fx2:
            b = x2
        if fx1 == fx2:
            a = x1
            b = x2
        #step 4
        k = k+1
    return [x1,x2]