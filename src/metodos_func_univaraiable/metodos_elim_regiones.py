import numpy as np



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




def interval_halving(a,b,e,funcion):
    #step 1
    Xm=(a+b)/2
    L=b-a

    while abs(L) > e:
        # step 2
        x1 = a + L /4 #a
        x2 = b - L /4 #b
        fx1 = funcion(x1)
        fx2 = funcion(x2)
        Fxm = funcion(Xm)
        #step 3
        if fx1 < Fxm:
            b = Xm
            Xm = x1
        else:
            if fx2 < Fxm:
                a = Xm
                Xm=x2
            else:
                a = x1
                b = x2
        
        L = b-a
        if abs(L) < e:
            return [x1,x2]