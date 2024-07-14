import numpy as np
import math


def hooke_jeeves(func, x0, step_size=0.5, step_reduction=0.5, tolerance=1e-6, max_iterations=1000):
    n = len(x0)
    x = np.array(x0)
    best = np.copy(x)
    step = np.full(n, step_size)

    def explore(base_point, step_size):
        new_point = np.copy(base_point)
        for i in range(n):
            for direction in [1, -1]:
                candidate = np.copy(new_point)
                candidate[i] += direction * step_size[i]
                if func(candidate) < func(new_point):
                    new_point = candidate
                    break
        return new_point

    iteration = 0
    while np.max(step) > tolerance and iteration < max_iterations:
        new_point = explore(x, step)
        if func(new_point) < func(x):
            best = new_point + (new_point - x)
            x = new_point
        else:
            step = step * step_reduction
        iteration += 1
        # print(f"Iteration {iteration}, x: {x}, f(x): {func(x)}")

    return x




def nelder_mead(func, x_start, tol=1e-6, max_iter=1000):
    # Parámetros del algoritmo
    alpha = 1.0
    gamma = 2.0
    rho = 0.5
    sigma = 0.5

    # Inicializar el simplex
    n = len(x_start)
    simplex = np.zeros((n + 1, n))
    simplex[0] = x_start
    for i in range(n):
        y = np.array(x_start, copy=True)
        y[i] += 0.05 if x_start[i] == 0 else 0.05 * x_start[i]
        simplex[i + 1] = y

    # Evaluar la función en los vértices del simplex
    f_values = np.apply_along_axis(func, 1, simplex)
    iter_count = 0
    
    while iter_count < max_iter:
        # Ordenar el simplex por los valores de la función
        indices = np.argsort(f_values)
        simplex = simplex[indices]
        f_values = f_values[indices]

        # Centroid de los mejores n puntos
        centroid = np.mean(simplex[:-1], axis=0)

        # Reflejar
        xr = centroid + alpha * (centroid - simplex[-1])
        fxr = func(xr)

        if fxr < f_values[0]:
            # Expansión
            xe = centroid + gamma * (xr - centroid)
            fxe = func(xe)
            if fxe < fxr:
                simplex[-1] = xe
                f_values[-1] = fxe
            else:
                simplex[-1] = xr
                f_values[-1] = fxr
        else:
            if fxr < f_values[-2]:
                simplex[-1] = xr
                f_values[-1] = fxr
            else:
                # Contracción
                xc = centroid + rho * (simplex[-1] - centroid)
                fxc = func(xc)
                if fxc < f_values[-1]:
                    simplex[-1] = xc
                    f_values[-1] = fxc
                else:
                    # Reducción
                    for i in range(1, len(simplex)):
                        simplex[i] = simplex[0] + sigma * (simplex[i] - simplex[0])
                    f_values = np.apply_along_axis(func, 1, simplex)
        
        iter_count += 1

        # Verificar convergencia
        if np.max(np.abs(simplex[0] - simplex[1:])) < tol:
            break
    
    return f_values[0]





def simplex_search_meth(x,func,gama=2.0,beta=0.2,epsilon=0.001):
    # step 1
    #no cero hipervolumen
    alpha=1
    N = len(x)
    d1 = ((math.sqrt(N+1)+N-1)/N*math.sqrt(2))*alpha
    d2 = ((math.sqrt(N+1)-1)/N*math.sqrt(2))*alpha
    simplex = np.zeros((N + 1,N))
    for i in range(len(simplex)):
        for j in range(N):
            if j == i:
                simplex[i,j] = x[j]+d1
            if j != i:
                simplex[i,j] = x[j]+d2
    i_max = 10
    i = 0

    # step 2
    f_values = np.apply_along_axis(func, 1, simplex)
    xi=0
    
    while i < i_max:
        val_orden = np.argsort(f_values)
        simplex = simplex[val_orden]
        xl,xg,xh = f_values[val_orden]
        #Xc
        xc = np.mean(simplex[:-1])
        i+=1
        #step 3
        xr = 2*xc - xh
        xnew = xr
        
        if func(xr) < func(xl):
            xnew = (1+gama)*xc - (gama*xh) 
        elif func(xr) >= func(xh):
            xnew = (1-beta)*xc+(beta*xh)
        elif func(xg) < func(xr) < func(xh):
            xnew = (1+beta)*xc-(beta*xh)
        xh = xnew
        #step 4
        xi= np.sum(func(simplex))
        term1=np.sum((xi-xc)**2/(N+1))
        if term1**0.5 < epsilon:
            break
    return xnew



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

