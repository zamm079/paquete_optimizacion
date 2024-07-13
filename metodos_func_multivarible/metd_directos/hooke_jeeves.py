import numpy as np

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


