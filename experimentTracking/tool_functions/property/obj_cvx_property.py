"""
## Check the convexity of the objective function
"""

import numpy as np
import scipy.linalg
import scipy.optimize as opt
import matplotlib.pyplot as plt
import jax.numpy as jnp
from jax import grad, hessian, jit


def generate_TV_prob(n, p, qy, qd, num_instance):
    A, B, B2, E, C, D = generate_system_matrices(n, p, qy, qd)
    rho = 0.5

    dx, dy = generate_disturbances(qd, num_instance)

    sensMat, sensMat2, dxMat = define_steady_state_map(A, B, B2, E, C)
    sensMatIA = set_inexact_sensitivities(sensMat)

    quadMat, linearMat = generate_coefficient_matrices(p, num_instance)

    lb = -2 * np.random.rand(p, 1)
    ub = 2 * np.random.rand(p, 1)

    uOpt, valOpt, hessians = obtain_optimal_solutions_TV(n, p, num_instance, sensMat, sensMat2, dxMat, dx, dy, D,
                                                         quadMat, linearMat, lb, ub, rho)

    plt.figure(1)
    plt.plot(range(1, num_instance + 1), np.linalg.norm(uOpt, axis=0))
    plt.xlabel('Index of the instance')
    plt.ylabel('Norm of the optimal solution')
    plt.grid(True, which='both')

    plt.figure(2)
    plt.plot(range(1, num_instance + 1), valOpt)
    plt.xlabel('Index of the instance')
    plt.ylabel('Optimal value')
    plt.grid(True, which='both')

    problem = {
        'A': A, 'B': B, 'B2': B2, 'E': E, 'C': C, 'D': D, 'rho': rho, 'dx': dx, 'dy': dy,
        'quadMat': quadMat, 'linearMat': linearMat,
        'sensMat': sensMat, 'sensMat2': sensMat2, 'sensMatIA': sensMatIA, 'dxMat': dxMat,
        'uOpt': uOpt, 'valOpt': valOpt, 'lb': lb, 'ub': ub, 'numInstance': num_instance
    }

    return problem


def generate_system_matrices(n, p, qy, qd):
    unscaledA = np.random.randn(n, n)
    unscaledA = np.triu(unscaledA) + np.triu(unscaledA, 1).T
    A = 0.05 / max(abs(np.linalg.eigvals(unscaledA))) * unscaledA
    B = 0.5 * np.random.randn(n, p)
    B2 = 0.5 * np.random.randn(n, p)
    E = 0.5 * np.random.randn(n, qd)
    C = 0.5 * np.random.randn(qy, n)
    D = 0.5 * np.random.randn(qy, qd)
    return A, B, B2, E, C, D


def generate_disturbances(qd, num_instance):
    dx = 2e-2 * np.random.randn(qd, num_instance)
    dy = 2e-2 * np.random.randn(qd, num_instance)
    return dx, dy


def define_steady_state_map(A, B, B2, E, C):
    xssMat = np.eye(A.shape[0]) - A
    sensMat = C @ np.linalg.solve(xssMat, B)
    sensMat2 = C @ np.linalg.solve(xssMat, B2)
    dxMat = C @ np.linalg.solve(xssMat, E)
    return sensMat, sensMat2, dxMat


def set_inexact_sensitivities(sensMat):
    sensErr = 0.15 * np.mean(sensMat)
    return sensMat + 2 * sensErr * np.random.rand(*sensMat.shape) - sensErr


def generate_coefficient_matrices(p, num_instance):
    quadMat = np.zeros((p, p, num_instance))
    linearMat = np.zeros((1, p, num_instance))

    P = 0.5 * np.random.randn(p, p)
    M1 = P @ P.T
    M2 = 0.5 * np.random.randn(1, p)

    for t in range(num_instance):
        P = 1e-1 * np.random.randn(p, p)
        perturb1 = P @ P.T
        perturb2 = 2e-2 * np.random.randn(1, p)
        quadMat[:, :, t] = M1 + perturb1
        linearMat[:, :, t] = M2 + perturb2

    return quadMat, linearMat


def obtain_optimal_solutions_TV(n, p, num_instance, sensMat, sensMat2, dxMat, dx, dy, D, quadMat, linearMat, lb, ub,
                                rho):
    uOpt = np.zeros((p, num_instance))
    valOpt = np.zeros(num_instance)
    hessians = []

    def PhiTilde(u, M1Cur, M2Cur, yss):
        return u.T @ M1Cur @ u + M2Cur @ u + 0.5*jnp.linalg.norm(yss(u)) ** 2

    for t in range(num_instance):
        M1Cur = quadMat[:, :, t]
        M2Cur = linearMat[:, :, t]

        yss = lambda u: sensMat @ u + rho * sensMat2 @ (jnp.sin(u) + u ** 2) + dxMat @ dx[:, t] + D @ dy[:, t]

        Phi_hessian = hessian(lambda u: PhiTilde(u, M1Cur, M2Cur, yss))

        result = opt.minimize(lambda u: PhiTilde(u, M1Cur, M2Cur, yss),
                              jnp.zeros(p), bounds=[(lb[i], ub[i]) for i in range(p)])

        uOpt[:, t] = result.x
        valOpt[t] = result.fun

        # Compute Hessian at the optimal solution
        H_optimal = Phi_hessian(result.x)
        H_lb = Phi_hessian(lb.flatten())
        H_ub = Phi_hessian(ub.flatten())

        eigvals_opt = check_hessians(H_optimal)
        eigvals_lb = check_hessians(H_lb)
        eigvals_ub = check_hessians(H_ub)
        hessians.append(H_optimal)

    print('Finished solving all problem instances')
    return uOpt, valOpt, hessians


def check_hessians(H):
    eigvals = np.linalg.eigvals(H)
    is_positive_definite = np.all(eigvals >= 0)
    # print("Hessian Matrix:\n", H)
    print("Eigenvalues of Hessian:\n", eigvals)
    print("Is Hessian Positive Semi-definite?:", is_positive_definite)
    return eigvals


if __name__ == '__main__':
    seed = 0
    np.random.seed(seed)
    n = 30
    p = 15
    qy = 10
    qd = 10
    num_instance = 1
    problem = generate_TV_prob(n, p, qy, qd, num_instance)
    print(problem['valOpt'])
